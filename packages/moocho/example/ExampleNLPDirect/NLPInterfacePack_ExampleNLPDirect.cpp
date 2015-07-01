// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include <assert.h>

#include <stdexcept>

#include "NLPInterfacePack_ExampleNLPDirect.hpp"
#include "ExampleNLPDirectRTOps.h"
#include "AbstractLinAlgPack_BasisSystemComposite.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_MatrixSymDiagStd.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "RTOpPack_RTOpC.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"

namespace {

static RTOpPack::RTOpC explnlp2_calc_py_D_op;

class init_rtop_server_t {
public:
  init_rtop_server_t() {
    TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_explnlp2_calc_py_D_construct(0,&explnlp2_calc_py_D_op.op()));
  }
}; 
init_rtop_server_t  init_rtop_server;

} // end namespace

namespace NLPInterfacePack {

ExampleNLPDirect::ExampleNLPDirect(
  const VectorSpace::space_ptr_t&  vec_space
  ,value_type                      xo
  ,bool                            has_bounds
  ,bool                            dep_bounded
  )
  :ExampleNLPObjGrad(vec_space,xo,has_bounds,dep_bounded),initialized_(false)
{
  namespace rcp = MemMngPack;

  // Create the factory object for D
  factory_D_ = Teuchos::rcp(new Teuchos::AbstractFactoryStd<MatrixOp,MatrixSymDiagStd>());
  NLPDirect::set_factories(
    Teuchos::rcp(
      new Teuchos::AbstractFactoryStd<MatrixSymOp,MatrixSymDiagStd>())               // D'*D
    ,Teuchos::rcp(
      new Teuchos::AbstractFactoryStd<MatrixSymOpNonsing,MatrixSymDiagStd>())    // S
    );
}

// Overridden public members from NLP

void ExampleNLPDirect::initialize(bool test_setup)
{

  if( initialized_ ) {
    NLPDirect::initialize(test_setup);
    ExampleNLPObjGrad::initialize(test_setup);
    return;
  }

  NLPDirect::initialize(test_setup);
  ExampleNLPObjGrad::initialize(test_setup);

  initialized_ = true;
}

bool ExampleNLPDirect::is_initialized() const
{
  return initialized_;
}

// Overridden public members from NLPDirect

Range1D ExampleNLPDirect::var_dep() const
{
  return ExampleNLPObjGrad::var_dep();
}

Range1D ExampleNLPDirect::var_indep() const
{
  return ExampleNLPObjGrad::var_indep();
}

const NLPDirect::mat_fcty_ptr_t
ExampleNLPDirect::factory_D() const
{
  return factory_D_;
}

void ExampleNLPDirect::calc_point(
  const Vector     &x
  ,value_type      *f
  ,VectorMutable   *c
  ,bool            recalc_c
  ,VectorMutable   *Gf
  ,VectorMutable   *py
  ,VectorMutable   *rGf
  ,MatrixOp        *GcU
  ,MatrixOp        *D
  ,MatrixOp        *Uz
  ) const
{
  using Teuchos::dyn_cast;
  using LinAlgOpPack::Vp_MtV;

  assert_is_initialized();

  const size_type
    n = this->n();

  // Validate the input

#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    x.dim() != n, std::invalid_argument
    ,"ExampleNLPDirect::calc_point(...), Error x.dim() = " << x.dim()
    << " != n = " << n );
  TEUCHOS_TEST_FOR_EXCEPTION(
    c && !this->space_c()->is_compatible(c->space()), std::invalid_argument
    ,"ExampleNLPDirect::calc_point(...), Error c is not compatible" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    Gf && !this->space_x()->is_compatible(Gf->space()), std::invalid_argument
    ,"ExampleNLPDirect::calc_point(...), Error, Gf is not compatible" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    py && !this->space_x()->sub_space(this->var_dep())->is_compatible(py->space()), std::invalid_argument
    ,"ExampleNLPDirect::calc_point(...), Error, py is not compatible" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    rGf && !this->space_x()->sub_space(this->var_dep())->is_compatible(rGf->space()), std::invalid_argument
    ,"ExampleNLPDirect::calc_point(...), Error, py is not compatible" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    GcU, std::invalid_argument
    ,"ExampleNLPDirect::calc_point(...), Error, there are no undecomposed equalities" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    D && !dynamic_cast<MatrixSymDiagStd*>(D), std::invalid_argument
    ,"ExampleNLPDirect::calc_point(...), Error, D is not compatible" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    py!=NULL && c==NULL, std::invalid_argument
    ,"ExampleNLPDirect::calc_point(...) : "
    "Error, if py!=NULL then c!=NULL must also be true" );
#endif

  // ///////////////////////////////////
  // Compute f(x), c(x) and Gf(x)

  typedef ExampleNLPDirect  this_t;

  // Make temp Gf if needed
  VectorSpace::vec_mut_ptr_t  Gf_ptr;
  if( rGf && !Gf ) {
    Gf_ptr = this->space_x()->create_member();
    Gf = Gf_ptr.get();
  }

  // Make temp D if needed
  mat_fcty_ptr_t::element_type::obj_ptr_t  D_ptr;
  if( rGf && !D ) {
    D_ptr = this->factory_D()->create();
    D = D_ptr.get();
  }

  // Remember what references are already set
  value_type             *f_saved  = const_cast<this_t*>(this)->get_f();
  VectorMutable    *c_saved  = const_cast<this_t*>(this)->get_c();
  VectorMutable    *Gf_saved = const_cast<this_t*>(this)->get_Gf();
  // Set and compute the quantities
  const_cast<this_t*>(this)->set_f(f);
  const_cast<this_t*>(this)->set_c(c);
  const_cast<this_t*>(this)->set_Gf(Gf);
  if(Gf)
    this->calc_Gf(x,true);
  if(c && recalc_c)
    this->calc_c(x,false);
  if(f)
    this->calc_f(x,false);
  // Reset the remembered references
  const_cast<this_t*>(this)->set_f(f_saved);
  const_cast<this_t*>(this)->set_c(c_saved);
  const_cast<this_t*>(this)->set_Gf(Gf_saved);

  // ////////////////////////////////////////////////////////////////////////
  // Compute py = -inv(C)*c and/or D = -inv(C)*N at the same time
  // 
  //				[ 1/(1-x(m+1))											]
  //				[				1/(1-x(m+2))							]
  //	-inv(C)	= 	[								.						]
  //				[									.					]
  //				[										1/(1-x(m+m))	]
  //
  //
  //				[ x(1) - 10												]
  //				[				x(2) - 10								]
  //	N 		= 	[								.						]
  //				[									.					]
  //				[										x(m) - 10		]

  MatrixSymDiagStd
    *D_diag = dynamic_cast<MatrixSymDiagStd*>(D);

  Vector::vec_ptr_t
    xD= x.sub_view(this->var_dep()),
    xI = x.sub_view(this->var_indep());

  int task;
  if     ( py  && !D )  task = 0;
  else if( !py &&  D )  task = 1;
  else if( py  &&  D )  task = 2;
  
  TEUCHOS_TEST_FOR_EXCEPT(0!=RTOp_TOp_explnlp2_calc_py_D_set_task(task,&explnlp2_calc_py_D_op.op()));

  const int              num_vecs = task < 2 ? 2 : 3;
  const Vector*          vecs[3] = { NULL, NULL, NULL };
  const int              num_targ_vecs = task < 2 ? 1 : 2;
  VectorMutable*         targ_vecs[2] = { NULL, NULL };

  // targ_vecs[0] will have apply_op(...) called on it.
  if(D) {
    D_diag->init_identity( *this->space_c(), 0.0 );
    targ_vecs[0]= &D_diag->diag();
  }
  else if(py)
    targ_vecs[0] = py;
  else
    TEUCHOS_TEST_FOR_EXCEPT(true); // Only local error?
  // targ_vecs[1] will be passed to apply_op(...)
  if(py && D)
    targ_vecs[1] = py;
  
  // vecs[...]
  int k = 0;
  vecs[k] = xD.get(); ++k;
  if(D)  { vecs[k] = xI.get(); ++k; }
  if(py) { vecs[k] = c;        ++k; }

  AbstractLinAlgPack::apply_op(
    explnlp2_calc_py_D_op, num_vecs, vecs, num_targ_vecs, num_targ_vecs?targ_vecs:NULL
    ,NULL
    );

  // rGf = Gf(var_indep) + D' * Gf(var_dep)
  if(rGf) {
    *rGf = *Gf->sub_view(this->var_indep());
    Vp_MtV( rGf, *D, BLAS_Cpp::trans, *Gf->sub_view(this->var_dep()) );
  }
}

void ExampleNLPDirect::calc_semi_newton_step(
  const Vector    &x
  ,VectorMutable  *c
  ,bool           recalc_c
  ,VectorMutable  *py
  ) const
{
  // In this case just call calc_point(...).
  // In a more specialized application, this could be much cheaper!
  calc_point(x,NULL,c,recalc_c,NULL,py,NULL,NULL,NULL,NULL);
}

}	// end namespace NLPInterfacePack
