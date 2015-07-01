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

#include "NLPInterfacePack_ExampleNLPFirstOrder.hpp"
#include "NLPInterfacePack_ExampleBasisSystem.hpp"
#include "AbstractLinAlgPack_BasisSystemComposite.hpp"
#include "AbstractLinAlgPack_MatrixSymDiagStd.hpp"
#include "AbstractLinAlgPack_MatrixComposite.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "ReleaseResource_ref_count_ptr.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Assert.hpp"

namespace NLPInterfacePack {

ExampleNLPFirstOrder::ExampleNLPFirstOrder(
  const VectorSpace::space_ptr_t&  vec_space
  ,value_type                      xo
  ,bool                            has_bounds
  ,bool                            dep_bounded
  )
  :ExampleNLPObjGrad(vec_space,xo,has_bounds,dep_bounded)
  ,initialized_(false)
{
  basis_sys_ = Teuchos::rcp(
    new ExampleBasisSystem(
      this->space_x()
      ,this->var_dep()
      ,this->var_indep()
      )
    );
}

// Overridden public members from NLPFirstOrder

void ExampleNLPFirstOrder::set_Gc(MatrixOp* Gc)
{
  if(Gc) // Throw an exception if this matrix is the wrong type!
    Teuchos::dyn_cast<AbstractLinAlgPack::MatrixComposite>(*Gc);
  NLPFirstOrder::set_Gc(Gc);
}

const NLPFirstOrder::mat_fcty_ptr_t
ExampleNLPFirstOrder::factory_Gc() const
{
  return factory_Gc_;
}

const NLPFirstOrder::basis_sys_ptr_t
ExampleNLPFirstOrder::basis_sys() const
{
  return basis_sys_;
}

// Overridden public members from NLP

void ExampleNLPFirstOrder::initialize(bool test_setup)
{
  namespace rcp = MemMngPack;

  ExampleNLPObjGrad::initialize(test_setup);
  NLPFirstOrder::initialize(test_setup);

  factory_Gc_ = BasisSystemComposite::factory_Gc();
  
  initialized_ = true;
}

bool ExampleNLPFirstOrder::is_initialized() const
{
  return initialized_;
}

// Overridden protected members from NLPFirstOrder

void ExampleNLPFirstOrder::imp_calc_Gc(
  const Vector& x, bool newx, const FirstOrderInfo& first_order_info) const
{
  namespace rcp = MemMngPack;
  using Teuchos::dyn_cast;
  using AbstractLinAlgPack::Vp_S; // Should not have to do this!

  const index_type
    n = this->n(),
    m = this->m();
  const Range1D
    var_dep   = this->var_dep(),
    var_indep = this->var_indep();

  // Get references to aggregate C and N matrices (if allocated)
  MatrixOpNonsing
    *C_aggr = NULL;
  MatrixOp
    *N_aggr = NULL;
  BasisSystemComposite::get_C_N(
    first_order_info.Gc, &C_aggr, &N_aggr ); // Will return NULLs if Gc is not initialized

  // Allocate C and N matrix objects if not done yet!
  Teuchos::RCP<MatrixOpNonsing>
    C_ptr = Teuchos::null;
  Teuchos::RCP<MatrixOp>
    N_ptr = Teuchos::null;
  if( C_aggr == NULL ) {
    const VectorSpace::space_ptr_t
      space_x  = this->space_x(),
      space_xD = space_x->sub_space(var_dep);
    C_ptr  = Teuchos::rcp(new MatrixSymDiagStd(space_xD->create_member()));
    N_ptr  = Teuchos::rcp(new MatrixSymDiagStd(space_xD->create_member()));
    C_aggr = C_ptr.get();
    N_aggr = N_ptr.get();
  }

  // Get references to concreate C and N matrices
  MatrixSymDiagStd
    &C = dyn_cast<MatrixSymDiagStd>(*C_aggr);
  MatrixSymDiagStd
    &N = dyn_cast<MatrixSymDiagStd>(*N_aggr);
  // Get x = [ x_D' x_I ]
  Vector::vec_ptr_t
    x_D = x.sub_view(var_dep),
    x_I = x.sub_view(var_indep);
  // Set the diagonals of C and N (this is the only computation going on here)
  C.diag() = *x_I;          // C.diag = x_I - 1.0
  Vp_S( &C.diag(),  -1.0 ); // ...
  N.diag() = *x_D;          // N.diag = x_D - 10.0
  Vp_S( &N.diag(), -10.0 ); // ...
  // Initialize the matrix object Gc if not done so yet
  if( C_ptr.get() != NULL ) {
    BasisSystemComposite::initialize_Gc(
      this->space_x(), var_dep, var_indep
      ,this->space_c()
      ,C_ptr, N_ptr
      ,first_order_info.Gc
      );
  }
}

}	// end namespace NLPInterfacePack
