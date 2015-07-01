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

#include "MoochoPack_EvalNewPointTailoredApproachOrthogonal_Step.hpp"
#include "ConstrainedOptPack_MatrixIdentConcatStd.hpp"
#include "NLPInterfacePack_NLPDirect.hpp"
#include "AbstractLinAlgPack_MatrixComposite.hpp"
#include "AbstractLinAlgPack_MatrixSymOpNonsing.hpp"
#include "AbstractLinAlgPack_MatrixSymInitDiag.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_MatrixOpOut.hpp"
#include "AbstractLinAlgPack_AssertOp.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Assert.hpp"

namespace MoochoPack {

EvalNewPointTailoredApproachOrthogonal_Step::EvalNewPointTailoredApproachOrthogonal_Step(
  const deriv_tester_ptr_t                &deriv_tester
  ,const bounds_tester_ptr_t              &bounds_tester
  ,EFDDerivTesting                        fd_deriv_testing
  )
  :EvalNewPointTailoredApproach_Step(deriv_tester,bounds_tester,fd_deriv_testing)
{}

// protected

void EvalNewPointTailoredApproachOrthogonal_Step::uninitialize_Y_Uy(
  MatrixOp         *Y
  ,MatrixOp        *Uy
  )
{
  using Teuchos::dyn_cast;

  MatrixIdentConcatStd
    *Y_orth = Y ? &dyn_cast<MatrixIdentConcatStd>(*Y)  : NULL;
  MatrixComposite
    *Uy_cpst = Uy ? &dyn_cast<MatrixComposite>(*Uy) : NULL;			

  if(Y_orth)
    Y_orth->set_uninitialized();
  TEUCHOS_TEST_FOR_EXCEPT( !( Uy_cpst == NULL ) ); // ToDo: Implement for undecomposed equalities
}

void EvalNewPointTailoredApproachOrthogonal_Step::calc_py_Y_Uy(
  const NLPDirect       &nlp
  ,const D_ptr_t        &D
  ,VectorMutable        *py
  ,MatrixOp             *Y
  ,MatrixOp             *Uy
  ,EJournalOutputLevel  olevel
  ,std::ostream         &out
  )
{
  namespace rcp = MemMngPack;
  using Teuchos::dyn_cast;
  using LinAlgOpPack::syrk;

  const size_type
    n = nlp.n(),
    r = nlp.r();
  const Range1D
    var_dep(1,r),
    var_indep(r+1,n),
    con_decomp   = nlp.con_decomp(),
    con_undecomp = nlp.con_undecomp();

  //
  // Get pointers to concreate matrices
  //
  
  MatrixIdentConcatStd
    *Y_orth = Y ? &dyn_cast<MatrixIdentConcatStd>(*Y)  : NULL;
  MatrixComposite
    *Uy_cpst = Uy ? &dyn_cast<MatrixComposite>(*Uy) : NULL;			

  //
  // Initialize the matrices
  //

  // Y
  if(Y_orth) {
    D_ptr_t  D_ptr = D;
//		if(mat_rel == MATRICES_INDEP_IMPS) {
//			D_ptr = D->clone();
//			TEUCHOS_TEST_FOR_EXCEPTION(
//				D_ptr.get() == NULL, std::logic_error
//				,"DecompositionSystemOrthogonal::update_decomp(...) : Error, "
//				"The matrix class used for the direct sensitivity matrix D = inv(C)*N of type \'"
//				<< typeName(*D) << "\' must return return.get() != NULL from the clone() method "
//				"since mat_rel == MATRICES_INDEP_IMPS!" );
//		}
    Y_orth->initialize(
      nlp.space_x()                                     // space_cols
      ,nlp.space_x()->sub_space(var_dep)->clone()       // space_rows
      ,MatrixIdentConcatStd::BOTTOM                     // top_or_bottom
      ,-1.0                                             // alpha
      ,D_ptr                                            // D_ptr
      ,BLAS_Cpp::trans                                  // D_trans
      );
  }

  // S
  if(S_ptr_.get() == NULL) {
    S_ptr_ = nlp.factory_S()->create();
  }
  // S = I + (D)'*(D')'
  dyn_cast<MatrixSymInitDiag>(*S_ptr_).init_identity(D->space_rows());
  syrk(*D,BLAS_Cpp::trans,1.0,1.0,S_ptr_.get());

  TEUCHOS_TEST_FOR_EXCEPT( !( Uy_cpst == NULL ) ); // ToDo: Implement for undecomposed equalities

  recalc_py(*D,py,olevel,out);

}

void EvalNewPointTailoredApproachOrthogonal_Step::recalc_py(
  const MatrixOp           &D
  ,VectorMutable           *py
  ,EJournalOutputLevel     olevel
  ,std::ostream            &out
  )
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using AbstractLinAlgPack::Vp_StMtV;
  using AbstractLinAlgPack::V_InvMtV;
  using LinAlgOpPack::V_MtV;

  const MatrixSymOpNonsing   &S = *S_ptr_;

  VectorSpace::vec_mut_ptr_t               // ToDo: make workspace!
    tIa = D.space_rows().create_member(),
    tIb = D.space_rows().create_member();
  //
  // py = -inv(R)*c
  // py = -((I - D*inv(S)*D')*inv(C))*c
  //    = -(I - D*inv(S)*D')*(-py)
  //    = py - D*inv(S)*D'*py
  //
  // =>
  //
  // tIa  = D'*py
  // tIb  = inv(S)*tIa
  // py   += -D*tIb
  //
  V_MtV( tIa.get(), D, trans, *py );              // tIa  = D'*py
  V_InvMtV( tIb.get(), S, no_trans, *tIa );       // tIb  = inv(S)*tIa
  Vp_StMtV( py, -1.0, D, no_trans, *tIb );        // py   += -D*tIb
}

void EvalNewPointTailoredApproachOrthogonal_Step::print_calc_py_Y_Uy(
  std::ostream& out, const std::string& L
  ) const
{
  out
    << L << "*** Orthogonal decomposition\n"
    << L << "py = inv(I + D*D') * py <: space_range\n"
    << L << "Y = [ I ; -D' ] <: space_x|space_range\n"
    << L << "Uy = ???\n"
    ;
}

}	// end namespace MoochoPack 
