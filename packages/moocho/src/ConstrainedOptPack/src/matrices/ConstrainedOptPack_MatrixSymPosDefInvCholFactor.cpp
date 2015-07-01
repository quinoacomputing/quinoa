#if 0

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

#include "ConstrainedOptPack_MatrixSymPosDefInvCholFactor.hpp"
#include "SymInvCholMatrixOp.hpp"
#include "AbstractLinAlgPack_SpVectorOp.hpp"
#include "DenseLinAlgPack_DVectorOp.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "DenseLinAlgPack_DMatrixOp.hpp"
#include "DenseLinAlgPack_DMatrixOut.hpp"

namespace LinAlgOpPack {

using AbstractLinAlgPack::Vp_StV;
using AbstractLinAlgPack::Vp_StMtV;
using AbstractLinAlgPack::Mp_StM;
using ConstrainedOptPack::Vp_StMtV;

}	// end namespace LinAlgOpPack

namespace ConstrainedOptPack {

// Overridden from Matrix 

size_type MatrixSymPosDefInvCholFactor::cols() const
{
  return rows();
}

// Overridden from MatrixOp

MatrixOp& MatrixSymPosDefInvCholFactor::operator=(const MatrixOp& m)
{
  return MatrixWithOpConcreteEncap<SymInvCholMatrix>::operator=(m);
}

std::ostream& MatrixSymPosDefInvCholFactor::output(std::ostream& out) const
{	return out << "\n*** Inverse Cholesky factor:\n" << m().UInv(); }

// Level-2 BLAS

void MatrixSymPosDefInvCholFactor::Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
  , const DVectorSlice& vs_rhs2, value_type beta) const
{
  ConstrainedOptPack::Vp_StMtV(vs_lhs,alpha,m(),trans_rhs1,vs_rhs2,beta);
}

void MatrixSymPosDefInvCholFactor::Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
  , const SpVectorSlice& sv_rhs2, value_type beta) const
{
  using LinAlgOpPack::assign;
  DVector vs_rhs2;
  assign(&vs_rhs2,sv_rhs2);
  ConstrainedOptPack::Vp_StMtV(vs_lhs,alpha,m(),trans_rhs1,vs_rhs2,beta);
}

value_type MatrixSymPosDefInvCholFactor::transVtMtV(const DVectorSlice& vs_rhs1, BLAS_Cpp::Transp trans_rhs2
  , const DVectorSlice& vs_rhs3) const
{
  return ConstrainedOptPack::transVtMtV(vs_rhs1,m(),vs_rhs3);
}

value_type MatrixSymPosDefInvCholFactor::transVtMtV(const SpVectorSlice& sv_rhs1, BLAS_Cpp::Transp trans_rhs2
  , const SpVectorSlice& sv_rhs3) const
{
  using LinAlgOpPack::assign;
  DVector vs_rhs1, vs_rhs3;
  assign(&vs_rhs1,sv_rhs1);
  assign(&vs_rhs3,sv_rhs3);
  return ConstrainedOptPack::transVtMtV(vs_rhs1,m(),vs_rhs3);
}

// Overridden from MatrixFactorized

void MatrixSymPosDefInvCholFactor::V_InvMtV(DVector* v_lhs, BLAS_Cpp::Transp trans_rhs1
  , const DVectorSlice& vs_rhs2) const
{
  ConstrainedOptPack::V_InvMtV(v_lhs,m(),vs_rhs2);
}

void MatrixSymPosDefInvCholFactor::V_InvMtV(DVectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
  , const DVectorSlice& vs_rhs2) const
{
  ConstrainedOptPack::V_InvMtV(vs_lhs,m(),vs_rhs2);
}

void MatrixSymPosDefInvCholFactor::V_InvMtV(DVector* v_lhs, BLAS_Cpp::Transp trans_rhs1
  , const SpVectorSlice& sv_rhs2) const
{
  ConstrainedOptPack::V_InvMtV(v_lhs,m(),sv_rhs2);
}

void MatrixSymPosDefInvCholFactor::V_InvMtV(DVectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
  , const SpVectorSlice& sv_rhs2) const
{
  ConstrainedOptPack::V_InvMtV(vs_lhs,m(),sv_rhs2);
}

value_type MatrixSymPosDefInvCholFactor::transVtInvMtV(const DVectorSlice& vs_rhs1
  , BLAS_Cpp::Transp trans_rhs2, const DVectorSlice& vs_rhs3) const
{
  return ConstrainedOptPack::transVtInvMtV(vs_rhs1,m(),vs_rhs3);
}

value_type MatrixSymPosDefInvCholFactor::transVtInvMtV(const SpVectorSlice& sv_rhs1
  , BLAS_Cpp::Transp trans_rhs2, const SpVectorSlice& sv_rhs3) const
{
  return ConstrainedOptPack::transVtInvMtV(sv_rhs1,m(),sv_rhs3);}

// Overridden from MatrixSymFactorized

void MatrixSymPosDefInvCholFactor::M_StMtInvMtM(
    DMatrixSliceSym* S, value_type a, const MatrixOp& B
  , BLAS_Cpp::Transp B_trans, EMatrixDummyArg dummy_arg ) const
{
//	// Uncomment to use the defalut implementation (for debugging)
//	MatrixSymFactorized::M_StMtInvMtM(S,a,B,B_trans,dummy_arg); return;

  using BLAS_Cpp::trans;
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans_not;
  using BLAS_Cpp::upper;
  using BLAS_Cpp::nonunit;
  using AbstractLinAlgPack::M_StInvMtM;
  using DenseLinAlgPack::tri;
  using DenseLinAlgPack::syrk;
  using DenseLinAlgPack::M_StInvMtM;
  using LinAlgOpPack::M_StMtM;
  using LinAlgOpPack::assign;

  DenseLinAlgPack::MtM_assert_sizes( rows(), cols(), no_trans
    , B.rows(), B.cols(), trans_not(B_trans) );
  DenseLinAlgPack::Mp_MtM_assert_sizes( S->rows(), S->cols(), no_trans
    , B.rows(), B.cols(), B_trans
    , B.rows(), B.cols(), trans_not(B_trans) );
  //
  // S = a * op(B) * inv(M) * op(B)'
  // 
  // M = L * L'
  // inv(M) = inv(L * L') = inv(L') * inv(L) = UInv * UInv'
  // 
  // S = a * op(B) * UInv * UInv' * op(B)'
  // 
  // T = op(B)'
  // 
  // T = UInv' * T (inplace with BLAS)
  // 
  // S = a * T' * T
  // 

  // T = op(B)'
  DMatrix T;
  assign( &T, B, trans_not(B_trans) );
  // T = UInv' * T (inplace with BLAS)
  M_StMtM( &T(), 1.0, tri(m().UInv(),upper,nonunit), trans, T(), no_trans );
  // S = a * T' * T
  syrk( trans, a, T(), 0.0, S );
}

// Overridden from MatrixSymSecant

void MatrixSymPosDefInvCholFactor::init_identity(size_type n, value_type alpha)
{
  if( alpha <= 0.0 ) {
    std::ostringstream omsg;
    omsg	<< "MatrixSymPosDefInvCholFactor::init_identity(...) : Error, alpha = " << alpha
        << " <= 0.0 and therefore this is not a positive definite matrix.";
    throw UpdateSkippedException( omsg.str() );	
  }
  m().resize(n);
  m().UInv() = 0.0;
  m().UInv().diag() = 1.0 / ::sqrt( alpha );
}

void MatrixSymPosDefInvCholFactor::init_diagonal( const DVectorSlice& diag )
{
  DVectorSlice::const_iterator
    min_ele_ptr = std::min_element( diag.begin(), diag.end() );
  if( *min_ele_ptr <= 0.0 ) {
    std::ostringstream omsg;
    omsg	<< "MatrixSymPosDefInvCholFactor::init_diagonal(...) : Error, "
        << "diag(" << min_ele_ptr - diag.begin() + 1 << " ) = "
        << (*min_ele_ptr) << " <= 0.0.\n"
        << "Therefore this is not a positive definite matrix.";
    throw UpdateSkippedException( omsg.str() );	
  }
  const size_type n = diag.size();
  m().resize(n);
  m().UInv() = 0.0;

  DVectorSlice::const_iterator
    diag_itr = diag.begin();
  DVectorSlice::iterator
    inv_fact_diag_itr = m().UInv().diag().begin();

  while( diag_itr != diag.end() )
    *inv_fact_diag_itr++ = 1.0 / ::sqrt( *diag_itr++ );
}

void MatrixSymPosDefInvCholFactor::secant_update(DVectorSlice* s, DVectorSlice* y, DVectorSlice* _Bs)
{
  using LinAlgOpPack::V_MtV;
  try {
    if(!_Bs) {
      DVector Bs;
      V_MtV( &Bs, *this, BLAS_Cpp::no_trans, *s );
      ConstrainedOptPack::BFGS_update(&m(),s,y,&Bs());
    }
    else {
      ConstrainedOptPack::BFGS_update(&m(),s,y,_Bs);
    }
  }
  catch(const BFGSUpdateSkippedException& excpt) {
    throw UpdateSkippedException( excpt.what() );
  }
}

// Overridden from MatrixExtractInvCholFactor

void MatrixSymPosDefInvCholFactor::extract_inv_chol( DMatrixSliceTriEle* InvChol ) const
{
  DenseLinAlgPack::assign( InvChol, DenseLinAlgPack::tri_ele( m().UInv(), BLAS_Cpp::upper ) );
}

// Overridden from Serializable

void MatrixSymPosDefInvCholFactor::serialize( std::ostream &out ) const
{
  TEUCHOS_TEST_FOR_EXCEPT(true);
}

void MatrixSymPosDefInvCholFactor::unserialize( std::istream &in )
{
  TEUCHOS_TEST_FOR_EXCEPT(true);
}

}	// end namespace ConstrainedOptPack

#endif // 0
