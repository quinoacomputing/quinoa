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

#include "AbstractLinAlgPack_MatrixSymNonsingSerial.hpp"
#include "AbstractLinAlgPack_MatrixSymOpGetGMSSymMutable.hpp"
#include "AbstractLinAlgPack_MatrixOpSerial.hpp"
#include "AbstractLinAlgPack_EtaVector.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "DenseLinAlgPack_DMatrixOp.hpp"
#include "DenseLinAlgPack_DMatrixAsTriSym.hpp"
#include "AbstractLinAlgPack_LinAlgOpPackHack.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace LinAlgOpPack {
  using AbstractLinAlgPack::Vp_StMtV;
  using AbstractLinAlgPack::Mp_StMtM;
}

namespace AbstractLinAlgPack {

void MatrixSymNonsingSerial::M_StMtInvMtM(
    DMatrixSliceSym* S, value_type a, const MatrixOpSerial& B
  , BLAS_Cpp::Transp B_trans, EMatrixDummyArg ) const
{
  using BLAS_Cpp::trans;
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans_not;
  using AbstractLinAlgPack::M_StInvMtM;
  using DenseLinAlgPack::nonconst_tri_ele;
  using DenseLinAlgPack::tri_ele;
  using DenseLinAlgPack::assign;
  using LinAlgOpPack::M_StMtM;
  //
  // S = a * op(B) * inv(M) * op(B')
  //
  // We will form S won column at a time:
  //
  // S(:,j) = a * op(B) * inv(M) * op(B') * e(j)
  //
  // for j = 1 ... op(B').cols()
  //   t1 = op(B')*e(j)
  //   t2 = inv(M)*t1
  //   t3 = a*op(B)*t2
  //   S(:,j) = t3
  //
  // Above we only need to set the lower (lower triangle stored)
  // or upper (upper triangle stored) part of S(:,k)
  //
  DenseLinAlgPack::MtM_assert_sizes( rows(), cols(), no_trans
    , B.rows(), B.cols(), trans_not(B_trans) );
  DenseLinAlgPack::Mp_MtM_assert_sizes( S->rows(), S->cols(), no_trans
    , B.rows(), B.cols(), B_trans
    , B.rows(), B.cols(), trans_not(B_trans) );

  DVector t1, t2, t3; // ToDo: Use temp workspace!
  const size_type
    opBT_cols = BLAS_Cpp::cols( B.cols(), B.rows(), B_trans ),
    m         = S->rows();
  for( size_type j = 1; j <= m; ++j ) {
    EtaVector e_j(j,opBT_cols);                               // e(j)
    LinAlgOpPack::V_MtV( &t1, B, trans_not(B_trans), e_j() ); // t1 = op(B')*e(j)
    AbstractLinAlgPack::V_InvMtV( &t2, *this, no_trans, t1() ); // t2 = inv(M)*t1
    LinAlgOpPack::V_StMtV( &t3, a, B, B_trans, t2() );        // t3 = a*op(B)*t2
    Range1D
      rng = ( S->uplo() == BLAS_Cpp::upper ? Range1D(1,j) : Range1D(j,m) );
    S->gms().col(j)(rng) = t3(rng);
  }
}

// Overridden from MatrixSymNonsing

void MatrixSymNonsingSerial::M_StMtInvMtM(
  MatrixSymOp* symwo_lhs, value_type alpha
  ,const MatrixOp& mwo, BLAS_Cpp::Transp mwo_trans
  ,EMatrixDummyArg dummy
  ) const
{
  using Teuchos::dyn_cast;
  this->M_StMtInvMtM(
    &MatrixDenseSymMutableEncap(symwo_lhs)(), alpha
    ,dyn_cast<const MatrixOpSerial>(mwo), mwo_trans
    ,dummy );
}

} // end namespace AbstractLinAlgPack
