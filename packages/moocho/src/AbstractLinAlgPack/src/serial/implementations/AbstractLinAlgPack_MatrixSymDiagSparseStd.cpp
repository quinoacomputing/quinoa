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

#include "AbstractLinAlgPack_MatrixSymDiagSparseStd.hpp"
#include "AbstractLinAlgPack_SpVectorOut.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"
#include "Teuchos_Assert.hpp"

namespace AbstractLinAlgPack {

MatrixSymDiagSparseStd::MatrixSymDiagSparseStd( const SpVectorSlice& diag )
  : diag_(diag)
{}

void MatrixSymDiagSparseStd::initialize( const SpVectorSlice& diag )
{
  diag_ = diag;
}

// Overridden from MatrixOp

MatrixOp& MatrixSymDiagSparseStd::operator=(const MatrixOp& m)
{
  if(&m == this) return *this;	// assignment to self
  const MatrixSymDiagSparseStd
    *p_m = dynamic_cast<const MatrixSymDiagSparseStd*>(&m);
  if(p_m) {
    diag_ = p_m->diag_;
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::invalid_argument
      ,"MatrixSymDiagSparseStd::operator=(const MatrixOp& m) : Error!"
      "The concrete type of m = \'" << typeName(m) << "\' is not a subclass of "
      "MatrixSymDiagSparseStd as expected"
      );
  }
  return *this;
}

// Overridden from MatrixDiagonalSparse

const SpVectorSlice MatrixSymDiagSparseStd::diag() const
{
  return diag_();
}


}	// end namespace AbstractLinAlgPack
