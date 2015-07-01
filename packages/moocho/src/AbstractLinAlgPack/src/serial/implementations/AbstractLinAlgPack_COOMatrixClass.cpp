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

#include <sstream>

#include "AbstractLinAlgPack_COOMatrixClass.hpp"
#include "AbstractLinAlgPack_SparseCOOReadMatrix.hpp"
#include "DenseLinAlgPack_IVector.hpp"

// Junk, test compilation
//#include "MemMngPackDef.h"
//MemMngPack::RefCount<double> ref1;

//#include "SequentialAllocatorPack.hpp"
//template SequentialAllocatorPack::SequentialAllocator<double>;

// ///////////////////////////////////////////////////////////////////////////////////
// COOMatrix

AbstractLinAlgPack::COOMatrix& AbstractLinAlgPack::COOMatrix::operator=(const COOMatrix& coom)
{
  if(this == &coom) return *this;	// assignment to self

  val_.resize(coom.nz_);	// must resize, this is why you can't use default assignment.

  // Now assign the members(this is what the default would have done).
  rows_				= coom.rows_;
  cols_				= coom.cols_;
  nz_					= coom.nz_;
  val_				= coom.val_;
  ivect_ref_			= coom.ivect_ref_;
  jvect_ref_			= coom.jvect_ref_;

  return *this;
}

void AbstractLinAlgPack::COOMatrix::resize(size_type rows, size_type cols, size_type nz)
{
   // don't resize if you don't have to to presearve row and column access just in case.
  if(rows == rows_ && cols == cols_ && nz == nz_) return;
  
  rows_ = rows;
  cols_ = cols;
  nz_ = nz;
  val_.resize(nz);
  ivect_ref_.obj().resize(nz);
  jvect_ref_.obj().resize(nz);
}

void AbstractLinAlgPack::COOMatrix::initialize(std::istream& istrm) {
  // Read COO matrix into val, ivect, jvect
  read_coo_into_valarrays(istrm,rows_,cols_,nz_,val_,ivect_ref_.obj()
    ,jvect_ref_.obj());
}

#endif // 0
