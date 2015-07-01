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

#include "AbstractLinAlgPack_MatrixSymOpNonsing.hpp"

namespace AbstractLinAlgPack {

MatrixSymOpNonsing::mat_mswons_mut_ptr_t
MatrixSymOpNonsing::clone_mswons()
{
  return Teuchos::null;
}

MatrixSymOpNonsing::mat_mswons_ptr_t
MatrixSymOpNonsing::clone_mswons() const
{
  return Teuchos::null;
}

// Overridden from MatrixOp

MatrixSymOpNonsing::mat_mut_ptr_t
MatrixSymOpNonsing::clone()
{
  return clone_mswons();
}

MatrixSymOpNonsing::mat_ptr_t
MatrixSymOpNonsing::clone() const
{
  return clone_mswons();
}

// Overridden from MatrixNonsing

MatrixSymOpNonsing::mat_mns_mut_ptr_t
MatrixSymOpNonsing::clone_mns()
{
  return clone_mswons();
}

MatrixSymOpNonsing::mat_mns_ptr_t
MatrixSymOpNonsing::clone_mns() const
{
  return clone_mswons();
}

// Overridden from MatrixSymOp

MatrixSymOpNonsing::mat_mswo_mut_ptr_t
MatrixSymOpNonsing::clone_mswo()
{
  return clone_mswons();
}

MatrixSymOpNonsing::mat_mswo_ptr_t
MatrixSymOpNonsing::clone_mswo() const
{
  return clone_mswons();
}

// Overridden from MatrixSymNonsing

MatrixSymOpNonsing::mat_msns_mut_ptr_t
MatrixSymOpNonsing::clone_msns()
{
  return clone_mswons();
}

MatrixSymOpNonsing::mat_msns_ptr_t
MatrixSymOpNonsing::clone_msns() const
{
  return clone_mswons();
}

// Overridden from MatrixOpNonsing

MatrixSymOpNonsing::mat_mwons_mut_ptr_t
MatrixSymOpNonsing::clone_mwons()
{
  return clone_mswons();
}

MatrixSymOpNonsing::mat_mwons_ptr_t
MatrixSymOpNonsing::clone_mwons() const
{
  return clone_mswons();
}

}	// end namespace AbstractLinAlgPack
