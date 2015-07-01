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

#ifndef ETA_VECTOR_H
#define ETA_VECTOR_H

#include "AbstractLinAlgPack_SpVectorClass.hpp"

namespace AbstractLinAlgPack {

/** \brief Create an eta vector (scaled by alpha = default 1).
  *
  * The created vector is of size n and has the single nonzero
  * element of eta(i) = alpha.
  * 
  * The default constructor and assignment functions are not
  * allowed.
  */
class EtaVector {
public:
  
  typedef SpVectorSlice::element_type		ele_t;


  /** \brief . */
  EtaVector( ele_t::index_type i, size_type n, ele_t::value_type alpha = 1.0 )
    : ele_(i,alpha), sv_(&ele_,1,0,n,true)
  {}

  /// Implicit conversion to a SpVectorSlice object.
  operator const SpVectorSlice() const
  {
    return sv_;
  }

  /// Explicit conversion to a SpVectorSlice object.
  const SpVectorSlice& operator()() const
  {
    return sv_;
  }

private:
  ele_t			ele_;
  SpVectorSlice	sv_;

  // not defined and not to be called
  EtaVector();
  EtaVector& operator=(const EtaVector&);

};	// end class EtaVector


}	// namespace AbstractLinAlgPack

#endif	// ETA_VECTOR_H
