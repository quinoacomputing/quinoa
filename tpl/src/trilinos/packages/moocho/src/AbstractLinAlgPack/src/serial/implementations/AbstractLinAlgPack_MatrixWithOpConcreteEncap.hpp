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

#ifndef MATRIX_WITH_OP_CONCRETE_ENCAP_H
#define MATRIX_WITH_OP_CONCRETE_ENCAP_H

#include "AbstractLinAlgPack_MatrixOp.hpp"

namespace AbstractLinAlgPack {

/** \brief This template class defines the storage for a concrete matrix
  * class that operations are based on.
  *
  * The default copy constructor and assignment operator are allowed.
  */
template<class M>
class MatrixWithOpConcreteEncap : public virtual MatrixOp
{
public:

  // /////////////////////////////////////////////////////
  /** @name Representation access */
  //@{

  /// The compiler did not generate this default constructor
  MatrixWithOpConcreteEncap()
  {}

  /// This constructor will have to be overridden.
  MatrixWithOpConcreteEncap(const M& m) : m_(m)
  {}

  /// Get the underlying M object
  M& m() {
    return m_;
  }

  /** \brief . */
  const M& m() const {
    return m_;
  }

  //@}	// end Representation access

  // /////////////////////////////////////////////////////
  // Overridden from Matrix

  /** \brief . */
  size_type rows() const;

  /** \brief . */
  size_type cols() const;

  // /////////////////////////////////////////////////////
  // Overridden from MatrixOp

  /** \brief . */
  MatrixOp& operator=(const MatrixOp& m);

private:
  M m_;

};	// end class MatrixWithOpConcreteEncap<M>

// Template definitions

template<class M>
size_type MatrixWithOpConcreteEncap<M>::rows() const {
  return m().rows();
}

template<class M>
size_type MatrixWithOpConcreteEncap<M>::cols() const {
  return m().cols();
}

template<class M>
MatrixOp& MatrixWithOpConcreteEncap<M>::operator=(const MatrixOp& m) {
  if(&m == this) return *this;	// assignment to self
  const MatrixWithOpConcreteEncap<M> *p_m = dynamic_cast<const MatrixWithOpConcreteEncap<M>*>(&m);
  if(p_m) {
    m_ = p_m->m_;
  }
  else {
    throw std::invalid_argument("MatrixWithOpConcreteEncap<M>::operator=(const MatrixOp& m)"
      " : The concrete type of m is not a subclass of MatrixWithOpConcreteEncap<M> as expected" );
  }
  return *this;
}

}	// end namespace AbstractLinAlgPack 

#endif	// MATRIX_WITH_OP_CONCRETE_ENCAP_H
