/*
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
*/

// /////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_VectorDenseEncap.hpp

#ifndef SLAP_VECTOR_DENSE_ENCAP_H
#define SLAP_VECTOR_DENSE_ENCAP_H

#include "AbstractLinAlgPack_Types.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"

namespace AbstractLinAlgPack {

/** \brief Extract a constant <tt>DenseLinAlgPack::DVectorSlice</tt> view of a <tt>Vector</tt> object.
 *
 * This class is only to be used with <tt>Vector</tt> objects that store all of their
 * elements in the local address space or can easily access all of the vector elements in
 * this process (or thread).  It is generally to be used in serial applications but might
 * also find use in parallel appliations where a vector is replicated across processes.
 *
 * This utility class is defined purly in terms of the abstract interfaces.  It is only to
 * be used as an automatic variable on the stack.  For example, to extract a <tt>DVectorSlice</tt>
 * view of an abstract vector and use it to copy to another <tt>DVectorSlice</tt> object you could
 * write a function like:
 \code
 void copy(const Vector& vec_in, DVectorSlice* vs_out ) {
     VectorDenseEncap  vs_in(vec_in);
   *vs_out = vs_in();
 }
 \endcode
 * In the above code, if the underlying <tt>Vector</tt> object does not have to
 * perform any dynamic memory allocations and copy in the method <tt>Vector::get_sub_vector()</tt>
 * then the above code will only have a constant time overhead.  However, the above approach will work
 * for any <tt>Vector</tt> object (no matter how inefficient it may be).
 */
class VectorDenseEncap {
public:

  /// Calls <tt>vec.get_sub_vector(Range1D(),DENSE,&sub_vec)</tt> to get the view.  
  VectorDenseEncap( const Vector&  vec );
  /// Calls <tt>vec.free_sub_vector(&sub_vec)</tt> to release the view.  
  ~VectorDenseEncap();
  /// Returns a reference to a constant view of the dense vector.
  const DVectorSlice& operator()() const;

private:

  const Vector            &vec_;
  RTOpPack::SubVector     sub_vec_;
  DVectorSlice            vs_;
  VectorDenseEncap();                                     // Not defined and not to be called!
  VectorDenseEncap(const VectorDenseEncap&);              // ""
  VectorDenseEncap& operator=(const VectorDenseEncap&);   // ""

}; // end class VectorDenseEncap

/** \brief Extract a non-const <tt>DenseLinAlgPack::DVectorSlice</tt> view of a <tt>VectorMutable</tt> object.
 *
 * This utility class is defined purly in terms of the abstract interfaces.  It is only to
 * be used as an automatic variable on the stack.  Note that the underlying <tt>VectorMutable</tt>
 * object is  not guarrenteed to be modified until the destructor for \c this is called.
 */
class VectorDenseMutableEncap {
public:

  /// Calls <tt>vec.get_sub_vector(Range1D(),&sub_vec)</tt> to get the view.  
  VectorDenseMutableEncap( VectorMutable&  vec );
  /// Calls <tt>vec.commit_sub_vector(&sub_vec)</tt> to release the view.  
  ~VectorDenseMutableEncap();
  /// Returns a reference to a constant view of the dense vector.
  DVectorSlice& operator()();
  /// Returns a reference to a non-const view of the dense vector.
  const DVectorSlice& operator()() const;

private:

  VectorMutable                  &vec_;
  RTOpPack::MutableSubVector     sub_vec_;
  DVectorSlice                   vs_;
  VectorDenseMutableEncap();                                            // Not defined and not to be called!
  VectorDenseMutableEncap(const VectorDenseMutableEncap&);              // ""
  VectorDenseMutableEncap& operator=(const VectorDenseMutableEncap&);   // ""

}; // end class VectorDenseMutableEncap

// ///////////////////////////////////////////
// Inline members

// VectorDenseEncap

inline
VectorDenseEncap::VectorDenseEncap( const Vector&  vec )
  :vec_(vec)
{
  vec_.get_sub_vector(Range1D(),&sub_vec_);
  vs_.bind( DVectorSlice(
          const_cast<value_type*>(sub_vec_.values())
          ,sub_vec_.subDim()
          ,sub_vec_.stride()
          )
    );
}

inline
VectorDenseEncap::~VectorDenseEncap()
{
  vec_.free_sub_vector(&sub_vec_);
}

inline
const DVectorSlice& VectorDenseEncap::operator()() const
{
  return vs_;
}

// VectorDenseMutableEncap

inline
VectorDenseMutableEncap::VectorDenseMutableEncap( VectorMutable&  vec )
  :vec_(vec)
{
  vec_.get_sub_vector(Range1D(),&sub_vec_);
  vs_.bind( DVectorSlice(
          sub_vec_.values()
          ,sub_vec_.subDim()
          ,sub_vec_.stride()
          )
    );
}

inline
VectorDenseMutableEncap::~VectorDenseMutableEncap()
{
  vec_.commit_sub_vector(&sub_vec_);
}

inline
DVectorSlice& VectorDenseMutableEncap::operator()()
{
  return vs_;
}

inline
const DVectorSlice& VectorDenseMutableEncap::operator()() const
{
  return vs_;
}

} // end namespace SparseLinALgPack

#endif // SLAP_VECTOR_DENSE_ENCAP_H
