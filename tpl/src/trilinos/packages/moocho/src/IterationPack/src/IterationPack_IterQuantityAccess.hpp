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

#ifndef ITER_QUANTITY_ACCESS_H
#define ITER_QUANTITY_ACCESS_H

#include "IterationPack_IterQuantity.hpp"

namespace IterationPack {

/** \brief Interface to typed iteration quantities.
 *
 * Quantities are updated, read and queried given the
 * offset to the current iteration k.  For example, to set
 * a quantity for the <tt>k+1</tt> iteration you would call <tt>set_k(+1)</tt>.
 * The functions ending with <tt>prefix_k(offset)</tt> are meant to suggest
 * <tt>prefix k + offset</tt>.  For example:
 \verbatim

 has_storage_k(+1) => has storage k+1
 get_k(-1) => get k-1
 \endverbatim
 * Subclasses can implement this interface in a variety of ways.  But
 * they must follow a few simple rules: <ul>
 * <li>  Only forward transitions are allowed.  This effects the behavior of
 *	      \c IterQuantity::has_storage_k() and \c set_k().
 * </ul>
 *
 * Note that any type that is to be used as an iteration quantity must at
 * the bare minimum have the assignment operation defined for it.
 *
 * The client should not have to worry about how much memory is
 * available.  Instead, it is for the object that configures the client
 * to provide the appropriate subclass to meet the needs of the client.
 *
 * <b> Usage: </b>
 *
 * There are sevearl different techniques for using <tt>IterQuantityAccess<T_info></tt>.
 * An implementation subclass for <tt>IterQuantityAccess<T_info></tt> will maintain
 * one or more storage locations for iteration quantities that can be updateded
 * and accessed.  These storage locations need not be contiguous but this is the
 * most common scenario (see <tt>IterQuantityAccessContiguous</tt>).
 *
 * It is important to understand what a client is requesting and what is implied
 * by various calls to get and set methods.
 *
 * By calling the const version of <tt>get_k(offset)</tt>, the client is expecting
 * to gain read access to previously updated iteration quantity.  This quantity
 * should not be modified through this constant reference.
 *
 * The non-const version of <tt>get_k(offset)</tt> is called by a client
 * to gain read-write access to a reviously updated iteration quantity.  The quantity
 * can be further modified through this non-constant reference.
 *
 * The method <tt>set_k(offset)</tt> is called by a client so that it may
 * fully initialize the iteration quantity which may not have been previously
 * updated.  The object pointed to by the non-const reference returned from
 * <tt>set_k(offset)</tt> can not expected to have any more than a default
 * initialization.
 *
 * Finally, the method <tt>set_k(set_offset,get_offset)</tt> is called to
 * explicitly update the iteration quantitiy at <tt>k + set_offset</tt> to the
 * value at iteration <tt>k + get_offset</tt>.  This method is most useful
 * when the client simply wants to set the iteration quantity at the current
 * iteration (<tt>set_offset = 0</tt> to its value at the previous iteration
 * (<tt>get_offset = -1</tt>).  The outcome of this operation is insured no
 * matter how the storage is managed by the subclass.  The non-const reference
 * returned from this method may be used to futher modify the just updated
 * iteration quantity or the reference can just be discarded.
 */
template<class T_info>
class IterQuantityAccess : public IterQuantity {
public:

  /** \brief . */
  typedef	IterQuantity::NoStorageAvailable	NoStorageAvailable;
  /** \brief . */
  typedef IterQuantity::QuanityNotSet			QuanityNotSet;

  /** \brief Return a reference for the <tt>k + offset</tt> iteration quanity.
   *
   * Clients call this member function to access a quantity for a
   * given iteration or modify the quantity which has already been
   * set for that iteration.
   *
   * Preconditions:<ul>
   * <li> <tt>this->updated_k(offset) == true</tt> (throw QuanityNotSet)
   * </ul>
   */
  virtual T_info& get_k(int offset) = 0;

  /** \brief Return a const reference for the <tt>k + offset</tt> iteration quanity.
    *
    * Preconditions:<ul>
    * <li> <tt>this->updated_k(offset) == true</tt> (throw QuanityNotSet)
    * </ul>
    *
    * Clients call this member function to access a const quantity for a
    * given iteration.
    */
  virtual const T_info& get_k(int offset) const = 0;

  /** \brief Return a reference to the storage location for the <tt>k + offset</tt> iteration quanity.
   *
   * Precondtions:<ul>
   * <li> <tt>this->has_storage_k(offset) == true</tt> (throw <tt>NoStorageAvailable</tt>)
   * </ul>
   *
   * Postcondtions:<ul>
   * <li> <tt>this->updated_k(offset) == true</tt>
   * <li> <tt>this->updated_k(i) == false</tt> for <tt>i>/tt> in the set of
   *      <tt>{ i : <tt>this->will_loose_mem(i,offset) == true</tt> }</tt> before this call.
   * </ul> 
   *
   * This function will return a reference to the storage for the
   * <tt>k + offset</tt> iteration quantity.  Calling this function
   * may cause the loss of memory for a back iteration.  If
   * <tt>will_loose_mem(back_offset,offset)</tt> returns
   * <tt>true</tt> then <tt>updated_k(back_offset)</tt> will return
   * <tt>false</tt> after this function returns (assuming an
   * exception is not thrown).
   *
   * The client should expect nothing more than simple default
   * initialization of the object who's reference is returned from
   * this method.  The client is expected to use this reference to
   * initalize this object appropriately.  If the client does not
   * sufficiently update the object at <tt>k + offset</tt> before
   * the reference is let go, the object's reference can be required
   * with a call to the non-const version of <tt>get_k(offset)</tt>.
   */
  virtual T_info& set_k(int offset) = 0;

  /** \brief Set the iteration quantity for the <tt>k + set_offset</tt>
   * iteration to the <tt>k + get_offset</tt> iteration and return
   * the reference to the <tt>k + set_offset</tt> iteration
   * quantity.
   *
   * @param  set_offset  [in] The iteration offset to be set.
   * @param  get_offset  [in[ The iteration offset to copy into the
   *                     <tt>k + set_offset</tt> iteration.
   *
   * Precondtions:<ul>
   * <li> <tt>this->has_storage_k(set_offset) == true</tt> (throw <tt>NoStorageAvailable</tt>)
   * <li> <tt>this->updated_k(get_offset) == true</tt> (throw QuanityNotSet)
   * </ul>
   *
   * Postcondtions:<ul>
   * <li> <tt>this->updated_k(set_offset) == true</tt>
   * <li> <tt>this->updated_k(i) == false</tt> for <tt>i>/tt> in the set of
   *      <tt>{ i : <tt>this->will_loose_mem(i,offset) == true</tt> }</tt> before this call.
   * </ul> 
   *
   * This method blends the functionality of the <tt>get_k()</tt>
   * and the other <tt>set_k()</tt> methods.  This method ensures
   * that a quantity from one iteration (<tt>k + get_offset</tt>)
   * will be properly and efficienlty copied into the storage
   * location of another iteration (<tt>k + set_offset</tt>) not
   * matter how storage is handled by the implementing subclass.
   */
  virtual T_info& set_k(int set_offset, int get_offset) = 0;

};	// end class IterQuantityAccess 

}	// end namespace IterationPack

#endif	// ITER_QUANTITY_ACCESS_H
