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
//
// Change Log:
//	11/18/99:
//		* last_updated() Added
//		* set_not_updated(k) Added

#ifndef ITER_QUANTITY_H
#define ITER_QUANTITY_H

#include <stdexcept>
#include <string>
#include <iomanip>
#include <limits>

#include "IterationPack_Types.hpp"

namespace IterationPack {

/** \brief Iterface for information about Iteration Quantities.
  *
  * This class provides and interface to all concrete types of iteration quantities
  * and provides all of the services except storage access.  It is assumed that
  * every concrete iteration quantity subclass will will be derived from
  * the tempalted interface class \c IterQuantityAccess.
  *
  * ToDo: Finish the documentation and give examples.
  */
class IterQuantity {
public:

  /** @name Public types */
  //@{

  /// Constant for value returned when no iteration quantity has been updated.
  enum { NONE_UPDATED = INT_MIN };

  /// Thrown memory if attempted to be set that storage can not be allocated to.
  class NoStorageAvailable : public std::logic_error
  {public: NoStorageAvailable(const std::string& what_arg) : std::logic_error(what_arg) {}};

  /// Thrown when memory access is attempted when it has not yet been updated.
  class QuanityNotSet : public std::logic_error
  {public: QuanityNotSet(const std::string& what_arg) : std::logic_error(what_arg) {}};
  
  //@}

  /** @name Constructors, destructors and memory managment */
  //@{

  // Virtual destructor
  virtual ~IterQuantity() {}

  // Clone this iteration quantity
  virtual IterQuantity* clone() const = 0;

  //@}

  /** @name Misc query (const) methods */
  //@{

  /// Return the name (zero terminated string) of this quantity.
  virtual const char* name() const = 0; 
  
  /** \brief Determine if there is storage advailable for the k <tt>offset</tt> iteration quanity.
    *
    * If this function returns true then <tt>set_f(offset)</tt> can be called to set
    * the quanity for the kth iteration or <tt>get_k(offset)</tt> (see \c IterQuantityAccess)
    * can be called if <tt>updated_k(offset)</tt> is already true.  
    */
  virtual bool has_storage_k(int offset) const = 0;
  
  /** \brief Determine if the quanity for the k <tt>offset</tt> iteration has been accessed
    * by a call to <tt>set_k()</tt> (see \c IterQuantityAccess).
    *
    * This function does not confirm that the k <tt>offset</tt> quanity has been
    * set to a meaningfull value, only that <tt>set_k()</tt> was called
    * to get a reference to that quanity and <tt>get_k()</tt> can be called to
    * get that same reference.
    */
  virtual bool updated_k(int offset) const = 0;

  /** \brief Return the highest k such that <tt>updated_k(k)</tt> returns true.
    *
    * If <tt>updated_k(k) == false</tt> false for all \c k then this function
    * will return \c NONE_UPDATED.
    */
  virtual int last_updated() const = 0;

  /** \brief Determine if the memory for the k + <tt>offset</tt> quantityy will be lost if
    * <tt>set_k(set_offset)</tt> is called (see \c IterQuantityAccess).
    *
    * This member function allows clients to know a little about the
    * specific behavior of the subclass.  Clients can use this function
    * to determine if it is safe to call <tt>get_k(offset)</tt> after <tt>set_k(set_offset)</tt>
    * is called.  For example, imagine the case where you wanted to update a 
    * vector in iteration k+1 given the elements in the k iteraiton.  For
    * a subclass with only single storage (<tt>info.will_loose_mem(0,+1) == true</tt>)
    * the following code would not work:
    \code
    for(int i = 1; i <= n; ++i) {
        info.set_k(+1)(i) = info.get_k(0)(i);
    }
    \endcode
    * For <tt>i</tt> == 1, <tt>set_k(+1)</tt> would cause a state transition and for <tt>i</tt> == 2
    * <tt>info.get_k(0)</tt> would throw an exception.  Actually, the compiler may evaluate
    * <tt>info.set_k(+1)</tt> before <tt>info.get_k(0)</tt> so <tt>info.get_k(0)</tt> whould throw
    * an exception right away for <tt>i</tt> == 1.
    *
    * If the client knows that only single storage is needed then it could use
    * something like the following code:
    \code
    if(info.will_loose_mem(0,+1) {
        info.set_k(+1) = info.get_k(0);
        for(int i = 1; i <= n; ++i) info.set_k(+1)(i) = info.get_k(+1)(i) * 2;
    }
    else {
        for(int i = 1; i <= n; ++i) info.set_k(+1)(i) = info.get_k(0)(i);
    }
    \endcode
    * In an actually implemention one would use temporary references and would not
    * call \c set_k() and \c get_k() multiple times like this but you get
    * the basic idea.  The above code works for both single and multiple storage
    * and will not result in any unnecessary copying since assingment to self
    * should be detected.  In the above code <tt>info.set_k(+1) = info.get_k(0);</tt>
    * is called to effect the state transistion.
    *
    * On the other hand if you need dual storage you will need a temporary
    * copy in the event that <tt>will_loose_mem(offset, set_offset)</tt> returns true.
    * For example you need dual storage for the code:
    \code
    for(int i = 2; i <= n; ++i) info.set_k(+1)(i) = info.get_k(0)(i) * info.get_k(0)(i-1);
    \endcode
    * Even the above operation can be implemented without a temporary vector but you get the
    * idea, the (i-1) quanity is modifed and is not the original for i > 2.
    *
     * Preconditions:<ul>
    * <li> <tt>updated_k(offset) == true</tt> [throw <tt>QuanityNotSet</tt>]
    * </ul> 
    */
  virtual bool will_loose_mem(int offset, int set_offset) const = 0;

  //@}

  /** @name Misc modifier (non-const) methods */
  //@{

  /** \brief Causes <tt>updated_k(k)</tt> to return false.
    *
     * Preconditions:<ul>
    * <li> <tt>updated_k(offset) == true</tt> [throw <tt>QuanityNotSet</tt>]
    * </ul> 
    *
     * Postconditions:<ul>
    * <li> <tt>updated_k(offset) == false</tt>
    * </ul> 
    */
  virtual void set_not_updated_k(int offset) = 0;

  /** \brief Causes <tt>updated_k(k)</tt> to return false for all <tt>k</tt>.
    */
  virtual void set_all_not_updated() = 0;

  //@}

  /** @name Iteration incrementation */
  //@{

  /** \brief Shift the reference point from the k to the k+1 iteration.
    *
     * Postcondtions:<ul>
    * <li> <tt>updated_k(offset)</tt> before the call equals <tt>updated_k(offset-1)</tt> after return
    * <li> <tt>&this->get_k(offset)</tt> for <tt>this->updated_k(offset) == true</tt> before the call, equals
    *			<tt>&this->get_k(offset-1)</tt> for <tt>this->updated_k(offset-1) == true</tt> after return
    * </ul> 
    */
  virtual void next_iteration() = 0;

  //@}

  /** @name Runtime information */
  //@{

  /** \brief Print to an output stream a description of this iteration quantity.
   *
   * The purpose if this method is allow the client get information as to what the
   * type of the iteration quantity really is for debugging and informational purposes.
   * This should just include information on types and nothing else.
   *
   * The concrete type of \c this can be printed using <tt>typeName(*this)</tt>.
   */
  virtual void print_concrete_type( std::ostream& out ) const = 0;

  //@}

  /** @name Assert state */
  //@{

  /// Assert <tt>has_storage_k(offset) == true</tt> (throw <tt>NoStorageAvailable</tt>).
  void assert_has_storage_k(int offset) const;

  /// Assert updated_k(offset) == true<tt> (throw QuanityNotSet).
  void assert_updated_k(int offset) const;

  //@}

};	// end class IterQuantity 

}	// end namespace IterationPack

#endif	// ITER_QUANTITY_H
