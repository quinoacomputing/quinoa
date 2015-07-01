/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */

#ifndef SUNDANCEORDEREDHANDLE_HPP
#define SUNDANCEORDEREDHANDLE_HPP

#include "SundanceDefs.hpp"
#include "PlayaHandle.hpp"
#include <typeinfo>


#define ORDERED_HANDLE_CTORS(handle, contents) \
/** Empty ctor */ \
handle() : OrderedHandle<contents >() {;} \
/** Construct a #handle with a raw pointer to a #contents */ \
handle(Playa::Handleable<contents >* rawPtr) : OrderedHandle<contents >(rawPtr) {;} \
/** Construct a #handle with a smart pointer to a #contents */ \
handle(const RCP<contents >& smartPtr) : OrderedHandle<contents >(smartPtr){;}




namespace Sundance
{
  using namespace Teuchos;

  /**
   * Class OrderedHandle is an extension to Playa::Handle that 
   * includes a comparison operator ("<" operator) so that
   * the handle can be used in ordered containers such as STL maps and sets.
   */
  template <class PointerType>
  class OrderedHandle : public Playa::Handle<PointerType>
  {
  public:
    /** empty ctor */
    OrderedHandle() : Playa::Handle<PointerType>() {;}

    /** Construct from a raw ptr */
    OrderedHandle(Playa::Handleable<PointerType>* rawPtr) : Playa::Handle<PointerType>(rawPtr) {;}

    /** Construct from a smart ptr*/
    OrderedHandle(const RCP<PointerType>& smartPtr) 
      : Playa::Handle<PointerType>(smartPtr) {;}

    /** comparison operator */
    bool operator<(const OrderedHandle<PointerType>& other) const 
    {
      /* first compare types */
      const PointerType* me = this->ptr().get();
      const PointerType* you = other.ptr().get();
      if (me==0 && you==0) 
        {
          return false;
        }
      if (me==0) 
        {
          return true;
        }
      if (you==0) 
        {
          return false;
        }

      if (typeid(*me).before(typeid(*you))) 
        {
          return true;
        }

      if (typeid(*you).before(typeid(*me)))
        {
          return false;
        }
      
      /* if the types are equal, compare values of the contents using
       * the lessThan() method. */
      bool rtn = this->ptr()->lessThan(other.ptr().get());
      return rtn;
    }

    /** */
    bool operator!=(const OrderedHandle<PointerType>& other) const
      {
        return (*this < other) || (other < *this);
      }


    /** */
    bool operator==(const OrderedHandle<PointerType>& other) const
      {
        return !(*this != other);
      }

    
  };
}

namespace std
{
  /** \relates OrderedHandle */
  template <class P> inline 
  std::ostream& operator<<(std::ostream& os, 
                           const Sundance::OrderedHandle<P>& h)
  {
    h.print(os);
    return os;
  }
}




#endif

