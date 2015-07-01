/* @HEADER@ */
// ************************************************************************
// 
//                 Playa: Programmable Linear Algebra
//                 Copyright 2012 Sandia Corporation
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

#ifndef PLAYA_HANDLEABLE_HPP
#define PLAYA_HANDLEABLE_HPP

#include "PlayaDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"

#define GET_RCP(Base) \
  /** Handleable<##Base> interface */ \
  virtual Teuchos::RCP<Base > getRcp() {return rcp(this);}

namespace Playa
{
using namespace Teuchos;

/**
 * Class Handleable provides an abstract interface for polymorphic
 * conversion from raw pointers to smart pointers. Recall from the
 * Teuchos RefCountPtr documentation that one should never create
 * directly a smart pointer from a raw pointer; rather, smart pointers
 * should be created through a call to rcp(). The type of the argument
 * to rcp() must be known at compile time. This makes the syntax
 * \code
 * Handle h = new Derived();
 * \endcode
 * impossible with the straightforward implementation in which Handle takes
 * a raw pointer to a Base. In order to preserve this clean syntax, we
 * require any handles supporting this syntax to take a raw
 * pointer to a Handleable<Base>, where Handleable<Base> provides a 
 * getRcp() method which returns the result of a call to rcp() on this.
 */
template <class Base>
class Handleable
{
public:
  /** Virtual dtor */
  virtual ~Handleable(){;}

  /** Return a safely-created RefCountPtr to the base type */
  virtual RCP<Base> getRcp() = 0 ;
    
};
  
}




#endif
