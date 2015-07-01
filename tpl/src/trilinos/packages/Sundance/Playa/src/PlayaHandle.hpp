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

#ifndef PLAYA_HANDLE_HPP
#define PLAYA_HANDLE_HPP

#include "PlayaDefs.hpp"
#include "PlayaOut.hpp"
#include "PlayaPrintable.hpp"
#include "Teuchos_Describable.hpp"
#include "PlayaHandleable.hpp"
#include "PlayaExceptions.hpp"
#include "PlayaObjectWithVerbosity.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_TypeNameTraits.hpp"


 /** \def This helper macro defines boilerplate constructors for classes deriving
  * from Handle.
  *
  * If class MyHandle is a handle to a type MyType, simply 
  * put
  * \code
  * HANDLE_CTORS(MyHandle, MyType);
  * \endcode
  * in the class declaration of MyHandle and the macro will create 
  * an empty ctor, a ctor from a smart ptr, and a ctor from a raw pointer. 
  * The macro will also create appropriate doxygen for the handle ctors */
#define HANDLE_CTORS(handle, contents) \
  /** Empty ctor */ \
handle() : Playa::Handle<contents >() {;} \
  /** Construct a #handle with a raw pointer to a #contents */ \
  handle(Playa::Handleable<contents >* rawPtr) : Playa::Handle<contents >(rawPtr) {;} \
  /** Construct a #handle with a smart pointer to a #contents */ \
  handle(const Teuchos::RCP<contents >& smartPtr) : Playa::Handle<contents >(smartPtr){;}



namespace Playa
{
using Teuchos::RCP;
using Teuchos::rcp;


/** This traits class is used to extract the non-const version of
 * a template argument. The generic case returns the template argument. */
template <class X>
class ConstHandleTraits
{
public:
  typedef X NonconstType;
};

/** Specialization of CHT to types "const X". The nonconst type can
 * be extracted from the template argument. This lets us figure out
 * type X when the handle is to a "const X" */
template <class X>
class ConstHandleTraits<const X>
{
public:
  typedef X NonconstType;
};



/**
 * Class Playa::Handle provides a general implementation
 * of the common features of reference-counted handles.
 */
template <class PointerType>
class Handle 
{
public:
  /** Empty ctor  */
  Handle() : ptr_() {;}

  /** Construct from a smart pointer */
  Handle(const RCP<PointerType>& _ptr) : ptr_(_ptr) {;}

  /** Construct from a raw pointer to a Handleable.  */
  Handle(Handleable<PointerType>* rawPtr) : ptr_(rawPtr->getRcp()) {;}

  /** Read-only access to the underlying smart pointer. */
  const RCP<PointerType>& ptr() const {return ptr_;}

  /** Read-write access to the underlying smart pointer. */
  RCP<PointerType>& ptr() {return ptr_;}

  /** 
   * Print to a stream. This tries to use, in order, the 
   * Printable, Describable, and fallback forms of output. 
   */
  void print(std::ostream& os) const ;


  /** 
   * Return a short descriptive std::string using the Describable interface.
   * If the contents of the handle cannot be 
   * downcasted or crosscasted to a Describable*, an exception
   * will be thrown. 
   */
  std::string description() const ;

  /** Write a fallback description to be used in objects that are
   * neither named, printable, or describable */
  std::string fallbackDescription() const ;

  /** Return the verbosity */
  int verb() const ;

  /** Set the verbosity */
  void setVerb(int v);

private:
  RCP<PointerType> ptr_;
};

/* implementation of print() */
template <class PointerType> inline 
void Handle<PointerType>::print(std::ostream& os) const 
{
  const Printable* p = dynamic_cast<const Printable*>(ptr_.get());
  const Describable* d = dynamic_cast<const Describable*>(ptr_.get());

  if (p!=0) p->print(os);
  else if (d!=0) os << d->description();
  else os << fallbackDescription();
}

/* implementation of description() */
template <class PointerType> inline
std::string Handle<PointerType>::description() const 
{
  const Describable* d = dynamic_cast<const Describable*>(ptr_.get());
  TeuchosOStringStream oss;

  if (d!=0) oss << d->description();
  else oss << fallbackDescription();

  return oss.str();
}

template <class PointerType> inline
std::string Handle<PointerType>::fallbackDescription() const
{
  typedef typename ConstHandleTraits<PointerType>::NonconstType NC;
  TeuchosOStringStream oss;

  oss << "Handle[" << TypeNameTraits<NC>::name()
      << ", ptr=" << ptr_.get() << "]";
  return oss.str();
}



/* implementation of verb() */
template <class PointerType> inline
int Handle<PointerType>::verb() const 
{
  const ObjectWithVerbosity* v 
    = dynamic_cast<const ObjectWithVerbosity*>(ptr_.get());
  if (v) return v->verb();
  return 0;
}


/* implementation of verb() */
template <class PointerType> inline
void Handle<PointerType>::setVerb(int verbosity) 
{
  ObjectWithVerbosity* v 
    = dynamic_cast<ObjectWithVerbosity*>(ptr_.get());
  if (v) v->setVerb(verbosity);
  TEUCHOS_TEST_FOR_EXCEPTION(v==0, RuntimeError, 
    "attempt to set verbosity on a handle that doesn't wrap "
    "an ObjectWithVerbosity subtype");
}

template <class PointerType> inline
std::ostream& operator<<(std::ostream& os, const Playa::Handle<PointerType>& h)
{
  h.print(os);
  return os;
}

}

#define STREAM_OUT(handleType) \
                        inline std::ostream& operator<<(std::ostream& os, const handleType& h) \
                        {h.print(os); return os;}



#endif

