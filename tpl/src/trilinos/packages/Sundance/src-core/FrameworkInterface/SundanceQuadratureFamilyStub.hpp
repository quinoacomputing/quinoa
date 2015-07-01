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

#ifndef SUNDANCE_QUADRATUREFAMILYSTUB_H
#define SUNDANCE_QUADRATUREFAMILYSTUB_H

#include "SundanceDefs.hpp"
#include "SundanceNoncopyable.hpp"
#include "PlayaHandleable.hpp"
#include "PlayaPrintable.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_XMLObject.hpp"

namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;

class QuadratureFamilyStub 
  : public Playa::Handleable<QuadratureFamilyStub>,
    public Playa::Printable,
    public Teuchos::Describable,
    public Noncopyable
{
public:
  /** Empty ctor */
  QuadratureFamilyStub(int order);

  /** virtual dtor */
  virtual ~QuadratureFamilyStub(){;}

  /** Return the order of the quadrature rule */
  int order() const {return order_;}

  /** Write to XML */
  virtual XMLObject toXML() const ;

  /** Ordering for storage in STL maps */
  virtual bool lessThan(const QuadratureFamilyStub* other) const 
    {return order() < other->order();}

  /** \name Printable interface */
  //@{
  /** Print to a stream */
  virtual void print(std::ostream& os) const {os << toXML();}
  //@}

  /** \name Describable interface */
  //@{
  /** Print to a stream */
  virtual std::string description() const 
    {return "QuadratureFamilyStub[order=" + Teuchos::toString(order()) 
        +  "]";}
  //@}

  /** */
  virtual RCP<QuadratureFamilyStub> getRcp() {return rcp(this);}

  /** */
  static RCP<QuadratureFamilyStub>& defaultQuadrature()
    {static RCP<QuadratureFamilyStub> rtn
        = rcp(new QuadratureFamilyStub(1)); return rtn;}
      
private:
  int order_;
};
}


#endif
