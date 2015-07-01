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

#ifndef SUNDANCE_PARAMETER_H
#define SUNDANCE_PARAMETER_H


#include "SundanceDefs.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceSpatiallyConstantExpr.hpp"
#include "SundanceSymbolicFuncDescriptor.hpp"
#include "SundanceParameterData.hpp"


namespace Sundance
{
  using namespace Sundance;
  using namespace Teuchos;

  
  

  /** 
   * A Parameter is an expression subtype used to represent
   * a spatially-constant parameter that may change during the
   * course of a simulation, for example, the time in a transient
   * simulation or a continuation parameter when using a homotopy method.
   * While it is possible to use simple double-precision 
   * constants in expressions, their values are immutable once created.
   * When a constant's value may need to be changed, use a Parameter rather
   * than a simple double.
   */
  class Parameter : public virtual DiscreteFuncElement,
                    public SymbolicFuncDescriptor,
                    public virtual SpatiallyConstantExpr
  {
  public:
    /** */
    Parameter(const double& value, const std::string& name="");

    /** virtual destructor */
    virtual ~Parameter() {;}

    /** */
    virtual XMLObject toXML() const ;

    /** */
    virtual const double& value() const ;
    
    /** */
    virtual void setValue(const double& value) ;

    //@{
    /** \name Functions for internal use */
    /** */
    Set<MultipleDeriv> 
    internalFindW(int order, const EvalContext& context) const ;
          
    /** */
    Set<MultipleDeriv> 
    internalFindV(int order, const EvalContext& context) const ;
    
    /** */
    Set<MultipleDeriv> 
    internalFindC(int order, const EvalContext& context) const ;

    /** */
    RCP<Array<Set<MultipleDeriv> > > 
    internalDetermineR(const EvalContext& context,
                       const Array<Set<MultipleDeriv> >& RInput) const ;


    /** */
    virtual Evaluator* createEvaluator(const EvaluatableExpr* expr,
                                       const EvalContext& context) const ;

    /** */
    virtual RCP<ExprBase> getRcp() {return rcp(this);}
    //@}

    /** Write self in text form */
    virtual std::ostream& toText(std::ostream& os, bool paren) const 
    {os << "Parameter[" << name() << " = " << value() << "]"; return os;}

    /** */
    bool isParameter() const {return true;}


  protected:
    /** */
    const ParameterData* data() const ;

    /** */
    ParameterData* data();


  };
}
#endif
