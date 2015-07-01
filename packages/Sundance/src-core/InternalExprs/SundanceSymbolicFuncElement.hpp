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

#ifndef SUNDANCE_SYMBOLICFUNCELEMENT_H
#define SUNDANCE_SYMBOLICFUNCELEMENT_H


#include "SundanceDefs.hpp"
#include "SundanceFuncElementBase.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceSymbolicFuncEvaluator.hpp"
#include "SundanceSymbolicFuncDescriptor.hpp"
#include "SundanceCommonFuncDataStub.hpp"

namespace Sundance
{
  using namespace Sundance;
    class DiscreteFuncElement;
    using namespace Teuchos;

    
    

    /** 
     * SymbolicFuncElement represents a scalar-valued element of a (possibly)
     * list-valued SymbolicFunction. 
     */
    class SymbolicFuncElement : public FuncElementBase,
                                public SymbolicFuncDescriptor,
                                virtual public EvaluatableExpr,
                                public GenericEvaluatorFactory<SymbolicFuncElement, SymbolicFuncElementEvaluator>
    {
    public:
      /** */
      SymbolicFuncElement(const std::string& name, 
        const std::string& suffix,
        const FunctionIdentifier& fid,
        const RCP<const CommonFuncDataStub>& data);
      
      /** virtual destructor */
      virtual ~SymbolicFuncElement() {;}

      /** Append to the set of func IDs present in this expression. */
      void accumulateFuncSet(Set<int>& funcDofIDs, 
        const Set<int>& activeSet) const ;

      /** */
      virtual bool hasTestFunctions() const {return false;}


      /** Specify that expressions involving this function are to be evaluated
       * with this function set to zero. Test functions should always be
       * evaluated at zero. For unknown functions, 
       * substituting zero is appropriate for computing
       * the functional derivatives that arise in a linear problem.
       * */
      void substituteZero() const ;

      /** Specify that expressions involving this function are to be evaluated
       * with this function set to the discrete function (or constant parameter) \f$u_0\f$. 
       * This is appropriate for computing
       * the functional derivatives that arise in a nonlinear expression
       * being linearized about \f$u_0\f$. 
       */
      void substituteFunction(const RCP<DiscreteFuncElement>& u0) const ;

      /** Return the point in function space at which this symbolic 
       * function is to be evaluated. */
      const EvaluatableExpr* evalPt() const {return evalPt_.get();}

      /** Return the point in function space at which this symbolic 
       * function is to be evaluated. */
      EvaluatableExpr* evalPt() {return evalPt_.get();}


      /** */
      bool evalPtIsZero() const ;

      /** */
      const RCP<const CommonFuncDataStub>& commonData() const {return commonData_;}


      /** \name Preprocessing */
      //@{
      /** */
      virtual Set<MultipleDeriv> 
      internalFindW(int order, const EvalContext& context) const ;
          
      /** */
      virtual Set<MultipleDeriv> 
      internalFindV(int order, const EvalContext& context) const ;
          
      /** */
      virtual Set<MultipleDeriv> 
      internalFindC(int order, const EvalContext& context) const ;

      /** */
      virtual RCP<Array<Set<MultipleDeriv> > > 
      internalDetermineR(const EvalContext& context,
                         const Array<Set<MultipleDeriv> >& RInput) const ;
      /** */
      virtual void registerSpatialDerivs(const EvalContext& context, 
                                         const Set<MultiIndex>& miSet) const ;
      //@}
      

      /** Indicate whether the expression is independent of the given 
       * functions */
      virtual bool isIndependentOf(const Expr& u) const ;

      
      /** Indicate whether the expression is linear in the given 
       * functions */
      virtual bool isLinearForm(const Expr& u) const ;
      
      /** */
      virtual RCP<ExprBase> getRcp() {return rcp(this);}
      
    private:
      RCP<const CommonFuncDataStub> commonData_;

      mutable RCP<EvaluatableExpr> evalPt_;

      mutable Array<int> evalPtDerivSetIndices_;
    };
}

#endif
