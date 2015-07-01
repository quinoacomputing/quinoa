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



#ifndef SUNDANCE_SYMBPREPROCESSOR_H
#define SUNDANCE_SYMBPREPROCESSOR_H

#include "SundanceExpr.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceComputationType.hpp"



namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;




class Parameter;


/** */
class SymbPreprocessor 
{
public:
        
  /** */
  static DerivSet setupVariations(const Expr& expr, 
    const Expr& vars,
    const Expr& varEvalPts,
    const Expr& unks,
    const Expr& unkEvalPts,
    const Expr& unkParams,
    const Expr& unkParamEvalPts,
    const Expr& fixedFields,
    const Expr& fixedFieldEvalPts, 
    const Expr& fixedParams,
    const Expr& fixedParamEvalPts, 
    const EvalContext& context,
    const ComputationType& compType);

  /** */
  static DerivSet setupFunctional(const Expr& expr, 
    const Expr& fixedParams,
    const Expr& fixedParamEvalPts,
    const Expr& fixedFields,
    const Expr& fixedFieldEvalPts,
    const EvalContext& context,
  const ComputationType& compType);

  /** */
  static DerivSet setupGradient(const Expr& expr, 
    const Expr& vars,
    const Expr& varEvalPts,
    const Expr& fixedParams,
    const Expr& fixedParamEvalPts,
    const Expr& fixedFields,
    const Expr& fixedFieldEvalPts, 
    const EvalContext& contex,
  const ComputationType& compType);

  /** */
  static DerivSet setupSensitivities(const Expr& expr, 
    const Expr& tests,
    const Expr& unks,
    const Expr& unkEvalPts, 
    const Expr& unkParams,
    const Expr& unkParamEvalPts,
    const Expr& fixedParams,
    const Expr& fixedParamEvalPts,
    const Expr& fixedFields,
    const Expr& fixedFieldEvalPts,
    const EvalContext& context,
  const ComputationType& compType);

  /** */
  static DerivSet setupFwdProblem(const Expr& expr, 
    const Expr& tests,
    const Expr& unks,
    const Expr& unkEvalPts, 
    const Expr& unkParams,
    const Expr& unkParamEvalPts,
    const Expr& fixedParams,
    const Expr& fixedParamEvalPts,
    const Expr& fixedFields,
    const Expr& fixedFieldEvalPts,
    const EvalContext& context,
  const ComputationType& compType);

  /** check the input functions for redundancies and functions
   * of unexpected type. Set evaluation points and collect
   * a set of function IDs. This function is templated so we can use
   * the same code for processing unknown and variational functions */
  template <class T>
  static Set<int> processInputFuncs(const Expr& u, const Expr& u0)
    {
      Set<int> idSet;

      for (int i=0; i<u.size(); i++)
      {
        /* make sure all input functions have the correct type */
        const T* uPtr = dynamic_cast<const T*>(u[i].ptr().get());
        TEUCHOS_TEST_FOR_EXCEPTION(uPtr==0, std::runtime_error, 
          "Unexpected function type error: function " << u[i].toString()
          << " is of type=" << typeid(u).name() 
          << ", but we expected type=" << typeid(T).name());

        /* Add the function's ID to the ID set. While we're here, check
         * to ensure we have no duplicates in the input list. */
        int fid = uPtr->fid().dofID();
        TEUCHOS_TEST_FOR_EXCEPTION(idSet.contains(fid), std::runtime_error,
          "duplicate function in input list " << u.toString());
        idSet.put(fid);


        /* Check that the evaluation point is either a discrete function
         * or a zero expression. */
        RCP<DiscreteFuncElement> u0Ptr
          = rcp_dynamic_cast<DiscreteFuncElement>(u0[i].ptr());
        RCP<ZeroExpr> u0ZeroPtr
          = rcp_dynamic_cast<ZeroExpr>(u0[i].ptr());
        TEUCHOS_TEST_FOR_EXCEPTION(u0Ptr.get()==NULL && u0ZeroPtr.get()==NULL,
          std::runtime_error,
          "evaluation point " << u0[i].toString() << " for func=" << u[i]
          << " is neither a discrete function nor a zero expr");

        /* Set the evaluation point */
        if (u0Ptr.get()==NULL)
        {
          uPtr->substituteZero();
        }
        else
        {
          uPtr->substituteFunction(u0Ptr);
        }
      }

      return idSet;
    }

  /** check the input parameters for redundancies and type. Set 
   * evaluation points and collect
   * a set of parameter IDs. This function is templated so we can use
   * the same code for processing unknown and variational parameters */
  template <class T>
  static Set<int> processInputParams(const Expr& alpha, const Expr& alpha0)
    {
      Set<int> paramID;

      for (int i=0; i<alpha.size(); i++)
      {
        /* ensure everyone has the correct type */
        const T* aPtr = dynamic_cast<const T*>(alpha[i].ptr().get());
        TEUCHOS_TEST_FOR_EXCEPTION(aPtr==0, std::runtime_error,
          "list of purported parameters "
          "contains a function that is not an unknown parameter:"
          << alpha[i].toString());

        int fid = aPtr->fid().dofID();
        TEUCHOS_TEST_FOR_EXCEPTION(paramID.contains(fid), std::runtime_error,
          "duplicate input parameter in list "
          << alpha.toString());
        paramID.put(fid);

        RCP<Parameter> a0Ptr
          = rcp_dynamic_cast<Parameter>(alpha0[i].ptr());
        TEUCHOS_TEST_FOR_EXCEPTION(a0Ptr.get()==NULL,
          std::runtime_error,
          "parameter evaluation point " << alpha0[i].toString()
          << " is not a parameter");
        aPtr->substituteFunction(a0Ptr);
      }
      return paramID;
    }


  /** */
  TEUCHOS_TIMER(preprocTimer, "symbolic preprocessing");
};

/** */
Expr makeZeros(const Expr& e);


}

#endif
