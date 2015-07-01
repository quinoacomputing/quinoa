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

#ifndef SUNDANCE_STDMATHOPS_H
#define SUNDANCE_STDMATHOPS_H

#include "SundanceDefs.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceExpr.hpp"
#include "SundanceScalarExpr.hpp"
#include "SundanceStdMathFunctors.hpp"
#include "SundanceNonlinearUnaryOp.hpp"



#ifndef DOXYGEN_DEVELOPER_ONLY


using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;

#define SUNDANCE_UNARY_OP(opName, functorName, description) \
/** \relates Expr description */\
inline Expr opName(const Expr& expr) \
{\
RCP<ScalarExpr> arg = rcp_dynamic_cast<ScalarExpr>(expr[0].ptr());\
    TEUCHOS_TEST_FOR_EXCEPTION(arg.get()==0, std::runtime_error,\
                       "non-scalar argument in " #opName " function");\
    return new NonlinearUnaryOp(arg, rcp(new functorName()));\
}

namespace Sundance
{
  inline Expr pow(const Expr& expr, const double& p)
  {
    RCP<ScalarExpr> arg = rcp_dynamic_cast<ScalarExpr>(expr[0].ptr());
    TEUCHOS_TEST_FOR_EXCEPTION(arg.get()==0, std::runtime_error,
                       "non-scalar argument in pow function");
    return new NonlinearUnaryOp(arg, rcp(new PowerFunctor(p)));
  }

  /** \name Elementary math functions */
  //@{
  SUNDANCE_UNARY_OP(reciprocal, StdReciprocal, "reciprocal function")

  SUNDANCE_UNARY_OP(fabs, StdFabs, "absolute value")

  SUNDANCE_UNARY_OP(sign, StdSign, "sign function")

  SUNDANCE_UNARY_OP(exp, StdExp, "exponential function")

  SUNDANCE_UNARY_OP(log, StdLog, "logarithm")

  SUNDANCE_UNARY_OP(sqrt, StdSqrt, "square root"])

  SUNDANCE_UNARY_OP(sin, StdSin, "sine function")

  SUNDANCE_UNARY_OP(cos, StdCos, "cosine function")

  SUNDANCE_UNARY_OP(tan, StdTan, "tangent function")

  SUNDANCE_UNARY_OP(asin, StdASin, "inverse sine")

  SUNDANCE_UNARY_OP(acos, StdACos, "inverse cosine")

  SUNDANCE_UNARY_OP(atan, StdATan, "inverse tangent")

  SUNDANCE_UNARY_OP(sinh, StdSinh, "hyperbolic sine")

  SUNDANCE_UNARY_OP(cosh, StdCosh, "hyperbolic cosine")

  SUNDANCE_UNARY_OP(tanh, StdTanh, "hyperbolic tangent")

  SUNDANCE_UNARY_OP(asinh, StdASinh, "inverse hyperbolic sine")

  SUNDANCE_UNARY_OP(acosh, StdACosh, "inverse hyperbolic cosine")

  SUNDANCE_UNARY_OP(atanh, StdATanh, "inverse hyperbolic tangent")
//@}

}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
