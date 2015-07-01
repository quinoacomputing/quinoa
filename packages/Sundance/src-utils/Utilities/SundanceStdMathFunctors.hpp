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

#ifndef SUNDANCE_STDMATHFUNCTORS_H
#define SUNDANCE_STDMATHFUNCTORS_H

#include "SundanceDefs.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceUnaryFunctor.hpp"
#ifdef _MSC_VER
# include "winmath.h"
#endif


namespace Sundance
{
  using namespace Teuchos;

  /** */
  class PowerFunctor : public UnaryFunctor
  {
  public:
    /** */
    PowerFunctor(const double& p);
    
    /** Evaluate power function and deriv at an array of values */ 
    virtual void eval1(const double* const x, 
              int nx, 
              double* f, 
              double* df) const ;
    /** Evaluate power function at an array of values */ 
    virtual void eval0(const double* const x, int nx, double* f) const ;

    /** Evaluate power function and first two derivs at an array of values */
    virtual void eval2(const double* const x, 
                      int nx, 
                      double* f, 
                      double* df_dx,
                      double* d2f_dxx) const ;

    /** Evaluate power function and first three derivs at an array of values */
    virtual void eval3(const double* const x, 
                      int nx, 
                      double* f, 
                      double* df_dx,
                      double* d2f_dxx,
                      double* d3f_dxxx) const ;

    static RCP<FunctorDomain> powerDomain(const double& p);

  protected:
    bool acceptX(int diffOrder, const double& x) const
    {
      if (powerIsInteger_)
	{
	  return p_>=0.0 || x!=0.0;
	}
      else
	{
	  if (x<0.0) return false; 
	  if (x==0.0) return (p_ > diffOrder);
	}
      return true;      
    }
  private:
    double p_;
    bool powerIsInteger_;
  };


  
  

  SUNDANCE_UNARY_FUNCTOR(reciprocal, StdReciprocal, "reciprocal function", 
                         NonzeroDomain(), 1.0/x[i], -f[i]*f[i], -2.0*df[i]/x[i])

    SUNDANCE_UNARY_FUNCTOR(fabs, StdFabs, "absolute value", UnboundedDomain(), ::fabs(x[i]), ((x[i]>=0.0) ? x[i] : -x[i]), 0.0)

  SUNDANCE_UNARY_FUNCTOR(sign, StdSign, "sign function", UnboundedDomain(), 
                         ((x[i]>0.0) ? 1.0 : ( (x[i]<0.0) ? -1.0 : 0.0)), 
                         0.0, 0.0)

    SUNDANCE_UNARY_FUNCTOR3(exp, StdExp, "exponential function", UnboundedDomain(), ::exp(x[i]), f[i], f[i], f[i])

    SUNDANCE_UNARY_FUNCTOR3(log, StdLog, "logarithm", PositiveDomain(), ::log(x[i]), 1.0/x[i], -df[i]*df[i], -2.0*d2f[i]*df[i])

  SUNDANCE_UNARY_FUNCTOR(sqrt, StdSqrt, "square root", PositiveDomain(), ::sqrt(x[i]), 0.5/f[i], -0.5*df[i]/x[i])

    SUNDANCE_UNARY_FUNCTOR3(sin, StdSin, "sine function", UnboundedDomain(), ::sin(x[i]), ::cos(x[i]), -f[i], -df[i])

    SUNDANCE_UNARY_FUNCTOR3(cos, StdCos, "cosine function", UnboundedDomain(), ::cos(x[i]), -::sin(x[i]), -f[i], -df[i])

  SUNDANCE_UNARY_FUNCTOR(tan, StdTan, "tangent function", UnboundedDomain(),
                         ::tan(x[i]), 1.0 + f[i]*f[i], 2.0*f[i]*df[i])

  SUNDANCE_UNARY_FUNCTOR(asin, StdASin, "inverse sine", 
                         BoundedDomain(-1.0, 1.0),
                         ::asin(x[i]), 1.0/::sqrt(1.0-x[i]*x[i]),
                         x[i]*df[i]*df[i]*df[i])

  SUNDANCE_UNARY_FUNCTOR(acos, StdACos, "inverse cosine",
                         BoundedDomain(-1.0, 1.0), 
                         ::acos(x[i]), -1.0/::sqrt(1.0-x[i]*x[i]),
                         x[i]*df[i]*df[i]*df[i])

  SUNDANCE_UNARY_FUNCTOR(atan, StdATan, "inverse tangent", 
                         UnboundedDomain(),
                         ::atan(x[i]), 1.0/(1.0 + x[i]*x[i]),
                         -2.0*x[i]*df[i]*df[i])

  SUNDANCE_UNARY_FUNCTOR(sinh, StdSinh, "hyperbolic sine",
                         UnboundedDomain(),
                         ::sinh(x[i]), ::cosh(x[i]), f[i])

  SUNDANCE_UNARY_FUNCTOR(cosh, StdCosh, "hyperbolic cosine",
                         UnboundedDomain(),
                         ::cosh(x[i]), ::sinh(x[i]), f[i])

  SUNDANCE_UNARY_FUNCTOR(tanh, StdTanh, "hyperbolic tangent",
                         UnboundedDomain(),
                         ::tanh(x[i]), 1.0 - f[i]*f[i], -2.0*f[i]*df[i])

  SUNDANCE_UNARY_FUNCTOR(asinh, StdASinh, "inverse hyperbolic sine",
                         UnboundedDomain(),
                         ::asinh(x[i]), 1.0/::sqrt(1.0 + x[i]*x[i]),
                         -x[i]*df[i]*df[i]*df[i])

  SUNDANCE_UNARY_FUNCTOR(acosh, StdACosh, "inverse hyperbolic cosine",
                         LowerBoundedDomain(1.0),
                         ::acosh(x[i]), 1.0/::sqrt(x[i]*x[i]-1.0),
                         -x[i]*df[i]*df[i]*df[i])

  SUNDANCE_UNARY_FUNCTOR(atanh, StdATanh, "inverse hyperbolic tangent",
                         BoundedDomain(-1.0, 1.0), 
                         ::atanh(x[i]), 1.0/(1.0 - x[i]*x[i]),
                         2.0*x[i]*df[i]*df[i])


}

                  

#endif
