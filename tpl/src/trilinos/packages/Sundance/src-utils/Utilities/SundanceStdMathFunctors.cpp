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

#include <math.h>
#include "SundanceStdMathFunctors.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_Utils.hpp"

using namespace Sundance;
using namespace Teuchos;
using std::endl;

PowerFunctor::PowerFunctor(const double& p) 
  : UnaryFunctor("pow("+Teuchos::toString(p)+")",
		 powerDomain(p)), 
    p_(p),
    powerIsInteger_(p==floor(p))
{}

RCP<FunctorDomain> PowerFunctor::powerDomain(const double& p)
{
  /* There are four cases: 
   * (1) p is a positive integer, in which case the domain is all of R.
   * (2) p is a negative integer, in which case the domain is R\0.
   * (3) p is a positive real \notin Z, in which case the domain is [0,\infty).
   * (3) p is a negative real \notin Z, in which case the domain is (0,\infty).
   */
  bool isInZ = floor(p)==p;
  bool isNegative = p < 0.0;

  if (isInZ)
    {
      if (isNegative) return rcp(new NonzeroDomain());
    }
  else
    {
      if (isNegative) return rcp(new StrictlyPositiveDomain());
      return rcp(new PositiveDomain());
    }
  return rcp(new UnboundedDomain());
}

void PowerFunctor::eval1(const double* const x, 
                        int nx, 
                        double* f, 
                        double* df) const
{
  if (p_==2)
    {
      for (int i=0; i<nx; i++) 
        {
          df[i] = 2.0*x[i];
          f[i] = x[i]*x[i];
        }
    }
  else if (p_==1)
    {
      for (int i=0; i<nx; i++) 
        {
          df[i] = 1.0;
          f[i] = x[i];
        }
    }
  else if (p_==0)
    {
      for (int i=0; i<nx; i++) 
        {
          df[i] = 0.0;
          f[i] = 1.0;
        }
    }
  else
    {
      if (checkResults())
        {
          for (int i=0; i<nx; i++) 
            {
	      TEUCHOS_TEST_FOR_EXCEPTION(!acceptX(1,x[i]), std::runtime_error,
				 "first deriv of pow(" << x[i] 
				 << ", " << p_ << ") "
				 "is undefined");
	      
              double px = ::pow(x[i], p_-1);
              df[i] = p_*px;
              f[i] = x[i]*px;
              //bvbw tried to include math.h, without success
#ifdef REDDISH_PORT_PROBLEM
              TEUCHOS_TEST_FOR_EXCEPTION(fpclassify(f[i]) != FP_NORMAL 
                                 || fpclassify(df[i]) != FP_NORMAL,
                                 std::runtime_error,
                                 "Non-normal floating point result detected in "
                                 "evaluation of unary functor " << name());
#endif
            }
        }
      else 
        {
	  for (int i=0; i<nx; i++) 
	    {
	      TEUCHOS_TEST_FOR_EXCEPTION(!acceptX(1,x[i]), std::runtime_error,
				 "first deriv of pow(" << x[i] 
				 << ", " << p_ << ") "
				 "is undefined");
	      double px = ::pow(x[i], p_-1);
	      df[i] = p_*px;
	      f[i] = x[i]*px;
	    }
        }
    }
}




void PowerFunctor::eval3(const double* const x, 
                         int nx, 
                         double* f, 
                         double* df,
                         double* d2f_dxx,
                         double* d3f_dxxx) const
{
  if (p_==3)
    {
      for (int i=0; i<nx; i++) 
        {
          d3f_dxxx[i] = 6.0;
          d2f_dxx[i] = 6.0*x[i];
          df[i] = 3.0*x[i]*x[i];
          f[i] = x[i]*x[i]*x[i];
        }
    }
  else if (p_==2)
    {
      for (int i=0; i<nx; i++) 
        {
          d3f_dxxx[i] = 0.0;
          d2f_dxx[i] = 2.0;
          df[i] = 2.0*x[i];
          f[i] = x[i]*x[i];
        }
    }
  else if (p_==1)
    {
       for (int i=0; i<nx; i++) 
        {
          d3f_dxxx[i] = 0.0;
          d2f_dxx[i] = 0.0;
          df[i] = 1.0;
          f[i] = x[i];
        }
    }
  else if (p_==0)
    {
      for (int i=0; i<nx; i++) 
        {
          d3f_dxxx[i] = 0.0;
          d2f_dxx[i] = 0.0;
          df[i] = 0.0;
          f[i] = 1.0;
        }
    }
  else
    {
      if (checkResults())
        {
          for (int i=0; i<nx; i++) 
            {
              double px = ::pow(x[i], p_-3);
              d3f_dxxx[i] = p_ * (p_-1) * (p_-2) * px;
              d2f_dxx[i] = p_ * (p_-1) * x[i] * px;
              df[i] = p_*x[i]*x[i]*px;
              f[i] = x[i]*x[i]*x[i]*px;
	      TEUCHOS_TEST_FOR_EXCEPTION(!acceptX(3,x[i]), std::runtime_error,
				 "third deriv of pow(" << x[i] 
				 << ", " << p_ << ") "
				 "is undefined");


#ifdef REDDISH_PORT_PROBLEM
              TEUCHOS_TEST_FOR_EXCEPTION(fpclassify(f[i]) != FP_NORMAL 
                                 || fpclassify(df[i]) != FP_NORMAL,
                                 std::runtime_error,
                                 "Non-normal floating point result detected in "
                                 "evaluation of unary functor " << name());
#endif
            }
        }
      else
        {
          for (int i=0; i<nx; i++) 
            {
	      TEUCHOS_TEST_FOR_EXCEPTION(!acceptX(3,x[i]), std::runtime_error,
				 "third deriv of pow(" << x[i] 
				 << ", " << p_ << ") "
				 "is undefined");

              double px = ::pow(x[i], p_-3);
              d3f_dxxx[i] = p_ * (p_-1) * (p_-2) * px;
              d2f_dxx[i] = p_ * (p_-1) * x[i] * px;
              df[i] = p_*x[i]*x[i]*px;
              f[i] = x[i]*x[i]*x[i]*px;
            }
        }
    }
}


void PowerFunctor::eval2(const double* const x, 
                        int nx, 
                        double* f, 
                        double* df,
                        double* d2f_dxx) const
{
  if (p_==2)
    {
      for (int i=0; i<nx; i++) 
        {
          d2f_dxx[i] = 2.0;
          df[i] = 2.0*x[i];
          f[i] = x[i]*x[i];
        }
    }
  else if (p_==1)
    {
       for (int i=0; i<nx; i++) 
        {
          d2f_dxx[i] = 0.0;
          df[i] = 1.0;
          f[i] = x[i];
        }
    }
  else if (p_==0)
    {
      for (int i=0; i<nx; i++) 
        {
          d2f_dxx[i] = 0.0;
          df[i] = 0.0;
          f[i] = 1.0;
        }
    }
  else
    {
      if (checkResults())
        {
          for (int i=0; i<nx; i++) 
            {
	      TEUCHOS_TEST_FOR_EXCEPTION(!acceptX(2,x[i]), std::runtime_error,
				 "second deriv of pow(" << x[i] 
				 << ", " << p_ << ") "
				 "is undefined");


              double px = ::pow(x[i], p_-2);
              d2f_dxx[i] = p_ * (p_-1) * px;
              df[i] = p_*x[i]*px;
              f[i] = x[i]*x[i]*px;
#ifdef REDDISH_PORT_PROBLEM
              TEUCHOS_TEST_FOR_EXCEPTION(fpclassify(f[i]) != FP_NORMAL 
                                 || fpclassify(df[i]) != FP_NORMAL,
                                 std::runtime_error,
                                 "Non-normal floating point result detected in "
                                 "evaluation of unary functor " << name());
#endif
            }
        }
      else
        {
	  for (int i=0; i<nx; i++) 
	    {
	      TEUCHOS_TEST_FOR_EXCEPTION(!acceptX(2,x[i]), std::runtime_error,
				 "second deriv of pow(" << x[i] 
				 << ", " << p_ << ") "
				 "is undefined");
	      
	      double px = ::pow(x[i], p_-2);
	      
	      d2f_dxx[i] = p_ * (p_-1) * px;
	      df[i] = p_*x[i]*px;
	      f[i] = x[i]*x[i]*px;
	    }
	}
    }
}



void PowerFunctor::eval0(const double* const x, 
                        int nx, 
                        double* f) const
{
  if (checkResults())
    {
      for (int i=0; i<nx; i++) 
        {
	  TEUCHOS_TEST_FOR_EXCEPTION(!acceptX(0,x[i]), std::runtime_error,
			     "pow(" << x[i] 
			     << ", " << p_ << ") "
			     "is undefined");


          f[i] = ::pow(x[i], p_);
#ifdef REDDISH_PORT_PROBLEM
          TEUCHOS_TEST_FOR_EXCEPTION(fpclassify(f[i]) != FP_NORMAL, 
                             std::runtime_error,
                             "Non-normal floating point result detected in "
                             "evaluation of unary functor " << name());
#endif
	}
    }
  else
    {
      for (int i=0; i<nx; i++) 
	{
	  TEUCHOS_TEST_FOR_EXCEPTION(!acceptX(0,x[i]), std::runtime_error,
			     "pow(" << x[i] 
			     << ", " << p_ << ") "
			     "is undefined");

	  f[i] = ::pow(x[i], p_);
	}
    }
}
