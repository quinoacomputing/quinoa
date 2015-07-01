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

#ifndef SUNDANCE_ADFIELD_H
#define SUNDANCE_ADFIELD_H

#include "SundanceDefs.hpp"
#include "SundanceADBasis.hpp"

namespace SundanceTesting
{
  using namespace Sundance;
  using namespace Teuchos;
  using namespace Sundance;
  using namespace Sundance;

  /** 
   *
   */
  class ADField
  {
  public:
    /** */
    ADField(){;}

    /** */
    ADField(const ADBasis& basis, const double& coeff);

    /** */
    static Point& evalPoint() {static Point x(0.123, 0.456, 0.789); return x;}


    /** */
    void setCoeff(const double& c);

    /** */
    const ADBasis& basis() const {return basis_;}

    /** */
    ADReal evaluate() const ;

    /** */
    double coeff() const {return *coeff_;}

    double value() const {return evaluate().value();}

    ADReal operator+(const ADReal& x) const ;

    ADReal operator+(const double& x) const ;

    ADReal operator+(const ADField& x) const ;



    ADReal operator-(const ADReal& x) const ;

    ADReal operator-(const double& x) const ;

    ADReal operator-(const ADField& x) const ;



    ADReal operator*(const ADReal& x) const ;

    ADReal operator*(const double& x) const ;

    ADReal operator*(const ADField& x) const ;



    ADReal operator/(const ADReal& x) const ;

    ADReal operator/(const double& x) const ;

    ADReal operator/(const ADField& x) const ;



    ADReal operator-() const ;

    ADReal reciprocal() const ;
    

    
  private:
    ADBasis basis_;
    RCP<double> coeff_;
  };


  inline ADReal operator+(const ADReal& x, const ADField& y)
  {
    return y + x;
  }

  inline ADReal operator+(const double& x, const ADField& y)
  {
    return y + x;
  }

  inline ADReal operator-(const ADReal& x, const ADField& y)
  {
    return -y + x;
  }

  inline ADReal operator-(const double& x, const ADField& y)
  {
    return -y + x;
  }

  inline ADReal operator*(const ADReal& x, const ADField& y)
  {
    return y * x;
  }

  inline ADReal operator*(const double& x, const ADField& y)
  {
    return y * x;
  } 

  inline ADReal operator/(const ADReal& x, const ADField& y)
  {
    return x * y.reciprocal();
  }

  inline ADReal operator/(const double& x, const ADField& y)
  {
    return x * y.reciprocal();
  } 
  
  inline ADReal sin(const ADField& x)
  {
    return sin(x.evaluate());
  }
  
  inline ADReal cos(const ADField& x)
  {
    return cos(x.evaluate());
  }
  
  inline ADReal pow(const ADField& x, const double& y)
  {
    return pow(x.evaluate(), y);
  }
  
  inline ADReal exp(const ADField& x)
  {
    return exp(x.evaluate());
  }
  
  inline ADReal log(const ADField& x)
  {
    return log(x.evaluate());
  }
}



#endif
