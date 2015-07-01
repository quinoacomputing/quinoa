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

#ifndef SUNDANCE_FUNCTORDOMAIN_H
#define SUNDANCE_FUNCTORDOMAIN_H

#include "SundanceDefs.hpp"
#include "Teuchos_StrUtils.hpp"

namespace Sundance
{
  using namespace Teuchos;
  using std::string;

  class FunctorDomain
  {
  public:
    FunctorDomain();

    virtual ~FunctorDomain(){;}

    virtual bool hasLowerBound() const {return false;}

    virtual double lowerBound() const ;

    virtual bool hasUpperBound() const {return false;}

    virtual double upperBound() const ;

    virtual bool hasExcludedPoint() const {return false;}

    virtual double excludedPoint() const ;

    virtual string description() const = 0 ;
  };

  class UnboundedDomain : public FunctorDomain
  {
  public:
    UnboundedDomain();

    string description() const {return "UnboundedDomain()";}
  };


  class PositiveDomain : public FunctorDomain
  {
  public:
    PositiveDomain();

    bool hasLowerBound() const {return true;}
    
    double lowerBound() const {return 0.0;}
    
    string description() const {return "PositiveDomain()";}
  };

  class StrictlyPositiveDomain : public FunctorDomain
  {
  public:
    StrictlyPositiveDomain();

    bool hasLowerBound() const {return true;}
    
    double lowerBound() const {return 0.0;}
    
    bool hasExcludedPoint() const {return true;}
    
    double excludedPoint() const {return 0.0;}

    string description() const {return "StrictlyPositiveDomain()";}

  };


  class BoundedDomain : public FunctorDomain
  {
  public:
    BoundedDomain(const double& lower, const double& upper);

    bool hasLowerBound() const {return true;}
    
    double lowerBound() const {return lower_;}
    
    bool hasUpperBound() const {return true;}
    
    double upperBound() const {return upper_;}

    string description() const {return "BoundedDomain("
	+ Teuchos::toString(lowerBound()) + ", "
	+ Teuchos::toString(upperBound()) + ")";}

  private:
    double lower_;

    double upper_;
  };


  class LowerBoundedDomain : public FunctorDomain
  {
  public:
    LowerBoundedDomain(const double& lower);

     bool hasLowerBound() const {return true;}

     double lowerBound() const {return lower_;}

    string description() const {return "LowerBoundedDomain("
	+ Teuchos::toString(lowerBound()) + ")";}



  private:
    double lower_;
  };

class NonzeroDomain : public FunctorDomain
  {
  public:
    NonzeroDomain();

    bool hasExcludedPoint() const {return true;}
    
    double excludedPoint() const {return 0.0;}

    string description() const {return "NonzeroDomain()";}


  };

  inline std::ostream& operator<<(std::ostream& os, const FunctorDomain& f)
  {
    os << f.description();
    return os;
  }
}



#endif
