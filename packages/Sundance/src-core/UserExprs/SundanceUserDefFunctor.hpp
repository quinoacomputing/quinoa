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

#ifndef SUNDANCE_USERDEFFUNCTOR_H
#define SUNDANCE_USERDEFFUNCTOR_H

#include "SundanceDefs.hpp"
#include "SundanceMap.hpp"
#include "SundanceMultiSet.hpp"
#include "PlayaExceptions.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"



namespace Sundance
{
using namespace Sundance;
using namespace Teuchos;




class EvalVector;
class EvalManager;

/**
 * UserDefFunctor defines an interface for callbacks used to implement
 * user-defined nonlinear operators in the Sundance Expr system.
 */
class UserDefFunctor
{
public:
  /** ctor */
  UserDefFunctor(const std::string& name, int domainDim, int rangeDim) ;

  /** */
  virtual ~UserDefFunctor(){;}

  /** */
  const std::string& name(int elemIndex) const {return elemNames_[elemIndex];}

  /** */
  const std::string& name() const {return name_;}


  /** */
  virtual void evaluationCallback(int nPoints, int maxDiffOrder,
    const double** in,
    double** out) const = 0 ;

  /** */
  virtual void eval0(const Array<double>& in, double* outVals) const ;

  /**
   * Evaluate the expression and its derivative. The values should be put into
   * the outVals array. The derivatives should be put into the outDerivs array,
   * ordered with the domain index running fastest. That is, 
   * \f[
   * outDerivs[i*N_R + j] = \frac{\partial F_i}{\partial q_j}
   * \f]
   */
  virtual void eval1(const Array<double>& in, double* outVals, 
    double* outDerivs) const ;

    

  /** */
  int domainDim() const {return domainDim_;}

  /** */
  int rangeDim() const {return rangeDim_;}

  /** */
  virtual int maxOrder() const = 0 ;

  /** */
  void reset() const ;

private:
  const std::string name_;
  Array<string> elemNames_;
  const int domainDim_;
  const int rangeDim_;
};


}


#endif
