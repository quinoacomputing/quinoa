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

#ifndef SUNDANCE_POINTWISEUSERDEFFUNCTOR_H
#define SUNDANCE_POINTWISEUSERDEFFUNCTOR_H

#include "SundanceDefs.hpp"
#include "SundanceMap.hpp"
#include "SundanceUserDefFunctor.hpp"
#include "SundanceMultiSet.hpp"
#include "PlayaExceptions.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"



namespace Sundance
{
  using namespace Sundance;
  using namespace Teuchos;

  
  

    /**
     * PointwiseUserDefFunctor0 is an implementation of UserDefFunctor for which
     * the user writes code to evaluate the function at a single quadrature point.
     * Looping over quadrature points is done by the this class.
     */
  class PointwiseUserDefFunctor0 : public UserDefFunctor
    {
    public:
      /** ctor */
      PointwiseUserDefFunctor0(const std::string& name, int domainDim, int rangeDim) ;

      /** */
      virtual ~PointwiseUserDefFunctor0(){;}

      /** */
      void evaluationCallback(int nPoints, int maxDiffOrder,
                              const double** in,
                              double** out) const ;

      /** */
      virtual void eval0(const double* in, double* out) const = 0 ;

      /** */
      virtual int maxOrder() const {return 0;}

    private:
    };


    /**
     * PointwiseUserDefFunctor1 is an implementation of UserDefFunctor for which
     * the user writes code to evaluate the function at a single quadrature point.
     * Looping over quadrature points is done by the this class.
     */
  class PointwiseUserDefFunctor1 : public PointwiseUserDefFunctor0
    {
    public:
      /** ctor */
      PointwiseUserDefFunctor1(const std::string& name, int domainDim, int rangeDim) ;

      /** */
      virtual ~PointwiseUserDefFunctor1(){;}

      /** */
      void evaluationCallback(int nPoints, int maxDiffOrder,
                              const double** in,
                              double** out) const ;

      /** */
      virtual void eval0(const double* in, double* out) const ;

      /** */
      virtual void eval1(const double* in, double* outVals, double* outDerivs) const = 0 ;

      /** */
      virtual int maxOrder() const {return 1;}

    private:
    };


    /**
     * PointwiseUserDefFunctor2 is an implementation of UserDefFunctor for which
     * the user writes code to evaluate the function at a single quadrature point.
     * Looping over quadrature points is done by the this class.
     */
  class PointwiseUserDefFunctor2 : public PointwiseUserDefFunctor1
    {
    public:
      /** ctor */
      PointwiseUserDefFunctor2(const std::string& name, int domainDim, int rangeDim) ;

      /** */
      virtual ~PointwiseUserDefFunctor2(){;}

      /** */
      void evaluationCallback(int nPoints, int maxDiffOrder,
                              const double** in,
                              double** out) const ;

      /** */
      virtual void eval0(const double* in, double* out) const ;

      /** */
      virtual void eval1(const double* in, double* outVals, double* outDerivs) const ;

      virtual void eval2(const double* in, double* outVals, double* outDerivs,
                         double* outDerivs2) const = 0 ;

      /** */
      virtual int maxOrder() const {return 2;}

    private:
    };


}


#endif
