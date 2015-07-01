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

#ifndef SUNDANCEDUMMYPARAMETRIZEDCURVE_H_
#define SUNDANCEDUMMYPARAMETRIZEDCURVE_H_

#include "SundanceCurveBase.hpp"

namespace Sundance
{

/** Dummy curve */
class DummyParametrizedCurve : public CurveBase
{
public:

	DummyParametrizedCurve(){;}

	virtual ~DummyParametrizedCurve(){;}

	/** */
	virtual Expr getParams() const { return Expr();}

	/** */
	virtual void getIntegrationParams(Array<double>& alphas) const {
		alphas.resize(2);
		alphas[0] = 1.0;
		alphas[1] = 1.0;
	}

	/** */
	virtual double integrationParameter(const Point& evaluationPoint) const{ return 1.0;}

protected:
	virtual double curveEquation_intern(const Point& evaluationPoint) const { return 1.0; }
public:

	virtual void returnIntersectPoints(const Point& startEdgePoint, const Point& endEdgePoint,
			                      int& nrPoints ,Array<Point>& result) const {;}

	virtual void returnIntersect(const Point& startEdgePoint, const Point& endEdgePoint,
			                      int& nrPoints ,Array<double>& result) const {;}

	/** Return a ref count pointer to self */
	virtual RCP<CurveBase> getRcp() {return rcp(this);}

	virtual bool isCurveValid() const { return false; }


};

}
#endif /* SUNDANCEDUMMYPARAMETRIZEDCURVE_H_ */
