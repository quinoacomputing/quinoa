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


/*
 * SundanceCurveExpr.hpp
 *
 *  Created on: Jun 7, 2011
 *      Author: benk
 */

#ifndef SUNDANCECURVEEXPR_HPP_
#define SUNDANCECURVEEXPR_HPP_

#include "SundancePointwiseUserDefFunctor.hpp"
#include "SundanceParametrizedCurve.hpp"

namespace Sundance{

/** class which serves as expression for interface (curve or surface) values */
class CurveExpr : public PointwiseUserDefFunctor0
{
public:
    /** Ctor the input is a ParamCurve and an index for the scalar field */
	CurveExpr(const ParametrizedCurve& interface , int scalarFieldIndex)
    : PointwiseUserDefFunctor0("CurveExpr", interface.getCurveDim()+1 , 1)
      , interface_(interface) , scalarFieldIndex_(scalarFieldIndex){;}

    /** to evaluate the expression just pass the control to the interface object */
    void eval0(const double* vars, double* f) const {
    	interface_.eval0( vars, f , scalarFieldIndex_ ) ;
    }

private:
  const ParametrizedCurve& interface_;
  int scalarFieldIndex_;
};

//Expr CurveExpression(const ParametrizedCurve& interface , int scalarFieldIndex, const Expr& x, const Expr& y)
//{
//  return new UserDefOp( List(x,y), rcp(new CurveExpr(interface, scalarFieldIndex)));
//}

}

#endif /* SUNDANCECURVEEXPR_HPP_ */
