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

#ifndef SUNDANCE_POLYGONQUADRATURE_H
#define SUNDANCE_POLYGONQUADRATURE_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceQuadratureFamilyBase.hpp"
#include "SundanceQuadratureFamily.hpp"

namespace Sundance
{
using namespace Teuchos;


/** 
 * Quadrature special for the polygon curve integration <br>.
 * There we need for each line segment inside one cell a line integral.
 * So for one cell we might need several line quadrature points, and this class
 * gives us per cell(line) a constant number of quadrature points for a line,
 * so we can integrate along the curve exactly. (each line segment exactly inside one cell) <br>
 * As an input we need a quadrature class, and the points of that quadrature rule for a line
 * will be used for each line segment of the polygon inside one cell (2D).
 */
class PolygonQuadrature : public QuadratureFamilyBase
{
public:
  /** */
  PolygonQuadrature( const QuadratureFamily& quad );

  /** */
  virtual ~PolygonQuadrature(){;}


  /** */
  virtual XMLObject toXML() const ;

  /** Describable interface */
  virtual std::string description() const 
    {return "PolygonQuadrature[order=" + Teuchos::toString(order())
        +  "]";}

  /* handleable boilerplate */
  GET_RCP(QuadratureFamilyStub);

  /** method for the user to change this value
   * too much line nr is an overhead
   * too less will cause error thrown in the code
   * @param maxNrLine */
  static void setNrMaxLinePerCell( int maxNrLine ) { nrMaxLinePerCell_ = maxNrLine; }

  /** return the maximal number of line segments inside one cell */
  static int getNrMaxLinePerCell() { return nrMaxLinePerCell_; }

protected:

  /** for polygon curve integrals only this method should be used */
  virtual void getLineRule(Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;

private:

  /** The quadrature which will be used for the polygon lines integration */
  const QuadratureFamily& quad_;

  /** maximum number of lines inside a cell, the default value is 6 but this can be changed
   * by the user */
  static int nrMaxLinePerCell_;

};
}


#endif
