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

#ifndef SUNDANCE_QUADRATUREFAMILYBASE_H
#define SUNDANCE_QUADRATUREFAMILYBASE_H

#include "SundanceDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Array.hpp"
#include "SundanceQuadratureFamilyStub.hpp"
#include "SundanceCellType.hpp"
#include "SundancePoint.hpp"
#include "SundanceParametrizedCurve.hpp"
#include "SundanceMesh.hpp"

namespace Sundance
{
using namespace Teuchos;

/** 
 * QuadratureFamilyBase extends QuadratureFamilyStub to provide
 * an interface for getting quadrature points for a given cell type. 
 */
class QuadratureFamilyBase : public QuadratureFamilyStub 
{
public:
  /** */
  QuadratureFamilyBase(int order) : QuadratureFamilyStub(order) {;}

  /** */
  virtual ~QuadratureFamilyBase(){;}

  /** Gets number of points associated with a particular cell type:
      WARNING: this is slow.  Call it once and store the result. 
      TODO: make it pure virtual and override with queries in
      the derived classes, making them supply the information.  */
  virtual int getNumPoints( const CellType &cellType ) const 
    {
      Array<Point> qp;
      Array<double> qw;
      this->getPoints(cellType,qp,qw);
      return qp.size();
    }

  /** Get the quadrature points and weights for the given cell type */
  virtual void getPoints(const CellType& cellType, 
    Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;
      
  /** This method is used for the Adaptive Cell Integration, which returns
   * special quadrature weights for cells which are cut by the curve <br>
   * The quadPoints and the quadWeights are set to the default values.
   * (quadrature without the curve , no curve)
   * quadWeights will be changed depending on the curve, if the curve cuts the cell
   * isCut should be set by this method, if it is true then the quadWeights will be
   * used for quadrature of the cell */
  virtual void getAdaptedWeights(const CellType& cellType ,
	int cellDim,
	int celLID ,
    int facetIndex ,
    const Mesh& mesh ,
    const ParametrizedCurve& globalCurve ,
    Array<Point>& quadPoints ,
    Array<double>& quadWeights,
    bool& isCut) const ;

protected:

  /** compute a rule for the reference line cell */
  virtual void getLineRule(Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;

  /** compute a rule for the reference triangle cell */
  virtual void getTriangleRule(Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;

  /** compute a rule for the reference quad cell */
  virtual void getQuadRule(Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;

  /** compute a rule for the reference tet cell */
  virtual void getTetRule(Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;


  /** compute a rule for the reference brick cell */
  virtual void getBrickRule(Array<Point>& quadPoints,
    Array<double>& quadWeights) const ;
private:
};

}

#endif

