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

#ifndef SUNDANCE_CURVEEVALMEDIATOR_H
#define SUNDANCE_CURVEEVALMEDIATOR_H

#include "SundanceDefs.hpp"
#include "SundanceMap.hpp"
#include "SundanceStdFwkEvalMediator.hpp"
#include "SundanceQuadratureFamily.hpp"
#include "SundanceBasisFamily.hpp"
#include "SundanceOrderedTuple.hpp"
#include "SundanceParametrizedCurve.hpp"

namespace Sundance
{
using namespace Teuchos;

/**
 * 
 */
class CurveEvalMediator : public StdFwkEvalMediator
{
public:


  /** */
  CurveEvalMediator(const Mesh& mesh,
	const ParametrizedCurve& paramcurve ,
    int cellDim,
    const QuadratureFamily& quad);


  /** */
  virtual ~CurveEvalMediator(){;}

      
  /** Evaluate the given coordinate expression, putting
   * its numerical values in the given EvalVector. */
  virtual void evalCoordExpr(const CoordExpr* expr,
    RCP<EvalVector>& vec) const ;
      
  /** Evaluate the given discrete function, putting
   * its numerical values in the given EvalVector. */
  virtual void evalDiscreteFuncElement(const DiscreteFuncElement* expr,
    const Array<MultiIndex>& mi,
    Array<RCP<EvalVector> >& vec) const ;

  /** Evaluate the given cell diameter expression, putting
   * its numerical values in the given EvalVector. */
  virtual void evalCellDiameterExpr(const CellDiameterExpr* expr,
    RCP<EvalVector>& vec) const ;

  /** Evaluates one component of the normal vector to a given parameterized curve
   * i.e. x,y or z component of that vector in 3D <br>
   * , this method is only in the CurveEvalMediator class implemented */
  virtual void evalCurveNormExpr(const CurveNormExpr* expr,
    RCP<EvalVector>& vec) const ;

  /** Evaluate the given cell vector expression, putting
   * its numerical values in the given EvalVector. */
  virtual void evalCellVectorExpr(const CellVectorExpr* expr,
    RCP<EvalVector>& vec) const ;
            
  /** */
  virtual void setCellType(const CellType& cellType,
    const CellType& maxCellType,
    bool isInternalBdry) ;

  /** */
  virtual void print(std::ostream& os) const ;

  /** */
  RCP<Array<Array<Array<double> > > > getFacetRefBasisVals(const BasisFamily& basis, int diffOrder) const ;

  /** */
  int numQuadPts(const CellType& cellType) const ;

  static double& totalFlops() {static double rtn = 0; return rtn;}

      

  static void addFlops(const double& flops) {totalFlops() += flops;}

  /**
   * Return the number of different cases for which reference
   * basis functions must be evaluated. If we're on maximal cells,
   * this will be one. If we're on lower-dimensional cells, this will
   * be equal to the number of cellDim-dimensional facets of the maximal
   * cells.
   */
  int numEvaluationCases() const {return numEvaluationCases_;}


  /** */
  static Time& coordEvaluationTimer() ;

private:

  /** */
  int numEvaluationCases_;

  /** */
  QuadratureFamily quad_;

  /** */
  int numQuadPtsForMaxCell_;

  /** */
  CellType maxCellType_;

  /** */
  CellType curveCellType_;

  /** */
  ParametrizedCurve paramcurve_;
      
};
}


#endif /* SUNDANCE_CURVEEVALMEDIATOR_H */
