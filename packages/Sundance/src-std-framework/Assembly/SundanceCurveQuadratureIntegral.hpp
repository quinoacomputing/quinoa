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

#ifndef SUNDANCE_CURVEQUADRATUREINTEGRAL_H
#define SUNDANCE_CURVEQUADRATUREINTEGRAL_H

#include "SundanceDefs.hpp"
#include "SundanceQuadratureIntegralBase.hpp"

namespace Sundance
{

using namespace Teuchos;

/** 
 * Integration by quadrature over a maximal cell  
 */
class CurveQuadratureIntegral : public ElementIntegral
{
public:
  /** Construct a zero-form to be computed by quadrature */
  CurveQuadratureIntegral(
    const CellType& maxCellType,
    const bool isConstantIntegral,
    const QuadratureFamily& quad,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh,
    int verb);

  /** Construct a one form to be computed by quadrature */
  CurveQuadratureIntegral(
    const CellType& maxCellType,
    const bool isConstantIntegral,
    const BasisFamily& testBasis,
    int alpha,
    int testDerivOrder,
    const QuadratureFamily& quad,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh,
    int verb);

  /** Construct a two-form to be computed by quadrature */
  CurveQuadratureIntegral(
    const CellType& maxCellType,
    const bool isConstantIntegral,
    const BasisFamily& testBasis,
    int alpha,
    int testDerivOrder,
    const BasisFamily& unkBasis,
    int beta,
    int unkDerivOrder,
    const QuadratureFamily& quad,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh,
    int verb);

  /** virtual dtor */
  virtual ~CurveQuadratureIntegral(){;}
      
     /** */
  virtual void transform(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& isLocalFlag,
    const Array<int>& facetNum,
    const RCP<Array<int> >& cellLIDs,
    const double constCoeff,
    const double* const coeff,
    RCP<Array<double> >& A) const 
    {
      if (order()==2) transformTwoForm(JTrans, JVol, facetNum, cellLIDs, constCoeff, coeff, A);
      else if (order()==1) transformOneForm(JTrans, JVol, facetNum, cellLIDs, constCoeff, coeff, A);
      else transformZeroForm(JTrans, JVol, isLocalFlag, facetNum, cellLIDs, constCoeff, coeff, A);
    }

  /** */
  virtual void transformZeroForm(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& isLocalFlag,
    const Array<int>& facetIndex,
    const RCP<Array<int> >& cellLIDs,
    const double constCoeff,
    const double* const coeff,
    RCP<Array<double> >& A) const ;
  
  /** */
  virtual void transformTwoForm(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& facetIndex,
    const RCP<Array<int> >& cellLIDs,
    const double constCoeff,
    const double* const coeff,
    RCP<Array<double> >& A) const ;
  
  /** */
  void transformOneForm(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& facetIndex,
    const RCP<Array<int> >& cellLIDs,
    const double constCoeff,
    const double* const coeff,
    RCP<Array<double> >& A) const ;
  
private:

  /** updates the reference cell information (quadPoints,derivatives,normals)*/
  void updateRefCellInformation(int maxCellLID , const ParametrizedCurve& curve) const;

  /** updates the W_ for the given cell for one form Integral */
  void updateRefCellIntegralOneForm(int maxCellLID , int cellInBatch) const ;

  /** updates the W_ for the given cell for two form Integral */
  void updateRefCellIntegralTwoForm(int maxCellLID , int cellInBatch) const ;

  /** Do the integration by summing reference quantities over quadrature
   * points and then transforming the sum to physical quantities.  */
  void transformSummingFirst(int nCells,
	const CellJacobianBatch& JVol,
    const Array<int>& facetIndex,
    const RCP<Array<int> >& cellLIDs,
    const double constCoeff,
    const double* const GPtr,
    const double* const coeff,
    RCP<Array<double> >& A) const ;

  /** Do the integration by transforming to physical coordinates 
   * at each quadrature point, and then summing */
  void transformSummingLast(int nCells,
    const Array<int>& facetIndex,
    const RCP<Array<int> >& cellLIDs,
    const double* const GPtr,
    const double* const coeff,
    RCP<Array<double> >& A) const ;

  /** Determine whether to do this batch of integrals using the
   * sum-first method or the sum-last method */
  bool useSumFirstMethod() const {return useSumFirstMethod_;}
      
  /** */
  inline double& wValue(int q, int testDerivDir, int testNode,
    int unkDerivDir, int unkNode) const
    {return W_[unkNode
        + nNodesUnk()
        *(testNode + nNodesTest()
          *(unkDerivDir + nRefDerivUnk()
            *(testDerivDir + nRefDerivTest()*q)))];}

      

  /** */
  inline const double& wValue(int facetCase, 
    int q, 
    int testDerivDir, int testNode,
    int unkDerivDir, int unkNode) const 
    {
      return W_[unkNode
        + nNodesUnk()
        *(testNode + nNodesTest()
          *(unkDerivDir + nRefDerivUnk()
            *(testDerivDir + nRefDerivTest()*q)))];
    }
      
  /** */
  inline double& wValue(int q, int testDerivDir, int testNode) const
    {return W_[testNode + nNodesTest()*(testDerivDir + nRefDerivTest()*q)];}

  /** The quadrature family needed for curve integration */
  QuadratureFamily quad_;

  /** The quadrature points*/
  mutable Array<Point> quadPts_;

  /** The standard weights */
  Array<double> quadWeights_;

  /* this must be changeable, because each cell will have different reference cell values <br>
   * here we store the result from the reference basis evaluation*/
  mutable Array<double> W_;

  /* */
  bool useSumFirstMethod_;

  /** weather to use constant coefficients */
  bool useConstCoeff_;

  /** The derivative of the curve at the curve */
  mutable Array<Point> quadCurveDerivs_;

  /** The derivative of the curve at the curve */
  mutable Array<Point> quadCurveNormals_;

};
}


#endif /* SUNDANCE_CURVEQUADRATUREINTEGRAL_H */
