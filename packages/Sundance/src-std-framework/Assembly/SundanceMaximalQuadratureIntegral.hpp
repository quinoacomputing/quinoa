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

#ifndef SUNDANCE_MAXIMALQUADRATUREINTEGRAL_H
#define SUNDANCE_MAXIMALQUADRATUREINTEGRAL_H

#include "SundanceDefs.hpp"
#include "SundanceQuadratureIntegralBase.hpp"

namespace Sundance
{

using namespace Teuchos;

/** 
 * Integration by quadrature over a maximal cell  
 */
class MaximalQuadratureIntegral : public ElementIntegral
{
public:
  /** Construct a zero-form to be computed by quadrature */
  MaximalQuadratureIntegral(
    const CellType& maxCellType,
    const QuadratureFamily& quad,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh,
    int verb);

  /** Construct a one form to be computed by quadrature */
  MaximalQuadratureIntegral(
    const CellType& maxCellType,
    const BasisFamily& testBasis,
    int alpha,
    int testDerivOrder,
    const QuadratureFamily& quad,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh,
    int verb);

  /** Construct a two-form to be computed by quadrature */
  MaximalQuadratureIntegral(
    const CellType& maxCellType,
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
  virtual ~MaximalQuadratureIntegral(){;}
      
     /** */
  virtual void transform(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& isLocalFlag,
    const Array<int>& facetNum,
    const RCP<Array<int> >& cellLIDs,
    const double* const coeff,
    RCP<Array<double> >& A) const 
    {
      if (order()==2) transformTwoForm(JTrans, JVol, facetNum, cellLIDs,coeff, A);
      else if (order()==1) transformOneForm(JTrans, JVol, facetNum, cellLIDs,coeff, A);
      else transformZeroForm(JTrans, JVol, isLocalFlag, facetNum, cellLIDs,coeff, A);
    }

  /** */
  virtual void transformZeroForm(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& isLocalFlag,
    const Array<int>& facetIndex,
    const RCP<Array<int> >& cellLIDs,
    const double* const coeff,
    RCP<Array<double> >& A) const ;
  
  /** */
  virtual void transformTwoForm(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& facetIndex,
    const RCP<Array<int> >& cellLIDs,
    const double* const coeff,
    RCP<Array<double> >& A) const ;
  
  /** */
  void transformOneForm(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& facetIndex,
    const RCP<Array<int> >& cellLIDs,
    const double* const coeff,
    RCP<Array<double> >& A) const ;
  
private:

  /** Do the integration by summing reference quantities over quadrature
   * points and then transforming the sum to physical quantities.  */
  void transformSummingFirst(int nCells,
    const Array<int>& facetIndex,
    const RCP<Array<int> >& cellLIDs,
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
    int unkDerivDir, int unkNode)
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
  inline double& wValue(int q, int testDerivDir, int testNode)
    {return W_[testNode + nNodesTest()*(testDerivDir + nRefDerivTest()*q)];}


  /** */
  inline const double& wValue(int q, int testDerivDir, int testNode) const 
    {return W_[testNode + nNodesTest()*(testDerivDir + nRefDerivTest()*q)];}

  /* */
  Array<double> W_;

  /* */
  bool useSumFirstMethod_;
      
  /** For ACI (ACI = Adaptive Cell Integration), store the reference integral values for one form
   * The indexes facet, quadPoints, nRefDerivTest , nNodesTest */
  Array<Array<Array<double> > >  W_ACI_F1_;

  /** For ACI (ACI = Adaptive Cell Integration), store the reference integral values for two form
   * The indexes quadPoints, nRefDerivTest , nNodesTest , nRefDerivUnk , nNodesUnk */
  Array<Array<Array<Array<Array<double> > > > > W_ACI_F2_;

  /** The quadrature family needed for special integration (ACI)*/
  QuadratureFamily quad_;

  /** The quadrature points*/
  Array<Point> quadPts_;

  /** The standard weights (in case of ACI these might change)*/
  Array<double> quadWeights_;
};
}


#endif
