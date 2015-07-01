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

#ifndef SUNDANCE_REDUCED_INTEGRAL_H
#define SUNDANCE_REDUCED_INTEGRAL_H

#include "SundanceDefs.hpp"
#include "SundanceElementIntegral.hpp"

namespace Sundance
{

using namespace Teuchos;

/** 
 * 
 */
class ReducedIntegral : public ElementIntegral
{
public:
  /** Construct a zero-form to be computed by reference integration
   * with coefficients that are piecewise constant */
  ReducedIntegral(int spatialDim,
    const CellType& maxCellType,
    int dim, 
    const CellType& cellType,
    const QuadratureFamily& quad,
    bool isInternalBdry,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh,
    int verb);

  /** Construct a one form to be computed by reference integration
   * with coefficients that are piecewise constant  */
  ReducedIntegral(int spatialDim,
    const CellType& maxCellType,
    int dim, 
    const CellType& cellType,
    const BasisFamily& testBasis,
    int alpha,
    int testDerivOrder,
    const QuadratureFamily& quad,
    bool isInternalBdry,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh,
    int verb);

  /** Construct a two-form to be computed by reference integration
   * with coefficients that are piecewise constant */
  ReducedIntegral(int spatialDim,
    const CellType& maxCellType,
    int dim,
    const CellType& cellType,
    const BasisFamily& testBasis,
    int alpha,
    int testDerivOrder,
    const BasisFamily& unkBasis,
    int beta,
    int unkDerivOrder,
    const QuadratureFamily& quad,
    bool isInternalBdry,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh,
    int verb);

  /** virtual dtor */
  virtual ~ReducedIntegral(){;}

  /** */
  void transform(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& isLocalFlag,
    const Array<int>& facetNum,
    const RCP<Array<int> >& cellLIDs,
    const double* const coeffs,
    RCP<Array<double> >& A) const
    {
      if (order()==2) transformTwoForm(JTrans, JVol, facetNum, cellLIDs, coeffs, A);
      else if (order()==1) transformOneForm(JTrans, JVol, facetNum, cellLIDs, coeffs, A);
      else transformZeroForm(JTrans, JVol, isLocalFlag, facetNum,
        cellLIDs, coeffs, A);
    }

  /** */
  virtual void transformZeroForm(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& isLocalFlag,
    const Array<int>& facetIndex,
    const RCP<Array<int> >& cellLIDs,
    const double* const coeffs,
    RCP<Array<double> >& A) const ;
      
  /** */
  virtual void transformTwoForm(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& facetIndex,
    const RCP<Array<int> >& cellLIDs,
    const double* const coeffs,
    RCP<Array<double> >& A) const ;
      
  /** */
  void transformOneForm(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& facetIndex,
    const RCP<Array<int> >& cellLIDs,
    const double* const coeffs,
    RCP<Array<double> >& A) const ;

private:

  

  /** */
  inline double& value(int facetCase, int testDerivDir, int testNode,
    int unkDerivDir, int unkNode)
    {return W_[facetCase][unkNode + nNodesUnk()*testNode 
        + nNodes()*(unkDerivDir 
          + nRefDerivUnk()*testDerivDir)];}

  /** */
  inline const double& value(int facetCase, 
    int testDerivDir, int testNode,
    int unkDerivDir, int unkNode) const 
    {
      return W_[facetCase][unkNode + nNodesUnk()*testNode 
        + nNodes()*(unkDerivDir 
          + nRefDerivUnk()*testDerivDir)];
    }
      
  /** */
  inline double& value(int facetCase, int testDerivDir, int testNode)
    {return W_[facetCase][nNodesTest()*testDerivDir + testNode];}

  /** */
  inline const double& value(int facetCase, 
    int testDerivDir, int testNode) const 
    {return W_[facetCase][nNodesTest()*testDerivDir + testNode];}

  static double& totalFlops() {static double rtn = 0; return rtn;}

protected:

  static void addFlops(const double& flops) {totalFlops() += flops;}
      
private:

  Array<Array<double> > W_;

};
}


#endif
