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

#ifndef SUNDANCE_QUADRATUREINTEGRALBASE_H
#define SUNDANCE_QUADRATUREINTEGRALBASE_H

#include "SundanceDefs.hpp"
#include "SundanceElementIntegral.hpp"

namespace Sundance
{
using namespace Teuchos;
    
/** 
 *  
 *
 */
class QuadratureIntegralBase
  : public ElementIntegral
{
public:

  /** Construct a zero form to be computed by quadrature */
  QuadratureIntegralBase(int spatialDim,
    const CellType& maxCellType,
    int dim, 
    const CellType& cellType,
    const QuadratureFamily& quad,
    bool isInternalBdry,
    const ParametrizedCurve& globalCurve,
    const Mesh& mesh,
    int verb);

  /** Construct a one form to be computed by quadrature */
  QuadratureIntegralBase(int spatialDim,
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
      
  /** Construct a two-form to be computed by quadrature */
  QuadratureIntegralBase(int spatialDim,
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
  virtual ~QuadratureIntegralBase(){;}
      
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
    RCP<Array<double> >& A) const = 0;
      
  /** */
  virtual void transformTwoForm(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& facetIndex,
    const RCP<Array<int> >& cellLIDs,
    const double* const coeff,
    RCP<Array<double> >& A) const = 0;
      
  /** */
  virtual void transformOneForm(const CellJacobianBatch& JTrans,
    const CellJacobianBatch& JVol,
    const Array<int>& facetIndex,
    const RCP<Array<int> >& cellLIDs,
    const double* const coeff,
    RCP<Array<double> >& A) const = 0;
      
      
  /** */
  virtual int nQuad() const {return nQuad_;}
      
  static double& totalFlops() {static double rtn = 0; return rtn;}

protected:
  static void addFlops(const double& flops) {totalFlops() += flops;}
      
  const QuadratureFamily& quad() const {return quad_;}
  /** */
  int nQuad_ ;

private:

  /* */
  QuadratureFamily quad_;
      
};

}


#endif
