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

#include "SundanceElementIntegral.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "Teuchos_TimeMonitor.hpp"

using namespace Sundance;
using namespace Teuchos;

using std::endl;

static Time& transCreationTimer() 
{
  static RCP<Time> rtn 
    = TimeMonitor::getNewTimer("building integral transformation matrices"); 
  return *rtn;
}


bool& ElementIntegral::alwaysUseCofacets()
{
  static bool rtn = true;
  return rtn;
}

ElementIntegral::ElementIntegral(int spatialDim,
  const CellType& maxCellType,
  int dim, 
  const CellType& cellType,
  bool isInternalBdry,
  const ParametrizedCurve& globalCurve,
  const Mesh& mesh,
  int verb)
  : setupVerb_(verb),
    integrationVerb_(0),
    transformVerb_(0),
    spatialDim_(spatialDim),
    dim_(dim),
    isInternalBdry_(isInternalBdry),
    nFacetCases_(1),
    testDerivOrder_(-1), 
    nRefDerivTest_(-1),
    nNodesTest_(-1),
    unkDerivOrder_(-1), 
    nRefDerivUnk_(-1),
    nNodesUnk_(-1),
    nNodes_(-1),
    order_(0),
    alpha_(),
    beta_(),
    cellType_(cellType),
    maxCellType_(maxCellType),
    evalCellType_(cellType),
    testBasis_(),
    unkBasis_(),
    globalCurve_(globalCurve),
    mesh_(mesh)
{
  Tabs tab0;
  SUNDANCE_MSG2(setupVerb(), tab0 << "constructing 0-form ElementIntegral");
  /* if we're integrating a derivative along a facet, we need to refer back
   * to the maximal cell. */
  if (alwaysUseCofacets() || (dim != spatialDim && !isInternalBdry))
  {
    evalCellType_ = maxCellType_;
    nFacetCases_ = numFacets(maxCellType, dim);
  }
}

ElementIntegral::ElementIntegral(int spatialDim,
  const CellType& maxCellType,
  int dim, 
  const CellType& cellType,
  const BasisFamily& testBasis,
  int alpha,
  int testDerivOrder,
  bool isInternalBdry,
  const ParametrizedCurve& globalCurve,
  const Mesh& mesh,
  int verb)
  : setupVerb_(verb),
    integrationVerb_(0),
    transformVerb_(0),
    spatialDim_(spatialDim),
    dim_(dim),
    isInternalBdry_(isInternalBdry),
    nFacetCases_(1),
    testDerivOrder_(testDerivOrder), 
    nRefDerivTest_(ipow(spatialDim, testDerivOrder)),
    nNodesTest_(testBasis.nReferenceDOFsWithFacets(maxCellType, cellType)),
    unkDerivOrder_(-1), 
    nRefDerivUnk_(-1),
    nNodesUnk_(-1),
    nNodes_(nNodesTest_),
    order_(1),
    alpha_(alpha),
    beta_(-1),
    cellType_(cellType),
    maxCellType_(maxCellType),
    evalCellType_(cellType),
    testBasis_(testBasis),
    unkBasis_(),
    globalCurve_(globalCurve),
    mesh_(mesh)
{
  Tabs tab0(0);
  SUNDANCE_MSG2(setupVerb(), tab0 << "constructing 1-form ElementIntegral");
  /* if we're integrating a derivative along a facet, we 
   * may need to refer back to the maximal cell. */
  bool okToRestrictTestToBdry = basisRestrictableToBoundary(testBasis);
    
  Tabs tab1;
  SUNDANCE_MSG2(setupVerb(), tab1 << "dim=" << dim << " spatialDim=" << spatialDim);
  if (dim != spatialDim)
  {
    if (isInternalBdry)
    {
      TEUCHOS_TEST_FOR_EXCEPT(!okToRestrictTestToBdry);
    }
    if (alwaysUseCofacets() || testDerivOrder>0)
    {
      Tabs tab2;
      evalCellType_ = maxCellType_;
      nFacetCases_ = numFacets(maxCellType, dim);
      nNodesTest_ = testBasis.nReferenceDOFsWithFacets(maxCellType, maxCellType);
      SUNDANCE_MSG2(setupVerb(), tab2 << "nNodesTest=" << nNodesTest_);
      nNodes_ = nNodesTest_;
      TEUCHOS_TEST_FOR_EXCEPT(nNodes_ == 0);
    }
    else
    {
      TEUCHOS_TEST_FOR_EXCEPT(!okToRestrictTestToBdry);
    }
  }

  SUNDANCE_MSG2(setupVerb(), tab1 << "nNodes=" << nNodes_);
}



ElementIntegral::ElementIntegral(int spatialDim,
  const CellType& maxCellType,
  int dim,
  const CellType& cellType,
  const BasisFamily& testBasis,
  int alpha,
  int testDerivOrder,
  const BasisFamily& unkBasis,
  int beta,
  int unkDerivOrder,
  bool isInternalBdry,
  const ParametrizedCurve& globalCurve,
  const Mesh& mesh,
  int verb)
  : setupVerb_(verb),
    integrationVerb_(0),
    transformVerb_(0),
    spatialDim_(spatialDim),
    dim_(dim),
    isInternalBdry_(isInternalBdry),
    nFacetCases_(1),
    testDerivOrder_(testDerivOrder), 
    nRefDerivTest_(ipow(spatialDim, testDerivOrder)),
    nNodesTest_(testBasis.nReferenceDOFsWithFacets(maxCellType, cellType)), 
    unkDerivOrder_(unkDerivOrder), 
    nRefDerivUnk_(ipow(spatialDim, unkDerivOrder)),
    nNodesUnk_(unkBasis.nReferenceDOFsWithFacets(maxCellType, cellType)), 
    nNodes_(nNodesTest_*nNodesUnk_),
    order_(2),
    alpha_(alpha),
    beta_(beta),
    cellType_(cellType),
    maxCellType_(maxCellType),
    evalCellType_(cellType),
    testBasis_(testBasis),
    unkBasis_(unkBasis),
    globalCurve_(globalCurve),
    mesh_(mesh)
{
  Tabs tab0(0);
  SUNDANCE_MSG2(setupVerb(), tab0 
    << "***** constructing 2-form ElementIntegral");
  /* if we're integrating a derivative along a facet, we may need to refer back
   * to the maximal cell. */
  bool okToRestrictTestToBdry = basisRestrictableToBoundary(testBasis);
  bool okToRestrictUnkToBdry = basisRestrictableToBoundary(unkBasis);

    
  Tabs tab1;
  SUNDANCE_MSG2(setupVerb(), tab1 << "dim=" << dim << " spatialDim=" << spatialDim);
  if (dim != spatialDim)
  {
    if (isInternalBdry)
    {
      TEUCHOS_TEST_FOR_EXCEPT(!(okToRestrictTestToBdry && okToRestrictUnkToBdry));   
    }
    if (alwaysUseCofacets() || testDerivOrder>0 || unkDerivOrder>0)
    {
      Tabs tab2;
      evalCellType_ = maxCellType_;
      nFacetCases_ = numFacets(maxCellType, dim);
      nNodesTest_ = testBasis.nReferenceDOFsWithFacets(maxCellType, maxCellType);
      SUNDANCE_MSG2(setupVerb(), tab2 << "nNodesTest=" << nNodesTest_);
      nNodesUnk_ = unkBasis.nReferenceDOFsWithFacets(maxCellType, maxCellType);
      SUNDANCE_MSG2(setupVerb(), tab2 << "nNodesUnk=" << nNodesUnk_);
      nNodes_ = nNodesTest_ * nNodesUnk_;
      TEUCHOS_TEST_FOR_EXCEPT(nNodes_ == 0);
    }
    else
    {
      TEUCHOS_TEST_FOR_EXCEPT(okToRestrictTestToBdry != okToRestrictUnkToBdry);
    }
  }

  SUNDANCE_MSG2(setupVerb(), tab1 << "nNodes=" << nNodes_);
}


void ElementIntegral::setVerb(
  int integrationVerb,
  int transformVerb)
{
  integrationVerb_ = integrationVerb;
  transformVerb_ = transformVerb;
}

void ElementIntegral::describe(std::ostream& os) const 
{
  Tabs tab(0);
  bool hasTest = testDerivOrder() >= 0;
  bool hasUnk = unkDerivOrder() >= 0;

  if (hasTest)
  {
    os << tab << "Test functions:" << std::endl;
    Tabs tab1;
    os << tab1 << "test basis = " << testBasis() << std::endl;
    os << tab1 << "test differentiation order = " << testDerivOrder() << std::endl;
    os << tab1 << "alpha = " << alpha() << std::endl;
    os << tab1 << "num test gradient components=" << nRefDerivTest() << std::endl;
  }
  if (hasUnk)
  {
    os << tab << "Unknown functions:" << std::endl;
    Tabs tab1;
    os << tab1 << "unk basis = " << unkBasis() << std::endl;
    os << tab1 << "unk differentiation order = " << unkDerivOrder() << std::endl;
    os << tab1 << "beta = " << beta() << std::endl;
    os << tab1  << "num unk gradient components=" << nRefDerivUnk() << std::endl;
  }
  os << tab << "Geometry:" << std::endl;
  Tabs tab1;
  os << tab1 << "cell type on which integral is defined: " << cellType()
     << std::endl;
  os << tab1 << "maximal cell type: " << maxCellType() << std::endl;
  os << tab1 << "cell type on which quad pts defined: " << evalCellType()
     << std::endl;
  os << tab << "Number of evaluation cases: " << nFacetCases() << std::endl;
}


void ElementIntegral::assertBilinearForm() const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(testDerivOrder() < 0 || testDerivOrder() > 1,
    std::logic_error,
    "Test function derivative order=" << testDerivOrder()
    << " must be 0 or 1");
  
  TEUCHOS_TEST_FOR_EXCEPTION(unkDerivOrder() < 0 || unkDerivOrder() > 1,
    std::logic_error,
    "Unknown function derivative order=" << unkDerivOrder()
    << " must be 0 or 1");
}

void ElementIntegral::assertLinearForm() const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(testDerivOrder() < 0 || testDerivOrder() > 1,
    std::logic_error,
    "Test function derivative order=" << testDerivOrder()
    << " must be 0 or 1");
}


void ElementIntegral::getQuad(const QuadratureFamily& quad,
  int evalCase, Array<Point>& quadPts, Array<double>& quadWeights) const 
{
  Tabs tab(0);

  SUNDANCE_MSG2(setupVerb(), tab << "getting quad points for rule "
    << quad);

  if (nFacetCases()==1) 
  {
    quad.getPoints(cellType(), quadPts, quadWeights);
  }
  else 
  {
    int dim = dimension(cellType());
    quad.getFacetPoints(maxCellType(), dim, evalCase, 
      quadPts, quadWeights);
  }

  if (setupVerb() >= 4)
  {
    Tabs tab1;
    SUNDANCE_MSG4(setupVerb(), 
      tab1 << "quadrature points on ref cell are:");
    printQuad(Out::os(), quadPts, quadWeights);
  }
}
  


Array<double>& ElementIntegral::G(int alpha)
{
  static Array<Array<double> > rtn(3);

  return rtn[alpha];
}

Array<double>& ElementIntegral::G(int alpha, int beta)
{
  static Array<Array<double> > rtn(9);

  return rtn[3*alpha + beta];
}

int& ElementIntegral::transformationMatrixIsValid(int alpha, int beta)
{
  static Array<int> rtn(9, false);
  return rtn[3*alpha + beta];
}

int& ElementIntegral::transformationMatrixIsValid(int alpha)
{
  static Array<int> rtn(3, false);
  return rtn[alpha];
}

void ElementIntegral::invalidateTransformationMatrices()
{
  for (int i=0; i<3; i++)
  {
    transformationMatrixIsValid(i) = false;
    for (int j=0; j<3; j++)
    {
      transformationMatrixIsValid(i, j) = false;
    }
  }
}




int ElementIntegral::ipow(int base, int power) 
{
  int rtn = 1;
  for (int i=0; i<power; i++) rtn *= base;
  return rtn;
}

void ElementIntegral
::createTwoFormTransformationMatrix(const CellJacobianBatch& JTrans,
  const CellJacobianBatch& JVol) const
{
  TimeMonitor timer(transCreationTimer());
  Tabs tab;

  int flops = 0;

  int maxDim = JTrans.cellDim();
  //int cellDim = JVol.cellDim();

  if (testDerivOrder() == 1 && unkDerivOrder() == 1)
  {
    Tabs tab2;
    if (transformationMatrixIsValid(alpha(), beta())) return;
    transformationMatrixIsValid(alpha(), beta()) = true;

    G(alpha(), beta()).resize(JTrans.numCells() * JTrans.cellDim() * JTrans.cellDim());

    double* GPtr = &(G(alpha(),beta())[0]);
    int k = 0;

    for (int c=0; c<JTrans.numCells(); c++)
    {
      static Array<double> invJ;
      JTrans.getInvJ(c, invJ);
      double detJ = fabs(JVol.detJ()[c]);
      for (int gamma=0; gamma<maxDim; gamma++)
      {
        for (int delta=0; delta<maxDim; delta++, k++)
        {
          GPtr[k] =  detJ*invJ[alpha() + gamma*maxDim]
            * invJ[beta() + maxDim*delta];
        }
      }
    }
    flops = 2 * JTrans.numCells() * maxDim * maxDim + JTrans.numCells();
  }

  else if (testDerivOrder() == 1 && unkDerivOrder() == 0)
  {
    if (transformationMatrixIsValid(alpha())) return;
    transformationMatrixIsValid(alpha()) = true;

    G(alpha()).resize(JTrans.numCells() * JTrans.cellDim());

    int k = 0;
    double* GPtr = &(G(alpha())[0]);

    for (int c=0; c<JTrans.numCells(); c++)
    {
      static Array<double> invJ;
      JTrans.getInvJ(c, invJ);
      double detJ = fabs(JVol.detJ()[c]);
      for (int gamma=0; gamma<maxDim; gamma++,k++)
      {
        GPtr[k] = detJ*invJ[alpha() + maxDim * gamma];
      }
    }
    flops = JTrans.numCells() * maxDim + JTrans.numCells();
  }

  else 
  {
    if (transformationMatrixIsValid(beta())) return;
    transformationMatrixIsValid(beta()) = true;

    G(beta()).resize(JTrans.numCells() * JTrans.cellDim());

    int k = 0;
    double* GPtr = &(G(beta())[0]);

    for (int c=0; c<JTrans.numCells(); c++)
    {
      static Array<double> invJ;
      JTrans.getInvJ(c, invJ);
      double detJ = fabs(JVol.detJ()[c]);
      for (int gamma=0; gamma<maxDim; gamma++,k++)
      {
        GPtr[k] = detJ*invJ[beta() + maxDim * gamma];
      }
    }
    flops = JTrans.numCells() * maxDim + JTrans.numCells();
  }

  addFlops(flops);
}


void ElementIntegral
::createOneFormTransformationMatrix(const CellJacobianBatch& JTrans,
  const CellJacobianBatch& JVol) const 
{
  TimeMonitor timer(transCreationTimer());
  Tabs tab;
  SUNDANCE_MSG2(transformVerb(), 
    tab << "ElementIntegral creating linear form trans matrices");

  int maxDim = JTrans.cellDim();

  if (transformationMatrixIsValid(alpha())) return;
  transformationMatrixIsValid(alpha()) = true;

  int flops = JTrans.numCells() * maxDim + JTrans.numCells();

  G(alpha()).resize(JTrans.numCells() * JTrans.cellDim());

  int k = 0;
  double* GPtr = &(G(alpha())[0]);

  for (int c=0; c<JTrans.numCells(); c++)
  {
    Array<double> invJ;
    JTrans.getInvJ(c, invJ);
    double detJ = fabs(JVol.detJ()[c]);
    for (int gamma=0; gamma<maxDim; gamma++, k++)
    {
      GPtr[k] = detJ*invJ[alpha() + maxDim * gamma]; 
    }
  }
  
  addFlops(flops);
}

