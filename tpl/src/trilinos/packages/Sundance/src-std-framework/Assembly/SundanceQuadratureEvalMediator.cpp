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

#include "SundanceQuadratureEvalMediator.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceTempStack.hpp"
#include "SundanceCellDiameterExpr.hpp"
#include "SundanceCellVectorExpr.hpp"
#include "SundanceSpatialDerivSpecifier.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceCellJacobianBatch.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "PlayaExceptions.hpp"

#include "Teuchos_BLAS.hpp"


using namespace Teuchos;
TEUCHOS_TIMER(coordEvalTimer, "Quad mediator: coord eval")

using namespace Sundance;
using namespace Playa;
using std::endl;
using std::setw;


QuadratureEvalMediator
::QuadratureEvalMediator(const Mesh& mesh, 
  int cellDim,
  const QuadratureFamily& quad)
  : StdFwkEvalMediator(mesh, cellDim),
    numEvaluationCases_(-1),
    quad_(quad),
    numQuadPtsForCellType_(),
    quadPtsForReferenceCell_(),
    quadPtsReferredToMaxCell_(),
    physQuadPts_(),
    refFacetBasisVals_(2)
{}

Time& QuadratureEvalMediator::coordEvaluationTimer()
{
  return coordEvalTimer();
}

void QuadratureEvalMediator::setCellType(const CellType& cellType,
  const CellType& maxCellType,
  bool isInternBdry) 
{
  StdFwkEvalMediator::setCellType(cellType, maxCellType, isInternBdry);

  Tabs tab;

  SUNDANCE_MSG2(verb(), tab <<  "QuadEvalMed::setCellType: cellType=" 
    << cellType << " cellDim=" << cellDim() << " maxCellType=" << maxCellType);
  SUNDANCE_MSG2(verb(), tab << "integration spec: =" 
    << integrationCellSpec());
  SUNDANCE_MSG2(verb(), tab << "forbid cofacet integrations =" 
    << forbidCofacetIntegrations());
  if (isInternalBdry()) 
  {
    SUNDANCE_MSG2(verb(), tab << "working on internal boundary");
  }
  
//  TEUCHOS_TEST_FOR_EXCEPT(isInternalBdry()
//    && integrationCellSpec() != NoTermsNeedCofacets);

  if (cellType != maxCellType)
  {
    Tabs tab1;
    SUNDANCE_MSG2(verb(), tab1 << "working out #facet cases"); 
    numEvaluationCases_ = numFacets(maxCellType, cellDim());
  }
  else 
  {
    Tabs tab1;
    SUNDANCE_MSG2(verb(), tab1 << "no need for facet cases; work on original cell"); 
    numEvaluationCases_ = 1;
  }
  SUNDANCE_MSG2(verb(), tab << "num eval cases =" << numEvaluationCases_); 

  if (!isInternalBdry() 
    && quadPtsReferredToMaxCell_.containsKey(cellType)) return;

  if (!quadPtsForReferenceCell_.containsKey(cellType))
  {
    SUNDANCE_MSG2(verb(), tab << "creating quad points for ref cell type=" 
      << cellType);
    RCP<Array<Point> > pts = rcp(new Array<Point>());
    RCP<Array<double> > wgts = rcp(new Array<double>()); 
    
    quad_.getPoints(cellType, *pts, *wgts);
    quadPtsForReferenceCell_.put(cellType, pts);
    
    numQuadPtsForCellType_.put(cellType, pts->size());
  }

  if (!quadPtsReferredToMaxCell_.containsKey(cellType))
  {
    SUNDANCE_MSG2(verb(), tab << "creating quad points for max cell type=" 
      << maxCellType);
    RCP<Array<Array<Point> > > facetPts 
      = rcp(new Array<Array<Point> >(numEvaluationCases()));
    RCP<Array<Array<double> > > facetWgts 
      = rcp(new Array<Array<double> >(numEvaluationCases()));

    for (int fc=0; fc<numEvaluationCases(); fc++)
    {
      if (cellType != maxCellType)
      {
        quad_.getFacetPoints(maxCellType, cellDim(), fc, 
          (*facetPts)[fc], (*facetWgts)[fc]);
      }
      else
      {
        quad_.getPoints(maxCellType, (*facetPts)[fc], (*facetWgts)[fc]);
      }
    }
    quadPtsReferredToMaxCell_.put(cellType, facetPts);
  }
}

int QuadratureEvalMediator::numQuadPts(const CellType& ct) const 
{
  TEUCHOS_TEST_FOR_EXCEPTION(!numQuadPtsForCellType_.containsKey(ct),
    std::runtime_error,
    "no quadrature points have been tabulated for cell type=" << ct);
  return numQuadPtsForCellType_.get(ct);
}

void QuadratureEvalMediator::evalCellDiameterExpr(const CellDiameterExpr* expr,
  RCP<EvalVector>& vec) const 
{
  Tabs tabs;
  SUNDANCE_MSG2(verb(),tabs 
    << "QuadratureEvalMediator evaluating cell diameter expr " 
    << expr->toString());

  int nQuad = numQuadPts(cellType());
  int nCells = cellLID()->size();

  SUNDANCE_MSG3(verb(),tabs << "number of quad pts=" << nQuad);
  Array<double> diameters;
  mesh().getCellDiameters(cellDim(), *cellLID(), diameters);

  vec->resize(nQuad*nCells);
  double * const xx = vec->start();
  int k=0;
  for (int c=0; c<nCells; c++)
  {
    double h = diameters[c];
    for (int q=0; q<nQuad; q++, k++) 
    {
      xx[k] = h;
    }
  }
}

void QuadratureEvalMediator::evalCellVectorExpr(const CellVectorExpr* expr,
  RCP<EvalVector>& vec) const 
{
  Tabs tabs;
  SUNDANCE_MSG2(verb(),tabs 
    << "QuadratureEvalMediator evaluating cell vector expr " 
    << expr->toString());

  int nQuad = numQuadPts(cellType());
  int nCells = cellLID()->size();

  vec->resize(nQuad*nCells);

  SUNDANCE_MSG3(verb(),tabs << "number of quad pts=" << nQuad);
  int dir = expr->componentIndex();

  Array<Point> vectors;
  if (expr->isNormal())
  { 
    mesh().outwardNormals(*cellLID(), vectors);
  }
  else
  {
    TEUCHOS_TEST_FOR_EXCEPTION(cellDim() != 1, std::runtime_error,
      "unable to compute tangent vectors for cell dim = " << cellDim());
    mesh().tangentsToEdges(*cellLID(), vectors);
  }
    
  double * const xx = vec->start();
  int k=0;
  for (int c=0; c<nCells; c++)
  {
    double n = vectors[c][dir];
    for (int q=0; q<nQuad; q++, k++) 
    {
      xx[k] = n;
    }
  }
}

void QuadratureEvalMediator::evalCoordExpr(const CoordExpr* expr,
  RCP<EvalVector>& vec) const 
{
  Tabs tabs;
  SUNDANCE_MSG2(verb(),tabs 
    << "QuadratureEvalMediator evaluating coord expr " 
    << expr->toString());

  TimeMonitor timer(coordEvalTimer());
  
  computePhysQuadPts();
  int nQuad = physQuadPts_.length();
  int d = expr->dir();
  
  SUNDANCE_MSG3(verb(),tabs << "number of quad pts=" << nQuad);

  vec->resize(nQuad);
  double * const xx = vec->start();
  for (int q=0; q<nQuad; q++) 
  {
    xx[q] = physQuadPts_[q][d];
  }
}

RCP<Array<Array<Array<double> > > > QuadratureEvalMediator
::getFacetRefBasisVals(const BasisFamily& basis, int diffOrder) const
{
  Tabs tab;
  RCP<Array<Array<Array<double> > > > rtn ;

  CellType evalCellType = cellType();
  if (cellDim() != maxCellDim())
  {
    if (verb() >= 2)
    {
      Out::os() << tab << "alwaysUseCofacets = " 
                << ElementIntegral::alwaysUseCofacets() << std::endl;
      Out::os() << tab << "diffOrder = " << diffOrder << std::endl;
    }
    if (ElementIntegral::alwaysUseCofacets() || diffOrder>0)
    {
      if (!cofacetCellsAreReady()) setupFacetTransformations();
      evalCellType = maxCellType();
    
      TEUCHOS_TEST_FOR_EXCEPTION(!cofacetCellsAreReady(), std::runtime_error, 
        "cofacet cells not ready in getFacetRefBasisVals()");
    }
  }

  SUNDANCE_MSG2(verb(), tab << "eval cell type = " << evalCellType);
  typedef OrderedPair<BasisFamily, CellType> key;

  int nDerivResults = 1;
  if (diffOrder==1) nDerivResults = maxCellDim();

  if (!refFacetBasisVals_[diffOrder].containsKey(key(basis, cellType())))
  {
    SUNDANCE_OUT(this->verb() > 2,
      tab << "computing basis values on facet quad pts");
    rtn = rcp(new Array<Array<Array<double> > >(numEvaluationCases()));

    Array<Array<Array<Array<double> > > > tmp(nDerivResults);
    
    if (verb() >= 2)
    {
      Out::os() << tab << "numEvalCases = " << numEvaluationCases()
                << std::endl;
      Out::os() << tab << "diff order = " << diffOrder << std::endl;
      Out::os() << tab << "cell type = " << cellType() << std::endl;
      Out::os() << tab << "quad pt map = ";
      if (evalCellType!=cellType())
      { 
        Out::os() << quadPtsReferredToMaxCell_ << std::endl;
      }
      else
      {
        Out::os() << quadPtsForReferenceCell_ << std::endl;
      }
    }

    if (evalCellType!=cellType())
    { 
      TEUCHOS_TEST_FOR_EXCEPTION(quadPtsReferredToMaxCell_.size() == 0,
        std::runtime_error,
        "empty quadrature point map (max cell)");
    }
    else
    {
      TEUCHOS_TEST_FOR_EXCEPTION(quadPtsForReferenceCell_.size() == 0,
        std::runtime_error,
        "empty quadrature point map (ref cell)");
    }

    for (int fc=0; fc<numEvaluationCases(); fc++)
    {
      Tabs tab1;
      (*rtn)[fc].resize(basis.dim());
      SUNDANCE_MSG2(verb(), tab1 << "fc = " << fc);

      for (int r=0; r<nDerivResults; r++)
      {
        tmp[r].resize(basis.dim());
        MultiIndex mi;
        if (diffOrder==1)
        {
          mi[r]=1;
        }
        SpatialDerivSpecifier deriv(mi);
        /* Here we evaluate the basis functions at specified quadrature points
         * on the reference cell */
        if (evalCellType != cellType())
        {
          SUNDANCE_MSG2(verb(), tab1 << "referring to max cell");
          basis.refEval(evalCellType, 
            (*(quadPtsReferredToMaxCell_.get(cellType())))[fc], 
            deriv, tmp[r], verb());
        }
        else
        {
          SUNDANCE_MSG2(verb(), tab1 << "computing on reference cell");
          basis.refEval(evalCellType, 
            (*(quadPtsForReferenceCell_.get(cellType()))), 
            deriv, tmp[r], verb());
        }
      }
      /* the tmp array contains values indexed as [quad][node]. 
       * We need to put this into fortran order with quad index running
       * fastest */
      int dim = maxCellDim();
      int nQuad = tmp[0][0].size();
      int nNodes = tmp[0][0][0].size();
      int nTot = dim * nQuad * nNodes;
      for (int d=0; d<basis.dim(); d++)
      {
        (*rtn)[fc][d].resize(nTot);
        for (int r=0; r<nDerivResults; r++)
        {
          for (int q=0; q<nQuad; q++)
          {
            for (int n=0; n<nNodes; n++)
            {
              (*rtn)[fc][d][(n*nQuad + q)*nDerivResults + r] = tmp[r][d][q][n];
            }
          }
        }
      }
    }
    refFacetBasisVals_[diffOrder].put(key(basis, cellType()), rtn);
  }
  else
  {
    SUNDANCE_OUT(this->verb() > 2,
      tab << "reusing facet basis values on quad pts");
    rtn = refFacetBasisVals_[diffOrder].get(key(basis, cellType()));
  }

  return rtn;
}

Array<Array<double> >* QuadratureEvalMediator
::getRefBasisVals(const BasisFamily& basis, int diffOrder) const
{
  Tabs tab;
  RCP<Array<Array<Array<double> > > > fRtn 
    = getFacetRefBasisVals(basis, diffOrder);
  return &((*fRtn)[0]);
}


void QuadratureEvalMediator
::evalDiscreteFuncElement(const DiscreteFuncElement* expr,
  const Array<MultiIndex>& multiIndices,
  Array<RCP<EvalVector> >& vec) const
{
  const DiscreteFunctionData* f = DiscreteFunctionData::getData(expr);
  TEUCHOS_TEST_FOR_EXCEPTION(f==0, std::logic_error,
    "QuadratureEvalMediator::evalDiscreteFuncElement() called "
    "with expr that is not a discrete function");
  Tabs tab;

  SUNDANCE_MSG1(verb(),tab << "QuadEvalMed evaluating DF " << expr->name());

  int nQuad = numQuadPts(cellType());
  int myIndex = expr->myIndex();

  for (int i=0; i<multiIndices.size(); i++)
  {
    Tabs tab1;
    const MultiIndex& mi = multiIndices[i];
    SUNDANCE_MSG2(dfVerb(),
      tab1 << "evaluating DF " << expr->name() 
      << " for multiindex " << mi << std::endl
      << tab1 << "num cells = " << cellLID()->size() << std::endl
      << tab1 << "num quad points = " << nQuad << std::endl
      << tab1 << "my index = " << expr->myIndex() << std::endl
      << tab1 << "num funcs = " << f->discreteSpace().nFunc());

    vec[i]->resize(cellLID()->size() * nQuad);
  
    if (mi.order() == 0)
    {
      Tabs tab2;
      SUNDANCE_MSG2(dfVerb(),tab2 << "in mi.order() == 0 branch");
      if (!fCache().containsKey(f) || !fCacheIsValid()[f])
      {
        fillFunctionCache(f, mi);
      }
      else
      {
        SUNDANCE_MSG2(dfVerb(),tab2 << "reusing function cache");
      }

      const RCP<const MapStructure>& mapStruct = mapStructCache()[f];
      int chunk = mapStruct->chunkForFuncID(myIndex);
      int funcIndex = mapStruct->indexForFuncID(myIndex);
      int nFuncs = mapStruct->numFuncs(chunk);

      SUNDANCE_MSG3(dfVerb(),tab2 << "chunk number = " << chunk << std::endl
        << tab2 << "function index=" << funcIndex << " of nFuncs=" 
        << nFuncs);

      const RCP<Array<Array<double> > >& cacheVals 
        = fCache()[f];

      SUNDANCE_MSG4(dfVerb(),tab2 << "cached function values=" << (*cacheVals)[chunk]);

      const double* cachePtr = &((*cacheVals)[chunk][0]);
      double* vecPtr = vec[i]->start();
          
      int cellSize = nQuad*nFuncs;
      int offset = funcIndex*nQuad;
      SUNDANCE_MSG3(dfVerb(),tab2 << "cell size=" << cellSize << ", offset=" 
        << offset);
      int k = 0;
      for (int c=0; c<cellLID()->size(); c++)
      {
        for (int q=0; q<nQuad; q++, k++)
        {
          vecPtr[k] = cachePtr[c*cellSize + offset + q];
        }
      }
      SUNDANCE_MSG4(dfVerb(),tab2 << "result vector=");
      if (dfVerb() >= 5)
      {
        vec[i]->print(Out::os());
        computePhysQuadPts();
        k=0;
        for (int c=0; c<cellLID()->size(); c++)
        {
          Out::os() << tab2 << "-------------------------------------------"
                    << std::endl;
          Out::os() << tab2 << "c=" << c << " cell LID=" << (*cellLID())[c]
                    << std::endl;
          Tabs tab3;
          for (int q=0; q<nQuad; q++, k++)
          {
            Out::os() << tab3 << "q=" << q << " " << physQuadPts_[k]
                      << " val=" << vecPtr[k] << std::endl;
          }
        }
      }
    }
    else
    {
      Tabs tab2;
      SUNDANCE_MSG2(dfVerb(),tab2 << "in mi.order() != 0 branch");
      if (!dfCache().containsKey(f) || !dfCacheIsValid()[f])
      {
        fillFunctionCache(f, mi);
      }
      else
      {
        SUNDANCE_MSG3(dfVerb(),tab2 << "reusing function cache");
      }

      RCP<const MapStructure> mapStruct = mapStructCache()[f];
      int chunk = mapStruct->chunkForFuncID(myIndex);
      int funcIndex = mapStruct->indexForFuncID(myIndex);
      int nFuncs = mapStruct->numFuncs(chunk);


      SUNDANCE_MSG3(dfVerb(),tab2 << "chunk number = " << chunk << std::endl
        << tab2 << "function index=" << funcIndex << " of nFuncs=" 
        << nFuncs);

      const RCP<Array<Array<double> > >& cacheVals 
        = dfCache()[f];

      SUNDANCE_MSG4(dfVerb(),tab2 << "cached function values=" << (*cacheVals)[chunk]);

      int dim = maxCellDim();
      int pDir = mi.firstOrderDirection();
      const double* cachePtr = &((*cacheVals)[chunk][0]);
      double* vecPtr = vec[i]->start();

      int cellSize = nQuad*nFuncs*dim;
      int offset = funcIndex * nQuad * dim;
      int k = 0;

      SUNDANCE_MSG2(dfVerb(),tab2 << "dim=" << dim << ", pDir=" << pDir
        << ", cell size=" << cellSize << ", offset=" 
        << offset);
      for (int c=0; c<cellLID()->size(); c++)
      {
        for (int q=0; q<nQuad; q++, k++)
        {
          vecPtr[k] = cachePtr[c*cellSize + offset + q*dim + pDir];
        }
      }
      SUNDANCE_MSG4(dfVerb(),tab2 << "result vector=");
      if (dfVerb() >= 5)
      {
        vec[i]->print(Out::os());
        computePhysQuadPts();
        k=0;
        for (int c=0; c<cellLID()->size(); c++)
        {
          Out::os() << tab2 << "-------------------------------------------"
                    << std::endl;
          Out::os() << tab2 << "c=" << c << " cell LID=" << (*cellLID())[c]
                    << std::endl;
          Tabs tab3;
          for (int q=0; q<nQuad; q++, k++)
          {
            Out::os() << tab3 << "q=" << q << " " << physQuadPts_[k]
                      << " val=" << vecPtr[k] << std::endl;
          }
        }
      }
    }
  }
}

void QuadratureEvalMediator::fillFunctionCache(const DiscreteFunctionData* f,
  const MultiIndex& mi) const 
{
  Tabs tab0(0);
  
  SUNDANCE_MSG2(dfVerb(), tab0 << "QuadratureEvalMediator::fillFunctionCache()");
  SUNDANCE_MSG2(dfVerb(), tab0 << "multiIndex=" << mi);
  
  
  int diffOrder = mi.order();
  CellType evalCellType = cellType();

  int funcSupportDim = f->map()->cellDim();

  TEUCHOS_TEST_FOR_EXCEPT(diffOrder > 0 && funcSupportDim < maxCellDim());

  int flops = 0;
  double jFlops = CellJacobianBatch::totalFlops();

  RCP<Array<Array<double> > > localValues;
  RCP<const MapStructure> mapStruct;

  Teuchos::BLAS<int,double> blas;

  {
    Tabs tab1;
    if (cellDim() != maxCellDim())
      {
        if (dfVerb() >= 2)
        {
          Out::os() << tab1 << "alwaysUseCofacets = " 
                    << ElementIntegral::alwaysUseCofacets() << std::endl;
          Out::os() << tab1 << "diffOrder = " << diffOrder << std::endl;
        }
        if (diffOrder==0 && funcSupportDim < maxCellDim())
        {
          evalCellType = cellType();
        }
        else
        {
          evalCellType = maxCellType();
        }
      }
  
    SUNDANCE_MSG2(dfVerb(), tab1 << "cell type=" << cellType());
    SUNDANCE_MSG2(dfVerb(), tab1 << "max cell type=" << maxCellType());
    SUNDANCE_MSG2(dfVerb(), tab1 << "eval cell type=" << evalCellType);

    SUNDANCE_MSG2(dfVerb(), tab1 << "packing local values");

    if (!localValueCacheIsValid().containsKey(f) 
      || !localValueCacheIsValid().get(f))
    {
      Tabs tab2;
      SUNDANCE_MSG2(dfVerb(), tab2 << "filling cache");
      localValues = rcp(new Array<Array<double> >());
      if (cellDim() != maxCellDim())
      {
        if (diffOrder==0 && funcSupportDim < maxCellDim())
        {
           mapStruct = f->getLocalValues(cellDim(), *cellLID(), 
            *localValues);
        }
        else
        {
          if (!cofacetCellsAreReady()) setupFacetTransformations();
          mapStruct = f->getLocalValues(maxCellDim(), *cofacetCellLID(), 
            *localValues);
        }
      }
      else
      {
        mapStruct = f->getLocalValues(maxCellDim(), *cellLID(), *localValues);
      }

      TEUCHOS_TEST_FOR_EXCEPT(mapStruct.get() == 0);
      localValueCache().put(f, localValues);
      mapStructCache().put(f, mapStruct);
      localValueCacheIsValid().put(f, true);
    }
    else
    {
      Tabs tab2;
      SUNDANCE_MSG2(dfVerb(), tab2 << "reusing local value cache");
      localValues = localValueCache().get(f);
      mapStruct = mapStructCache().get(f);
    }
  }

  RCP<Array<Array<double> > > cacheVals;

  if (mi.order()==0)
  {
    if (fCache().containsKey(f))
    {
      cacheVals = fCache().get(f);
    }
    else
    {
      cacheVals = rcp(new Array<Array<double> >(mapStruct->numBasisChunks()));
      fCache().put(f, cacheVals);
    }
    fCacheIsValid().put(f, true);
  }
  else
  {
    if (dfCache().containsKey(f))
    {
      cacheVals = dfCache().get(f);
    }
    else
    {
      cacheVals = rcp(new Array<Array<double> >(mapStruct->numBasisChunks()));
      dfCache().put(f, cacheVals);
    }
    dfCacheIsValid().put(f, true);
  }


  
  for (int chunk=0; chunk<mapStruct->numBasisChunks(); chunk++)
  {
    Tabs tab1;
    SUNDANCE_MSG2(dfVerb(), tab1 << "processing dof map chunk=" << chunk
      << " of " << mapStruct->numBasisChunks());
    Tabs tab2;
    BasisFamily basis = rcp_dynamic_cast<BasisFamilyBase>(mapStruct->basis(chunk));
    SUNDANCE_MSG4(dfVerb(), tab2 << "basis=" << basis);

    int nFuncs = mapStruct->numFuncs(chunk);
    SUNDANCE_MSG2(dfVerb(), tab2 << "num funcs in this chunk=" << nFuncs);
    

    Array<double>& cache = (*cacheVals)[chunk];

    int nQuad = numQuadPts(cellType());
    int nCells = cellLID()->size();
    SUNDANCE_MSG2(dfVerb(), tab2 << "num quad points=" << nQuad);
    SUNDANCE_MSG2(dfVerb(), tab2 << "num cells=" << nCells);


    int nDir;

    if (mi.order()==1)
    {
      nDir = maxCellDim();
    }
    else
    {
      nDir = 1;
    }
    cache.resize(cellLID()->size() * nQuad * nDir * nFuncs);

      
    /* 
     * Sum over nodal values, which we can do with a matrix-matrix multiply
     * between the ref basis values and the local function values.
     *
     * There are two cases: (1) When we are evaluating 
     * on a facet, we must use different sets of reference function values
     * on the different facets. We must therefore loop over the evaluation
     * cells, using a vector of reference values chosen according to the
     * facet number of the current cell.  Let A be the
     * (nQuad*nDir)-by-(nNode) matrix of reference basis values for the
     * current cell's facet index and B be the (nNode)-by-(nFuncs) matrix of
     * function coefficient values for the current cell. Then C = A * B is
     * the (nQuad*nDir)-by-(nFunc) matrix of the function's derivative
     * values at the quadrature points in the current cell.  Each
     * matrix-matrix multiplication is done with a call to dgemm.
     *
     * (2) In other cases, we're either evaluating spatial derivatives on a
     * maximal cell or evaluating 0-order derivatives on a submaximal
     * cell. In these cases, all cells in the workset have the same
     * reference values. This lets us reuse the same matrix A on all matrix
     * multiplications, so that we can assemble one big
     * (nNode)-by-(nFuncs*nCells) matrix B and do all cells with a single
     * dgemm call to multiply A*B. The result C is then a single
     * (nQuad*nDir)-by-(nFuncs*nCells) matrix.

    */
    if (cellType() != evalCellType)
    {
      Tabs tab2;
      SUNDANCE_MSG2(dfVerb(), 
        tab2 << "evaluating by reference to maximal cell");
      
      RCP<Array<Array<Array<double> > > > refFacetBasisValues 
        = getFacetRefBasisVals(basis, mi.order());
      /* Note: even though we're ultimately not evaluating on 
       * maxCellType() here, use maxCellType() for both arguments
       * to nReferenceDOFs() because derivatives need to be
       * evaluated using all DOFs from the maximal cell, not just
       * those on the facet.
       */
      int nNodes = basis.nReferenceDOFsWithFacets(maxCellType(), maxCellType());
      int nRowsA = nQuad*nDir;
      int nColsA = nNodes;
      int nColsB = nFuncs; 
      int lda = nRowsA;
      int ldb = nNodes;
      int ldc = lda;
      double alpha = 1.0;
      double beta = 0.0;

      SUNDANCE_MSG2(dfVerb(), tab2 << "building transformation matrices cell-by-cell");
      int vecComp = 0;
      for (int c=0; c<nCells; c++)
      {
        int facetIndex = facetIndices()[c];
        double* A = &((*refFacetBasisValues)[facetIndex][vecComp][0]);
        double* B = &((*localValues)[chunk][c*nNodes*nFuncs]);
        double* C = &((*cacheVals)[chunk][c*nRowsA*nColsB]);
        blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, nRowsA, nColsB, nColsA,
          alpha, A, lda, B, ldb, beta, C, ldc);
      }
    }
    else /* cellType() == evalCellType */
    {
      /* 
       * Sum over nodal values, which we can do with a matrix-matrix multiply
       * between the ref basis values and the local function values.
       */
      Tabs tab2;
      SUNDANCE_MSG2(dfVerb(), tab2 << "building batch transformation matrix");

      Array<Array<double> >* refBasisValues 
        = getRefBasisVals(basis, diffOrder);
      int nNodes = basis.nReferenceDOFsWithFacets(maxCellType(), cellType());
      int nRowsA = nQuad*nDir;
      int nColsA = nNodes;
      int nColsB = nFuncs*nCells; 
      int lda = nRowsA;
      int ldb = nNodes;
      int ldc = lda;
      double alpha = 1.0;
      double beta = 0.0;
      int vecComp = 0;
      if (dfVerb() >= 5)
      {
        Tabs tab3;
        Out::os() << tab2 << "Printing values at nodes" << std::endl;
        for (int c=0; c<nCells; c++)
        {
          Out::os() << tab3 << "-------------------------------------------"
                    << std::endl;
          Out::os() << tab3 << "c=" << c << " cell LID=" << (*cellLID())[c]
                    << std::endl;
          Tabs tab4;
          int offset = c*nNodes*nFuncs;
          for (int n=0; n<nNodes; n++)
          {
            Out::os() << tab4 << "n=" << n;
            for (int f=0; f<nFuncs; f++)
            {
              Out::os() << " " << (*localValues)[chunk][offset + n*nFuncs + f];
            }
            Out::os() << std::endl;
          }
        }
        Out::os() << tab2 << "-------------------------------------------";
        Out::os() << tab2 << "Printing reference basis at nodes" << std::endl;
        Out::os() << tab2 << "-------------------------------------------"
                  << std::endl;
        for (int n=0; n<nNodes; n++)
        {
          Out::os() << tab3 << "node=" << n << std::endl;
          for (int q=0; q<nQuad; q++)
          {
            Tabs tab4;
            Out::os() << tab4 << "q=" << q;
            for (int r=0; r<nDir; r++)
            {
              Out::os () << " " 
                         << (*refBasisValues)[vecComp][(n*nQuad + q)*nDir + r];
            }
            Out::os() << std::endl;
          }
        }
      }
      double* A = &((*refBasisValues)[vecComp][0]);
      double* B = &((*localValues)[chunk][0]);
      double* C = &((*cacheVals)[chunk][0]);
      blas.GEMM( Teuchos::NO_TRANS, Teuchos::NO_TRANS, nRowsA, nColsB, nColsA, alpha,
        A, lda, B, ldb, beta, C, ldc );
    }


    /* Transform derivatives to physical coordinates */
    const CellJacobianBatch& J = JTrans();
    double* C = &((*cacheVals)[chunk][0]);
    if (mi.order()==1)
    {
      SUNDANCE_MSG2(dfVerb(), tab1 << "doing transformations via dgemm");
      Tabs tab2;
      SUNDANCE_MSG2(dfVerb(), tab2 << "Jacobian batch nCells=" << J.numCells());
      SUNDANCE_MSG2(dfVerb(), tab2 << "Jacobian batch cell dim=" << J.cellDim());
      SUNDANCE_MSG2(dfVerb(), tab2 << "Jacobian batch spatial dim=" << J.spatialDim());
    
      int nRhs = nQuad * nFuncs;
      for (int c=0; c<cellLID()->size(); c++)
      {
        double* rhsPtr = &(C[(nRhs * nDir)*c]);
        J.applyInvJ(c, 0, rhsPtr, nRhs, false);
      }
    }
    else
    {
      SUNDANCE_MSG2(dfVerb(), tab1 << "no derivatives to transform");
    }
    SUNDANCE_MSG2(dfVerb(), tab1 << "done transformations");
  }

  jFlops = CellJacobianBatch::totalFlops() - jFlops;
  addFlops(flops + jFlops);

  SUNDANCE_MSG2(dfVerb(), 
    tab0 << "done QuadratureEvalMediator::fillFunctionCache()");
}

void QuadratureEvalMediator::computePhysQuadPts() const 
{
  if (cacheIsValid()) 
  {
    Tabs tab0(0);
    SUNDANCE_MSG2(verb(), tab0 << "reusing phys quad points");
  }
  else
  {
    Tabs tab0(0);
    double jFlops = CellJacobianBatch::totalFlops();
    SUNDANCE_MSG2(verb(), tab0 << "computing phys quad points");
    physQuadPts_.resize(0);
    if (cellDim() != maxCellDim() && ElementIntegral::alwaysUseCofacets()
      && !isInternalBdry())
    {
      Tabs tab1;
      SUNDANCE_MSG2(verb(), tab1 << "using cofacets");
      SUNDANCE_MSG3(verb(), tab1 << "num cofacets = " << cofacetCellLID()->size());
      Array<Point> tmp;
      Array<int> cell(1);
      for (int c=0; c<cellLID()->size(); c++)
      {
        cell[0] = (*cofacetCellLID())[c];
        int facetIndex = facetIndices()[c];
        const Array<Point>& refFacetPts 
          =  (*(quadPtsReferredToMaxCell_.get(cellType())))[facetIndex];
        tmp.resize(refFacetPts.size());
        mesh().pushForward(maxCellDim(), cell, refFacetPts, tmp);
        for (int q=0; q<tmp.size(); q++)
        {
          physQuadPts_.append(tmp[q]);
        }
      }
    }
    else
    {
      const Array<Point>& refPts 
        = *(quadPtsForReferenceCell_.get(cellType()));
      mesh().pushForward(cellDim(), *cellLID(), 
        refPts, physQuadPts_); 
    }
    addFlops(CellJacobianBatch::totalFlops() - jFlops);
    cacheIsValid() = true;
    SUNDANCE_OUT(this->verb() > 2, 
      "phys quad: " << physQuadPts_);
  }
}


void QuadratureEvalMediator::print(std::ostream& os) const 
{
  if (cacheIsValid())
  {
    Tabs tab0;
    os << tab0 << "Physical quadrature points" << std::endl;
    int i=0;
    for (int c=0; c<cellLID()->size(); c++)
    {
      Tabs tab1;
      os << tab1 << "cell " << c << std::endl;
      for (int q=0; q<physQuadPts_.size()/cellLID()->size(); q++, i++)
      {
        Tabs tab2;
        os << tab2 << "q=" << q << " " << physQuadPts_[i] << std::endl;
      }
    }
  }
}


void QuadratureEvalMediator
::showResults(std::ostream& os,
	      const RCP<SparsitySuperset>& sp,
	      const Array<RCP<EvalVector> >& vecResults,
	      const Array<double>& constantResults) const
{
  Tabs tabs(0);

  

  /* find the maximum size of the std::string reps of the derivatives.
   * We'll use this to set the field width for printing derivatives. */
  int maxlen = 25;
  for (int i=0; i<sp->numDerivs(); i++)
    {
      int s = sp->deriv(i).toString().length();
      if (s > maxlen) maxlen = s;
    }

  
  int vecIndex=0;
  int constIndex = 0;
  os << tabs << "Results Superset" << std::endl;
  for (int i=0; i<sp->numDerivs(); i++)
    {
      Tabs tab1;
      os << tab1 << i << " " ;
      os.width(maxlen);
      os.setf(std::ios_base::left, std::ios_base::adjustfield);
      os << sp->deriv(i).toString() ;
      Tabs tab2;
      switch(sp->state(i))
        {
        case ZeroDeriv:
          os  << "Zero" << std::endl;
          break;
        case ConstantDeriv:
          os << "const val=" << constantResults[constIndex++] << std::endl;
          break;
        case VectorDeriv:
          if (vecResults[vecIndex].get()==0)
            {
              os << "{Null}";
            }
          else
            {
	      const double* data = vecResults[vecIndex]->start();
	      if (EvalVector::shadowOps())
		{
		  TEUCHOS_TEST_FOR_EXCEPT(vecResults[vecIndex]->str().size()==0);
		  os << vecResults[vecIndex]->str();
		}
	      os << endl;
	      int nQuad = numQuadPts(cellType());
	      int k=0;
	      TEUCHOS_TEST_FOR_EXCEPT(cellLID()->size() * nQuad 
			      != vecResults[vecIndex]->length());
	      for (int c=0; c<cellLID()->size(); c++)
		{
		  Tabs tab3;
		  os << tab3 << "cell LID=" << (*cellLID())[c] << endl;
		  for (int q=0; q<nQuad; q++, k++)
		    {
		      Tabs tab4;
		      os << tab4 << "q=" << setw(5) << q 
			 << setw(10) << Utils::chop(data[k]) << endl;
		    }
		}
            }
          vecIndex++;
          os << std::endl;
          break;
        }
    }
}
