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

#include "SundanceDiscreteSpace.hpp"
#include "SundanceDOFMapBuilder.hpp"
#include "SundanceFunctionSupportResolver.hpp"
#include "SundanceOut.hpp"
#include "SundanceMaximalCellFilter.hpp"
#include "PlayaDefaultBlockVectorSpaceDecl.hpp"
#include "PlayaMPIContainerComm.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaDefaultBlockVectorImpl.hpp"
#endif


using namespace Sundance;
using namespace Teuchos;
using Playa::blockSpace;
using Playa::MPIComm;
using Playa::MPIContainerComm;
using Playa::MPIDataType;
using Playa::MPIOp;

const int* vecPtr(const Array<int>& x)
{
  static const int* dum = 0;
  if (x.size()==0) return dum;
  else return &(x[0]);
}

DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisFamily& basis,
  const VectorType<double>& vecType,
  int setupVerb)
  : setupVerb_(setupVerb), 
    map_(),
    mesh_(mesh), 
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
  ,transformationBuilder_(0)
{
  init(maximalRegions(1), BasisArray(tuple(basis)));
}

DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
  const VectorType<double>& vecType,
  int setupVerb)
  : setupVerb_(setupVerb),
    map_(), 
    mesh_(mesh), 
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
  ,transformationBuilder_(0)
{
  init(maximalRegions(basis.size()), basis);
}

DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
  const Array<CellFilter>& funcDomains,
  const VectorType<double>& vecType,
  int setupVerb)
  : setupVerb_(setupVerb),
    map_(), 
    mesh_(mesh), 
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
  ,transformationBuilder_(0)
{
  init(funcDomains, basis);
}


DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisFamily& basis,
  const CellFilter& funcDomains,
  const VectorType<double>& vecType,
  int setupVerb)
  : setupVerb_(setupVerb),
    map_(), 
    mesh_(mesh), 
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
  ,transformationBuilder_(0)
{
  init(tuple(funcDomains), BasisArray(tuple(basis)));
}


DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
  const CellFilter& funcDomains,
  const VectorType<double>& vecType,
  int setupVerb)
  : setupVerb_(setupVerb),
    map_(), 
    mesh_(mesh), 
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
  ,transformationBuilder_(0)
{
  init(Array<CellFilter>(basis.size(), funcDomains), basis);
}



DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
  const RCP<DOFMapBase>& map,
  const RCP<Array<int> >& bcIndices,
  const VectorType<double>& vecType,
  int setupVerb)
  : setupVerb_(setupVerb),
    map_(map), 
    mesh_(mesh),
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
  ,transformationBuilder_(0)
{
  init(map->funcDomains(), basis, bcIndices, true);
}

DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
  const RCP<DOFMapBase>& map,
  const VectorType<double>& vecType,
  int setupVerb)
  : map_(map), 
    mesh_(mesh),
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
  ,transformationBuilder_(0)
{
  init(map->funcDomains(), basis);
}


DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisFamily& basis,
  const SpectralBasis& spBasis,
  const VectorType<double>& vecType,
  int setupVerb)
  : setupVerb_(setupVerb),
    map_(),
    mesh_(mesh), 
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
  ,transformationBuilder_(0)
{
  init(maximalRegions(spBasis.nterms()), 
    replicate(basis, spBasis.nterms()));
}

DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
  const SpectralBasis& spBasis,
  const VectorType<double>& vecType,
  int setupVerb)
  : setupVerb_(setupVerb),
    map_(), 
    mesh_(mesh), 
    subdomains_(),
    basis_(),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
  ,transformationBuilder_(0)
{
  init(maximalRegions(basis.size() * spBasis.nterms()), 
    replicate(basis, spBasis.nterms()));
}

DiscreteSpace::DiscreteSpace(const Mesh& mesh, const BasisArray& basis,
  const RCP<FunctionSupportResolver>& fsr,
  const VectorType<double>& vecType,
  int setupVerb)
  : setupVerb_(setupVerb),
    map_(), 
    mesh_(mesh), 
    subdomains_(),
    basis_(basis),
    vecSpace_(), 
    vecType_(vecType),
    ghostImporter_()
  ,transformationBuilder_(new DiscreteSpaceTransfBuilder())
{
  bool partitionBCs = false;
  DOFMapBuilder builder(mesh, fsr, partitionBCs, setupVerb);

  map_ = builder.colMap()[0];
  Array<Set<CellFilter> > cf = builder.unkCellFilters()[0];

  for (int i=0; i<cf.size(); i++)
  {
    Array<Array<CellFilter> > dimCF(mesh.spatialDim()+1);
    for (Set<CellFilter>::const_iterator 
           iter=cf[i].begin(); iter != cf[i].end(); iter++)
    {
      CellFilter f = *iter;
      int dim = f.dimension(mesh);
      dimCF[dim].append(f);
    }
    for (int d=mesh.spatialDim(); d>=0; d--)
    {
      if (dimCF[d].size() == 0) continue;
      CellFilter f = dimCF[d][0];
      for (int j=1; j<dimCF[d].size(); j++)
      {
        f = f + dimCF[d][j];
      }
      subdomains_.append(f);
      break;
    }
  }
  RCP<Array<int> > dummyBCIndices;
  
  // set up the transformation
  transformationBuilder_ = rcp(new DiscreteSpaceTransfBuilder( mesh , basis , map_ ));

  initVectorSpace(dummyBCIndices, partitionBCs);
  initImporter();
}

void DiscreteSpace::init(
  const Array<CellFilter>& regions,
  const BasisArray& basis)
{
  RCP<Array<int> > dummyBCIndices;
  init(regions, basis, dummyBCIndices, false);
}

void DiscreteSpace::init(
  const Array<CellFilter>& regions,
  const BasisArray& basis,
  const RCP<Array<int> >& isBCIndex, 
  bool partitionBCs)
{
  basis_ = basis;
  subdomains_ = regions;
  Array<RCP<BasisDOFTopologyBase> > basisTop(basis.size());
  for (int b=0; b<basis.size(); b++)
  {
    basisTop[b] = rcp_dynamic_cast<BasisDOFTopologyBase>(basis[b].ptr());
  }

  if (map_.get()==0) 
  {
    Array<Set<CellFilter> > cf(regions.size());
    for (int i=0; i<regions.size(); i++) cf[i] = makeSet(regions[i]);
    DOFMapBuilder b(setupVerb_);
    map_ = b.makeMap(mesh_, basisTop, cf);
  }

  // set up the transformation
  transformationBuilder_ = rcp(new DiscreteSpaceTransfBuilder( mesh_ , basis , map_ ));

  initVectorSpace(isBCIndex, partitionBCs);

  if (!partitionBCs)
  {
    initImporter();
  }
}

void DiscreteSpace::initVectorSpace(
  const RCP<Array<int> >& isBCIndex, 
  bool partitionBCs)
{
  TEUCHOS_TEST_FOR_EXCEPTION(map_.get()==0, std::logic_error,
    "uninitialized map");

  int nDof = map_->numLocalDOFs();
  int lowDof = map_->lowestLocalDOF();

  if (partitionBCs)
  {
    TEUCHOS_TEST_FOR_EXCEPT(isBCIndex.get() == 0);

    int nBCDofs = 0;
    for (int i=0; i<nDof; i++)
    {
      if ((*isBCIndex)[i]) nBCDofs++;
    }
    
    int nTotalBCDofs = nBCDofs;
    mesh().comm().allReduce(&nBCDofs, &nTotalBCDofs, 1, MPIDataType::intType(), MPIOp::sumOp());
    int nTotalInteriorDofs = map_->numDOFs() - nTotalBCDofs;

    Array<int> interiorDofs(nDof - nBCDofs);
    Array<int> bcDofs(nBCDofs);
    int iBC = 0;
    int iIn = 0;
    for (int i=0; i<nDof; i++)
    {
      if ((*isBCIndex)[i]) bcDofs[iBC++] = lowDof+i;
      else interiorDofs[iIn++] = lowDof+i;
    }
    const int* bcDofPtr = vecPtr(bcDofs);
    const int* intDofPtr = vecPtr(interiorDofs);
    VectorSpace<double> bcSpace = vecType_.createSpace(nTotalBCDofs, nBCDofs,
      bcDofPtr, mesh().comm());
    VectorSpace<double> interiorSpace = vecType_.createSpace(nTotalInteriorDofs, nDof-nBCDofs,
      intDofPtr, mesh().comm());

    vecSpace_ = blockSpace<double>(interiorSpace, bcSpace);
  }
  else
  {
    Array<int> dofs(nDof);
    for (int i=0; i<nDof; i++) dofs[i] = lowDof + i;
    
    vecSpace_ = vecType_.createSpace(map_->numDOFs(),
      map_->numLocalDOFs(),
      &(dofs[0]), mesh().comm());
  }
}

void DiscreteSpace::initImporter()
{
  TEUCHOS_TEST_FOR_EXCEPTION(map_.get()==0, std::logic_error,
    "uninitialized map");
  TEUCHOS_TEST_FOR_EXCEPTION(vecSpace_.ptr().get()==0, std::logic_error,
    "uninitialized vector space");
  TEUCHOS_TEST_FOR_EXCEPTION(vecType_.ptr().get()==0, std::logic_error,
    "uninitialized vector type");
  
  RCP<Array<int> > ghostIndices = map_->ghostIndices();
  int nGhost = ghostIndices->size();
  int* ghosts = 0;
  if (nGhost!=0) ghosts = &((*ghostIndices)[0]);
  ghostImporter_ = vecType_.createGhostImporter(vecSpace_, nGhost, ghosts);
}


Array<CellFilter> DiscreteSpace::maximalRegions(int n) const
{
  CellFilter cf = new MaximalCellFilter();
  Array<CellFilter> rtn(n, cf);
  return rtn;
}



void DiscreteSpace::importGhosts(const Vector<double>& x,
  RCP<GhostView<double> >& ghostView) const
{
  ghostImporter_->importView(x, ghostView);
}
