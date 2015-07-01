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

#ifndef SUNDANCE_STDFWKEVALMEDIATOR_H
#define SUNDANCE_STDFWKEVALMEDIATOR_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceAbstractEvalMediator.hpp"
#include "SundanceCellJacobianBatch.hpp"
#include "SundanceIntegralGroup.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceDiscreteFunction.hpp"


namespace Sundance
{
using namespace Teuchos;

/** 
 * StdFwkEvalMediator evaluates mesh-dependent functions in the 
 * standard framework. A number of subtypes are supported: 
 * QuadratureEvalMediator, which does evaluation on quadrature points,  
 * and NodalEvalMediator, which does evaluation at nodal points.  */

class StdFwkEvalMediator : public AbstractEvalMediator,
                           public Playa::Printable
{
public:
  /** */
  StdFwkEvalMediator(const Mesh& mesh, int cellDim);

  /** */
  virtual ~StdFwkEvalMediator(){;}

  /** */
  void setCellBatch(const RCP<const Array<int> >& cellLID);

  /** */
  void setIntegrationSpec(IntegrationCellSpecifier intCellSpec);


  /** Update the cell type */
  virtual void setCellType(const CellType& cellType,
    const CellType& maxCellType,
    bool isInternalBdry) ;

  /** Return the Jacobian to be used in computing the volume of cells
      being integrated. This will not necessarily be the same as the
      Jacobian used for transformations of vectors: when integrating
      derivatives over boundaries, the volume is the volume of the facet,
      while the transformations are computed on the maximal cofacets. */
  const CellJacobianBatch& JVol() const {return *JVol_;}

  /** Return the Jacobian to be used in derivative transformations. */
  const CellJacobianBatch& JTrans() const ;

  /** When evaluating derivatives on boundaries, we evaluate basis
      functions on the maximal cofacets of the boundary cells. This function
      returns the facet index, relative to the maximal cofacet, of each
      boundary cell in the batch.  */
  const Array<int>& facetIndices() const {return *facetIndices_;}

  /** */
  const Array<int>& maxCellLIDs() const {return *maxCellLIDs_;}

  /** */
  int cellDim() const {return cellDim_;}

  /** */
  int maxCellDim() const {return mesh_.spatialDim();}

  /** */
  const CellType& cellType() const {return cellType_;}

  /** */
  const CellType& maxCellType() const {return maxCellType_;}

  /** */
  const RCP<const Array<int> >& cellLID() const {return cellLID_;}

  /** */
  const RCP<Array<int> >& cofacetCellLID() const {return maxCellLIDs_;}

  /** */
  IntegrationCellSpecifier integrationCellSpec() const {return intCellSpec_;}

  /** */
  bool cofacetCellsAreReady() const {return cofacetCellsAreReady_;}

  /** */
  bool isInternalBdry() const {return isInternalBdry_;}

  /** */
  bool forbidCofacetIntegrations() const 
    {return forbidCofacetIntegrations_;}


protected:
  const Mesh& mesh() const {return mesh_;}

  Mesh& mesh() {return mesh_;}

  bool& cacheIsValid() const {return cacheIsValid_;}

  /** */
  void setupFacetTransformations() const ;

  /** */
  Map<const DiscreteFunctionData*, RCP<Array<Array<double> > > >& fCache() const {return fCache_;}
  /** */
  Map<const DiscreteFunctionData*, RCP<Array<Array<double> > > >& dfCache() const {return dfCache_;}
  /** */
  Map<const DiscreteFunctionData*, RCP<Array<Array<double> > > >& localValueCache() const {return localValueCache_;}

  /** */
  Map<const DiscreteFunctionData*, RCP<const MapStructure> >& mapStructCache() const
    {return mapStructCache_;}

  /** */
  Map<const DiscreteFunctionData*, bool>& fCacheIsValid() const {return fCacheIsValid_;}
  /** */
  Map<const DiscreteFunctionData*, bool>& dfCacheIsValid() const {return dfCacheIsValid_;}
  /** */
  Map<const DiscreteFunctionData*, bool>& localValueCacheIsValid() const {return localValueCacheIsValid_;}
      
private:
  Mesh mesh_;

  int cellDim_;

  CellType cellType_;

  CellType maxCellType_;

  bool isInternalBdry_;

  bool forbidCofacetIntegrations_;

  RCP<const Array<int> > cellLID_;

  mutable IntegrationCellSpecifier intCellSpec_;

  mutable RCP<CellJacobianBatch> JVol_;

  mutable RCP<CellJacobianBatch> JTrans_;

  mutable RCP<Array<int> > facetIndices_;

  mutable RCP<Array<int> > maxCellLIDs_;

  mutable bool cofacetCellsAreReady_;

  mutable bool cacheIsValid_;

  mutable bool jCacheIsValid_;



  /** */
  mutable Map<const DiscreteFunctionData*, RCP<Array<Array<double> > > > fCache_;
  /** */
  mutable Map<const DiscreteFunctionData*, RCP<Array<Array<double> > > > dfCache_; 
  /** */
  mutable Map<const DiscreteFunctionData*, RCP<Array<Array<double> > > > localValueCache_;
  /** */
  mutable Map<const DiscreteFunctionData*, RCP<const MapStructure> > mapStructCache_;

  /** */
  mutable Map<const DiscreteFunctionData*, bool> fCacheIsValid_;
  /** */
  mutable Map<const DiscreteFunctionData*, bool> dfCacheIsValid_;
  /** */
  mutable Map<const DiscreteFunctionData*, bool> localValueCacheIsValid_;

};
}


#endif
