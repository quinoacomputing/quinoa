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

#ifndef SUNDANCE_PERIODIC_MESH1D_H
#define SUNDANCE_PERIODIC_MESH1D_H

#include "SundanceDefs.hpp"
#include "SundanceMeshBase.hpp"
#include "SundanceSet.hpp"
#include "SundancePoint.hpp"
#include "SundanceCellType.hpp"

#include "Teuchos_Array.hpp"
#include "Teuchos_Hashtable.hpp"

namespace Sundance
{
using namespace Teuchos;

/** 
 * Periodic mesh in 1D
 */

class PeriodicMesh1D : public MeshBase
{
public:
  /** */
	PeriodicMesh1D(double xMin, double xMax, int numElems);

	/** virtual dtor */
	virtual ~PeriodicMesh1D(){}

  /**
   * Get the number of cells having dimension dim
   */
  int numCells(int dim) const;

  /**
   * Return the position of the i-th node
   */
  Point nodePosition(int i) const;

  /**
   * Return a view of the i-th node's position
   */
  const double* nodePositionView(int i) const;



  /**
   * Compute the jacobians of a batch of cells, returning the
   * result via reference argument
   *
   * @param cellDim dimension of the cells whose Jacobians are to
   * be computed
   * @param cellLID local indices of the cells for which Jacobians
   * are to be computed
   * @param jBatch reference to the resulting Jacobian batch
   */
  virtual void getJacobians(int cellDim, const Array<int>& cellLID,
    CellJacobianBatch& jBatch) const;

  /**
   * Compute the diameters of a batch of cells,
   * result via reference argument
   *
   * @param cellDim dimension of the cells whose diameters are to
   * be computed
   * @param cellLID local indices of the cells for which diameters
   * are to be computed
   * @param diameters reference to the array of cell diameters
   */
  virtual void getCellDiameters(int cellDim, const Array<int>& cellLID,
    Array<double>& cellDiameters) const;


  /**
   * Map reference quadrature points to physical points on the
   * given cells.
   */
  virtual void pushForward(int cellDim, const Array<int>& cellLID,
    const Array<Point>& refQuadPts,
    Array<Point>& physQuadPts) const;

  /**
   * Return the rank of the processor that owns the given cell
   */
  virtual int ownerProcID(int cellDim, int cellLID) const {return 0;}

  /**
   * Return the number of facets of the given cell
   */
  virtual int numFacets(int cellDim, int cellLID,
    int facetDim) const;

  /**
   * Return the local ID of a facet cell
   * @param cellDim dimension of the cell whose facets are being obtained
   * @param cellLID local index of the cell whose
   * facets are being obtained
   * @param facetDim dimension of the desired facet
   * @param facetIndex index into the list of the cell's facets
   * @param facetOrientation orientation of the facet as seen from
   * the given cell (returned via reference)
   */
  virtual int facetLID(int cellDim, int cellLID,
    int facetDim, int facetIndex,
    int& facetOrientation) const;

  /**
   * Return by reference argument an array containing&(elemVerts_.value(cellLID, 0))
   * the LIDs of the facetDim-dimensional facets of the
   * given batch of cells
   */
  virtual void getFacetLIDs(int cellDim,
    const Array<int>& cellLID,
    int facetDim,
    Array<int>& facetLID,
    Array<int>& facetSign) const;

  /**
   * Return a view of an element's zero-dimensional facets,
   * @return an array of integers with the indexes of the points which for it
   */
  const int* elemZeroFacetView(int cellLID) const;

  /**
   * Return the number of maximal cofacets of the given cell
   */
  virtual int numMaxCofacets(int cellDim, int cellLID) const;

  /**
   * Return the local ID of a maximal cofacet cell
   * @param cellDim dimension of the cell whose cofacets are being obtained
   * @param cellLID local index of the cell whose
   * cofacets are being obtained
   * @param cofacetIndex which maximal cofacet to get
   * @param facetIndex index of the cell cellLID into the list of the
   * maximal cell's facets
   */
  virtual int maxCofacetLID(int cellDim, int cellLID,
    int cofacetIndex,
    int& facetIndex) const;
  /**
   * Get the LIDs of the maximal cofacets for a batch of codimension-one
   * cells.
   *
   * \param cellLIDs [in] array of LIDs of the cells whose cofacets are
   * being obtained
   * \param cofacets [out] the batch of cofacets
   */
  virtual void getMaxCofacetLIDs(const Array<int>& cellLIDs,
    MaximalCofacetBatch& cofacets) const;


  /**
   * Find the cofacets of the given cell
   * @param cellDim dimension of the cell whose cofacets are being obtained
   * @param cellLID local index of the cell whose
   * cofacets are being obtained
   * @param cofacetDim dimension of the cofacets to get
   * @param cofacetLIDs LIDs of the cofacet
   */
  void getCofacets(int cellDim, int cellLID,
    int cofacetDim, Array<int>& cofacetLIDs) const;

  /**
   * Find the local ID of a cell given its global index
   */
  virtual int mapGIDToLID(int cellDim, int globalIndex) const ;

  /**
   * Determine whether a given cell GID exists on this processor
   */
  virtual bool hasGID(int cellDim, int globalIndex) const;


  /**
   * Find the global ID of a cell given its local index
   */
  virtual int mapLIDToGID(int cellDim, int localIndex) const;

  /**
   * Get the type of the given cell
   */
  virtual CellType cellType(int cellDim) const;


  /** Get the label of the given cell */
  virtual int label(int cellDim, int cellLID) const;

  /** Get the labels for a batch of cells */
  virtual void getLabels(int cellDim, const Array<int>& cellLID,
    Array<int>& labels) const;



  /** Get the list of all labels defined for cells of the given dimension */
  virtual Set<int> getAllLabelsForDimension(int cellDim) const;

  /**
   * Get the cells associated with a specified label. The array
   * cellLID will be filled with those cells of dimension cellDim
   * having the given label.
   */
  virtual void getLIDsForLabel(int cellDim, int label, Array<int>& cellLIDs) const;


  /** Set the label of the given cell */
  virtual void setLabel(int cellDim, int cellLID, int label);

  /** Coordinate intermediate cell definitions across processors  */
  virtual void assignIntermediateCellGIDs(int cellDim) {}

  /** */
  virtual bool hasIntermediateGIDs(int dim) const {return true;}



private:

  int numElems_;
  double xMin_;
  double xMax_;
  Array<Point> x_;
  Array<Array<int> > verts_;
  Array<Array<int> > labels_;
};
}


#endif 
