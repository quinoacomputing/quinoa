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

#ifndef SUNDANCE_MESH_H
#define SUNDANCE_MESH_H

#include "SundanceDefs.hpp"
#include "SundanceMeshBase.hpp"
#include "SundanceMaximalCofacetBatch.hpp"
#include "SundanceIncrementallyCreatableMesh.hpp"
#include "SundanceIdentityReorderer.hpp"
#include "SundanceCellReorderer.hpp"
#include "PlayaHandle.hpp"

namespace Sundance
{
using namespace Teuchos;


/**
 * Mesh is the user-level object representing discrete geometry. 
 * The Mesh class is a handle to a MeshBase, which is an abstract interface
 * for meshes. 
 */
class Mesh : public Playa::Handle<MeshBase>
{
public:

  /* */
  HANDLE_CTORS(Mesh, MeshBase);

  /** */
  int id() const {return ptr()->id();}

  /** \brief Get the ordering convention used by this mesh */
  const MeshEntityOrder& meshOrder() const {return ptr()->meshOrder();}
    
  /** 
   * Get the spatial dimension of the mesh
   */
  int spatialDim() const {return ptr()->spatialDim();}

  /** 
   * Get the number of cells having dimension dim
   */
  int numCells(int dim) const {return ptr()->numCells(dim);}

  /** 
   * Return the position of the i-th node
   */
  Point nodePosition(int i) const {return ptr()->nodePosition(i);}

  /** 
   * Return a view of the i-th node's position
   */
  const double* nodePositionView(int i) const {return ptr()->nodePositionView(i);}

  /** 
   * Return the centroid position of the cellLID-th cell of dimension
   * cellDim.
   */
  Point centroid(int cellDim, int cellLID) const 
    {return ptr()->centroid(cellDim, cellLID);}


  /** 
   * Get the outward normals for the batch of cells of dimension
   * spatialDim()-1. If any cell in the batch is not on the boundary,
   * an exception is thrown. 
   *
   * \param cellLIDs [in] LIDs for the cells whose normals are to be
   * computed. 
   * \param outwardNormals [out] Outward normal unit vectors for each
   * cell in the batch.
   */
  void outwardNormals(
    const Array<int>& cellLIDs,
    Array<Point>& outwardNormals
    ) const 
    {ptr()->outwardNormals(cellLIDs, outwardNormals);}

  /** 
   * Get tangent vectors for a batch of edges
   *
   * \param cellLIDs [in] LIDs for the cells whose tangents are to be
   * computed. 
   * \param tangentVectors [out] Unit tangents for each cell
   */
  void tangentsToEdges(
    const Array<int>& cellLIDs,
    Array<Point>& tangentVectors
    ) const 
    {ptr()->tangentsToEdges(cellLIDs, tangentVectors);}
      
      
      

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
  void getJacobians(int cellDim, const Array<int>& cellLID,
    CellJacobianBatch& jBatch) const 
    {ptr()->getJacobians(cellDim, cellLID, jBatch);}

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
    Array<double>& diameters) const 
    {ptr()->getCellDiameters(cellDim, cellLID, diameters);}


  /**
   * Map reference quadrature points to physical points on the
   * given cells. 
   */
  void pushForward(int cellDim, const Array<int>& cellLID,
    const Array<Point>& refQuadPts,
    Array<Point>& physQuadPts) const 
    {ptr()->pushForward(cellDim, cellLID, refQuadPts, physQuadPts);}

      

  /** 
   * Return the rank of the processor that owns the given cell
   */
  int ownerProcID(int cellDim, int cellLID) const 
    {return ptr()->ownerProcID(cellDim, cellLID);}
    

  /** 
   * Return the number of facets of the given cell
   */
  int numFacets(int cellDim, int cellLID, 
    int facetDim) const 
    {return ptr()->numFacets(cellDim, cellLID, facetDim);}

  /** 
   * Return the local ID of a facet cell
   * @param cellDim dimension of the cell whose facets are being obtained
   * @param cellLID local index of the cell whose
   * facets are being obtained
   * @param facetDim dimension of the desired facet
   * @param facetIndex index into the list of the cell's facets
   */
  int facetLID(int cellDim, int cellLID,
    int facetDim, int facetIndex,
    int& facetOrientation) const 
    {return ptr()->facetLID(cellDim, cellLID, 
        facetDim, facetIndex,
        facetOrientation);}

  /**
   * Return by reference argument an
   *  array containing the LIDs of the facetDim-dimensional 
   * facets of the given cell
   */
  void getFacetArray(int cellDim, int cellLID, int facetDim, 
    Array<int>& facetLIDs,
    Array<int>& facetOrientations) const 
    {ptr()->getFacetArray(cellDim, cellLID, 
        facetDim, facetLIDs,
        facetOrientations);}

  /** 
   * Return a view of an element's zero-dimensional facets
   */
  const int* elemZeroFacetView(int cellLID) const 
    {return ptr()->elemZeroFacetView(cellLID);}

  /** 
   * Return by reference argument an array containing
   * the LIDs of the facetDim-dimensional facets of the
   * given batch of cells 
   */
  void getFacetLIDs(int cellDim, 
    const Array<int>& cellLID,
    int facetDim,
    Array<int>& facetLID,
    Array<int>& facetOrientations) const 
    {ptr()->getFacetLIDs(cellDim, cellLID, 
        facetDim, facetLID, facetOrientations);}


  /** 
   * Return the number of maximal cofacets of the given cell
   */
  int numMaxCofacets(int cellDim, int cellLID) const 
    {return ptr()->numMaxCofacets(cellDim, cellLID);}

  /** 
   * Return the local ID of a maximal cofacet of a cell
   * @param cellDim dimension of the cell whose cofacets are being obtained
   * @param cellLID local index of the cell whose
   * cofacets are being obtained
   * @param cofacetIndex [in] index of the maximal cell 
   * into the list of the cell's cofacets
   * @param facetIndex [out] index of the calling cell
   * into the list of the maximal cell's facets
   */
  int maxCofacetLID(int cellDim, int cellLID,
    int cofacetIndex,
    int& facetIndex) const 
    {return ptr()->maxCofacetLID(cellDim, cellLID, cofacetIndex, facetIndex);}

  /** 
   * Get the LIDs of the maximal cofacets for a batch of 
   * codimension-one cells. 
   *
   * \param cellLIDs [in] array of LIDs of the cells whose cofacets are 
   * being obtained
   * \param cofacets [out] the batch of cofacets
   */
  void getMaxCofacetLIDs(const Array<int>& cellLIDs,
    MaximalCofacetBatch& cofacets) const 
    {ptr()->getMaxCofacetLIDs(cellLIDs, cofacets);}



  /** 
   * Find the cofacets of the given cell
   * @param cellDim dimension of the cell whose cofacets are being obtained
   * @param cellLID local index of the cell whose
   * cofacets are being obtained
   * @param cofacetDim dimension of the cofacets to get
   * @param cofacetLIDs LIDs of the cofacet
   */
  void getCofacets(int cellDim, int cellLID,
    int cofacetDim, Array<int>& cofacetLIDs) const 
    {ptr()->getCofacets(cellDim, cellLID, cofacetDim, cofacetLIDs);}

  /** 
   * Find the local ID of a cell given its global index
   */
  int mapGIDToLID(int cellDim, int globalIndex) const 
    {return ptr()->mapGIDToLID(cellDim, globalIndex);}

  /** 
   * Determine whether a given cell GID exists on this processor
   */
  bool hasGID(int cellDim, int globalIndex) const 
    {return ptr()->hasGID(cellDim, globalIndex);}

    

  /** 
   * Find the global ID of a cell given its local index
   */
  int mapLIDToGID(int cellDim, int localIndex) const 
    {return ptr()->mapLIDToGID(cellDim, localIndex);}

  /**
   * Get the type of the given cell
   */
  CellType cellType(int cellDim) const 
    {return ptr()->cellType(cellDim);}

  /** Get the label of the given cell */
  int label(int cellDim, int cellLID) const 
    {return ptr()->label(cellDim, cellLID);}

  /** Get the labels for a batch of cells */
  void getLabels(int cellDim, const Array<int>& cellLID, Array<int>& labels) const 
    {ptr()->getLabels(cellDim, cellLID, labels);}

  /** Set the label for the given cell */
  void setLabel(int cellDim, int cellLID, int label)
    {ptr()->setLabel(cellDim, cellLID, label);}

  /** Get the list of all labels defined for cells of the given dimension */
  Set<int> getAllLabelsForDimension(int cellDim) const 
    {return ptr()->getAllLabelsForDimension(cellDim);}

  /** 
   * Return the number of labels associated with the given dimension.
   */
  virtual int numLabels(int cellDim) const 
    {return getAllLabelsForDimension(cellDim).size();}

  /** 
   * Get the cells associated with a specified label. The array 
   * cellLID will be filled with those cells of dimension cellDim
   * having the given label.
   */
  void getLIDsForLabel(int cellDim, int label, Array<int>& cellLIDs) const 
    {ptr()->getLIDsForLabel(cellDim, label, cellLIDs);}

  /** Get the MPI communicator over which the mesh is distributed */
  const MPIComm& comm() const {return ptr()->comm();}

  /** \name Mesh creation methods */
  //@{
  /** Allocate space for an estimated number of vertices */
  void estimateNumVertices(int nPts) 
    {creatableMesh()->estimateNumVertices(nPts);}
    
  /** Allocate space for an estimated number of elements */
  void estimateNumElements(int nElems) 
    {creatableMesh()->estimateNumElements(nElems);}

    
  /** Add a new vertex to the mesh */
  int addVertex(int globalIndex, const Point& x,
    int ownerProcID, int label)
    {return creatableMesh()->addVertex(globalIndex, x, ownerProcID, label);}

  /** Add a new vertex to the mesh */
  int addElement(int globalIndex, const Array<int>& vertLID,
    int ownerProcID, int label)
    {return creatableMesh()->addElement(globalIndex, vertLID, 
        ownerProcID, label);}

  /** */
  void freezeTopology() 
    {creatableMesh()->freezeTopology();}

    
  /** \name Reordering */
  //@{
  /** Set the reordering method to be used with this mesh */
  void setReorderer(const CellReorderer& reorderer) 
    {ptr()->setReorderer(reorderer);}

  /** Set the reordering method to be used with this mesh */
  const CellReordererImplemBase* reorderer() const 
    {return ptr()->reorderer();}
  //@}
    
  /** */
  static CellReorderer& defaultReorderer()
    {
      static CellReorderer rtn = new IdentityReorderer();
      return rtn;
    }

  /** Work out global numberings for the cells of dimension cellDim */
  void assignIntermediateCellGIDs(int cellDim) 
    {
      if (!hasIntermediateGIDs(cellDim))
        ptr()->assignIntermediateCellGIDs(cellDim);
    }

  /** */
  bool hasIntermediateGIDs(int cellDim) const 
    {
      return ptr()->hasIntermediateGIDs(cellDim);
    }

  /** */
  void dump(const std::string& filename) const ;

  /** Test the consistency of the mesh numbering scheme 
   * across processors. This is meant as a check on Sundance's internal
   * logic rather than as a check on the validity of a user's mesh. */
  bool checkConsistency(const std::string& filename) const ;

  /** Test the consistency of the mesh numbering scheme 
   * across processors. This is meant as a check on Sundance's internal
   * logic rather than as a check on the validity of a user's mesh. */
  bool checkConsistency(std::ostream& os) const ;

  /** \name Functions for Mesh with hanging nodes */
    //@{
    /** Function returns true if the mesh allows hanging nodes (by refinement), false otherwise */
    bool allowsHangingHodes() const { return ptr()->allowsHangingHodes(); }

    /** Function returns true if the specified element is a "hanging" element
     * false otherwise */
    bool isElementHangingNode(int cellDim , int cellLID) const
    		{ return ptr()->isElementHangingNode(cellDim , cellLID); }

   /** Returns the index in the parent maxdim Cell of the refinement tree
    * @param maxCellLID [in] the LID of the cell */
    int indexInParent(int maxCellLID) const
    		{ return ptr()->indexInParent(maxCellLID); }

   /** How many children has a refined element. <br>
    * This function provides information of either we have bi or trisection */
   int maxChildren() const { return ptr()->maxChildren();}

   /** Function returns the facets of the maxdim parent cell (needed for HN treatment) */
   void returnParentFacets( int childCellLID , int dimFacets ,
    		                         Array<int> &facetsLIDs , int &parentCellLIDs ) const {
    	ptr()->returnParentFacets( childCellLID , dimFacets , facetsLIDs , parentCellLIDs );
   }
  //@}


  /** \name Special Weights Storage for Adaptive Cell Integration */
   //@{
  /** returns the status of the special weights if they are valid <br>
    *  These weights are usually computed for one setting of the curve (Adaptive Cell Integration)*/
  bool IsSpecialWeightValid() const {return ptr()->IsSpecialWeightValid();}

  /** specifies if the special weights are valid <br>
   *  if this is false then usually the special weights have to be recomputed */
  void setSpecialWeightValid(bool val) const { ptr()->setSpecialWeightValid(val);}

  /** deletes all the special weights*/
  void flushSpecialWeights() const { ptr()->flushSpecialWeights(); }

  /** verifies if the specified cell with the given dimension has special weights */
  bool hasSpecialWeight(int dim, int cellLID) const {return ptr()->hasSpecialWeight( dim, cellLID); }

  /** Sets the special weights */
  void setSpecialWeight(int dim, int cellLID, Array<double>& w) const {ptr()->setSpecialWeight(dim, cellLID, w);}

  /** Returns the special weights */
  void getSpecialWeight(int dim, int cellLID, Array<double>& w) const {ptr()->getSpecialWeight(dim, cellLID, w);}
   //@}



  /** \name Store the intersection/quadrature points for the curve/surf integral <br>
   *  for a curve or surf integral we need some quadrature points along the curve in one curve <br>
   *  These */
    //@{

  /** */
  bool IsCurvePointsValid() const {return ptr()->IsCurvePointsValid();}

  /**  */
  void setCurvePointsValid(bool val)  const {ptr()->setCurvePointsValid(val); }

  /** deletes all the curve points */
  void flushCurvePoints() const { ptr()->flushCurvePoints(); }

  /** verifies if the specified maxCell has already precalculated quadrature point for one curve */
  bool hasCurvePoints(int maxCellLID , int curveID) const { return ptr()->hasCurvePoints( maxCellLID , curveID); }

  /** Sets the points, curve derivatives and curve normals for one maxCell needed for curve/surf integral*/
  void setCurvePoints(int maxCellLID, int curveID , Array<Point>& points , Array<Point>& derivs , Array<Point>& normals) const
		{ptr()->setCurvePoints( maxCellLID, curveID , points , derivs , normals); }

  /** Gets the points, curve derivatives and curve normals for one maxCell needed for curve/surf integral*/
  void getCurvePoints(int maxCellLID,  int curveID , Array<Point>& points , Array<Point>& derivs , Array<Point>& normals) const
		{ptr()->getCurvePoints( maxCellLID,  curveID ,  points , derivs , normals); }

  //@}

private:
  /** */
  IncrementallyCreatableMesh* creatableMesh();


  /** */
  bool checkVertexConsistency(std::ostream& os) const ;
  /** */
  bool checkCellConsistency(std::ostream& os, int dim) const ;

  /** */
  bool checkRemoteEntity(std::ostream& os, int p, int dim, int gid, 
    int owner, bool mustExist, int& lid) const ;

  /** */
  bool testIdentity(std::ostream& os, int a, int b, const std::string& msg) const ;

  /** */
  bool testIdentity(std::ostream& os, 
    const Array<int>& a,
    const Array<int>& b, const std::string& msg) const ;
};
}

#endif
