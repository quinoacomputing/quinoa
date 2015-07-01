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
/*
 * SundanceHNMesh3D.h
 *
 *  Created on: May 30, 2010
 *      Author: benk
 */

#ifndef SUNDANCEHNMESH3D_H_
#define SUNDANCEHNMESH3D_H_

#include "SundanceDefs.hpp"
#include "SundanceMeshBase.hpp"
#include "SundanceSet.hpp"
#include "SundancePoint.hpp"
#include "SundanceCellType.hpp"

#include "Teuchos_Array.hpp"
#include "Teuchos_Hashtable.hpp"

#include "SundanceRefinementBase.hpp"
#include "SundanceRefinementClass.hpp"

#include "SundanceDomainDefinition.hpp"


namespace Sundance
{
    using namespace Sundance;
    using namespace Teuchos;

/**
 * Class for 3D hierarchical structured quad Mesh <br>
 * The main internal idea of the mesh ist that there are different numberings of the elements
 * ID -> each element has an ID (also those which are completely outside the mesh domain, and those
 * which are not leaf) <br>
 * In the code sometimes internaly LID is used instead of ID due to legacy !! :-P <br>
 * GID -> all leaf elements which are visible to Sundance (only the leaf elements should be visible) <br>
 * LID -> all the elements which have GID and belong and either belong to the local processor or are in the ghost cells <br>*/

class HNMesh3D : public MeshBase{

public:
	/**
	 * The Ctor for the dummy grid with hanging nodes */
	HNMesh3D(int dim,
    const MPIComm& comm,
    const MeshEntityOrder& order);

  /**
   * The Ctor for the HNMesh3D grid in 3D*/
	void createMesh(
    double position_x,
    double position_y,
    double position_z,
    double offset_x,
    double offset_y,
    double offset_z,
    int resolution_x,
    int resolution_y,
    int resolution_z,
    const RefinementClass& refineClass,
    const MeshDomainDef& meshDomain );

	/** Dtor */
	virtual ~HNMesh3D(){;}

  /**
   * Get the number of cells having dimension dim
   */
  virtual int numCells(int dim) const;

  /**
   * Return the position of the i-th node
   */
  virtual Point nodePosition(int i) const;

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
  virtual int ownerProcID(int cellDim, int cellLID) const;

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
  virtual int mapGIDToLID(int cellDim, int globalIndex) const;

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
  virtual void assignIntermediateCellGIDs(int cellDim);

  /** */
  virtual bool hasIntermediateGIDs(int dim) const;


  /** Function returns true if the mesh allows hanging nodes (by refinement),
   * false otherwise */
    virtual bool allowsHangingHodes() const { return true; }

  /** Function returns true if the specified element is a "hanging" element
   * false otherwise <br>
   * @param dim [in] should be between 0 , D-1
   * @param cellLID [in] the local ID of the element */
    virtual bool isElementHangingNode(int cellDim , int cellLID) const ;

    /** Returns the index in the parent maxdim Cell of the refinement tree
     * @param maxCellLID [in] the LID of the cell */
    virtual int indexInParent(int maxCellLID) const ;

   /** How many children has a refined element. <br>
    * This function provides information of either we have bi or trisection */
    virtual int maxChildren() const {return 27;}

  /** Function returns the facets of the maxdim parent cell (needed for HN treatment) <br>
   * @param childCellLID [in] the LID of the maxdim cell, whos parents facets we want
   * @param dimFacets [in] the dimension of the facets which we want to have
   * @param facetsLIDs [out] the LID of the parents facets (all) in the defined order
   * @param parentCellLIDs [out] the maxdim parent cell LID */
    virtual void returnParentFacets( int childCellLID , int dimFacets ,
           Array<int> &facetsLIDs , int &parentCellLIDs ) const;

private:

   /** For HN , returns parent facets, if the facet is not leaf, then return -1 at that place */
   int facetLID_tree(int cellDim, int cellLID,
                        int facetDim, int facetIndex) const;

   /** adds one vertex to the mesh */
   void addVertex(int vertexLID , int ownerProc , bool isHanging ,
		         double coordx , double coordy , double coordz , const Array<int> &maxCoF);

   /** adds one edge to the mesh */
   void addEdge(int edgeLID , int ownerProc , bool isHanging , int edgeOrientation ,
		        const Array<int> &vertexLIDs , const Array<int> &maxCoF);

   /** adds one edge to the mesh */
   void addFace(int faceLID , int ownerProc , bool isHanging , int faceOrientation ,
		        const Array<int> &vertexLIDs , const Array<int> &edgeLIDs ,
		        const Array<int> &maxCoF);

   /** adds one cell(3D) to the mesh <br>
    * cell must be always leaf*/
   void addCell(int cellLID , int ownerProc ,
           int indexInParent , int parentCellLID , int level ,
           const Array<int> &faceLIDs , const Array<int> &edgeLIDs ,
           const Array<int> &vertexLIDs) ;

   /** creates the mesh on the coarsest level as it is specified */
   void createCoarseMesh();

   /** Does one single refinement iteration. <br>
    *  Iterates trough all cells which are owned by the processor and refines if necessary */
   bool oneRefinementIteration();

   /** refine the given cell by cellID, (method assumes that this cell can be refined)*/
   void refineCell(int cellLID);

   /** Create Leaf numbering */
   void createLeafNumbering();

   /** Create Leaf numbering, but with a better load balancing for parallel case */
   void createLeafNumbering_sophisticated();

   /** estimate the load of one cell*/
   int estimateCellLoad(int ID);

   /** marks the cells recursivly in the tree (and facets) owned by one processor*/
   void markCellsAndFacets(int cellID , int procID);

   /** this updates the array with the local index of vertex, edge and face in the static array <br>
    * this method is only used in the coarse mesh creation */
   void updateLocalCoarseNumbering(int ix , int iy , int iz , int Nx , int Ny);

   /** Used for refinemement
    * @param cellDim the requested elements dimension (if exists)
    * @param usedge , whould we look for that element in edge or in face
    * @param parentID , ID of the parent edge or face
    * @param parentOffset , the offset in the parent */
   int getHangingElement(int cellDim, bool useEdge, int parentID , int parentOffset);

   /** Used for refinemement
    * @param cellDim the requested elements dimension (if exists)
    * @param cellID , ID of the new hanging cell
    * @param usedge , whould we look for that element in edge or in face
    * @param parentID , ID of the parent edge or face
    * @param parentOffset , the offset in the parent */
   void addHangingElement(int cellDim, int cellID ,bool useEdge, int parentID , int parentOffset);

   /**
    * @param cellDim
    * @param cellID */
   int numMaxCofacets_ID(int cellDim, int cellID);

  /** Number of processors */
  int nrProc_;
  int myRank_;
  /** The communicator */
  const MPIComm& _comm;
  /** */
  double _pos_x;
  /** */
  double _pos_y;
  /** */
  double _pos_z;
  /** */
  double _ofs_x;
  /** */
  double _ofs_y;
  /** */
  double _ofs_z;
  /**  */
  int _res_x;
  /**  */
  int _res_y;
  /**  */
  int _res_z;
  /** */
  mutable RefinementClass refineClass_;
  /** */
  mutable MeshDomainDef meshDomain_;

//------ Point storage ----
  /** all the global points index is [ID] */
  Array<Point> points_;
  /** [4] the nr of ID per dim*/
  Array<int> nrElem_;
  /** [4] the nr of owned elements per dim*/
  Array<int> nrElemOwned_;

//----- Facets -----  -1 -> MeshBoundary

  /** [cellID][8]*/
  Array< Array<int> > cellsPoints_;
  /** [cellID][12]*/
  Array< Array<int> > cellsEdges_;
  /** [cellID][6]*/
  Array< Array<int> > cellsFaces_;
  /** [faceID][4]*/
  Array< Array<int> > faceEdges_;
  /** [faceID][4]*/
  Array< Array<int> > facePoints_;
  /** [edgeID][2]*/
  Array< Array<int> > edgePoints_;

  /** stores the edge orientation {0,1,2}
   *  [edgeID]*/
  Array<short int> edgeOrientation_;
  /** stores the face orientation {0,1,2}
   *  [faceID]*/
  Array<short int> faceOrientation_;

// ----- MaxCofacets ----;  -1 -> MeshBoundary , -2 -> ProcBoundary

  /** [faceID][2] */
  Array< Array<int> > faceMaxCoF_;
  /** [edgeID][4] */
  Array< Array<int> > edgeMaxCoF_;
  /** [pointID][8]*/
  Array< Array<int> > pointMaxCoF_;

//------ Element (processor) ownership -----

  /** contains the ownership of the local elements [dim][ID] */
  Array< Array<short int> > elementOwner_;

//---- hierarchical storage -----

  /** [cellID] , the child index in the parent */
  Array<short int> indexInParent_;
  /** [cellID] , the LID of the parent cell */
  Array<int> parentCellLID_;
  /** [cellID] , actual level of the cell*/
  Array<short int> cellLevel_;
  /** [cellID] , if the element is leaf */
  Array<bool> isCellLeaf_;
  /** [cellID] , if the cell is complete outside the user defined mesh domain*/
  Array<bool> isCellOut_;
  /** [cellID] , children of the cell*/
  Array< Array<int> > cellsChildren_;

  // ---- "hanging" info storage ---
  /** [pointID] , true if the node is hanging , false otherwise */
  Array<bool> isPointHanging_;
  /** [edgeID] , true if the edge is hanging , false otherwise*/
  Array<bool> isEdgeHanging_;
  /** [faceID] , true if the face is hanging , false otherwise*/
  Array<bool> isFaceHanging_;

// ---- hanging element and refinement (temporary) storage ---

  /** [edgeID] - > { h P0 LID , h P1 LID , h E0 LID , h E1 LID , h E2 LID } <br>
   * These elements will only be put on not hanging when we access it 3 times <br>
   * together 5 elements*/
  Hashtable< int, Array<int> > edgeHangingElmStore_;

  /** the counter for each edge which stores hanging node information, when the counter reaches 2 (3 read access)<br>
   * then the hanging elements can be marked as not hanging <br>
   * At the beginning this should be equal with the numMaxCofacets */
  Array<short int> hangingAccessCount_;

  /** [faceID] - > { h P0 LID , h P1 LID , h P2 LID , h P3 LID , <br>
   *                 h E0 LID , h E1 LID , h E2 LID , h E3 LID , h E4 LID , h E5 LID , (4+)<br>
   *                 h E6 LID , h E7 LID , h E8 LID , h E9 LID , h E10 LID, h E11 LID , (16+) <br>
   *                 h F0 LID , h F1 LID , h F2 LID , h F3 LID , h F4 LID , h F5 LID , <br>
   *                 h F6 LID , h F7 LID , h F8 LID } <br>
   *                 together 16+9 = 25 elements*/
  Hashtable< int, Array<int> > faceHangingElmStore_;

  /** Neighbor Cell can mark the cell to provoke refinement */
  Array<short int> refineCell_;

// ---- leaf mapping , GID and LID --- (points need this because)

  /** [leaf_vertexLID] , the value must be a positive number */
  Array<int> vertexLeafToLIDMapping_;
  /** [leaf_edgeLID] , the value must be a positive number */
  Array<int> edgeLeafToLIDMapping_;
  /** [leaf_faceLID] , the value must be a positive number */
  Array<int> faceLeafToLIDMapping_;
  /** [leaf_cellLID] , the value must be a positive number */
  Array<int> cellLeafToLIDMapping_;
  /** [vertexLID] if vertex is inside the domain then > 0 , -1 otherwise */
  Array<int> vertexLIDToLeafMapping_;
  /** [edgeLID] if edge is leaf(or inside the domain) then > 0 , -1 otherwise */
  Array<int> edgeLIDToLeafMapping_;
  /** [faceLID] if face is leaf(or inside the domain) then > 0 , -1 otherwise */
  Array<int> faceLIDToLeafMapping_;
  /** [cellLID] if cell is leaf(or inside the domain) then > 0 , -1 otherwise */
  Array<int> cellLIDToLeafMapping_;

  /** leaf LID numbering */
  int nrVertexLeafLID_;
  int nrEdgeLeafLID_;
  int nrFaceLeafLID_;
  int nrCellLeafLID_;

  /** [leaf_vertexGID] , the value must be a positive number */
  Array<int> vertexLeafToGIDMapping_;
  /** [leaf_edgeGID] , the value must be a positive number */
  Array<int> edgeLeafToGIDMapping_;
  /** [leaf_faceGID] , the value must be a positive number */
  Array<int> faceLeafToGIDMapping_;
  /** [leaf_cellGID] , the value must be a positive number */
  Array<int> cellLeafToGIDMapping_;
  /** [vertexGID] if vertex is inside the domain then > 0 , -1 otherwise */
  Array<int> vertexGIDToLeafMapping_;
  /** [edgeGID] if edge is leaf(or inside the domain) then > 0 , -1 otherwise */
  Array<int> edgeGIDToLeafMapping_;
  /** [faceGID] if edge is leaf(or inside the domain) then > 0 , -1 otherwise */
  Array<int> faceGIDToLeafMapping_;
  /** [cellGID] if cell is leaf(or inside the domain) then > 0 , -1 otherwise */
  Array<int> cellGIDToLeafMapping_;

  /** leaf GID numbering */
  int nrVertexLeafGID_;
  int nrEdgeLeafGID_;
  int nrFaceLeafGID_;
  int nrCellLeafGID_;

// ------------- static data -------------------

  /** the offset in the X coordinate on the reference cell*/
  static int offs_Points_x_[8];

  /** the offset in the Y coordinate on the reference cell*/
  static int offs_Points_y_[8];

  /** the offset in the Z coordinate on the reference cell*/
  static int offs_Points_z_[8];

  /** stores the facet information on the reference Cell*/
  static int edge_Points_localIndex[12][2];

  /** the edge orientation (local orientation)*/
  static int edge_Orientation[12];
  static int edge_MaxCofacetIndex[3][4];
  static int edge_MaxCof[12];

  /** face edge-facet information */
  static int face_Edges_localIndex[6][4];

  /** face point-facet information */
  static int face_Points_localIndex[6][4];

  /** face orientation (local orientation)*/
  static int face_Orientation[6];
  static int face_MaxCofacetIndex[3][2];
  static int face_MaxCof[6];

  /** used for coarse mesh creation for the vertex indexing */
  static int vInd[8];

  /** used for coarse mesh creation for the edge indexing */
  static int eInd[12];

  /** used for coarse mesh creation for the face indexing */
  static int fInd[6];

  // the X and the Y coordinates of the newly
  static double vertex_X[64];

  static double vertex_Y[64];

  static double vertex_Z[64];

  // face index is above 20
  static int vertexToParentEdge[64];
  //
  static int vertexInParentIndex[64];
  //
  static int edgeToParentEdge[144];
  //
  static int edgeInParentIndex[144];
  //
  static int faceToParentFace[108];
  //
  static int faceInParentIndex[108];

};
}


#endif /* SUNDANCEHNMESH3D_H_ */
