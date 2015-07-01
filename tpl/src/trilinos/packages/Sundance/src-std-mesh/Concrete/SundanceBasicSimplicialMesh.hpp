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

#ifndef SUNDANCE_BASICSIMPLICIALMESH_H
#define SUNDANCE_BASICSIMPLICIALMESH_H



#include "SundanceDefs.hpp"
#include "SundanceMeshBase.hpp"
#include "SundanceSet.hpp"
#include "SundanceBasicVertexView.hpp"
#include "SundanceArrayOfTuples.hpp"
#include "SundanceIncrementallyCreatableMesh.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Hashtable.hpp"

namespace Sundance
{

/**
 * A no-frills parallel simplicial mesh
 */
class BasicSimplicialMesh : public IncrementallyCreatableMesh
{
public:
  /** */
  BasicSimplicialMesh(int dim, const MPIComm& comm, 
    const MeshEntityOrder& order);

  /** */
  virtual ~BasicSimplicialMesh(){;}

  /** 
   * Get the number of cells having dimension dim
   */
  virtual int numCells(int dim) const  ;

  /** 
   * Return the position of the i-th node
   */
  virtual Point nodePosition(int i) const {return points_[i];}

  /** 
   * Return a view of the i-th node's position
   */
  const double* nodePositionView(int i) const {return &(points_[i][0]);}

     

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
    CellJacobianBatch& jBatch) const  ;

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
    Array<double>& diameters) const ;


  /**
   * Map reference quadrature points to physical points on the
   * given cells. 
   */
  virtual void pushForward(int cellDim, const Array<int>& cellLID,
    const Array<Point>& refQuadPts,
    Array<Point>& physQuadPts) const ;

      

  /** 
   * Return the rank of the processor that owns the given cell
   */
  virtual int ownerProcID(int cellDim, int cellLID) const  ;

      

  /** 
   * Return the number of facets of the given cell
   */
  virtual int numFacets(int cellDim, int cellLID, 
    int facetDim) const  ;

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
    int& facetOrientation) const  ;


  /** 
   * Return by reference argument an array containing
   * the LIDs of the facetDim-dimensional facets of the
   * given batch of cells 
   */
  virtual void getFacetLIDs(int cellDim, 
    const Array<int>& cellLID,
    int facetDim,
    Array<int>& facetLID,
    Array<int>& facetOrientations) const ;

  /** 
   * Return a view of an element's zero-dimensional facets
   */
  const int* elemZeroFacetView(int cellLID) const 
    {return &(elemVerts_.value(cellLID, 0));}

  /** 
   * Return the number of maximal cofacets of the given cell
   */
  virtual int numMaxCofacets(int cellDim, int cellLID) const  ;

  /** 
   * Return the local ID of a maximal cofacet cell
   * @param cellDim dimension of the cell whose cofacets are being obtained
   * @param cellLID local index of the cell whose
   * cofacets are being obtained
   * @param cofacetIndex which maximal cofacet to get
   * @param cofacetIndex index of the cell cellLID into the list of the 
   * maximal cell's facets
   */
  virtual int maxCofacetLID(int cellDim, int cellLID,
    int cofacetIndex,
    int& facetIndex) const  ;

#ifdef MOVED_TO_BASE_CLASS
  /** 
   * Get the LIDs of the maximal cofacets for a batch of codimension-one
   * cells.
   *
   * \param cellLIDs [in] array of LIDs of the cells whose cofacets are 
   * being obtained
   * \param cofacets [out] the batch of cofacets
   */
  virtual void getMaxCofacetLIDs(const Array<int>& cellLIDs,
    MaximalCofacetBatch& cofacets) const ;
#endif
      
  /** 
   * Find the cofacets of the given cell
   * @param cellDim dimension of the cell whose cofacets are being obtained
   * @param cellLID local index of the cell whose
   * cofacets are being obtained
   * @param cofacetDim dimension of the cofacets to get
   * @param cofacetLIDs LIDs of the cofacet
   */
  void getCofacets(int cellDim, int cellLID,
    int cofacetDim, Array<int>& cofacetLIDs) const ;
      
  /** 
   * Find the local ID of a cell given its global index
   */
  virtual int mapGIDToLID(int cellDim, int globalIndex) const  ;

  /** 
   * Determine whether a given cell GID exists on this processor
   */
  virtual bool hasGID(int cellDim, int globalIndex) const ;

    

  /** 
   * Find the global ID of a cell given its local index
   */
  virtual int mapLIDToGID(int cellDim, int localIndex) const  ;

  /**
   * Get the type of the given cell
   */
  virtual CellType cellType(int cellDim) const  ;


  /** Get the label of the given cell */
  virtual int label(int cellDim, int cellLID) const ;

  /** Get the labels for a batch of cells */
  virtual void getLabels(int cellDim, const Array<int>& cellLID, 
    Array<int>& labels) const ;

      

  /** Get the list of all labels defined for cells of the given dimension */
  virtual Set<int> getAllLabelsForDimension(int cellDim) const ;

  /** 
   * Get the cells associated with a specified label. The array 
   * cellLID will be filled with those cells of dimension cellDim
   * having the given label.
   */
  virtual void getLIDsForLabel(int cellDim, int label, Array<int>& cellLIDs) const ;

  /** \name Incremental creation methods */
  //@{
  /** 
   * Add new new vertex to the mesh.
   * \param globalIndex the GID of the new vertex
   * \param x the spatial position of the new vertex
   * \param ownerProcID the processor that "owns" this vertex 
   * \param label a label for this vertex (optionally used in setting 
   * loads, boundary conditions, etc)
   * \return the LID of the vertex.
   */
  virtual int addVertex(int globalIndex, const Point& x,
    int ownerProcID, int label);

  /** 
   * Add a new element to the mesh.
   * \param globalIndex the GID of the new element
   * \param vertexGIDs tuple of GIDs for the vertices defining this element. 
   * \param ownerProcID the processor that "owns" this element
   * \param label a label for this element (optionally used in setting loads, 
   * material properties, etc)
   * \return the LID of the element
   */
  virtual int addElement(int globalIndex, const Array<int>& vertexGIDs,
    int ownerProcID, int label);

  /** Set the label of the given cell */
  virtual void setLabel(int cellDim, int cellLID, int label)
    {
      labels_[cellDim][cellLID] = label;
    }

  /** Optional preallocation of space for an estimated number of vertices */
  virtual void estimateNumVertices(int numVertices);

  /** Optional preallocation of space for an estimated number of elements */
  virtual void estimateNumElements(int numElements);

  /** Coordinate intermediate cell definitions across processors  */
  virtual void assignIntermediateCellGIDs(int cellDim) ;

  /** */
  virtual bool hasIntermediateGIDs(int dim) const 
    {
      if (dim==1) return hasEdgeGIDs_;
      return hasFaceGIDs_;
    }
      
  //@}
private:

  /** 
   * Add a new face, first checking to see if it already exists. 
   * This function is called from within addElement(),
   * not by the user, and is therefore private.
   *
   * \param vertGID The sorted GIDs for the three vertices of the face 
   * \param vertLID The LIDs for the three vertices of the face 
   * \param edgeLID{1,2,3} The LIDs for the three edges of the face 
   * \return LID of this face
   */
  int addFace(const Array<int>& vertLID, 
    const Array<int>& vertGID,
    const Array<int>& edgeGID,
    int elemLID,
    int elemGID);

  /**
   * Add a new edge, first checking to see if it already exists. 
   * This function is called from within addElement(),
   * not by the user, and is therefore private.
   *
   * \param vertLID{1,2} 
   * \param elemLID LID of the element that is adding the edge
   * \param elemGID GID of the element that is adding the edge
   * \param myFacetNumber facet number of the edge within the element
   * \return LID of this edge
   */
  int addEdge(int vertLID1, int vertLID2, int elemLID, int elemGID, 
    int myFacetNumber);

  /** 
   * Check for the presence of the edge (vertLID1, vertLID2) in the mesh.
   * \return the edge LID if the edge exists, -1 otherwise
   */
  int checkForExistingEdge(int vertLID1, int vertLID2);

  /** 
   * Check whether the face defined by the given vertices exists
   * in the mesh. Returns -1 if the face does not exist.
   * Called during the synchronization of intermediate cell GIDs.
   * \param vertGID{1,2,3} the global indices of the vertices defining the face
   * \return the LID of the face
   */
  int lookupFace(const Array<int>& vertGID) ;

  /** */
  void synchronizeGIDNumbering(int dim, int localCount) ;
      
  /** */
  void resolveEdgeOwnership(int cellDim);

  /** */
  std::string cellStr(int dim, const int* verts) const ;

  /** */
  std::string cellToStr(int dim, int cellLID) const ;

  /** */
  std::string printCells(int dim) const ;

  /** */
  void synchronizeNeighborLists();
      
      

  /** Number of cells of each dimension. */
  Array<int> numCells_;
    
  /** coordinates of vertices. The index into the array is the vertex LID .*/
  Array<Point> points_;

  /** pairs of local vertex indices for the edges, each pair ordered
   * from lower to higher <i>global</i> vertex index in order to define
   * an absolute edge orientation. Because global vertex indices are used, all
   * processors will agree on this orientation, regardless of the orientation
   * of the edge as seen by the element first defining it. 
   * The first index into this 2D array is the edge LID, the second the vertex 
   * number within the edge. 
   * */
  ArrayOfTuples<int> edgeVerts_;

  /** Tuples of local vertex indices for the faces, with each tuple ordered from
   * lowest to highest <i>global</i> index in order to define an absolute edge
   * orientation. Because global vertex indices are used, all
   * processors will agree on this orientation, regardless of the orientation
   * of the face as seen by the element first defining it. 
   * The first index into this 2D array is the face LID, the second the
   * vertex number within the face. */
  ArrayOfTuples<int> faceVertLIDs_;

  /** Tuples of global vertex indices for the faces, with each tuple ordered from
   * lowest to highest <i>global</i> index in order to define an absolute edge
   * orientation. Because global vertex indices are used, all
   * processors will agree on this orientation, regardless of the orientation
   * of the face as seen by the element first defining it. 
   * The first index into this 2D array is the face LID, the second the
   * vertex number within the face. 
   *
   * Notice that we duplicate the face vertex storage, storing both the
   * vertex LIDs and vertex GIDs for each face. This lets us do quick comparison
   * with the sorted GID array in order to identify pre-existing faces, while
   * also making it possible to retrieve face vertex LID information without
   * doing hashtable lookups.
   *
   */
  ArrayOfTuples<int> faceVertGIDs_;
    
  /** Tuples of local indices for the edges of all faces. The first index
   * into this 2D array is the face LID, the second the edge number. */
  ArrayOfTuples<int> faceEdges_;
    
  /** Tuples of edge signs for the faces. The edge sign indicates 
   * whether the orientation of the edge as given by moving around the face 
   * is parallel or antiparallel to the absolute orientation of the edge. */
  ArrayOfTuples<int> faceEdgeSigns_;
    
  /** tuples of local vertex indices for the elements. The first index into this
   * 2D array is the element LID, the second is the vertex number.  */
  ArrayOfTuples<int> elemVerts_;

  /** tuples of local edge indices for the elements. The first index into
   * this 2D array is the element LID, the second is the edge number. */
  ArrayOfTuples<int> elemEdges_;

  /** tuples of edge orientations for the elements, indicating whether 
   * the orientation of the edge as given by moving around the element 
   * is parallel or antiparallel to the absolute orientation of the edge. 
   * The first index into this 2D array is the element LID, the second 
   * the edge number. 
   * */
  ArrayOfTuples<int> elemEdgeSigns_;

  /** tuples of face LIDs for the elements. The first index is the
   * element LID, the second the face number. */
  ArrayOfTuples<int> elemFaces_;

  /** tuples of face rotations for the elements, defined relative to the 
   * absolute orientation of the face. */
  ArrayOfTuples<int> elemFaceRotations_;

  /** table for mapping vertex set -> face index */
  Hashtable<VertexView, int> vertexSetToFaceIndexMap_;

  /** array of face cofacets for the edges. The first index
   * is the edge LID, the second the cofacet number. */
  Array<Array<int> > edgeFaces_;

  /** array of element cofacets for the edges. The first index
   * is the edge LID, the second the cofacet number. */
  Array<Array<int> > edgeCofacets_;

  /** array of element cofacets for the faces. The first index is the
   * face LID, the second the cofacet number. */
  Array<Array<int> > faceCofacets_;

  /** array of edge cofacets for the vertices. The first index is the 
   * vertex LID, the second the edge cofacet number. */
  Array<Array<int> > vertEdges_;

  /** array of face cofacet LIDs for the vertices. The first index is the 
   * vertex LID, the second the cofacet number. */
  Array<Array<int> > vertFaces_;

  /** array of maximal cofacets for the vertices. The first index is the
   * vertex LID, the second the cafacet number. */
  Array<Array<int> > vertCofacets_;

  /** array of edge partners for the vertices. The partners are other
   * vertices sharing an edge with the specified vertex. */
  Array<Array<int> > vertEdgePartners_;

  /** map from local to global cell indices. The first index into this
   * 2D array is the cell dimension, the second the cell LID. */
  Array<Array<int> > LIDToGIDMap_;

  /** map from global to local cell indices. The array index is the 
   * cell dimension. The hashtable key is the cell GID, the value the
   * cell LID. */
  Array<Hashtable<int, int> > GIDToLIDMap_;
    
  /** Array of labels for the cells */
  Array<Array<int> > labels_;

  /** Array of owning processor IDs for the cells. Each cell is owned by
   * a single processor that is responsible for assigning global indices,
   * DOF numbers, and so on. */
  Array<Array<int> > ownerProcID_;
      
  /** 
   * Pointer to the pointer at the base of the face vertex GID array. This is
   * used to get small-array "views" of the face vertex GID array without
   * making copies, resulting in a significant performance improvement
   * in the vertex set hashtable lookups to identify pre-existing faces. 
   *
   * We use double rather than single indirection here because as elements
   * are added, the face vertex GID array will often be resized, thus changing
   * the base pointer. Each vertex set "view" keeps a pointer to the base pointer,
   * so that it always remains synchronized with the array of face vertices. 
   * 
   * IMPORTANT: any time faceVertGIDs_ is resized, faceVertGIDBase_[0] must
   * be reset to the base of the faceVertGIDs_ array so that the vertex
   * sets are pointing to the right place.
   */
  Array<int*> faceVertGIDBase_;


  /** flag indicating whether the edge GIDs have been synchronized */
  bool hasEdgeGIDs_;

  /** flag indicating whether the face GIDs have been synchronized */
  bool hasFaceGIDs_;

  /** Set of all neighbor processors sharing data with this one */
  Set<int> neighbors_;

  /** Whether the neighbor lists have been synchronized across 
   * processors */
  bool neighborsAreSynchronized_;


};

}




#endif
