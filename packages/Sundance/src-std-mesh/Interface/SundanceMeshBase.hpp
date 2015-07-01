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

#ifndef SUNDANCE_MESHBASE_H
#define SUNDANCE_MESHBASE_H

#include "SundanceDefs.hpp"
#include "SundancePoint.hpp"
#include "SundanceSet.hpp"
#include "SundanceMap.hpp"
#include "SundanceCellType.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceObjectWithInstanceID.hpp"
#include "PlayaMPIComm.hpp"
#include "Teuchos_Array.hpp"
#include "SundanceCellReorderer.hpp"
#include "SundanceCellReordererImplemBase.hpp"

namespace Sundance {


/** Identifier for ordering convention */
enum MeshEntityOrder {UFCMeshOrder, ExodusMeshOrder};

class MaximalCofacetBatch;

using namespace Teuchos;
using Playa::MPIComm;

class CellJacobianBatch;


/** \brief Abstract interface to a mesh.
 *
 * \section Sundance_MeshBase_Outline_sec Outline
 *
 * <ul>
 * <li>\ref Sundance_MeshBase_Introduction_sec
 * <li>\ref Sundance_MeshBase_Definitions_sec
 * <li>\ref Sundance_MeshBase_Common_Function_Arguments_sec
 * <li>\ref Sundance_MeshBase_Refactorings_sec
 * </ul>
 *
 * \section Sundance_MeshBase_Introduction_sec Introduction
 *
 * A mesh is a collection of connected cells (also called elements).  This
 * mesh interface is designed to allow for the representation of 1D, 2D, and 3D (and
 * even higher dimensional) meshes.  A mesh contains both topological
 * information and geometric information about the cells and their subparts
 * (i.e. facets).  Topological information refers to information about how
 * cells of various types and dimensions are connected to each other within
 * the mesh.  Geometric information refers to information about the position,
 * size, and shape of cells within the mesh.
 *
 * Currently, this interface only handles meshes where all cells of a given
 * cell dimension have the same cell type.  For example, this only allows
 * meshes with all triangles or all quads in 2D or all tetrahedrals in 3D.
 *
 * \todo Add some diagrams to show different cell types to help define cells,
 * facets, co-facets, vertices, etc.
 *
 * \section Sundance_MeshBase_Definitions_sec Definitions
 *
 * The following definitions are used in the specification of this mesh
 * interface.
 *
 *<ul>
 *
 * <li><b>Topological information</b>
 *
 *   <ul>
 *
 *   <li><b>Spatial Dimension (spatialDim)</b>: The maximum cell dimension of
 *   any cell type in the mesh.  For example, a 3D mesh will have
 *   spatialDim=3, a 2D mesh will have spatialDim=2 and a 1D mesh will have
 *   spatialDim=1.  Typically, when we use the short hand 1D, 2D or 3D, we are
 *   referring to the spatial dimension of the mesh.  Note that 4D and
 *   higher-dimensional meshes are also representable by this interface in
 *   principle with some small changes (e.g. like making the size of a point
 *   more flexible)
 *
 *   <li><b>Cell</b>: A mesh object of some dimension.  For example, in a 3D
 *   mesh, cells are things like elements, faces, edges, and nodes.
 *
 *   <li><b>Cell Type (cellType)</b>: An enumeration of type <tt>CellType</tt>
 *   that lists some common cell types.  For example, common cell types are
 *   things like vertexes (0D), lines or edges (1D), triangles (2D), and
 *   quadrilaterals (2D).
 *   
 *   <li><b>Cell Dimension (cellDim)</b>: The dimension of a cell.  More
 *   formally, the cell dimension is equal to the number of lower level facet
 *   types in the cell.  For example, A 3D tetrahedral cell (or element) has
 *   cellDim=3 and it has the three lower-level facet types faces (cellDim=2),
 *   edges (cellDim=1), and vertexes (cellDim=0).
 *
 *   <li><b>Relative Cell</b>: The term "relative cell" is used below to refer
 *   to the cell object for which we consider other cells called "facets" and
 *   "co-facets".
 *
 *   <li><b>Facet</b>: A lower dimensional cell directly connected to and
 *   contained within a relative <em>parent</em> cell. For example, a 3D cell
 *   (element) has three different facet types: faces (facetDim=2), edges
 *   (facetDim=1), and vertexes (facetDim=0).  A face (cellDim=2) relative
 *   cell (or a 2D element) has two different facet types: edges (facetDim=1),
 *   and vertexes (facetDim=0).  A particular relative cell
 *   <tt>(cellDim,cellLID)</tt> has
 *   <tt>nf=numFacets(cellDim,cellLID,facetDim)</tt> facets of a particular
 *   dimension <tt>facetDim</tt>.
 *
 *   <li><b>Facet Dim (facetDim)</b>: The cell dimension of a given facet.
 *   For example, the dimension of a face facet of a 3D element is facetDim=2.
 *
 *   <li><b>Facet Index (facetIndex)</b>: This is the local index of the facet
 *   with respect to its relative parent cell.  For example, a 3D brick
 *   relative cell/element will have 6 faces which can be accessed using the
 *   facet indexes 0, 1, 2, 3, 4, and 5.
 *
 *   <li><b>Co-facet</b>: A higher dimensional cell of which the given
 *   relative cell is a facet.  For example, the co-facets of a relative
 *   face cell in a 3D mesh are the one or two 3D element cells that share
 *   this face.  In the case of a boundary face in a 3D mesh, only one
 *   co-facet is connected to the boundary cell and that is its boundary
 *   element.  Note that in 3D, edges can have many different face and element
 *   co-facets.  There are some important properties of co-facets.  Relative
 *   cells of dimension cellDim=spatialDim-1 will have either one or two
 *   co-facets of dimension cofacetDim=spatialDim.  Therefore, relative cells
 *   of the dimension cellDim=spatialDim-1 that have only one co-facet must be
 *   a boundary cell (i.e. a boundary face in 3D, a boundary edge in 2D, or a
 *   boundary vertex in 1D).  Note that relative cells with dimension
 *   cellDim=spatialDim-n, where n >=2, can have more than two co-facets in
 *   general.
 *
 *   <li><b>Co-facet Dimension (cofacetDim)</b>: The cell dimension of a
 *   co-facet cell.
 *
 *   <li><b>Maximal co-facet</b>: A co-facet of a relative cell whose
 *   dimension is equal to the spatial dimension of the mesh.
 *
 *   <em>Assertion</em>: <tt>maxCofacetDim == spatialDim</tt>.
 *
 *   <li><b>Vertex</b>: A zero dimensional cell (cellDim=0).
 *
 *   <em>Assertion</em>: <tt>vertexCellDim == 0</tt>.
 *
 *   <li><b>Element</b>: The term "element" is used to refer to cells whose
 *   dimension is equal to the spatial dimension of the mesh.  Warning, the
 *   term "element" is an overloaded term and is also used to refer to an
 *   entry in an array of objects.
 *
 *   <li><b>Facet orientation (facetOrientation)</b>: The orientation of a
 *   facet with respect to its parent cell.  This is given as an integer value
 *   that requires agreement of interpretation.  In 2D, an edge for instance
 *   will have just a positive and a negative orientation.  In 3D, a face can
 *   have more than two orientations and each orientation is a permutation of
 *   the different nodal listings on the face.  The facet orientation integer
 *   value is typically interpreted as an index into an array of possible
 *   orientations.  The mesh reader/creator and the assembly code have to
 *   agree on how these orientations are interpreted.
 *
 *   <li><b>Co-facet index (cofacetIndex)</b>: The relative index of a
 *   co-facet of some relative cell.  For example, an internal face in a 3D
 *   mesh will have two co-facet element objects which can be accessed using
 *   the co-facet indexes 0 and 1.
 *
 *   <li><b>Label</b>: Every cell can be given a single integer value that
 *   can be used in a variety of ways.  For example, boundary faces in 3D can
 *   be given an integer label to specify different types of boundary
 *   conditions.  This integer could also be used with bit manipulation to
 *   allow for multiple classifications for a single cell.
 *
 *   <li><b>Intermediate Cell</b>: An intermediate cell is a cell that has a
 *   dimension cellDim in the range 0 < cellDim < spatialDim.  For example, in
 *   a 3D mesh, faces and edges are intermediate cell types.
 *
 *   </ul>
 *
 * <li><b>Geometric information</b>
 *
 *   <ul>
 *
 *   <li><b>Point</b>: A point is an array of coordinate positions.  Currently
 *   this is stored as a array with three doubles (i.e. x=0, y=1, and z=2).
 *   For 2D meshes only the x=0 and y=1 positions are used.  For 1D meshes,
 *   only the x=0 position is used.
 *
 *   <li><b>Cell Centroid</b>: A point that represents the geometric center of
 *   a cell.
 *
 *   <li> <b>Cell Diameter</b>: The cell diameter is defined as the maximum
 *   line segment that can be drawn in the cell.  This gives a measure of how
 *   big the cell is and, for instance, can be used in some types of
 *   finite-element stabilization methods.
 *
 *   <li><b>Cell Curvature</b>: The curvature of cell refers to the shape of a
 *   cell.  This is typically used to represent the shapes of faces in 3D or
 *   edges in 2D which is critical in the implementation of higher-order
 *   higher-order discrimination methods.  Cell curvature is currently not
 *   supported by this mesh interface but could be in the future.
 *
 *   <li><b>Reference Cell</b>: The term "reference cell" is used below to
 *   refer that a cell type's reference cell which is a geometric description
 *   of the cell in a simple coordinate system.  By using a standard
 *   definition of a "reference cell", we can decouple different types of
 *   mathematical operations defined of a cell of a certain basic type (e.g.
 *   triangles and quadrilaterals in 2D, and tetrahedral and bricks in 3D).
 *   Typically, the nodes of a reference cell will have coordinates that
 *   integer multiples of one (including zero).
 *
 *   <li><b>Cell Jacobian</b>: The cell Jacobian is a transformation matrix
 *   that maps cell vertex coordinate values (i.e. points) from the reference
 *   cell to the physical cell.  The inverse of the cell Jacobian maps points
 *   from the physical cell to the reference cell.  A cell with the cell
 *   dimension cellDim has a cell Jacobian matrix of dimension cellDim x
 *   cellDim.
 *
 *   </ul>
 *
 * <li><b>Parallel Information</b>
 *
 *   <ul>
 *   
 *   <li><b>Cell Ownership</b>: The ownership of a cell refers to which
 *   process owns the cell.  All other processes that need to access that cell
 *   get a "ghost" copy of the cell locally.
 *
 *   <li><b>Ghost Cell</b>: A ghost cell is a cell that is not owned by the
 *   local process that is copied to a process in order to perform certain
 *   computations.
 *
 *   <li><b>Cell Local ID (cellLID)</b>: The LID of a cell is an index to
 *   the cell within the local process and these IDs are typically ordered
 *   from 0 to numLocalCells-1.  Local cells can include both owned cells and
 *   ghosted cells.
 *
 *   <li><b>Cell Global ID (cellGID)</b>: The global ID of a cell is an index
 *   that is unique across all processes and is used as a process-independent
 *   identifier for the cell.
 *
 *   </ul>
 *
 * </ul>
 *
 * \section Sundance_MeshBase_Common_Function_Arguments_sec Common Function Arguments and Pre/Post-Conditions
 *
 * Below, a common set of Mesh arguments are described along with relevant
 * pre- and post conditions.
 *
 * <ul>
 *
 * <li><b>cellDim</b> [in] The dimension of a  cell
 *
 * <em>Precondition</em>: <tt>0 <= cellDim <= spatialDim</tt>
 *
 * <li><b>cellLID(s)</b> [in] Local ID (LID) of a cell(s).
 *
 * <em>Precondition</em>: <tt>0 <= cellLID < numCells(cellDim)</tt>, where
 * <tt>cellDim</tt> is the dimension of the given cell.
 *
 * <li><b>cellGID</b> [in] Global ID (GID) of a cell which is unique across
 * processes.
 *
 * <em>Precondition</em>: <tt>??? <= cellGID < ???</tt> (How do we write
 * this?).
 *
 * <li><b>facetDim</b> [in] The dimension of a given facet cell.
 *
 * <em>Precondition</em>: <tt>0 <= facetDim < cellDim</tt>, where
 * <tt>cellDim</tt> is the dimension of the relative cell
 *
 * <li><b>facetIndex</b> [in] Local index of the facet with respect to its
 * parent relative cell.
 *
 * <em>Precondition</em>: <tt>0 <= facetIndex <
 * numFacets(cellDim,cellLID,facetDim)</tt>, where <tt>facetDim</tt> is the
 * dimension of the facet of the relative parent cell
 * <tt>(cellDim,cellLID)</tt>
 *
 * <li><b>facetOrientation(s)</b> [out] The orientation of a facet with
 * respect to its parent cell.  This is given as an integer value that
 * requires agreement of interpretation (see above).
 *
 * <li><b>cofacetDim</b> [in] Cell dimension of a co-facet cell.
 *
 * <em>Precondition</em>: <tt>0 <= cellDim < cofacetDim <= spatialDim</tt>,
 * where <tt>cellDim</tt> is the dimension of the relative cell
 *
 * <li><b>cofacetIndex</b> [in] Relative index of a co-facet of some relative
 * cell.
 *
 * <em>Precondition</em>: <tt>0 <= cofacetIndex <
 * numMaxCofacets(cellDim,cellLID)</tt>, where <tt>(cellDim,cellLID)</tt> is
 * the relative cell.
 *
 * </ul>
 *
 * \section Sundance_MeshBase_Refactorings_sec Possible Refactorings
 *
 * <ul>
 *
 * <li>Copy ObjectWithInstanceID base class into a special Nihilo tools
 * package or into Teuchos.
 *
 * <li>Inherit MeshBase from Teuchos::VerboseObject instead of
 * ObjectWithVerbosity.
 *
 * <li>Remove any concrete implementations and data from this base interface
 * to provide a clean abstract interface.  This would include removing the
 * constructor.
 *
 * <li>Either remove the non-batch access functions all together, or make the
 * namings of batch and non-batch access functions more consistent.
 *
 * <li>Have all of the functions return just a single return value or when
 * more than one object is returned, return all values from the argument list.
 * Or, when the value that is being returned from the argument list is only
 * occasionally returned, change it to a pointer type and give it a default
 * value of NULL so that it can be safely ignored!  This will actually speed
 * up performance in some cases.  For example, <tt>facetOrientation</tt> is a
 * return argument that is often ignored and would benefit from making it an
 * optional argument.
 *
 * <li>Remove the staggered output function since the Teuchos::FancyOStream
 * should handle parallel output from every process well.
 *
 * <li>Add a new cell type called a "general" cell for very irregular cells
 * that defy fixed classification.  This could be used to represent many of
 * the irregular cells that remain after non-conforming mesh refinement.
 *
 * <li>Make sure to tell the DOF map how to interpret the facet orientation
 * integers.  Currently the interpretation is hard coded in Sundance.  These
 * permutations need to be enumerated in a C++ enum and the DOF map must be
 * given the array of these enums to interpret the orientation integer values.
 * To be 100 percent safe, we need a way of allowing the mesh to directly
 * return the interpretation of these orientations.  What you can do is to
 * require that the int value be an index into an an array of enumeration
 * values that corresponds to the cell type.  This would require an enum type
 * and an array of enum values for every facet cell type.  For the current
 * cell types for regular elements, this would be no problem.  However, for
 * the facets of a "general" cell type, the facets may also be classified as
 * "general" cell types and therefore may not allow such an orientation array
 * to be defined.
 *
 * <li>Clearly define the coordinate system of the reference cell of each of
 * the standard cell types in order to clearly define interpretation of the
 * reference points passed into the <tt>pushForward()</tt> function.  Without
 * a clear definition of the reference cell, these reference points are
 * meaningless.  Actually, for maximal flexibility, each cell type should
 * define, at runtime, the coordinate system for the reference cell.  Or, we
 * must just all agree on a single set of standard definitions of reference cells
 * for different cell types.
 *
 * </ul>
 */
class MeshBase
  : public ObjectWithClassVerbosity<MeshBase>
  , public ObjectWithInstanceID<MeshBase>
{
public:

  /** \name Constructors/Destructor */
  //@{

  /** \brief . */
  MeshBase(int dim, const MPIComm& comm, 
    const MeshEntityOrder& meshOrder);
  
  /** \brief . */
  virtual ~MeshBase(){;}

  //@}

  /** \name Topological Information */
  //@{

  /** \brief Get the ordering convention used by this mesh */
  const MeshEntityOrder& meshOrder() const {return order_;}
  
  /** \brief Get the spatial dimension of the mesh
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>0 < returnVal <= 3</tt>
   * </ul>
   */
  int spatialDim() const {return dim_;}
  
  /** \brief Get the number of local cells having dimension cellDim.
   *
   * \todo Change this to <tt>numLocalCells(cellDim)</tt>.
   */
  virtual int numCells(int cellDim) const = 0 ;
  
  /** \brief Return the number of facets of the given relative cell.
   *
   * \param  cellDim
   *           [in] The dimension of the relative cell
   * \param  cellLID
   *           [in] The LID of the relative cell
   * \param  facetDim
   *           [in] The dimension of the facets for the
   *           relative cell
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>returnVal >= 2</tt>: Every cell has two or more
   *     facets of a given dimension.
   * </ul>
   */
  virtual int numFacets(int cellDim, int cellLID, 
                        int facetDim) const = 0 ;
  
  /** \brief Return the LID of a single facet cell (and optionally its
   * orientation) with respect to a relative parent cell.
   *
   * \param  cellDim
   *           [in] Dimension of the parent cell whose facets
   *           are being obtained.
   * \param  cellLID
   *           [in] Local ID of the parent cell whose facets
   *           are being obtained
   * \param  facetDim
   *           [in] Dimension of the desired facet.
   * \param  facetIndex
   *           [in] The relative index into the list of the
   *           cell's facets.
   * \param  facetOrientation
   *           [out] The orientation of the facet w.r.t the
   *           parent cell.
   *
   * \todo Change the facetOrientation argument to be a raw pointer that can
   * be NULL and therefore easily ignored.
   */
  virtual int facetLID(int cellDim, int cellLID,
                       int facetDim, int facetIndex,
                       int& facetOrientation) const = 0 ;
  
  /** \brief Return an array containing the LIDs of facets of dimension
   * facetDim for the given batch of relative parent cells.
   *
   * \param  cellDim
   *           [in] Dimension of the relative parent cells whose facets
   *           are being obtained
   * \param  cellLIDs
   *           [in] Array of LIDs for the relative parent cells
   * \param  facetDim
   *           [in] Dimension of the desired facets
   * \param  facetLIDs
   *           [out] On exit this array gives the local facet IDs
   *           for all of the cells given in cellLIDs in one flat array
   *           with size = <tt>cellLIDs.size()*nf</tt>, where
   *           <tt>nf=numFacets(cellDim,cellLIDs[j],facetDim)</tt>) where
   *           <tt>j</tt> can be any index <tt>0 <= j < numCells(cellDim)</tt>.
   *           Specifically, the local facet IDs for <tt>cellLIDs[k]</tt>, where
   *           <tt>k=0...cellLIDs.size()-1</tt>, is given
   *           by <tt>facetLIDs[k*nf+j]</tt>, where <tt>j=0...nf-1</tt>.
   * \param  facetOrientations
   *           [out] On exist, if <tt>facetDim > 0</tt>, this array gives
   *           the integer orientation of the facet with respect to
   *           its parent cell (see above definition of "Facet Orientation").
   *           Oon exist this array will be
   *           resized to size = <tt>cellLIDs.size()*nf</tt>, where
   *           <tt>nf=numFacets(cellDim,cellLID,facetDim)</tt>).
   *           Specifically, the local facet orientation for the cell cellLIDs[k], where
   *           <tt>k=0...cellLIDs.size()-1</tt>, is given
   *           by <tt>facetOrientations[k*nf+j]</tt>, where <tt>j=0...nf-1</tt>.
   *           If <tt>facetDim==0</tt> then this array is ignored.
   *
   * <b>Warning!</b> This function will only work for homogeneous cell types!
   *
   * \todo Change the facetOrientation argument to an <tt>Array<int>*</tt>
   * type so that it can be NULL and therefore ignored (which is a common use
   * case).
   */
  virtual void getFacetLIDs(int cellDim, 
                            const Array<int>& cellLIDs,
                            int facetDim,
                            Array<int>& facetLIDs,
                            Array<int>& facetOrientations) const = 0 ;
  
  /** \brief Return an array containing the LIDs of the facets of
   * dimension facetDim of the single relative parent cell (optionally also
   * returning the facet orientations).
   *
   * \param  cellDim
   *           [in] Dimension of the parent cell whose facets
   *           are being obtained.
   * \param  cellLID
   *           [in] Local ID of the parent cell whose facets
   *           are being obtained
   * \param  facetDim
   *           [in] Dimension of the desired facets
   * \param  facetLIDs
   *           [out] On exit this array gives the local facet IDs
   *           for the parent cell
   *           with size = <tt>nf</tt>, where
   *           <tt>nf=numFacets(cellDim,cellLID,facetDim)</tt>.
   * \param  facetOrientations
   *           [out] On exist, if <tt>facetDim > 0</tt>, this array gives
   *           the integer orientation of the facet with respect to
   *           its parent cell (see above definition of "Facet Orientation").
   *           On exist this array will be
   *           resized to size = <tt>nf</tt>, where
   *           <tt>nf=numFacets(cellDim,cellLID,facetDim)</tt>).
   *           If <tt>facetDim==0</tt> this this array argument is ignored!
   *
   * The default implementation loops over calls to
   * <tt>facetLID()</tt>. Subclasses can provide a more efficient
   * implementation if desired.
   *
   * \todo Rename to something like getSingleCellFacetLIDs(...).
   *
   * \todo Change the facetOrientation argument to an <tt>Array<int>*</tt>
   * type so that it can be NULL and therefore ignored (which is a common use
   * case).
   */
  void getFacetArray(int cellDim, int cellLID, int facetDim, 
                     Array<int>& facetLIDs,
                     Array<int>& facetOrientations) const ;
  
  /** \brief Return a view of an array of LIDs for a maximum-dimensional
   * cell's zero-dimensional facets (i.e. vertexes).
   *
   * \param  maxCellLID
   *           [in] Local ID of the maximum dimensional cell (i.e. element).
   *
   * <b>Preconditions:</b><ul>
   * <li>The cell <tt>maxCellLID</tt> must be a maximum-dimensional cell!
   * </ul>
   *
   * \returns A raw pointer into an array which stores the local process IDs
   * of the vertices.  These IDs are the local process IDs, not the facet
   * indexes relative to the relative parent cell.  The dimension of this
   * array is determined by the cell type given by
   * <tt>cellType(maxCellLID)</tt>.
   *
   * \todo Return this array as a Teuchos::ArrayView<int> object which will
   * have the dimension embedded in it and will have full range checking!
   *
   * \todo Rename to something like getMaxCellZeroFacetsLIDsView(...).
   */
  virtual const int* elemZeroFacetView(int maxCellLID) const ;

  /** \brief Return the number of maximal co-facets of the given cell.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>0 <= cellDim < spatialDim()</tt>
   * </ul>
   *                                            
   * Note that if <tt>cellDim==spatialDim()-1</tt> and <tt>returnVal==1</tt> then
   * this cell must be a boundary cell (i.e. a boundary face in 3D, a boundary
   * edge in 1D, or a boundary node in 1D)!
   */
  virtual int numMaxCofacets(int cellDim, int cellLID) const = 0;

  /** \brief Return the LID of a maximal co-facet of a cell.
   *
   * \param  cellDim
   *           [in] Dimension of the cell whose co-facets are being obtained
   * \param  cellLID
   *           [in] Local ID of the cell whose co-facets are being obtained
   * \param  cofacetIndex
   *           [in] Index into the list of the cell's co-facets.
   * \param  facetIndex
   *           [out] Returns the local index into the facet array
   *           for the relative cell (cellDim,cellLID) with respect
   *           to the maximal co-facet.  In other words
   *           <tt>cellLID==facetLID(spatialDim(),returnVal,cellDim,facetIndex)</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>0 <= cellDim < spatialDim()</tt>
   * </ul>
   *
   * \returns The LID of the requested maximal co-facet.
   *
   * \todo Make the <tt>facetIndex</tt> an <tt>int*</tt> argument and give it a
   * default value of <tt>NUL</tt> so that it can be easily ignored!
   */
  virtual int maxCofacetLID(int cellDim, int cellLID,
                            int cofacetIndex, 
                            int& facetIndex) const = 0 ;
    /** 
     * Get the LIDs of the maximal cofacets for a batch of codimension-one
     * cells. The default implementation simply loops over the cells in the
     * cellLID array, taking no advantage of any internal data structures.
     *
     * \param cellLIDs [in] array of LIDs of the cells whose cofacets are 
     * being obtained
     * \param cofacets [out] 
     * \param facetIndex [out] index of each calling cell
     * into the list of its maximal cofacet's facets 
     */
  virtual void getMaxCofacetLIDs(const Array<int>& cellLIDs,
    MaximalCofacetBatch& cofacets) const ;

  /** \brief Return an array of the LIDs of all of the co-facets for a
   * given relative cell.
   *
   * \param  cellDim
   *           [in] Dimension of the relative cell whose co-facets are being obtained
   * \param  cellLID
   *           [in] Local index of the relative cell whose co-facets are being obtained
   * \param  cofacetDim
   *           [in] Dimension of the co-facets to get
   * \param  cofacetLIDs
   *           [out] Array containing the LIDs of the co-facets
   *           for the given relative cell (cellDim,cellLID).
   *
   * \todo Change name to <tt>getCofacetArray()</tt> to be consistent with
   * <tt>getFacetArray()</tt>!
   */
  virtual void getCofacets(int cellDim, int cellLID,
                           int cofacetDim, Array<int>& cofacetLIDs) const = 0 ;

  /** \brief Get the cell type of the given cell dimension.
   *
   * Note: This function by its very definition assumes that all cells of a
   * given dimension have the same cell type within a mesh!
   *
   * \todo This function must be changed in order to deal with mixed cell
   * types with the same cellDim!
   */
  virtual CellType cellType(int cellDim) const = 0 ;

  /** \brief Get the label of the given cell.
   *
   * \todo Change to getLabel(...)?
   */
  virtual int label(int cellDim, int cellLID) const = 0 ;

  /** \brief Get the labels for a batch of cells.
   *
   * \param  cellDim
   *           [in] Dimension of the parent cell whose facets
   *           are being obtained
   * \param  cellLIDs
   *           [in] Array of cell LIDs
   * \param  labels
   *           [out] On exit, contains an array (<tt>size=cellLIDs.size()</tt>)
   *           of the labels for each of the given cells.
   */
  void getLabels(int cellDim, const Array<int>& cellLIDs, 
    Array<int>& labels) const ;

  /** \brief Set the label for the given cell.
   *
   * \todo Move this out of this base interface and into a mesh loading
   * interface?
   */
  virtual void setLabel(int cellDim, int cellLID, int label) = 0 ;

  /** Get the list of all labels defined for cells of the given dimension */
  virtual Set<int> getAllLabelsForDimension(int cellDim) const ;

  /** 
   * Get the cells associated with a specified label. The array 
   * cellLID will be filled with those cells of dimension cellDim
   * having the given label.
   */
  virtual void getLIDsForLabel(int cellDim, int label, Array<int>& cellLIDs) const  ;
  //@}

  /** \name Geometric Information */
  //@{
  
  /** \brief Return the position of a local vertex.
   *
   * \param  vertexLID
   *           [in] The LID of the vertex.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>0 <= vertexLID < this->numCells(0)</tt>
   * </ul>
   *
   * \todo Change the name of this function to
   * <tt>getVertexPosition(...)</tt>.
   */
  virtual Point nodePosition(int vertexLID) const = 0 ;
  
  /** \brief Return a const view into an raw array for the position of a local
   * vertex.
   *
   * \param  vertexLID [in] The LID of the vertex
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>0 <= vertexLID < this->numCells(0)</tt>
   * </ul>
   *
   * \returns a raw pointer into an array of doubles of length
   * <tt>spatialDim</tt>.
   *
   * \todo Return this array as a Teuchos::ArrayView<double> object which will
   * have the dimension embedded in it and will have full range checking!
   *
   * \todo Change this function name to <tt>getVertexPositionView()</tt>.
   */
  virtual const double* nodePositionView(int vertexLID) const = 0 ;
  
  /** \brief Return the centroid of a cell.
   *
   * The default implementation just averages the positions of the
   * zero-dimensional facets (i.e. vertexes).
   */
  virtual Point centroid(int cellDim, int cellLID) const ;
  
  /** \brief Compute (or get) the Jacobians for a batch of cells.
   *
   * \param  cellDim
   *           [in] Dimension of the cells whose Jacobians are to be computed
   * \param  cellLID
   *           [in] Local IDs of the cells for which Jacobians are to be computed
   * \param  jBatch
   *           [out] Batch of cell Jacobians.
   *
   * Warning! The default implementation returns an empty batch of cell
   * Jacobians!
   *
   * \todo Add a query function to tell if this feature is supported and then
   * add a precondition based on this query function!
   */

  virtual void getJacobians(int cellDim, const Array<int>& cellLID,
    CellJacobianBatch& jBatch) const { ; }

  //bvbw  virtual void getJacobians(
  //  int cellDim, const Array<int>& cellLID,
  //  CellJacobianBatch& jBatch
  //  ) const
  //  { jBatch.resize(0,0,0); }
  
  /** \brief Compute the diameters of a batch of cells.
   *
   * \param  cellDim
   *           [in] Dimension of the cells whose diameters are to be computed
   * \param  cellLIDs
   *           [in] Local IDs of the cells for which diameters are to be computed
   * \param  diameters
   *           [out] Array (<tt>size = cellLIDs.size()</tt>) of cell diameters.
   *
   * Warning! The default implementation returns an empty array of cell
   * diameters!
   * ToDo: Change the default implementation to compute diameters based on
   * calls to the node position accessors. Going through the Mesh interface in
   * that way will be less efficient than a low-level implementation, but
   * would be a reasonable intermediate step for mesh developers. 
   * - KL 5 Aug 2007. 
   */
  virtual void getCellDiameters(
    int cellDim, const Array<int>& cellLIDs,
    Array<double>& diameters
    ) const
    { diameters.resize(0); }


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
    virtual void outwardNormals(
      const Array<int>& cellLIDs,
      Array<Point>& outwardNormals
      ) const ;

  /** 
   * Get tangent vectors for a batch of edges
   *
   * \param cellLIDs [in] LIDs for the cells whose tangents are to be
   * computed. 
   * \param tangentVectors [out] Unit tangents for each cell
   */
  virtual void tangentsToEdges(
    const Array<int>& cellLIDs,
    Array<Point>& tangenVectors
    ) const ;
      

  /** \brief Map points from a reference cell to physical points for a batch
   * of cells.
   *
   * \param  cellDim
   *           [in] Dimension of the cells
   * \param  cellLIDs
   *           [in] Local IDs of a batch of cells
   * \param  refPts
   *           [in] Array of points on the single reference cell with respect
   *           to the reference cell's coordinate system.  Note that the
   *           interpretation of these reference points is strictly determined
   *           by the coordinate system of the cell type
   *           <tt>cellType(cellDim)</tt> and must be clearly defined by this
   *           interface.
   * \param  physPts
   *           [out] Array (<tt>size = cellLIDs.size()*refPts.size()</tt>) of
   *           the physical points given in a flat array for the batch of
   *           cells.  Specifically, the physical points for each cell
   *           <tt>cellLIDs[k]</tt>, for <tt>k=0...cellLIDs.size()-1</tt>, is
   *           given by <tt>physPts[k*nrp+j]</tt>, for <tt>j=0...nrp-1</tt>
   *           and <tt>nrp=refPts.size()</tt>.
   *
   * Warning! The default implementation returns an empty array of physical
   * points!
   *
   * \todo Add a query function to tell if this feature is supported and then
   * add a precondition based on this query function!
   */
  virtual void pushForward(
    int cellDim, const Array<int>& cellLIDs,
    const Array<Point>& refPts,
    Array<Point>& physPts
    ) const
    { physPts.resize(0); }

  //@}

  /** \name Parallel Information */
  //@{

  /** \brief Return the MPI communicator over which this mesh is distributed. */
  const MPIComm& comm() const {return comm_;}
  
  /** \brief Return the rank of the processor that owns the given cell.
   *
   * If <tt>returnVal==comm().getRank()</tt> then this cell is owned by this
   * process.
   */
  virtual int ownerProcID(int cellDim, int cellLID) const = 0 ;
    
  /** \brief Determine whether a given cell GID exists on this processor. */
  virtual bool hasGID(int cellDim, int cellGID) const = 0 ;

  /** \brief Find the LID of a cell given its GID.
   *
   * \param  cellDim
   *           [in] Dimension of the cell
   * \param  cellGID
   *           [in] Global ID of the cell
   *
   * <b>Preconditions:</b><ul>
   * <li>hasGID(cellDim,cellGID)==true</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li>0 <= returnVal < numCells(cellDim)</tt>
   * </ul>
   */
  virtual int mapGIDToLID(int cellDim, int cellGID) const = 0 ;

  /** \brief Find the global ID of a cell given its LID.
   *
   * \param  cellDim
   *           [in] Dimension of the cell
   * \param  cellLID
   *           [in] Local ID of the cell
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>returnVal >= 0 </tt>
   * </ul>
   */
  virtual int mapLIDToGID(int cellDim, int cellLID) const = 0 ;

  /** \brief Return if cells of dimension cellDim have been assigned global
   * IDs or not.
   *
   * \param  cellDim
   *           [in] Dimension of the cell
   */
  virtual bool hasIntermediateGIDs(int cellDim) const = 0 ;

  /** \brief Assign global IDs to cells of dimension cellDim.
   *
   * \param  cellDim
   *           [in] Dimension of the cell
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>hasIntermediateGIDs(cellDim)==true</tt>
   * </ul>
   */
  virtual void assignIntermediateCellGIDs(int cellDim) = 0 ;

  //@}

  /** \name Reordering */
  //@{

  /** \brief Set the reordering strategy to be used with this mesh. */
  void setReorderer(const CellReorderer& reorderer) 
    {reorderer_ = reorderer.createInstance(this);}

  /** \brief Get the reordering strategy to be used with this mesh. */
  const CellReordererImplemBase* reorderer() const 
    {return reorderer_.get();}

  //@}



  /** \name Functions for Mesh with hanging nodes */
    //@{
    /** Function returns true if the mesh allows hanging nodes (by refinement),
     * false otherwise */
    virtual bool allowsHangingHodes() const { return false; }

    /** Function returns true if the specified element is a "hanging" element
     * false otherwise <br>
     * @param cellDim [in] should be between 0 , D-1
     * @param cellLID [in] the local ID of the element */
    virtual bool isElementHangingNode(int cellDim , int cellLID) const { return false; }

    /** Returns the index in the parent maxdim Cell of the refinement tree
     * @param maxCellLID [in] the LID of the cell */
    virtual int indexInParent(int maxCellLID) const { return 0; }

     /** How many children has a refined element. <br>
      * This function provides information of either we have bi or trisection */
     virtual int maxChildren() const { return 0;}

    /** Function returns the facets of the parent cell (needed for HN treatment) <br>
     * @param childCellLID [in] the LID of the maxdim cell, whos parents facets we want
     * @param dimFacets [in] the dimension of the facets which we want to have
     * @param facetsLIDs [out] the LID of the parents facets (all) in the defined order
     * @param parentCellLIDs [out] the maxdim parent cell LID */
    virtual void returnParentFacets( int childCellLID , int dimFacets ,
    		                         Array<int> &facetsLIDs , int &parentCellLIDs ) const { }
  //@}



  /** \name Store special weights in the mesh (for Adaptive Cell Integration) */
    //@{

    /** returns the status of the special weights if they are valid <br>
     *  These weights are usually computed for one setting of the curve (Adaptive Cell Integration)*/
    virtual bool IsSpecialWeightValid() const {return validWeights_;}

    /** specifies if the special weights are valid <br>
     *  if this is false then usually the special weights have to be recomputed */
    virtual void setSpecialWeightValid(bool val) const { validWeights_ = val;}

    /** deletes all special weights so those have to be recreated*/
    virtual void flushSpecialWeights() const;

    /** verifies if the specified cell with the given dimension has special weights */
    virtual bool hasSpecialWeight(int dim, int cellLID) const;

    /** Sets the special weights */
    virtual void setSpecialWeight(int dim, int cellLID, Array<double>& w) const;

    /** Returns the special weights */
    virtual void getSpecialWeight(int dim, int cellLID, Array<double>& w) const;
    //@}


    /** \name Store the intersection/quadrature points for the curve/surf integral <br>
     *  for a curve or surf integral we need some quadrature points along the curve in one curve <br>
     *  These */
      //@{

      /** */
      virtual bool IsCurvePointsValid() const {return curvePoints_Are_Valid_;}

      /**  */
      virtual void setCurvePointsValid(bool val) const { curvePoints_Are_Valid_ = val;}

      /** detletes all the points and its normals which have been stored */
      virtual void flushCurvePoints() const;

      /** verifies if the specified maxCell has already precalculated quadrature point for one curve */
      virtual bool hasCurvePoints(int maxCellLID , int curveID) const;

      /** Sets the points, curve derivatives and curve normals for one maxCell needed for curve/surf integral*/
      virtual void setCurvePoints(int maxCellLID, int curveID , Array<Point>& points , Array<Point>& derivs , Array<Point>& normals) const;

      /** Gets the points, curve derivatives and curve normals for one maxCell needed for curve/surf integral*/
      virtual void getCurvePoints(int maxCellLID,  int curveID , Array<Point>& points , Array<Point>& derivs , Array<Point>& normals) const;

private:
      virtual int mapCurveID_to_Index(int curveID) const;
public:
      //@}

  /** \name Output */
  //@{

  /** \brief Set whether to stagger output in parallel. Set to true for readable
   * output in parallel debugging sessions.
   *
   * This should be normally, as it causes one synchronization point per
   * process.
   *
   * \todo Get rid of this once we update to use Teuchos::FancyOStream since
   * parallel outputting will be done in a readable way automatically!
   */
  static bool& staggerOutput() {static bool rtn=false; return rtn;}

  //@}
    
private:

  int dim_;

  MPIComm comm_;

  MeshEntityOrder order_;

  RCP<CellReordererImplemBase> reorderer_;



  /** flag to indicate if the weights stored are valid */
  mutable bool validWeights_;

  /** Object to store the special weights for integration */
  mutable Array < Sundance::Map< int , Array<double> > > specialWeights_;



  /** true if the curve did not moved, false if those points are not reusable*/
  mutable bool curvePoints_Are_Valid_;

  /** how many curves are participating in curve integrals*/
  mutable int nrCurvesForIntegral_;

  /** store intersection informations to calculate the curve integral*/
  mutable Array < Sundance::Map< int , Array<Point> > > curvePoints_;

  /** store directional derivative informations to calculate the curve integral*/
  mutable Array < Sundance::Map< int , Array<Point> > > curveDerivative_;

  /** store normal directional used in the curve or in the surf integral <br>
   * in case of the surf integral it is the cross product from the integral */
  mutable Array < Sundance::Map< int , Array<Point> > > curveNormal_;

  /** map the curve ID to index in the previous arrays */
  mutable Sundance::Map< int , int > curveID_to_ArrayIndex_;
};


} // namespace Sundance

#endif
