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

#ifndef SUNDANCE_DOFMAPBASE_H
#define SUNDANCE_DOFMAPBASE_H

#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceCellSet.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceMapStructure.hpp"
#include "SundanceObjectWithVerbosity.hpp"

namespace Teuchos {class Time;}

namespace Sundance
{
using namespace Teuchos;

/** \brief Base interface for implementations of a degree of freedom map.
 *
 * A degree of freedom (DOF) map is a parallel-aware object that takes
 * takes DOFs on individual cells in the whole mesh, and creates global
 * IDs for them on the whole mesh across processes.
 *
 * A DOF map is constructed out of a mesh and assignment of various
 * discrete functions with various basis-es assigned to the different
 * cells in the mesh;
 *
 * This interface assumes that the global DOFs owned in this process are
 * ordered sequentially so the DOFs owned in this process are given by
 * <tt>this->lowestLocalDOF() + k</tt>, for <tt>k =
 * 0...this->numLocalDOFs()-1()</tt>, where <tt>this->numLocalDOFs()</tt>
 * is the number of DOFs owned by this process.  Therefore, any DOF with
 * value less than <tt>this->numLocalDOFs()</tt> and greater than or equal
 * to <tt>this->lowestLocalDOF()+this->numLocalDOFs()</tt> are necessarily
 * ghosted DOFs.  ??? ToDo: I don't think the above is correct! ???
 * 
 * ToDo: Finish documentation!
 *
 * \section Sundance_DOFMapBase_Defintions_sec Definitions
 *
 * <ul>
 *
 * <li><b>Degree of Freedom (DOF)</b>: ???
 *
 * <li><b>Homogeneous DOF Map</b>: ???
 *
 * </ul>
 *
 * \section Sundance_DOFMapBase_ToDo_sec Possible Refactorings
 *
 * <ul>
 *
 * <li>Inherit this base class from Teuchos::Describable and remove the
 * print() function, replacing it with the override to
 * Teuchos::Describable::describe()?
 *
 * <li>Break off the default implementation in this class into another
 * subclass (called something like <tt>DOFMapDefaultBase</tt>) and then
 * make this interface a true abstract interface?  There are lots of
 * advantages to having pure interface classes (e.g. less documenation,
 * standard delegation subclasses etc.).
 *
 * <li>Add a public function to return the <tt>MeshBase</tt> object?  Is
 * there any reason not to give the client access to the <tt>MeshBase</tt>
 * object?  If you do this, then you can remove the <tt>isRemote()</tt>
 * function and let the client just query the map object directly.
 *
 * <li>Refactor this interface and all objects accessed to only use
 * absrract interfaces and not handle classes?  This would involve the
 * same principle as is used in Thyra.  Is this workable?
 *
 * <li>Add some function somewhere to return the total number of functions
 * that is defined on the mesh.  This is needed to write precise
 * preconditions and postconditions for many of the functions.  For
 * example, a function <tt>this->numTotalFunctions()</tt> would be very
 * helpful in this regard.
 *
 * <li>???
 *
 * </ul>
 */
class DOFMapBase : public Playa::Printable
{
public:

  /** \brief . */
  DOFMapBase(const Mesh& mesh, int setupVerb);
      
  /** \brief .
   *
   * ToDo: Remove this virtual destructor since this interface already
   * inherits from a base interface that has a virtual destructor.
   */
  virtual ~DOFMapBase(){;}

  /** \brief . */
  const Mesh& mesh() const {return mesh_;}
  
  /** \brief Return <tt>true</tt> if the given cell is really owned by
   * another process and is only ghosted in this process (and optionally
   * return the owning process ID).
   *
   * \param  cellDim
   *           [in] The dimension of the cell.  
   * \param  cellLID
   *           [in] The LID of the cell in this process.  See
   *           <tt>MeshBase</tt> for more details, preconditions, etc.
   * \param  ownerProcID
   *           [out] The process rank ID which owns the cell
   *           <tt>(cellDim,cellLID)</tt>.
   *
   * <b>Preconditions:</b> See <tt>MeshBase</tt>
   *
   * ToDo: Change the <tt>ownerProcID</tt> argument to <tt>int*</tt> and
   * given it a default value of <tt>NULL</tt> so that it can be ignored
   * by the client.
   *
   * ToDo: Consider removing this function and just calling on the mesh
   * object directly.
   *
   * ToDo: Change name to <tt>isRemoteCell()</tt>? 
   */
  bool isRemote(int cellDim, int cellLID, int& ownerProcID) const 
    {return (ownerProcID=mesh_.ownerProcID(cellDim, cellLID)) != localProcID_;}
      
  /** \brief Get the global DOFs for a single function on a single
   * cell.
   *
   * \param  cellDim
   *           [in] The dimension of the cell
   * \param  cellLID
   *           [in] Local ID (LID) of the cell
   * \param  funcID
   *           [in] Function ID for which DOFs are requested
   * \param  dofs
   *           [out] Global IDs for DOFs of the requested function
   *           on the requested cell.
   *
   * 
   *
   */
  virtual void getDOFsForCell(int cellDim, int cellLID,
    int funcID,
    Array<int>& dofs) const;

  /** \brief Return the global DOFs for a batch of cells for a given set
   * of functions on those cells.
   *
   * \param  cellDim
   *           [in] The dimension of the cells in the batch of cells
   * \param  cellLIDs
   *           [in] Local IDs (LIDs) for the batch of cells
   * \param  requestedFuncSet
   *           [in] Set of function IDs for which DOFs are requested.
   *           Note that this must be equal to the allowed set of
   *           requested functions
   *           <tt>this->allowedFuncsOnCellBatch(cellDim,cellLIDs)</tt>.
   * \param  dofs
   *           [out] Global IDs for DOFs of the requested functions on the
   *           batch of cells.  The size of this array on output is \code
   *           dofs.size()==mapStructure.numBasisChunks() \endcode. The
   *           global DOFs for cell <tt>cellLIDs[c]</tt> (where <tt>0 <= c
   *           < cellLIDs.size()</tt>) for the cell for the basis chunk
   *           <tt>basisChunk</tt> (where <tt>0 <= basisChunk <
   *           mapStructure.numBasisChunks()</tt>) are given by \code
   *           dofs[c*nDOFsPerCell[c]+k] \endcode, where \code
   *           k=0...nDOFsPerCell[c]-1 \endcode and \code nDOFsPerCell[c]
   *           = mapStructure.basis(basisChunk).nNodes( spatialDim,
   *           this->mesh()->cellType(cellDim) ) *
   *           mapStructure.numFuncs(basisChunk) \endcode.  Note that the
   *           array <tt>nDOFsPerCell</tt> used above is not actually
   *           computed here or returned by this interface.  It is only
   *           used abstractly to define the <tt>dofs</tt> array above.
   * \param  nNodes
   *           [out] Array giving the number of coefficients for each type
   *           of basis family in <tt>mapStructure</tt> for each function.
   *           The size of this array on return is
   *           <tt>nNodes.size()==mapStructure.numBasisChunks()</tt> and
   *           <tt>nNodes[b]==mapStructure.basis(b).nNodes(spatialDim,cellType)</tt>,
   *           for <tt>b=0...mapStructure.numBasisChunks()-1</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>requestedFuncSet.setDifference(
   *           this->allowedFuncsOnCellBatch(cellDim,cellLIDs) ).size() == 0</tt>
   * <li>???Others???
   * </ul>
   *
   * \returns The map structure that groups sets of different functions
   * according with respect to those that have the same basis family.
   * Above, we will refer to this returned object as <tt>mapStructure</tt>
   * where <tt>mapStructure = *returnVal</tt>.
   *
   * ToDo: Remove the argument requestedFuncSet since all current use
   * cases and implemenations asume that all of the functions are
   * requested and returned.
   *
   * ToDo: Remove the nNodes return argument since this information can be
   * gotten just from the return <tt>mapStructure</tt> object.
   * Specifically, <tt>nNodes[basisChunk] =
   * mapStructure.basis(basisChunk).nNodes(spatialDim,cellType)</tt> so
   * what is the point in returning this array?  Since this is needed for
   * a whole batch of cells, it is cheap to just grab this from the
   * returned <tt>mapStructure</tt> object as shown above.
   */
  virtual RCP<const MapStructure> 
  getDOFsForCellBatch(
    int cellDim, const Array<int>& cellLIDs,
    const Set<int>& requestedFuncSet,
    Array<Array<int> >& dofs,
    Array<int>& nNodes,
    int verb
    ) const = 0 ;

  /** \brief Return the set of a function IDs for a batch of cells for
   * which DOFs can be obtained.
   *
   * \param  cellDim
   *           [in] The dimension of the cells in the batch of cells
   * \param  cellLIDs
   *           [in] Local IDs (LIDs) for the batch of cells
   */
  virtual RCP<const Set<int> >
  allowedFuncsOnCellBatch(int cellDim,
    const Array<int>& cellLIDs) const = 0 ;

  /** \brief Return an array of cell filters that gives, for each function
   * ID, the cells that the function lives on.
   *
   * \returns A reference to an array of size <tt>returnVal.size() ==
   * numTotalFunctcions</tt> where <tt>returnVal[funcID]</tt> gives a
   * handle to a <tt>CellFilterBase</tt> object where the function
   * <tt>funcID</tt> lives, where <tt>0 <= funcID <
   * numTotalFunctions</tt>.
   */
  virtual const Array<CellFilter>& funcDomains() const = 0 ;

  /** \brief Return the lowest DOF for DOFs owned in this process. */
  int lowestLocalDOF() const {return lowestLocalDOF_;}

  /** \brief Returns <tt>true</tt> if the given global DOF is owned in
   * this process.
   */
  bool isLocalDOF(int dof) const 
    {return (lowestLocalDOF_ <= dof && dof < lowestLocalDOF_+numLocalDOFs_);}

  /** \brief Return the number of DOFs owned in this process.
   *
   * ToDo: Is this the number of owned + ghosted DOFs or is it just owned
   * DOFs?
   */
  int numLocalDOFs() const {return numLocalDOFs_;}

  /** \brief Return the global number of DOFs over all processes.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>returnVal >= this->numLocalDOFs()</tt>
   * </ul>
   */
  int numDOFs() const {return numDOFs_;}

  /** \brief Return an array of the global DOFs for the DOFs that are
   * locally accessed in this process.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>returnVal->size() == ( this->numDOFs() - this->numLocalDOFs() )</tt>
   * <li><tt>returnVal[k] < this->lowsetLocalDOF() &&
   *     this->lowsetLocalDOF() + this->numLocalDOFs() <= returnVal[k]</tt>,
   *     for <tt>k=0...(this->numDOFs()-this->numLocalDOFs())-1</tt>.
   * </ul>
   *
   * ToDo: Change the return type to <tt>RCP<const Array<int> ></tt> since
   * I doubt that the client is allowed to change this array through this
   * function!
   *
   * ToDo: Change the name of this function to something like
   * <tt>getGhostDOFs()</tt>?
   */
  const RCP<Array<int> >& ghostIndices() const 
    {return ghostIndices_;}

  /** \brief Print the DOF map.
   *
   * ToDo: Replace this with override of Teuchos::Describable::describe().
   */
  virtual void print(std::ostream& os) const = 0 ;

  /** \brief Returns <tt>true</tt> if the map is homogeneous.
   *
   * See above defintion for a "Homogeneous" map.
   */
  virtual bool isHomogeneous() const {return false;}

  /**
   * \brief The largest dimension cell supported by this DOF map. Usually, this will
   * be the spatial dimension of the mesh. However, for functions defined only on a surface,
   * curve, or point set it may be lower. Such maps should override the default. 
   */
  virtual int cellDim() const {return mesh_.spatialDim();}


  int setupVerb() const {return setupVerb_;}
protected:

  void setLowestLocalDOF(int low) {lowestLocalDOF_ = low;}

  void setNumLocalDOFs(int numDOFs) {numLocalDOFs_ = numDOFs;}

  void setTotalNumDOFs(int numDOFs) {numDOFs_ = numDOFs;}

  const MPIComm& comm() const {return mesh().comm();}

  void addGhostIndex(int dof) {ghostIndices_->append(dof);}

  static Teuchos::Time& dofLookupTimer() ;

  static Teuchos::Time& batchedDofLookupTimer() ;



private:

  int setupVerb_;

  int localProcID_;

  Mesh mesh_;

  int lowestLocalDOF_;

  int numLocalDOFs_;

  int numDOFs_;

  RCP<Array<int> > ghostIndices_;

};
}


#endif
