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

#ifndef SUNDANCE_MESHSOURCEBASE_H
#define SUNDANCE_MESHSOURCEBASE_H


#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceMeshType.hpp"
#include "PlayaHandleable.hpp"
#include "Teuchos_Describable.hpp"
#include "PlayaPrintable.hpp"
#include "SundanceNoncopyable.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include "SundanceIncrementallyCreatableMesh.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Sundance
{
/**
 * MeshSourceBase provides the internal interface for mesh sources, i.e.,
 * objects such as mesh generators and file readers that can produce
 * a mesh object. The action of a mesh source should be independent
 * of the internal mesh representation used. To allow user-level
 * specification of the type of internal mesh representation to be
 * used, a MeshSourceBase is constructed with a MeshType object
 * which acts as a factory to produce empty meshes.
 *
 * Mesh sources are constructed with a MPI communicator. If the
 * communicator has more than one processor, the mesh created will
 * be distributed.
 *
 * <h4> Writing your own MeshSourceBase subtype </h4>
 *
 * The only method you will need to override is
 * <ul>
 * <li> <tt>virtual Mesh fillMesh() const </tt> 
 * </ul>
 * which is where you fill in a mesh object and return it. This method
 * should usually physically create the mesh with a call to createMesh(),
 * ensuring that the correct mesh representation type is created
 * using the MeshType factory with which the source was constructed.
 *
 * See the PartitionedLineMesher source code for a very simple
 * example of how to write a mesh source subtype. 
 *
 * Optionally, you can override the description() method to 
 * provide more informative descriptive output than the std::string
 * <tt>"MeshSourceBase[unknown subtype]".</tt>
 */
class MeshSourceBase : public Playa::Handleable<MeshSourceBase>,
                       public Playa::Printable,
                       public Teuchos::Describable,
                       public Noncopyable,
                       public Playa::ObjectWithVerbosity
{
public:
  /** Construct with a mesh type and MPI communicator */
  MeshSourceBase(const MeshType& meshType,
    int verbosity,
    const MPIComm& comm);

  /** Construct from a parameter list */
  MeshSourceBase(const ParameterList& params);

  /** virtual dtor */
  virtual ~MeshSourceBase(){;}

      
  /** Get a mesh from the source. If a mesh has already been created,
   * this method will return the cached mesh, otherwise it will 
   * build on with a call to fillMesh() */
  Mesh getMesh() const ;

  /** \name Extraction of attributes */
  //@{
  void getAttributes(RCP<Array<Array<double> > >& nodeAttributes,
    RCP<Array<Array<double> > >& elemAttributes) const ;
  //@}

  /** \name Printable interface */
  //@{
  /** Print to a stream */
  virtual void print(std::ostream& os) const {os << description();}
  //@}

  /** \name Describable interface */
  //@{
  /** Print to a stream */
  virtual std::string description() const 
    {return "MeshSourceBase[unknown subtype]";}

  /** access to the MPI communicator */
  const MPIComm& comm() const {return comm_;}
  //@}
      
protected:


  /** Get processor rank */
  int myRank() const {return comm().getRank();}
  /** Get number of processors */
  int nProc() const {return comm().getNProc();}

  /** Fill a mesh object with data from the source. Subclass
   * implementors will need to provide a fillMesh() for their
   * subclass. Implementors should use the createMesh() method
   * for the allocation of the new mesh. */
  virtual Mesh fillMesh() const = 0 ;

  /** createMesh() physically allocates the mesh object */
  Mesh createMesh(int dim) const ;

  /** internal access to the node attributes */
  RCP<Array<Array<double> > >& nodeAttributes() const 
    {return nodeAttributes_;}

  /** internal access to the element attributes */
  RCP<Array<Array<double> > >& elemAttributes() const 
    {return elemAttributes_;}



private:

  /** */
  mutable Mesh cachedMesh_;

  /** */
  mutable bool hasCachedMesh_;

  /** */
  MeshType meshType_;

  /** */
  MPIComm comm_;

  /** */
  mutable RCP<Array<Array<double> > > elemAttributes_;

  /** */
  mutable RCP<Array<Array<double> > > nodeAttributes_;

};
}


#endif
