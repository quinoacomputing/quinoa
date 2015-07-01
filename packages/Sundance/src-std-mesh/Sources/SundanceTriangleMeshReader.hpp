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

#ifndef SUNDANCE_TRIANGLEMESHREADER_H
#define SUNDANCE_TRIANGLEMESHREADER_H

#include "SundanceDefs.hpp"
#include "SundanceMeshReaderBase.hpp"

namespace Sundance
{
  using namespace Teuchos;
  
  /**
   * TriangleMeshReader reads a mesh stored in Shewchuk's Triangle format.
   * This format is documented at 
   * <A HREF="http://www-2.cs.cmu.edu/~quake/triangle.html"> 
   * the Triangle homepage. 
   * </A>
   * This reader expects to find node information in <tt>.node</tt> files
   * and element information in <tt>.ele</tt> files. The <tt> filename </tt>
   * constructor argument is the stem of the filenames, and so that 
   * a reader constructed with filename <tt>joe</tt> will look for node and
   * element data in <tt>joe.node</tt> and <tt>joe.ele</tt> respectively.
   * Node and element
   * attributes are read if present, and can be accessed with the 
   * <tt>getAttributes()</tt> method of <tt>MeshSource.</tt>
   * 
   * <h4> Parallel extensions </h4>
   * We have extended the Triangle format to deal with distributed meshes.
   * A TriangleMeshReader is constructed with an MPIComm object, and if
   * that communicator has more than one processor the mesh is assumed
   * to be split into files, one for each processor. Data
   * on mesh "joe" for the <i>nnn</i>-th processor will be found in the files
   * <ul>
   * <li> <tt>joe.node.</tt><i>nnn</i>
   * <li> <tt>joe.ele.</tt><i>nnn</i>
   * <li> <tt>joe.par.</tt><i>nnn</i>
   * </ul>
   * The <tt>.node.</tt><i>nnn</i> and <tt>.ele.</tt><i>nnn</i> files contain the
   * node and element data for the part of the mesh seen 
   * by the <i>nnn</i>-th
   * processor. The node and element 
   * numberings given in those two files are <b>local</b> indexes.
   * The <tt>.par.</tt><i>nnn</i> contains node and element 
   * ownership information for the part of the mesh seen 
   * by the <i>nnn</i>-th
   * processor. 
   *
   * <br> 
   *
   * A <tt>.par</tt> file is formatted as follows:
   * <ul>
   * <li> First line: <tt> rank numProcs </tt>
   * <li> Second line: <tt> numPoints </tt>
   * <li> Next <i> nPoints </i> lines: <tt> ptLID ptGID ptOwnerRank </tt>
   * <li> Next line: <tt> numElems </tt>
   * <li> Next <i> nElems </i> lines: <tt> elemLID elemGID elemOwnerRank </tt>
   * </ul>
   * 
   */
  class TriangleMeshReader : public MeshReaderBase
  {
  public:
    /** */
    TriangleMeshReader(const std::string& filename, 
      const MeshType& meshType,
      int verbosity=0,
      const MPIComm& comm = MPIComm::world());

    /** Construct from a ParameterList */
    TriangleMeshReader(const ParameterList& params);

    /** virtual dtor */
    virtual ~TriangleMeshReader(){;}


    /** Create a mesh */
    virtual Mesh fillMesh() const ;

    /** Print a short descriptive std::string */
    virtual std::string description() const 
    {return "TriangleMeshReader[file=" + filename() + "]";}
      

    /** Return a ref count pointer to self */
    virtual RCP<MeshSourceBase> getRcp() {return rcp(this);}

  private:
    /** */
    void readParallelInfo(Array<int>& ptGID, Array<int>& ptOwner,
                          Array<int>& elemGID, Array<int>& elemOwner) const ;

    /** */
    Mesh readNodes(Array<int>& ptGID,
                   Array<int>& ptOwner) const ;

    /** */
    void readSides(Mesh& mesh) const ;

    /** */
    void readEdges(Mesh& mesh) const ;

    /** */
    void readElems(Mesh& mesh,
                   const Array<int>& nodeGID,
                   Array<int>& elemGID,
                   Array<int>& elemOwner) const ;
    

    /** */
    std::string nodeFilename_;

    /** */
    std::string elemFilename_;

    /** */
    std::string parFilename_;

    /** */
    std::string sideFilename_;

    /** */
    std::string edgeFilename_;

    /** */
    mutable int offset_;
  };
}

#endif
