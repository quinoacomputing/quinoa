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

#ifndef SUNDANCE_PARTITIONEDRECTANGLEMESHER_H
#define SUNDANCE_PARTITIONEDRECTANGLEMESHER_H

#include "SundanceDefs.hpp"
#include "SundanceMeshSourceBase.hpp"

namespace Sundance
{
using namespace Teuchos;



/**
 * PartitionedRectangleMesher meshes the rectangle 
 * \f$ \left[a_x, b_x\right] \otimes \left[a_y, b_y\right] \f$
 * with \f$ n_x \otimes n_y \f$ elements per processor. The 
 * rectangle is partitioned among processors, with \f$np_x\f$
 * equal sized 
 * subdomains in the \f$x\f$ direction and \f$np_y\f$ in the \f$y\f$
 * direction.
 */
class PartitionedRectangleMesher : public MeshSourceBase
{
public:
  /** 
   * Set up meshing of the rectangle 
   * \f$ \left[a_x, b_x\right] \otimes \left[a_y, b_y\right] \f$
   * with \f$ n_x \otimes n_y \f$ elements per processor. The 
   * rectangle is partitioned among processors, with \f$np_x\f$
   * equal sized 
   * subdomains in the \f$x\f$ direction and \f$np_y\f$ in the \f$y\f$
   * direction.
   */
  PartitionedRectangleMesher(double ax, double bx, int nx, int npx,
    double ay, double by, int ny, int npy,
    const MeshType& meshType,
    int verbosity=0,
    const MPIComm& comm = MPIComm::world())
    : 
    MeshSourceBase(meshType, verbosity, comm),
    ax_(ax), bx_(bx), nx_(nx), npx_(npx),
    ay_(ay), by_(by), ny_(ny), npy_(npy) {;}

  /** 
   * Set up meshing of the rectangle 
   * \f$ \left[a_x, b_x\right] \otimes \left[a_y, b_y\right] \f$
   * with \f$ n_x \otimes n_y \f$ elements per processor. The
   * balance() function is used to choose npx and npy.
   */
  PartitionedRectangleMesher(double ax, double bx, int nx,
    double ay, double by, int ny,
    const MeshType& meshType,
    int verbosity=0,
    const MPIComm& comm = MPIComm::world())
    : 
    MeshSourceBase(meshType, verbosity, comm),
    ax_(ax), bx_(bx), nx_(nx), npx_(-1),
    ay_(ay), by_(by), ny_(ny), npy_(-1) 
    {
      balanceXY(comm.getNProc(), &npx_, &npy_);
    }

    
  /** Create a rectangle mesher from a ParameterList */
  PartitionedRectangleMesher(const ParameterList& params);
    
  /** */
  virtual ~PartitionedRectangleMesher() {;}

  /** Print a short descriptive std::string */
  virtual std::string description() const 
    {return "PartitionedRectangleMesher[ax=" + Teuchos::toString(ax_)
        + ", bx=" + Teuchos::toString(bx_)
        + ", nx=" + Teuchos::toString(nx_) +
        + ", ay=" + Teuchos::toString(ay_)
        + ", by=" + Teuchos::toString(by_)
        + ", ny=" + Teuchos::toString(ny_) + "]";}



  /** Find a nearly equal balance between X and Y partitions */
  static void balanceXY(int n, int* npx, int* npy);
    

  /** Return a ref count pointer to self */
  virtual RCP<MeshSourceBase> getRcp() {return rcp(this);}



protected:

  /** */
  virtual Mesh fillMesh() const ;

private:

  /** */
  double ax_;
  /** */
  double bx_;
  /** */
  int nx_;
  /** */
  int npx_;

  /** */
  double ay_;
  /** */
  double by_;
  /** */
  int ny_;
  /** */
  int npy_;


};
}

#endif
