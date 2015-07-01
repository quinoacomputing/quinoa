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

#ifndef SUNDANCE_PARTITIONEDLINEMESHER_H
#define SUNDANCE_PARTITIONEDLINEMESHER_H

#include "SundanceDefs.hpp"
#include "SundanceMeshSourceBase.hpp"

namespace Sundance
{
using namespace Teuchos;
  
/**
 * PartitionedLineMesher meshes the one-dimensional interval 
 * \f$\left[a_x, b_x\right]\f$
 * with \f$n_x\f$ elements per processor. 
 */
class PartitionedLineMesher : public MeshSourceBase
{
public:
  /** 
   * Set up a mesher for the interval \f$\left[a_x, b_x\right]\f$
   * with \f$n_x\f$ elements per processor. 
   */
  PartitionedLineMesher(double ax, double bx, int nx,
    const MeshType& meshType,
    int verbosity=0,
    const MPIComm& comm = MPIComm::world())
    : 
    MeshSourceBase(meshType, verbosity, comm),
    ax_(ax), bx_(bx), nx_(nx) {;}

  /** Create a line mesher from a ParameterList */
  PartitionedLineMesher(const ParameterList& params);
    
  /** */
  virtual ~PartitionedLineMesher() {;}

  /** Print a short descriptive std::string */
  virtual std::string description() const 
    {return "PartitionedLineMesher[ax=" + Teuchos::toString(ax_)
        + ", bx=" + Teuchos::toString(bx_)
        + ", nx=" + Teuchos::toString(nx_) + "]";}
      

  /** Return a ref count pointer to self */
  virtual RCP<MeshSourceBase> getRcp() {return rcp(this);}

private:

  /** */
  virtual Mesh fillMesh() const ;

  /** */
  double ax_;
  /** */
  double bx_;
  /** */
  int nx_;

};
}

#endif
