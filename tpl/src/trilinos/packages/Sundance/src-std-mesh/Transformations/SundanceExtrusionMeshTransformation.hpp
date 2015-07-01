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

#ifndef SUNDANCE_EXTRUSIONMESHFILTER_H
#define SUNDANCE_EXTRUSIONMESHFILTER_H

#include "SundanceDefs.hpp"
#include "SundanceMeshTransformationBase.hpp"

namespace Sundance
{
  
  /**
   * ExtrusionMeshTransformation extrudes a 2D mesh to 3D. 
   */
  class ExtrusionMeshTransformation : public MeshTransformationBase
  {
  public:
    /** Construct a filter to extrude a 2D mesh from
     * the plane \f$z=z_0\f$ to the plane \f$z=z_1\f$ in 
     * \f$n_z\f$ steps. */
    ExtrusionMeshTransformation(double z0, double z1, int nzLevels,
                        const MeshType& meshType)
      : MeshTransformationBase(meshType),
        z0_(z0), z1_(z1), nzLevels_(nzLevels) {;}

    /** virtual dtor */
    virtual ~ExtrusionMeshTransformation() {;}

    /** Apply the filter to an input mesh, returning an output mesh */
    virtual Mesh apply(const Mesh& inputMesh) const ;

    /** Print a short descriptive std::string */
    virtual std::string description() const 
    {return "ExtrusionMeshTransformation[z0=" + Teuchos::toString(z0_)
       + ", z1=" + Teuchos::toString(z1_)
       + ", nz=" + Teuchos::toString(nzLevels_) + "]";}
      

    /** Return a ref count pointer to self */
    virtual RCP<MeshTransformationBase> getRcp() {return rcp(this);}

  private:
    
    /** */
    double z0_;
    
    /** */
    double z1_;

    /** */
    int nzLevels_;
    
  };
}

#endif
