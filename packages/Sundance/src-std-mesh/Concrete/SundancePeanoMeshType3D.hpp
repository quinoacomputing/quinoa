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
 * SundancePeanoMeshType3D.hpp
 *
 *  Created on: Sep 16, 2009
 *      Author: benk
 */

#ifndef SUNDANCE_PEANOMESHTYPE3D_HPP_
#define SUNDANCE_PEANOMESHTYPE3D_HPP_

#include "SundanceDefs.hpp"
#include "SundanceMeshBase.hpp"
#include "SundancePeanoMesh3D.hpp"

namespace Sundance
{
  using namespace Teuchos;

  /**
   * PeanoMeshType is used to create
   * PeanoMesh objects.
   */
  class PeanoMeshType3D : public MeshTypeBase
  {
  public:
    /** Empty ctor */
	  PeanoMeshType3D(const MeshEntityOrder& order=ExodusMeshOrder)
	    : order_(order) {;}

    /** virtual dtor */
    virtual ~PeanoMeshType3D(){;}

    /** Create a mesh of the given dimension */
    virtual RCP<MeshBase> createEmptyMesh(int dim,
                                                  const MPIComm& comm) const
    // this line is never used since we create directly the mesh at the Mesher
    {return rcp(new PeanoMesh3D(dim, comm, order_));}

    /** */
    std::string description() const {return "PeanoMeshType3D";}

#ifndef DOXYGEN_DEVELOPER_ONLY
    /** Return a ref count pointer to self */
    virtual RCP<MeshTypeBase> getRcp() {return rcp(this);}
#endif  /* DOXYGEN_DEVELOPER_ONLY */

  private:
    MeshEntityOrder order_;

  };
}

#endif /* SUNDANCEPEANOMESHTYPE3D_HPP_ */
