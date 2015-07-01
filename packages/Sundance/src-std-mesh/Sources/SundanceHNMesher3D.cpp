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


#include "SundanceHNMesher3D.hpp"

using namespace Sundance;
using namespace Teuchos;

REFINE_MESH_ESTIMATE( DummyNoRefineClass , { return false; } , {return 0;} )
MESH_DOMAIN( DummyAllMeshDomain , {return true;})

const RefinementClass HNMesher3D::dummyRefineClass_ = (new DummyNoRefineClass());

const MeshDomainDef HNMesher3D::dummyMeshDomain_ = (new DummyAllMeshDomain());

HNMesher3D::HNMesher3D(const ParameterList& params)
: MeshSourceBase(params),
  _position_x(params.get<double>("position_x")),
  _position_y(params.get<double>("position_y")),
  _position_z(params.get<double>("position_z")),
  _offset_x(params.get<double>("offset_x")),
  _offset_y(params.get<double>("offset_y")),
  _offset_z(params.get<double>("offset_z")),
  _resolution_x(params.get<int>("resolution_x")),
  _resolution_y(params.get<int>("resolution_y")),
  _resolution_z(params.get<int>("resolution_z")),
  refineClass_(HNMesher3D::dummyRefineClass_),
  meshDomain_(HNMesher3D::dummyMeshDomain_)
{
	// nothing to do here
}

Mesh HNMesher3D::fillMesh() const
{
	// here we create the mesh and return to the Sundance
	HNMesh3D *hnodegrid;
	Mesh mesh = createMesh(3);
	// get the pointer
	hnodegrid = (HNMesh3D*) mesh.ptr().get();
	hnodegrid->createMesh( _position_x , _position_y ,_position_z ,
			               _offset_x , _offset_y , _offset_z ,
			               _resolution_x , _resolution_y , _resolution_z ,
			                refineClass_ , meshDomain_);
	return (mesh);
}

