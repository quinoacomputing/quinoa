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
 * SundancePeanoMesher2D.cpp
 *
 *  Created on: Sep 16, 2009
 *      Author: benk
 */

#include "SundancePeanoMesher2D.hpp"

using namespace Sundance;
using namespace Sundance;

using namespace Teuchos;
using namespace Sundance;


PeanoMesher2D::PeanoMesher2D(const ParameterList& params)
: MeshSourceBase(params),
  _position_x(params.get<double>("position_x")),
  _position_y(params.get<double>("position_y")),
  _offset_x(params.get<int>("offset_x")),
  _offset_y(params.get<double>("offset_y")),
  _resolution(params.get<double>("resolution"))
{
	// nothing to do here
}

Mesh PeanoMesher2D::fillMesh() const
{
	// here we create the Peano grid and return to the Sundance
	PeanoMesh2D *peanogrid;
	Mesh mesh = createMesh(2);
	// get the pointer to the Peano grid, and then create (in a complicated manner 2Xsmart pointers)
	peanogrid = (PeanoMesh2D*) mesh.ptr().get();
	peanogrid->createMesh( _position_x , _position_y , _offset_x , _offset_y , _resolution );
	return (mesh); // at this stage the Peano grid is completly created
}

