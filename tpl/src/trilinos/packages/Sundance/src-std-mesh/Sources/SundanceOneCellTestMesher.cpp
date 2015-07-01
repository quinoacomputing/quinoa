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


#include "SundanceOneCellTestMesher.hpp"

using namespace Sundance;
using namespace Teuchos;


OneTriangleMesher::OneTriangleMesher(
  const Point& A,
  const Point& B,
  const Point& C,
  const MeshType& meshType)
  : MeshSourceBase(meshType, 0, MPIComm::world()),
    A_(A), B_(B), C_(C)
{}


Mesh OneTriangleMesher::fillMesh() const
{
  Mesh mesh = createMesh(2);
  
  int a = mesh.addVertex(0, A_, 0, 0);
  int b = mesh.addVertex(1, B_, 0, 0);
  int c = mesh.addVertex(2, C_, 0, 0);

  mesh.addElement(0, tuple(a,b,c), 0, 0);

  mesh.freezeTopology();

  return mesh;
}
  


OneTetMesher::OneTetMesher(
  const Point& A,
  const Point& B,
  const Point& C,
  const Point& D,
  const MeshType& meshType)
  : MeshSourceBase(meshType, 0, MPIComm::world()),
    A_(A), B_(B), C_(C), D_(D)
{}


Mesh OneTetMesher::fillMesh() const
{
  Mesh mesh = createMesh(3);
  
  int a = mesh.addVertex(0, A_, 0, 0);
  int b = mesh.addVertex(1, B_, 0, 0);
  int c = mesh.addVertex(2, C_, 0, 0);
  int d = mesh.addVertex(3, D_, 0, 0);

  mesh.addElement(0, tuple(a,b,c,d), 0, 0);

  return mesh;
}
  
