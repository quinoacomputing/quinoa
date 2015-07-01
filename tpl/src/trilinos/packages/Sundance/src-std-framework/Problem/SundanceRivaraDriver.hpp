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


#ifndef SUNDANCE_RIVARADRIVER_HPP
#define SUNDANCE_RIVARADRIVER_HPP


#include "SundanceMeshType.hpp"
#include "SundanceMeshTransformationBase.hpp"
#include "SundanceRivaraMesh.hpp"
#include "SundanceExpr.hpp"

namespace Sundance
{

class Mesh;
using Sundance::Expr;

class RefinementTransformation : public MeshTransformationBase
{
public:
  /** */
  RefinementTransformation(const MeshType& meshType, const Expr& errExpr,
    const double& reqErr, const double& minArea)
    : MeshTransformationBase(meshType), meshType_(meshType),
      errExpr_(errExpr), reqErr_(reqErr), minArea_(minArea),
      numRefined_(-1) {}

  /** */
  Mesh apply(const Mesh& inputMesh) const ;

  /** */
  int numRefined() const {return numRefined_;}

  /* */
  GET_RCP(MeshTransformationBase);
private:
  /** */
  void meshToRivara(
    const Mesh& mesh,
    Array<int>& lidMap,
    RCP<Rivara::RivaraMesh>& rivMesh) const ;

  /** */
  Mesh rivaraToMesh(const RCP<Rivara::RivaraMesh>& rivMesh,
    const MPIComm& comm) const ;

  MeshType meshType_;
  Expr errExpr_;
  double reqErr_;
  double minArea_;
  mutable int numRefined_;
};


}


#endif
