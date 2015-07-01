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

#ifndef SUNDANCE_ONE_CELL_TEST_MESHER_H
#define SUNDANCE_ONE_CELL_TEST_MESHER_H

#include "SundanceDefs.hpp"
#include "SundanceMeshReaderBase.hpp"

namespace Sundance
{
using namespace Teuchos;
  
/**
 * Create a mesh having one triangle
 */
class OneTriangleMesher : public MeshSourceBase
{
public:
  /** */
  OneTriangleMesher(
    const Point& A, 
    const Point& B,
    const Point& C,
    const MeshType& meshType);

  /** virtual dtor */
  virtual ~OneTriangleMesher(){;}


  /** Create a mesh */
  virtual Mesh fillMesh() const ;

  /** Print a short descriptive std::string */
  virtual std::string description() const 
    {return "OneTriangleMesher";}
      

  /** Return a ref count pointer to self */
  virtual RCP<MeshSourceBase> getRcp() {return rcp(this);}

private:
  Point A_;
  Point B_;
  Point C_;
};

/**
 * Create a mesh having one tet
 */
class OneTetMesher : public MeshSourceBase
{
public:
  /** */
  OneTetMesher(
    const Point& A, 
    const Point& B,
    const Point& C,
    const Point& D,
    const MeshType& meshType);

  /** virtual dtor */
  virtual ~OneTetMesher(){;}


  /** Create a mesh */
  virtual Mesh fillMesh() const ;

  /** Print a short descriptive std::string */
  virtual std::string description() const 
    {return "OneTetMesher";}
      

  /** Return a ref count pointer to self */
  virtual RCP<MeshSourceBase> getRcp() {return rcp(this);}

private:
  Point A_;
  Point B_;
  Point C_;
  Point D_;
};
}

#endif
