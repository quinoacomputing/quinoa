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

#include "SundanceFieldWriterBase.hpp"
#include "PlayaExceptions.hpp"
#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"


using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;


FieldWriterBase::FieldWriterBase(const std::string& filename) 
  : filename_(filename),
    mesh_(),
    nProc_(0), 
    myRank_(-1),
    meshID_(-1),
    comments_(),
    pointScalarFields_(),
    cellScalarFields_(),
    pointVectorFields_(),
    cellVectorFields_(),
    pointScalarNames_(),
    cellScalarNames_(),
    pointVectorNames_(),
    cellVectorNames_(),
    undefinedValue_(0.0)
{;}


void FieldWriterBase::impersonateParallelProc(int nProc, int rank)
{
  nProc_ = nProc;
  myRank_ = rank;
}

int FieldWriterBase::nProc() const
{
  if (nProc_ < 1) return mesh().comm().getNProc(); 
  return nProc_;
}

int FieldWriterBase::myRank() const
{
  if (myRank_ < 0) return mesh().comm().getRank(); 
  return myRank_;
}




void FieldWriterBase::addMesh(const Mesh& mesh) 
{
  if (meshID_ < 0)
    {
      mesh_ = mesh;
      meshID_ = mesh.id();
    }
                     
  TEUCHOS_TEST_FOR_EXCEPTION(meshID_ != mesh.id(), std::runtime_error,
                     "FieldWriterBase::setMesh(): inconsistent meshes: "
                     "existing mesh has meshID=" << meshID_ << ", newly "
                     "added mesh has meshID=" << mesh.id());
}

void FieldWriterBase::addField(const std::string& name, 
                               const RCP<FieldBase>& expr) 
{

  std::string fieldName = name;

  if (expr->numElems() > 1)
    {
      //TEUCHOS_TEST_FOR_EXCEPTION(expr->numElems() > 1, std::runtime_error,
      //                   "FieldWriterBase::addField not ready for vector fields");

	  std::cout << "WARNING! : expr->numElems() > 1 , FieldWriterBase::addField only VTK can plot vector field " << std::endl;
	  std::cout << "WARNING! : All expressions(in the list of the expressions) must be of the same kind!!! " << std::endl;
	  /* Vector field plotting should be implemented for VTK files */
	  /* We assume that all the expressions are the same !!!! */
	  if (expr->isPointData()){
		  pointVectorFields_.append(expr);
		  pointVectorNames_.append(fieldName);
	  }else{
		  cellVectorFields_.append(expr);
		  cellVectorNames_.append(fieldName);
	  }

    } 
  else if (expr->isPointData()) 
    {
      /* expr is a single scalar field defined at points */
      pointScalarFields_.append(expr);
      pointScalarNames_.append(fieldName);
    }
  else if (expr->isCellData())
    {
      /* expr is a single scalar field defined on cells */
      cellScalarFields_.append(expr);
      cellScalarNames_.append(fieldName);
    }
}

void FieldWriterBase::addCommentLine(const std::string& line) 
{
  comments_.append(line);
}






