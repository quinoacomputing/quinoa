/* @HEADER@ */
// ************************************************************************
// 
//                 Playa: Programmable Linear Algebra
//                 Copyright 2012 Sandia Corporation
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

#include "PlayaEpetraGhostView.hpp"
#include "Epetra_Import.h"


namespace Playa
{

using namespace Teuchos;

const double& EpetraGhostView::getElement(int globalIndex) const 
{
  const Epetra_BlockMap& myMap = ghostView_->Map();
  return (*ghostView_)[myMap.LID(globalIndex)];
}

void EpetraGhostView::getElements(const int* globalIndices, int numElems,
  Array<double>& elems) const
{
  elems.resize(numElems);
  const Epetra_BlockMap& myMap = ghostView_->Map();

  for (int i=0; i<numElems; i++)
  {
    elems[i] = (*ghostView_)[myMap.LID(globalIndices[i])];
  }
}

void  EpetraGhostView::import(const Epetra_Import& importer,
  const Epetra_Vector& srcObject)
{
  /* If my vector does not yet exist, create it using the target map of the
   * importer */
  if (ghostView_.get()==0)
  {
    ghostView_ = rcp(new Epetra_Vector(importer.TargetMap()));
  }

  /* do the import */
  int ierr = ghostView_->Import(srcObject, importer, Insert);

  if (ierr < 0)
  {
    Out::os() << "target map=" << std::endl;
    importer.TargetMap().Print(Out::os());
    Out::os() << "source map=" << std::endl;
    srcObject.Map().Print(Out::os());
  }
  TEUCHOS_TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error, "ierr=" << ierr << " in EpetraGhostView::import()");
}

void EpetraGhostView::print(std::ostream& os) const
{
  if (ghostView_.get()==0) 
  {
    os << "[null Epetra ghost view]" << std::endl;
  }
  else
  {
    ghostView_->Print(os);
  }
}

}
