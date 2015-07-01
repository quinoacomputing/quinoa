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

#include "PlayaEpetraGhostImporter.hpp"
#include "PlayaEpetraGhostView.hpp"
#include "PlayaEpetraVector.hpp"
#include "PlayaTabs.hpp"
#include "PlayaOut.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#endif


namespace Playa
{
using namespace Teuchos;


EpetraGhostImporter
::EpetraGhostImporter(const RCP<const Epetra_Map>& localMap,
  int nGhost,
  const int* ghostElements)
  : localMap_(localMap),
    ghostMap_(),
    importer_()
{
  if (false && nGhost==0)
  {
    ghostMap_ = localMap_;

    importer_ = rcp(new Epetra_Import(*ghostMap_, *localMap_));
  }
  else
  {
    //bvbw not used      int nGlobal = localMap_->NumGlobalElements();
    int nLocal = localMap_->NumMyElements();
    int nGhostView = nLocal+nGhost;
    std::vector<int> globalIndices(nGhostView);
    for (int i=0; i<nLocal; i++) globalIndices[i] = localMap_->GID(i);
    for (int i=0; i<nGhost; i++) globalIndices[i+nLocal] = ghostElements[i];

    const Epetra_Comm& comm = localMap_->Comm();

    ghostMap_ = rcp(new Epetra_Map(-1, nGhostView, 
        &(globalIndices[0]), 0, comm));

    importer_ = rcp(new Epetra_Import(*ghostMap_, *localMap_));
  }
}

void EpetraGhostImporter
::importView(const Vector<double>& x,
  RCP<GhostView<double> >& ghostView) const
{
  Tabs tab;

  /* If given an uninitialized ghost view, create a EpetraGhostView */
  if (ghostView.get()==0) 
  {
    ghostView = rcp(new EpetraGhostView());
  }

  /* Ensure that the ghost view contains an EpetraGhostView */
  EpetraGhostView* epgv 
    = dynamic_cast<EpetraGhostView*>(ghostView.get());

  TEUCHOS_TEST_FOR_EXCEPTION(epgv==0, std::runtime_error,
    "argument ghostView to EpetraGhostImporter::importView() "
    "could not be cast to a EpetraGhostView pointer");

  const Epetra_Vector& xVec = EpetraVector::getConcrete(x);

  /* Do the import */
  epgv->import(*importer_, xVec);
}
    
}
