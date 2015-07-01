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

#ifndef PLAYA_EPETRAVECTORSPACE_HPP
#define PLAYA_EPETRAVECTORSPACE_HPP

#include "PlayaDefs.hpp"
#include "Epetra_Map.h"
#include "Teuchos_RefCountPtr.hpp"
#include "PlayaMPIComm.hpp"
#include "PlayaVectorSpaceBaseDecl.hpp"

namespace Playa
{
using namespace Teuchos;


/**
 * Adaptor wrapping Epetra map in the Thyra vector space system.
 */
class EpetraVectorSpace : public VectorSpaceBase<double>
{
public:
  /** */
  EpetraVectorSpace(const RCP<const Epetra_Map>& map);

  /** */
  RCP<VectorBase<double> > 
  createMember(const VectorSpace<double>& self) const  ;
    
  /** @name Overridden form Teuchos::Describable */
  //@{
  /** \brief . */
  std::string description() const;
  //@}

  /** @name Public overridden from VectorSpace */
  //@{

  /** */
  int dim() const {return globalDim_;}
  
  /** */
  int numLocalElements() const {return numLocalElements_;}
  
  /** */
  int baseGlobalNaturalIndex() const {return baseGlobalNaturalIndex_;}

  /** */
  bool isCompatible(const VectorSpaceBase<double>* other) const ;

  int numBlocks() const {return 1;}

  //@}

 /** */
  const RCP<const Epetra_Map>& epetraMap() const 
    {return epetraMap_;}


  /** */
  virtual const MPIComm& comm() const {return comm_;}

protected:

  MPIComm epetraCommToTeuchosMPIComm(const Epetra_Comm& epComm) ;
  
private:
  RCP<const Epetra_Map> epetraMap_;

  MPIComm comm_;

  int globalDim_;

  int baseGlobalNaturalIndex_;

  int numLocalElements_;
};
  
}

#endif
