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

#ifndef PLAYA_EPETRAVECTORTYPE_HPP
#define PLAYA_EPETRAVECTORTYPE_HPP

#include "PlayaEpetraVectorSpace.hpp"
#include "PlayaHandleable.hpp"
#include "PlayaPrintable.hpp"
#include "Teuchos_Describable.hpp"
#include "PlayaVectorTypeBase.hpp"


namespace Playa
{
using namespace Teuchos;
  
/**
 * \!brief Epetra vector type is a factory for epetra vector spaces
 */
class EpetraVectorType : public VectorTypeBase<double>,
                         public Playa::Handleable<VectorTypeBase<double> >,
                         public Printable,
                         public Describable
{
public:
  /** Construct a vector type */
  EpetraVectorType();
      
  /** virtual dtor */
  virtual ~EpetraVectorType() {;}

  /** create a distributed vector space.
   * @param dimension the dimension of the space 
   * @param nLocal number of indices owned by the local processor
   * @param locallyOwnedIndices array of indices owned by this processor  
   */
  RCP<const VectorSpaceBase<double> > 
  createSpace(int dimension, 
    int nLocal,
    const int* locallyOwnedIndices,
    const MPIComm& comm) const ;

  /**  
   * Create an importer for accessing ghost elements.
   * @param space the distributed vector space on which ghost elements
   * are to be shared
   * @param nGhost number of ghost elements needed by this processor
   * @param ghostIndices read-only C array of off-processor indices needed
   * by this processor.
   * @return A RCP to a GhostImporter object.
   */
  RCP<GhostImporter<double> > 
  createGhostImporter(const VectorSpace<double>& space,
    int nGhost,
    const int* ghostIndices) const ;


  /**
   * Create a matrix factory of type compatible with this vector type,
   * sized according to the given domain and range spaces.
   */
  RCP<MatrixFactory<double> >
  createMatrixFactory(const VectorSpace<double>& domain,
    const VectorSpace<double>& range) const ;

    

  /** \name Printable interface */
  //@{
  /** Print to stream */
  void print(std::ostream& os) const {os << description();}
  //@}

  GET_RCP(VectorTypeBase<double>);
};
  
}

#endif
