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

#ifndef PLAYA_EPETRAMATRIXFACTORY_HPP
#define PLAYA_EPETRAMATRIXFACTORY_HPP

#include "PlayaEpetraVectorSpace.hpp"
#include "PlayaIncrementallyConfigurableMatrixFactory.hpp"
#include "PlayaCollectivelyConfigurableMatrixFactory.hpp"
#include "PlayaMatrixFactory.hpp"
#include "PlayaHandleable.hpp"
#include "PlayaPrintable.hpp"
#include "Epetra_CrsGraph.h"

namespace Playa
{
using namespace Teuchos;
  

/** */
class EpetraMatrixFactory : public MatrixFactory<double>,
                            public IncrementallyConfigurableMatrixFactory,
                            public CollectivelyConfigurableMatrixFactory
{
public:

  /** Construct an uninitialized EpetraMatrixFactory */
  EpetraMatrixFactory(const RCP<const EpetraVectorSpace>& domain,
    const RCP<const EpetraVectorSpace>& range);

  /** */
  const RCP<const EpetraVectorSpace>& epRange() const {return range_;}

  /** */
  const RCP<const EpetraVectorSpace>& epDomain() const {return domain_;}


  /** Initialize a set of nonzero elements in the matrix's graph.
   * @param globalRowIndex the global index of the row to which these
   * elements belong.
   * @param nElemsToInsert the number of elements being inserted in this
   * step
   * @param globalColumnIndices array of column indices. Must 
   * be nElemsToInsert in length. 
   */
  virtual void initializeNonzerosInRow(int globalRowIndex,
    int nElemsToInsert,
    const int* globalColumnIndices) ;

  /** 
   * Initialize nonzeros in a batch of rows. 
   */
  virtual void initializeNonzeroBatch(int numRows, 
    int rowBlockSize,
    const int* globalRowIndices,
    int numColumnsPerRow,
    const int* globalColumnIndices,
    const int* skipRow);

  /** Configure all rows at once */
  virtual void configure(int lowestRow,
    const std::vector<int>& rowPtrs,
    const std::vector<int>& nnzPerRow,
    const std::vector<int>& data);

  /** */
  void finalize();

  /** */
  const Epetra_CrsGraph& graph() const ;

  /** */
  virtual LinearOperator<double> createMatrix() const ;

protected:

private:

  /** */
  RCP<Epetra_CrsGraph> graph_;

  /** */
  RCP<const EpetraVectorSpace> range_;

  /** */
  RCP<const EpetraVectorSpace> domain_;
};
}

#endif
