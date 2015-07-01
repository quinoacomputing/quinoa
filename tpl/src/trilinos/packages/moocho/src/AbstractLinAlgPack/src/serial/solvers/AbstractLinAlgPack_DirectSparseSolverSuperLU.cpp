// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifdef SPARSE_SOLVER_PACK_USE_SUPERLU

#include <assert.h>

#include <fstream>
#include <algorithm>

#include "AbstractLinAlgPack_DirectSparseSolverSuperLU.hpp"
#include "AbstractLinAlgPack_VectorDenseEncap.hpp"
#include "DenseLinAlgPack_PermVecMat.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace {

// Convert to compressed sparse row
void convet_to_csr(
  int                                   n
  ,int                                  m
  ,int                                  nz
  ,const DenseLinAlgPack::value_type    a_val[]
  ,const DenseLinAlgPack::index_type    a_row_i[]
  ,const DenseLinAlgPack::index_type    a_col_j[]
  ,double                               acsr_val[]
  ,int                                  acsr_col_j[]
  ,int                                  acsr_row_ptr[]
  )
{
  // Count the number of entries per row and put in acsr_row_ptr[1...m+1]
  std::fill_n( &acsr_row_ptr[0], m+1, 0 );
  {for( int k = 0; k < nz; ++k ) {
    ++acsr_row_ptr[a_row_i[k]]; // a_row_i[] is 1-based so this works out.
  }}

  // Transform the counts of entries per row into the start pointers for the rows.
  // We will make acsr_row_ptr[0] = 0 and then add form there.  We will then
  // shift this data so that acsr_row_ptr[1] = 0. This data
  // structure will be used to fill the entries per row.
  acsr_row_ptr[0] = 0;
  {for( int i = 2; i < m + 1; ++i ) {
    acsr_row_ptr[i] += acsr_row_ptr[i-1];
  }}
  {for( int i = m; i > 0; --i ) {
    acsr_row_ptr[i] = acsr_row_ptr[i-1];
  }}

  // Now copy into the compressed sparse row data structure
  {for( int k = 0; k < nz; ++k ) {
    const int row_i   = a_row_i[k];            // one-based
    const int row_ptr = acsr_row_ptr[row_i];  // returned value is zero-based
    acsr_val[row_ptr]   = a_val[k];
    acsr_col_j[row_ptr] = a_col_j[row_ptr] - 1; // from one-based to zero-based
    ++acsr_row_ptr[row_i];
  }}
  TEUCHOS_TEST_FOR_EXCEPT( !(  acsr_row_ptr[m] == nz  ) );

}

} // end namespace

namespace AbstractLinAlgPack {

//
// Implementation of DirectSparseSolver(Imp) interface using SuperLU.
//
// Here are some little bits of knowledge about SuperLU that I need
// to record after many hours of hard work.
//
// ToDo: Finish this!
//

// ToDo:
// a) Add an option for printing the values of the common
//    block parameters to out or to a file.  This can
//    be used to get a feel for the performance of
//    ma28
// b) Add provisions for an external client to change
//    the control options of SuperLU.  Most of these are
//    stored as common block variables.

// //////////////////////////////////////////////////
// DirectSparseSolverSuperLU::BasisMatrixSuperLU

// Overridden from BasisMatrixImp

Teuchos::RCP<DirectSparseSolverImp::BasisMatrixImp>
DirectSparseSolverSuperLU::BasisMatrixSuperLU::create_matrix() const
{
  return Teuchos::rcp(new BasisMatrixSuperLU);
}

void DirectSparseSolverSuperLU::BasisMatrixSuperLU::V_InvMtV(
  VectorMutable* y, BLAS_Cpp::Transp M_trans, const Vector& x
  ) const 
{
  using Teuchos::dyn_cast;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  size_type k;

  // Get concrete objects
  const FactorizationStructureSuperLU
    &fs = dyn_cast<const FactorizationStructureSuperLU>(*this->get_fact_struc());
  const FactorizationNonzerosSuperLU
    &fn = dyn_cast<const FactorizationNonzerosSuperLU>(*this->get_fact_nonzeros());

  VectorDenseMutableEncap  yd(*y);
  VectorDenseEncap         xd(x);

  yd() = xd(); // Copy rhs into lhs for SuperLU

  fs.superlu_solver_->solve(
    *fs.fact_struct_
    ,*fn.fact_nonzeros_
    ,M_trans == BLAS_Cpp::no_trans
    ,yd().dim()
    ,1
    ,&yd()[0]
    ,yd().dim()
    );

}

// //////////////////////////////////////////////////
// DirectSparseSolverSuperLU::FactorizationStructureSuperLU

DirectSparseSolverSuperLU::FactorizationStructureSuperLU::FactorizationStructureSuperLU()
  :superlu_solver_(SuperLUPack::SuperLUSolver::create_solver())
{}

// //////////////////////////////////////////////////
// DirectSparseSolverSuperLU

// Constructors/initializers

DirectSparseSolverSuperLU::DirectSparseSolverSuperLU()
{}

// Overridden from DirectSparseSolver

const DirectSparseSolver::basis_matrix_factory_ptr_t
DirectSparseSolverSuperLU::basis_matrix_factory() const
{
  namespace mmp = MemMngPack;
  return Teuchos::rcp(new Teuchos::AbstractFactoryStd<BasisMatrix,BasisMatrixSuperLU>());
}

void DirectSparseSolverSuperLU::estimated_fillin_ratio(
  value_type estimated_fillin_ratio
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(true);
}

// Overridden from DirectSparseSolverImp

const Teuchos::RCP<DirectSparseSolver::FactorizationStructure>
DirectSparseSolverSuperLU::create_fact_struc() const
{
  return Teuchos::rcp(new FactorizationStructureSuperLU);
}

const Teuchos::RCP<DirectSparseSolverImp::FactorizationNonzeros>
DirectSparseSolverSuperLU::create_fact_nonzeros() const
{
  return Teuchos::rcp(new FactorizationNonzerosSuperLU);
}

void DirectSparseSolverSuperLU::imp_analyze_and_factor(
  const AbstractLinAlgPack::MatrixConvertToSparse   &A
  ,FactorizationStructure                         *fact_struc
  ,FactorizationNonzeros                          *fact_nonzeros
  ,DenseLinAlgPack::IVector                            *row_perm
  ,DenseLinAlgPack::IVector                            *col_perm
  ,size_type                                      *rank
  ,std::ostream                                   *out
  )
{
  namespace mmp = MemMngPack;
  using Teuchos::dyn_cast;
  typedef MatrixConvertToSparse MCTS;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  if(out)
    *out << "\nUsing SuperLU to analyze and factor a new matrix ...\n";

  // Get the concrete factorization and nonzeros objects
  FactorizationStructureSuperLU
    &fs = dyn_cast<FactorizationStructureSuperLU>(*fact_struc);
  FactorizationNonzerosSuperLU
    &fn = dyn_cast<FactorizationNonzerosSuperLU>(*fact_nonzeros);

  // Allocate new storage if not done so already
  if(!fs.fact_struct_.get())
    fs.fact_struct_ = SuperLUPack::SuperLUSolver::create_fact_struct();
  if(!fn.fact_nonzeros_.get())
    fn.fact_nonzeros_ = SuperLUPack::SuperLUSolver::create_fact_nonzeros();

  // Get the dimensions of things.
  const index_type
    m = A.rows(),
    n = A.cols(),
    nz = A.num_nonzeros( MCTS::EXTRACT_FULL_MATRIX ,MCTS::ELEMENTS_ALLOW_DUPLICATES_SUM );

  // Validate input
  TEUCHOS_TEST_FOR_EXCEPTION(
    n <= 0 || m <= 0 || m > n, std::invalid_argument
    ,"DirectSparseSolverSuperLU::imp_analyze_and_factor(...) : Error!" );

  // Extract the matrix in coordinate format
  Workspace<value_type>   a_val(wss,nz);
  Workspace<index_type>   a_row_i(wss,nz);
  Workspace<index_type>   a_col_j(wss,nz);
  A.coor_extract_nonzeros(
    MCTS::EXTRACT_FULL_MATRIX
    ,MCTS::ELEMENTS_ALLOW_DUPLICATES_SUM
    ,nz
    ,&a_val[0]
    ,nz
    ,&a_row_i[0]
    ,&a_col_j[0]
    );

  //
  // Convert to compressed sparse row format (which is compressed sparse
  // column of the transposed matrix which will be passed to
  // SuperLU for factorization).
  //

  Workspace<double>   acsr_val(wss,nz);
  Workspace<int>      acsr_col_j(wss,nz);    // Zero-based for SuperLU
  Workspace<int>      acsr_row_ptr(wss,m+1);
  
  convet_to_csr(
    n,m,nz,&a_val[0],&a_row_i[0],&a_col_j[0]
    ,&acsr_val[0],&acsr_col_j[0],&acsr_row_ptr[0]
    );

  //
  // Have SuperLU factor this matrix.
  //
  // SuperLU works with the transpose of the matrix
  // That DirectSparseSolver works with.
  //
  
  Workspace<int>    perm_r(wss,m); // Zero-based for SuperLU
  Workspace<int>    perm_c(wss,n); // Zero-based for SuperLU
  int                    slu_rank = 0;

  fs.superlu_solver_->analyze_and_factor(
    n                         // m
    ,m                        // n
    ,nz                       // nz
    ,&acsr_val[0]             // a_val[]
    ,&acsr_col_j[0]           // a_row_i[]
    ,&acsr_row_ptr[0]         // a_col_ptr[]
    ,fs.fact_struct_.get()    // fact_struct
    ,fn.fact_nonzeros_.get()  // fact_nonzeros
    ,&perm_c[0]               // perm_r[]
    ,&perm_r[0]               // perm_c[]
    ,&slu_rank                // rank
    );

  // Copy the data to the output
  row_perm->resize(m);
  for( int i = 0; i < m; ++i )
    (*row_perm)[i] = perm_r[i] + 1; // Convert from zero based to one based
  col_perm->resize(n);
  for( int j = 0; j < n; ++j )
    (*col_perm)[j] = perm_c[j] + 1; // Convert from zero based to one based
  *rank = slu_rank;

  // Sort partitions into assending order (required!)
  std::sort(&(*row_perm)[0] , &(*row_perm)[0] + (*rank) );
  std::sort(&(*col_perm)[0] , &(*col_perm)[0] + (*rank) );
  if( *rank < m )
    std::sort(&(*row_perm)[0] + (*rank)	, &(*row_perm)[0] + m );
  if( *rank < n )
    std::sort(&(*col_perm)[0] + (*rank)	, &(*col_perm)[0] + n );

}

void DirectSparseSolverSuperLU::imp_factor(
  const AbstractLinAlgPack::MatrixConvertToSparse   &A
  ,const FactorizationStructure                   &fact_struc
  ,FactorizationNonzeros                          *fact_nonzeros
  ,std::ostream                                   *out
  )
{
  namespace mmp = MemMngPack;
  using Teuchos::dyn_cast;
  typedef MatrixConvertToSparse MCTS;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  if(out)
    *out << "\nUsing SuperLU to refactor the basis matrix ...\n";

  // Get the concrete factorization and nonzeros objects
  const FactorizationStructureSuperLU
    &fs = dyn_cast<const FactorizationStructureSuperLU>(fact_struc);
  FactorizationNonzerosSuperLU
    &fn = dyn_cast<FactorizationNonzerosSuperLU>(*fact_nonzeros);

  // Allocate new storage if not done so already
  TEUCHOS_TEST_FOR_EXCEPTION(
    !fs.fact_struct_.get(), std::logic_error
    ,"DirectSparseSolverSuperLU::imp_factor(...): Error, the factorization sturcture must "
    "have already been computed!"
    );
  if(!fn.fact_nonzeros_.get())
    fn.fact_nonzeros_ = SuperLUPack::SuperLUSolver::create_fact_nonzeros();

  // Get the dimensions of things.
  const index_type
    m = A.rows(),
    n = A.cols(),
    nz = A.num_nonzeros( MCTS::EXTRACT_FULL_MATRIX ,MCTS::ELEMENTS_ALLOW_DUPLICATES_SUM );

  // Validate input
  TEUCHOS_TEST_FOR_EXCEPTION(
    n <= 0 || m <= 0 || m > n, std::invalid_argument
    ,"DirectSparseSolverSuperLU::imp_factor(...) : Error!" );

  // Extract the matrix in coordinate format
  Workspace<value_type>   a_val(wss,nz);
  Workspace<index_type>   a_row_i(wss,nz);
  Workspace<index_type>   a_col_j(wss,nz);
  A.coor_extract_nonzeros(
    MCTS::EXTRACT_FULL_MATRIX
    ,MCTS::ELEMENTS_ALLOW_DUPLICATES_SUM
    ,nz
    ,&a_val[0]
    ,nz
    ,&a_row_i[0]
    ,&a_col_j[0]
    );

  //
  // Convert to compressed sparse row format (which is compressed sparse
  // column of the transposed matrix which will be passed to
  // SuperLU for factorization).
  //

  Workspace<double>   acsr_val(wss,nz);
  Workspace<int>      acsr_col_j(wss,nz);    // Zero-based for SuperLU
  Workspace<int>      acsr_row_ptr(wss,m+1);
  
  convet_to_csr(
    n,m,nz,&a_val[0],&a_row_i[0],&a_col_j[0]
    ,&acsr_val[0],&acsr_col_j[0],&acsr_row_ptr[0]
    );

  //
  // Have SuperLU factor this matrix.
  //
  // SuperLU works with the transpose of the matrix
  // That DirectSparseSolver works with.
  //

  fs.superlu_solver_->factor(
    n                         // m
    ,m                        // n
    ,nz                       // nz
    ,&acsr_val[0]             // a_val[]
    ,&acsr_col_j[0]           // a_row_i[]
    ,&acsr_row_ptr[0]         // a_col_ptr[]
    ,*fs.fact_struct_         // fact_struct
    ,fn.fact_nonzeros_.get()  // fact_nonzeros
    );

}

// private

}	// end namespace AbstractLinAlgPack 

#endif // SPARSE_SOLVER_PACK_USE_SUPERLU
