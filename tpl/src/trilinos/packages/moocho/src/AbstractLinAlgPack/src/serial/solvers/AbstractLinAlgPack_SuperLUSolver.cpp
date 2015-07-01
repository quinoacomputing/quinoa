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
#include <valarray>
#include <stdexcept>

#include "AbstractLinAlgPack_SuperLUSolver.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_Assert.hpp"

// SuperLU
#include "dsp_defs.h"
#include "util.h"

namespace {

// Static SuperLU stuff

int local_panel_size  = 0;
int local_relax       = 0;

class StaticSuperLUInit {
public:
  StaticSuperLUInit()
    {
      local_panel_size = sp_ienv(1);
      local_relax      = sp_ienv(2);
      StatInit(local_panel_size,local_relax);
    }
  ~StaticSuperLUInit()
    {
      StatFree();
    }
};

StaticSuperLUInit static_super_lu_init; // Will be created early and destroyed late!

// ToDo: RAB: 2002/10/14: We must find a better way to work with
// SuperLU than this.  It should not be too hard
// to do better in the future.

// A cast to const is needed because the standard does not return a reference from
// valarray<>::operator[]() const.
template <class T>
std::valarray<T>& cva(const std::valarray<T>& va )
{
  return const_cast<std::valarray<T>&>(va);
}

} // end namespace

namespace SuperLUPack {

class SuperLUSolverImpl;

/** \brief Implementation of SuperLUSolver.
 *
 * ToDo: Finish documentation!
 */
class SuperLUSolverImpl : public SuperLUSolver {
public:

  /** @name Public Types */
  //@{

  /** \brief . */
  class FactorizationStructureImpl : public FactorizationStructure {
  public:
    friend class SuperLUSolverImpl;
  private:
    int                   rank_;        // For square basis
    int                   nz_;          // ...
    std::valarray<int>    perm_r_;      // ...
    std::valarray<int>    perm_c_;      // ...
    std::valarray<int>    etree_;       // ...
    int                   m_orig_;      // For original rectangular matrix
    int                   n_orig_;      // ...
    int                   nz_orig_;     // ...
    std::valarray<int>    perm_r_orig_; // ...
    std::valarray<int>    perm_c_orig_; // ...
  };

  /** \brief . */
  class FactorizationNonzerosImpl : public FactorizationNonzeros {
  public:
    friend class SuperLUSolverImpl;
  private:
    SuperMatrix   L_;
    SuperMatrix   U_;
  };

  //@}

  /** @name Overridden from SuperLUSolver */
  //@{

  /** \brief . */
  void analyze_and_factor(
    int                         m
    ,int                        n
    ,int                        nz
    ,const double               a_val[]
    ,const int                  a_row_i[]
    ,const int                  a_col_ptr[]
    ,FactorizationStructure     *fact_struct
    ,FactorizationNonzeros      *fact_nonzeros
    ,int                        row_perm[]
    ,int                        col_perm[]
    ,int                        *rank
    );
  /** \brief . */
  void factor(
    int                             m
    ,int                            n
    ,int                            nz
    ,const double                   a_val[]
    ,const int                      a_row_i[]
    ,const int                      a_col_ptr[]
    ,const FactorizationStructure   &fact_struct
    ,FactorizationNonzeros          *fact_nonzeros
    );
  /** \brief . */
  void solve(
    const FactorizationStructure    &fact_struct
    ,const FactorizationNonzeros    &fact_nonzeros
    ,bool                           transp
    ,int                            n
    ,int                            nrhs
    ,double                         rhs[]
    ,int                            ldrhs
    ) const;

  //@}

private:

  /** \brief . */
  void copy_basis_nonzeros(
    int                             m_orig
    ,int                            n_orig
    ,int                            nz_orig
    ,const double                   a_orig_val[]
    ,const int                      a_orig_row_i[]
    ,const int                      a_orig_col_ptr[]
    ,const int                      a_orig_perm_r[]
    ,const int                      a_orig_perm_c[]
    ,const int                      rank
    ,double                         b_val[]
    ,int                            b_row_i[]
    ,int                            b_col_ptr[]
    ,int                            *b_nz
    ) const;

}; // end class SuperLUSolver

//
// SuperLUSolver
//

Teuchos::RCP<SuperLUSolver>
SuperLUSolver::create_solver()
{
  return Teuchos::rcp(new SuperLUSolverImpl());
}

Teuchos::RCP<SuperLUSolver::FactorizationStructure>
SuperLUSolver::create_fact_struct()
{
  return Teuchos::rcp(new SuperLUSolverImpl::FactorizationStructureImpl());
}

Teuchos::RCP<SuperLUSolver::FactorizationNonzeros>
SuperLUSolver::create_fact_nonzeros()
{
  return Teuchos::rcp(new SuperLUSolverImpl::FactorizationNonzerosImpl());
}

//
// SuperLUSolverImp
//

// Overridden from SuperLUSolver

void SuperLUSolverImpl::analyze_and_factor(
  int                         m
  ,int                        n
  ,int                        nz
  ,const double               a_val[]
  ,const int                  a_row_i[]
  ,const int                  a_col_ptr[]
  ,FactorizationStructure     *fact_struct
  ,FactorizationNonzeros      *fact_nonzeros
  ,int                        perm_r[]
  ,int                        perm_c[]
  ,int                        *rank
  )
{
  using Teuchos::dyn_cast;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  FactorizationStructureImpl
    &fs = dyn_cast<FactorizationStructureImpl>(*fact_struct);
  FactorizationNonzerosImpl
    &fn = dyn_cast<FactorizationNonzerosImpl>(*fact_nonzeros);

  char refact[] = "N";

  // Resize storage.
  // Note: if this function was called recursively for m>n on the last call
  // then m_orig, n_orig etc. will already be set and should not be
  // disturbed.
  fs.rank_ = n; // Assume this for now
  fs.nz_   = nz;
  fs.perm_r_.resize(m);
  fs.perm_c_.resize(n);
  fs.etree_.resize(n);

    // Create matrix A in the format expected by SuperLU
  SuperMatrix A;
  dCreate_CompCol_Matrix(
    &A, m, n, nz
    ,const_cast<double*>(a_val)
    ,const_cast<int*>(a_row_i)
    ,const_cast<int*>(a_col_ptr)
    ,NC, D_, GE
    );

  // Get the columm permutations
  int permc_spec = 0; // ToDo: Make this an external parameter
  get_perm_c(permc_spec, &A, &fs.perm_c_[0]);

  // Permute the columns of the matrix
  SuperMatrix AC;
  sp_preorder(refact,&A,&fs.perm_c_[0],&fs.etree_[0],&AC);

  int info = -1;
  dgstrf(
    refact
    ,&AC  
    ,1.0    /* diag_pivot_thresh */
    ,0.0    /* drop_tol */
    ,local_relax
    ,local_panel_size
    ,&fs.etree_[0]
    ,NULL   /* work */
    ,0      /* lwork */
    ,&fs.perm_r_[0]
    ,&fs.perm_c_[0]
    ,&fn.L_
    ,&fn.U_
    ,&info
    );

  TEUCHOS_TEST_FOR_EXCEPTION(
    info != 0, std::runtime_error
    ,"SuperLUSolverImpl::analyze_and_factor(...): Error, dgstrf(...) returned info = " << info
    );

  std::copy( &fs.perm_r_[0], &fs.perm_r_[0] + m, perm_r );
  std::copy( &fs.perm_c_[0], &fs.perm_c_[0] + n, perm_c );
  *rank = n; // We must assume this until I can figure out a way to do better!

  if(m > n) {
    // Now we must refactor the basis by only passing in the elements for the basis
    // determined by SuperLU.  This is wasteful but it is the easiest thing to do
    // for now.
    fs.rank_        = *rank;
    fs.m_orig_      = m;
    fs.n_orig_      = n;
    fs.nz_orig_     = nz;
    fs.perm_r_orig_ = fs.perm_r_;
    fs.perm_c_orig_ = fs.perm_c_;
    // Copy the nonzeros for the sqare factor into new storage
    Workspace<double>       b_val(wss,nz);
    Workspace<int>          b_row_i(wss,nz);
    Workspace<int>          b_col_ptr(wss,n+1);
    copy_basis_nonzeros(
      m,n,nz,a_val,a_row_i,a_col_ptr
      ,&fs.perm_r_orig_[0],&fs.perm_c_orig_[0],fs.rank_
      ,&b_val[0],&b_row_i[0],&b_col_ptr[0]
      ,&fs.nz_
      );
    // Analyze and factor the new matrix
    int b_rank = -1;
    analyze_and_factor(
      fs.rank_, fs.rank_, fs.nz_, &b_val[0], &b_row_i[0], &b_col_ptr[0]
      ,fact_struct, fact_nonzeros
      ,&fs.perm_r_[0], &fs.perm_c_[0], &b_rank
      );
    TEUCHOS_TEST_FOR_EXCEPTION(
      (b_rank != *rank), std::runtime_error
      ,"SuperLUSolverImpl::analyze_and_factor(...): Error, the rank determined by "
      "the factorization of the rectangular " << m << " x " << n << " matrix of "
      << (*rank) << " is not the same as the refactorization of the basis returned as "
      << b_rank << "!"
      );
  }
  else {
    fs.m_orig_  = m;
    fs.n_orig_  = n;
    fs.nz_orig_ = nz;
  }
}

void SuperLUSolverImpl::factor(
  int                             m
  ,int                            n
  ,int                            nz
  ,const double                   a_val[]
  ,const int                      a_row_i[]
  ,const int                      a_col_ptr[]
  ,const FactorizationStructure   &fact_struct
  ,FactorizationNonzeros          *fact_nonzeros
  )
{
  using Teuchos::dyn_cast;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  const FactorizationStructureImpl
    &fs = dyn_cast<const FactorizationStructureImpl>(fact_struct);
  FactorizationNonzerosImpl
    &fn = dyn_cast<FactorizationNonzerosImpl>(*fact_nonzeros);

  char refact[] = "Y";

  // Copy the nonzeros for the sqare factor into new storage
  Workspace<double>       b_val(wss,fs.nz_);
  Workspace<int>          b_row_i(wss,fs.nz_);
  Workspace<int>          b_col_ptr(wss,fs.rank_+1);
  if(fs.m_orig_ > fs.n_orig_) {
    int b_nz = -1;
    copy_basis_nonzeros(
      m,n,nz,a_val,a_row_i,a_col_ptr
      ,&cva(fs.perm_r_orig_)[0],&cva(fs.perm_c_orig_)[0],fs.rank_
      ,&b_val[0],&b_row_i[0],&b_col_ptr[0]
      ,&b_nz
      );
    TEUCHOS_TEST_FOR_EXCEPTION(
      (b_nz != fs.nz_), std::runtime_error
      ,"SuperLUSolverImpl::factor(...): Error!"
      );
  }
  else {
    std::copy( a_val,     a_val     + nz,  &b_val[0]     );
    std::copy( a_row_i,   a_row_i   + nz,  &b_row_i[0]   );
    std::copy( a_col_ptr, a_col_ptr + n+1, &b_col_ptr[0] );
  }

    // Create matrix A in the format expected by SuperLU
  SuperMatrix A;
  dCreate_CompCol_Matrix(
    &A, fs.rank_, fs.rank_, fs.nz_
    ,&b_val[0]
    ,&b_row_i[0]
    ,&b_col_ptr[0]
    ,NC, D_, GE
    );

  // Permute the columns
  SuperMatrix AC;
  sp_preorder(
    refact,&A
    ,&cva(fs.perm_c_)[0]
    ,&cva(fs.etree_)[0]
    ,&AC
    );

  int info = -1;
  dgstrf(
    refact
    ,&AC  
    ,1.0    /* diag_pivot_thresh */
    ,0.0    /* drop_tol */
    ,local_relax
    ,local_panel_size
    ,const_cast<int*>(&cva(fs.etree_)[0])
    ,NULL   /* work */
    ,0      /* lwork */
    ,&cva(fs.perm_r_)[0]
    ,&cva(fs.perm_c_)[0]
    ,&fn.L_
    ,&fn.U_
    ,&info
    );

  TEUCHOS_TEST_FOR_EXCEPTION(
    info != 0, std::runtime_error
    ,"SuperLUSolverImpl::factor(...): Error, dgstrf(...) returned info = " << info
    );

}

void SuperLUSolverImpl::solve(
  const FactorizationStructure    &fact_struct
  ,const FactorizationNonzeros    &fact_nonzeros
  ,bool                           transp
  ,int                            n
  ,int                            nrhs
  ,double                         rhs[]
  ,int                            ldrhs
  ) const
{

  using Teuchos::dyn_cast;

  const FactorizationStructureImpl
    &fs = dyn_cast<const FactorizationStructureImpl>(fact_struct);
  const FactorizationNonzerosImpl
    &fn = dyn_cast<const FactorizationNonzerosImpl>(fact_nonzeros);

  TEUCHOS_TEST_FOR_EXCEPTION(
    n != fs.rank_, std::runtime_error
    ,"SuperLUSolverImpl::solve(...): Error, the dimmensions n = " << n << " and fs.rank = " << fs.rank_
    << " do not match up!"
    );

  SuperMatrix B;
    dCreate_Dense_Matrix(&B, n, nrhs, rhs, ldrhs, DN, D_, GE);

  char transc[1];
  transc[0] = ( transp ? 'T' : 'N' );

  int info = -1;
    dgstrs(
    transc
    ,const_cast<SuperMatrix*>(&fn.L_)
    ,const_cast<SuperMatrix*>(&fn.U_)
    ,&cva(fs.perm_r_)[0]
    ,&cva(fs.perm_c_)[0]
    ,&B, &info
    );

  TEUCHOS_TEST_FOR_EXCEPTION(
    info != 0, std::runtime_error
    ,"SuperLUSolverImpl::solve(...): Error, dgssv(...) returned info = " << info
    );

}

// private

void SuperLUSolverImpl::copy_basis_nonzeros(
  int                             m_orig
  ,int                            n_orig
  ,int                            nz_orig
  ,const double                   a_orig_val[]
  ,const int                      a_orig_row_i[]
  ,const int                      a_orig_col_ptr[]
  ,const int                      a_orig_perm_r[]
  ,const int                      a_orig_perm_c[]
  ,const int                      rank
  ,double                         b_val[]
  ,int                            b_row_i[]
  ,int                            b_col_ptr[]
  ,int                            *b_nz
  ) const
{
  *b_nz = 0;
  b_col_ptr[0] = *b_nz;
  for( int j = 0; j < rank; ++j ) {
    const int col_start_k = a_orig_col_ptr[j];
    const int col_end_k   = a_orig_col_ptr[j+1];
    for( int k = col_start_k; k < col_end_k; ++k ) {
      const int i_orig = a_orig_row_i[k];
      if(i_orig < rank) {
        b_val[*b_nz]     = a_orig_val[k];
        b_row_i[*b_nz]   = a_orig_row_i[k];
        ++(*b_nz);
      }
    }
    b_col_ptr[j+1] = *b_nz;
  }
}

} // end namespace SuperLUPack

#endif // SPARSE_SOLVER_PACK_USE_SUPERLU
