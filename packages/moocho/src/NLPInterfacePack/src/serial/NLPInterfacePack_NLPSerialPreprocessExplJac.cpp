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

#include <assert.h>

#include <typeinfo>
#include <algorithm>

#include "NLPInterfacePack_NLPSerialPreprocessExplJac.hpp"
#include "AbstractLinAlgPack_MatrixPermAggr.hpp"
#include "AbstractLinAlgPack_BasisSystemFactory.hpp"
#include "AbstractLinAlgPack_MatrixComposite.hpp"
#include "AbstractLinAlgPack_MatrixSparseCOORSerial.hpp"
#include "AbstractLinAlgPack_PermutationSerial.hpp"
#include "AbstractLinAlgPack_VectorDenseEncap.hpp"
#include "DenseLinAlgPack_DVectorOp.hpp"
#include "DenseLinAlgPack_IVector.hpp"
#include "DenseLinAlgPack_PermVecMat.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "OptionsFromStreamPack_OptionsFromStream.hpp"

namespace NLPInterfacePack {

// NLPSerialPreprocessExplJac

// Constructors / initializers

NLPSerialPreprocessExplJac::NLPSerialPreprocessExplJac(
  const basis_sys_fcty_ptr_t  &basis_sys_fcty
  ,const factory_mat_ptr_t    &factory_Gc_full
  )
  :initialized_(false),test_setup_(false)
{
  this->set_basis_sys_fcty(basis_sys_fcty);
  this->set_factory_Gc_full(factory_Gc_full);
}

void NLPSerialPreprocessExplJac::set_factory_Gc_full(
  const factory_mat_ptr_t     &factory_Gc_full
  )
{
  if(factory_Gc_full.get())
    factory_Gc_full_ = factory_Gc_full;
  else 
    factory_Gc_full_ = Teuchos::rcp(
      new Teuchos::AbstractFactoryStd<MatrixOp,MatrixSparseCOORSerial>() );
  factory_Gc_ = Teuchos::rcp( new Teuchos::AbstractFactoryStd<MatrixOp,MatrixPermAggr>() );
}

// Overridden public members from NLP

void NLPSerialPreprocessExplJac::set_options( const options_ptr_t& options )
{
  options_ = options;
}

const NLP::options_ptr_t&
NLPSerialPreprocessExplJac::get_options() const
{
  return options_;
}

void NLPSerialPreprocessExplJac::initialize(bool test_setup)
{
  namespace mmp = MemMngPack;

  test_setup_ = test_setup;

  if( initialized_  && !imp_nlp_has_changed() ) {
    // The subclass NLP has not changed so we can just
    // slip this preprocessing.
    NLPFirstOrder::initialize(test_setup);
    NLPSerialPreprocess::initialize(test_setup);  // Some duplication but who cares!
    return;
  }

  // Initialize the base object first
  NLPFirstOrder::initialize(test_setup);
  NLPSerialPreprocess::initialize(test_setup);  // Some duplication but who cares!

  // Initialize the storage for the intermediate quanities
  Gc_nz_orig_ = imp_Gc_nz_orig();         // Get the estimated number of nonzeros in Gc
  Gc_val_orig_.resize(Gc_nz_orig_);
  Gc_ivect_orig_.resize(Gc_nz_orig_);
  Gc_jvect_orig_.resize(Gc_nz_orig_);
  Gh_nz_orig_ = imp_Gh_nz_orig();			// Get the estimated number of nonzeros in Gh
  Gh_val_orig_.resize(Gh_nz_orig_);
  Gh_ivect_orig_.resize(Gh_nz_orig_);
  Gh_jvect_orig_.resize(Gh_nz_orig_);

  Gc_perm_new_basis_updated_ = false;

  // If you get here then the initialization went Ok.
  initialized_ = true;
}

bool NLPSerialPreprocessExplJac::is_initialized() const {
  return initialized_;
}

// Overridden public members from NLPFirstOrder

const NLPFirstOrder::mat_fcty_ptr_t
NLPSerialPreprocessExplJac::factory_Gc() const
{
  return factory_Gc_;
}

const NLPFirstOrder::basis_sys_ptr_t
NLPSerialPreprocessExplJac::basis_sys() const
{
  BasisSystemFactory &fcty = const_cast<NLPSerialPreprocessExplJac*>(this)->basis_sys_fcty();
  fcty.set_options(options_);
  return fcty.create();
}

void NLPSerialPreprocessExplJac::set_Gc(MatrixOp* Gc)
{
  using Teuchos::dyn_cast;
  assert_initialized();
  if( Gc != NULL ) {
    dyn_cast<MatrixPermAggr>(*Gc); // With throw exception if not correct type!
  }
  NLPFirstOrder::set_Gc(Gc);
}

// Overridden public members from NLPVarReductPerm

bool NLPSerialPreprocessExplJac::get_next_basis(
  Permutation*  P_var,   Range1D* var_dep
  ,Permutation* P_equ,   Range1D* equ_decomp
  )
{
  const bool new_basis = NLPSerialPreprocess::get_next_basis(
    P_var,var_dep,P_equ,equ_decomp
    );
  if( new_basis ) {
    Gc_perm_new_basis_updated_ = false;
  }
  return new_basis;
}

void NLPSerialPreprocessExplJac::set_basis(
  const Permutation   &P_var,   const Range1D  &var_dep
  ,const Permutation  *P_equ,   const Range1D  *equ_decomp
  )
{
  NLPSerialPreprocess::set_basis(
    P_var,var_dep,P_equ,equ_decomp
    );
  Gc_perm_new_basis_updated_ = false;
}

// Overridden protected members from NLPFirstOrder

void NLPSerialPreprocessExplJac::imp_calc_Gc(
  const Vector& x, bool newx
  ,const FirstOrderInfo& first_order_info
  ) const
{
  namespace mmp = MemMngPack;
  using Teuchos::dyn_cast;

  assert_initialized();

  const Range1D
    var_dep      = this->var_dep(),
    equ_decomp   = this->equ_decomp();
  // Get the dimensions of the original NLP
  const size_type
    n       = this->n(),
    n_orig  = this->imp_n_orig(),
    m_orig  = this->imp_m_orig(),
    mI_orig = this->imp_mI_orig();
  // Get the dimensions of the full matrices
  const size_type
    n_full  = n_orig + mI_orig,
    m_full  = m_orig + mI_orig;
  // Get the number of columns in the matrix being constructed here
  const size_type
    num_cols = m_full;

  //
  // Get references to the constituent objects
  //

  // Get the concrete type for the Jacobian matrix
  MatrixPermAggr
    &G_aggr = dyn_cast<MatrixPermAggr>( *first_order_info.Gc );
  // Get smart pointers to the constituent members
  Teuchos::RCP<MatrixOp>
    G_full = Teuchos::rcp_const_cast<MatrixOp>( G_aggr.mat_orig() );
  Teuchos::RCP<PermutationSerial>
    P_row = Teuchos::rcp_dynamic_cast<PermutationSerial>(
      Teuchos::rcp_const_cast<Permutation>( G_aggr.row_perm() ) );  // variable permutation
  Teuchos::RCP<PermutationSerial>
    P_col = Teuchos::rcp_dynamic_cast<PermutationSerial>(
      Teuchos::rcp_const_cast<Permutation>( G_aggr.col_perm() ) );  // constraint permutation
  Teuchos::RCP<const MatrixOp>
    G_perm = G_aggr.mat_perm();
  // Remove references to G_full, G_perm, P_row and P_col.
  G_aggr.set_uninitialized();
  // Allocate the original matrix object if not done so yet
  if( G_full.get() == NULL || G_full.total_count() > 1 )
    G_full = factory_Gc_full_->create();
  // Get reference to the MatrixLoadSparseElements interface
  MatrixLoadSparseElements
    &G_lse = dyn_cast<MatrixLoadSparseElements>(*G_full);

  //
  // Calcuate the full explicit Jacobian
  //

  set_x_full( VectorDenseEncap(x)(), newx, &x_full() );
  if( m_orig )
    imp_calc_Gc_orig( x_full(), newx, first_order_expl_info() );
  if( mI_orig )
    imp_calc_Gh_orig( x_full(), newx, first_order_expl_info() );

  // Now get the actual number of nonzeros
  const size_type nz_full
    = Gc_nz_orig_ + Gh_nz_orig_ + mI_orig;  // Gc_orig, Gh_orig, -I

  // Determine if we need to set the structure and the nonzeros or just the nonzero values
  const bool load_struct = (G_lse.nz() == 0);

  size_type G_nz_previous;
  if( load_struct ) {
    G_lse.reinitialize(n,num_cols,nz_full); // The actual number of nonzeros will be minus the fixed variables
  }
  else {
    G_nz_previous = G_lse.nz();
    G_lse.reset_to_load_values();           // Use row and column indexes already set (better be same insert order!)
  }

  //
  // Load the matrix entries where we remove variables fixed by bounds
  //

  // Get pointers to buffers to fill with nonzero entries
  value_type			*val    = NULL;
   index_type			*ivect  = NULL,
                    *jvect  = NULL;
  G_lse.get_load_nonzeros_buffers(
    nz_full // We may actually load less
    ,&val
    ,load_struct ? &ivect : NULL
    ,load_struct ? &jvect : NULL
    );
  // Pointers to the full COOR matrix just updated
  const value_type      *val_orig     = NULL;
  const value_type      *val_orig_end = NULL;
   const index_type      *ivect_orig   = NULL;
  const index_type      *jvect_orig   = NULL;

  index_type nz = 0;
  if( m_orig ) {
    // Load entries for Gc_orig
    val_orig		= &Gc_val_orig_[0];
    val_orig_end	= val_orig + Gc_nz_orig_;
    ivect_orig		= &Gc_ivect_orig_[0];
    jvect_orig		= &Gc_jvect_orig_[0];
    imp_fill_jacobian_entries(
      n, n_full, load_struct
      ,0 // column offset
      ,val_orig, val_orig_end, ivect_orig, jvect_orig
      ,&nz // This will be incremented
      ,val, ivect, jvect
      );
  }
  if( mI_orig > 0 ) {
    // Load entires for Gc_orig and -I
    val_orig		= &Gh_val_orig_[0];
    val_orig_end	= val_orig + Gh_nz_orig_;
    ivect_orig		= &Gh_ivect_orig_[0];
    jvect_orig		= &Gh_jvect_orig_[0];
    imp_fill_jacobian_entries(
      n, n_full, load_struct
      ,m_orig // column offset (i.e. [ Gc_orig, Gh_orig ] )
      ,val_orig, val_orig_end, ivect_orig, jvect_orig
      ,&nz // This will be incremented
      ,val + nz, ivect + nz, jvect + nz
      );
    // -I
    value_type         *val_itr    = val   + nz;
    index_type         *ivect_itr  = ivect + nz;
    index_type         *jvect_itr  = jvect + nz;
    const IVector& var_full_to_remove_fixed = this->var_full_to_remove_fixed();
    if( load_struct ) {
      // Fill values and i and j
      for( index_type k = 1; k <= mI_orig; ++k ) {
        size_type var_idx = var_full_to_remove_fixed(n_orig+k); // Knows about slacks
#ifdef TEUCHOS_DEBUG
        TEUCHOS_TEST_FOR_EXCEPT( !(  0 < var_idx && var_idx <= n_full  ) );
#endif
        if(var_idx <= n) {
          // This is not a fixed variable
          *val_itr++ = -1.0;
          *ivect_itr++ = var_idx;
          *jvect_itr++ = m_orig + k; // (i.e. [ 0,  -I ] )
          ++nz;
        }
      }
    }
    else {
      // Just fill values
      for( index_type k = 1; k <= mI_orig; ++k ) {
        size_type var_idx = var_full_to_remove_fixed(n_orig+k); // Knows about slacks
#ifdef TEUCHOS_DEBUG
        TEUCHOS_TEST_FOR_EXCEPT( !(  0 < var_idx && var_idx <= n_full  ) );
#endif
        if(var_idx <= n) {
          // This is not a fixed variable
          *val_itr++ = -1.0;
          ++nz;
        }
      }
    }
  }

  if( !load_struct ) {
    // Check that the number of nonzeros added matches the number of nonzeros in G
    TEUCHOS_TEST_FOR_EXCEPTION(
      G_nz_previous != nz, std::runtime_error
      ,"NLPSerialPreprocessExplJac::imp_calc_Gc(...): Error, "
      "The number of added nonzeros does not match the number of nonzeros "
      "in the previous matrix load!." );
  }
  
  // Commit the nonzeros
  G_lse.commit_load_nonzeros_buffers(
    nz  // The actual number of nonzeros to set
    ,&val
    ,load_struct ? &ivect : NULL
    ,load_struct ? &jvect : NULL
    );
  G_lse.finish_construction(test_setup_);

  //
  // Setup permuted view
  //

  // Setup row (variable) permutation
  if( P_row.get() == NULL || P_col.total_count() > 1 )
      P_row = Teuchos::rcp(new PermutationSerial());
  Teuchos::RCP<IVector>        var_perm;
  if( P_row->perm().get() == NULL )  var_perm = Teuchos::rcp(new IVector(n_full));
  else                               var_perm = Teuchos::rcp_const_cast<IVector>(P_row->perm());
  *var_perm = this->var_perm();
  P_row->initialize(var_perm,Teuchos::null);
  // Setup column (constraint) permutation
  if( P_col.get() == NULL || P_col.total_count() > 1 )
      P_col = Teuchos::rcp(new PermutationSerial());
  Teuchos::RCP<IVector>        con_perm;
  if( P_col->perm().get() == NULL )  con_perm = Teuchos::rcp(new IVector(m_full));
  else                               con_perm = Teuchos::rcp_const_cast<IVector>(P_col->perm());
  *con_perm = this->equ_perm();
  P_col->initialize(con_perm,Teuchos::null);
  // Setup G_perm
  int num_row_part, num_col_part;
  index_type row_part[3], col_part[3];
  if(var_dep.size()) {
    num_row_part = 2;
    row_part[0] = 1;
    row_part[1] = (var_dep.lbound() == 1 ? var_dep.ubound()+1 : var_dep.lbound());
    row_part[2] = n+1;
  }
  else {
    num_row_part = 1;
    row_part[0] = 1;
    row_part[1] = n+1;
  }
  if(equ_decomp.size()) {
    num_col_part = 2;
    col_part[0] = 1;
    col_part[1] = (equ_decomp.lbound() == 1 ? equ_decomp.ubound()+1 : equ_decomp.lbound());
    col_part[2] = m_full+1;
  }
  else {
    num_col_part = 1;
    col_part[0] = 1;
    col_part[1] = m_full+1;
  }
  if( G_perm.get() == NULL || !Gc_perm_new_basis_updated_ ) {
    G_perm = G_full->perm_view(
      P_row.get(),row_part,num_row_part
      ,P_col.get(),col_part,num_col_part
      );
  }
  else {
    G_perm = G_full->perm_view_update(
      P_row.get(),row_part,num_row_part
      ,P_col.get(),col_part,num_col_part
      ,G_perm
      );
  }
  Gc_perm_new_basis_updated_ = true;

  //
  // Reinitialize the aggregate matrix object.
  //

  G_aggr.initialize(G_full,P_row,P_col,G_perm);
}

// protected members

void NLPSerialPreprocessExplJac::assert_initialized() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !initialized_, UnInitialized
    ,"NLPSerialPreprocessExplJac : The nlp has not been initialized yet" );
}

// private

void NLPSerialPreprocessExplJac::imp_fill_jacobian_entries(
  size_type           n
  ,size_type          n_full
  ,bool               load_struct
  ,const index_type   col_offset
  ,const value_type   *val_orig
  ,const value_type   *val_orig_end
  ,const index_type   *ivect_orig
  ,const index_type   *jvect_orig
  ,index_type         *nz
  ,value_type         *val_itr
  ,index_type         *ivect_itr
  ,index_type         *jvect_itr
  ) const
{
  const IVector& var_full_to_remove_fixed = this->var_full_to_remove_fixed();
  if( load_struct ) {
    // Fill values and i and j
    for( ; val_orig != val_orig_end ; ++val_orig, ++ivect_orig, ++jvect_orig) {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPT( !(  0 <= *ivect_orig && *ivect_orig <= n_full  ) );
#endif
      size_type var_idx = var_full_to_remove_fixed(*ivect_orig);
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPT( !(  0 < var_idx && var_idx <= n_full  ) );
#endif
      if(var_idx <= n) {
        // This is not a fixed variable
        *val_itr++ = *val_orig;
        // Also fill the row and column indices
        *ivect_itr++ = var_idx;
        *jvect_itr++ = col_offset + (*jvect_orig);
        ++(*nz);
      }
    }
  }
  else {
    // Just fill values
    for( ; val_orig != val_orig_end ; ++val_orig, ++ivect_orig) {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPT( !(  0 <= *ivect_orig && *ivect_orig <= n_full  ) );
#endif
      size_type var_idx = var_full_to_remove_fixed(*ivect_orig);
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPT( !(  0 < var_idx && var_idx <= n_full  ) );
#endif
      if(var_idx <= n) {
        // This is not a fixed variable
        *val_itr++ = *val_orig;
        ++(*nz);
      }
    }
  }
}

} // end namespace NLPInterfacePack
