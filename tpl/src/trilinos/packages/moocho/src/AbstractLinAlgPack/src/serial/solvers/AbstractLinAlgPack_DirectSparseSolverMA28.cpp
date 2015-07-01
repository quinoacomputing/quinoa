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

#include <fstream>
#include <algorithm>

#include "Moocho_ConfigDefs.hpp"


#ifdef HAVE_MOOCHO_MA28


#include "AbstractLinAlgPack_DirectSparseSolverMA28.hpp"
#include "AbstractLinAlgPack_MatrixScaling_Strategy.hpp"
#include "AbstractLinAlgPack_VectorDenseEncap.hpp"
#include "DenseLinAlgPack_PermVecMat.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "FortranTypes_f_open_file.hpp"

namespace {
//
template< class T >
inline
T my_max( const T& v1, const T& v2 ) { return v1 > v2 ? v1 : v2; }
// A cast to const is needed because the standard does not return a reference from
// valarray<>::operator[]() const.
template <class T>
std::valarray<T>& cva(const std::valarray<T>& va )
{
  return const_cast<std::valarray<T>&>(va);
}
}

namespace AbstractLinAlgPack {

//
// Implementation of DirectSparseSolver(Imp) interface using MA28.
//
// Here are some little bits of knowledge about MA28 that I need
// to record after many hours of hard work.
//
// * The 1979 paper in ACM TOMS (Vol. 5, No. 1, pages 27), seems
// to suggest that MA28 pivots by column for numerical stability
// but I am not sure about this.
//
// * When factoring a rectangular matrix, you must set 
// LBLOCK = .FALSE. or the row and column permutations
// extracted from IKEEP(:,2) and IKEEP(:,3) respectivly
// are meaningless.
//
// ToDo: Finish this discussion!
//

// ToDo:
// a) Add an option for printing the values of the common
//    block parameters to out or to a file.  This can
//    be used to get a feel for the performance of
//    ma28
// b) Add provisions for an external client to change
//    the control options of MA28.  Most of these are
//    stored as common block variables.

// //////////////////////////////////////////////////
// DirectSparseSolverMA28::FactorizationStructureMA28

DirectSparseSolverMA28::FactorizationStructureMA28::FactorizationStructureMA28()
  :m_(0),n_(0),max_n_(0),nz_(0),licn_(0),lirn_(0)
   ,u_(0.1),scaling_(NO_SCALING)
{}

// //////////////////////////////////////////////////
// DirectSparseSolverMA28::BasisMatrixMA28

// Overridden from BasisMatrixImp

Teuchos::RCP<DirectSparseSolverImp::BasisMatrixImp>
DirectSparseSolverMA28::BasisMatrixMA28::create_matrix() const
{
  return Teuchos::rcp(new BasisMatrixMA28);
}

void DirectSparseSolverMA28::BasisMatrixMA28::V_InvMtV(
  VectorMutable* y, BLAS_Cpp::Transp M_trans, const Vector& x
  ) const 
{
  using Teuchos::dyn_cast;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  size_type k;

  // Get concrete objects
  const FactorizationStructureMA28
    &fs = dyn_cast<const FactorizationStructureMA28>(*this->get_fact_struc());
  const FactorizationNonzerosMA28
    &fn = dyn_cast<const FactorizationNonzerosMA28>(*this->get_fact_nonzeros());

  // Validate input
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    y == NULL, std::invalid_argument
    ,"DirectSparseSolverMA28::BasisMatrixMA28::V_InvMtV(...) : Error! " );
#endif
  const size_type y_dim = y->dim(), x_dim = x.dim();
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    fs.rank_ != y_dim, std::invalid_argument
    ,"DirectSparseSolverMA28::BasisMatrixMA28::V_InvMtV(...) : Error! " );
  TEUCHOS_TEST_FOR_EXCEPTION(
    fs.rank_ != x_dim, std::invalid_argument
    ,"DirectSparseSolverMA28::BasisMatrixMA28::V_InvMtV(...) : Error! " );
#endif

  VectorDenseMutableEncap  yd(*y);
  VectorDenseEncap         xd(x);

  // Allocate workspace memory
  Workspace<value_type>  xfull_s(wss,fs.max_n_,false);
  DVectorSlice                 xfull(&xfull_s[0],xfull_s.size());
  Workspace<value_type>  ws(wss,fs.max_n_,false);
  DVectorSlice                 w(&ws[0],ws.size());

  // Get a context for transpose or no transpose
  const IVector
    &row_perm = (M_trans == BLAS_Cpp::no_trans) ? fs.row_perm_ : fs.col_perm_,
    &col_perm = (M_trans == BLAS_Cpp::no_trans) ? fs.col_perm_ : fs.row_perm_;

  // Copy x into positions in full w
  // Here, the rhs vector is set with only those equations that
  // are part of the nonsingular set.  It is important that the
  // ordering be the same as the original ordering sent to
  // MA28AD().
  xfull = 0.0;
  for( k = 1; k <= x_dim; ++k ) 
    xfull(row_perm(k)) = xd()(k);
  
  // Scale the rhs
  if( fs.matrix_scaling_.get() )
    fs.matrix_scaling_->scale_rhs( M_trans, xfull.raw_ptr() );

  // Solve for the rhs
  FortranTypes::f_int mtype = ( (M_trans == BLAS_Cpp::no_trans) ? 1 : 0 );
  fs.ma28_.ma28cd(
    fs.max_n_, &cva(fn.a_)[0], fs.licn_, &cva(fs.icn_)[0], &cva(fs.ikeep_)[0]
    ,xfull.raw_ptr(), w.raw_ptr(), mtype );

  // Scale the lhs
  if( fs.matrix_scaling_.get() )
    fs.matrix_scaling_->scale_rhs( M_trans, xfull.raw_ptr() );

  // Copy the solution into y
  // Here, the solution vector is set with only those variables that
  // are in the basis.  It is important that the
  // ordering be the same as the original ordering sent to
  // MA28AD().
  for( k = 1; k <= y_dim; ++k )
    yd()(k) = xfull(col_perm(k));
  
}

// //////////////////////////////////////////////////
// DirectSparseSolverMA28

// Constructors/initializers

DirectSparseSolverMA28::DirectSparseSolverMA28(
  value_type          estimated_fillin_ratio
  ,value_type         u
  ,bool               grow
  ,value_type         tol
  ,index_type         nsrch
  ,bool               lbig
  ,bool               print_ma28_outputs
  ,const std::string& output_file_name
  )
  :estimated_fillin_ratio_(estimated_fillin_ratio)
  ,u_(u)
  ,grow_(grow)
  ,tol_(tol)
  ,nsrch_(nsrch)
  ,lbig_(lbig)
  ,print_ma28_outputs_(print_ma28_outputs)
  ,output_file_name_(output_file_name)
  ,file_output_num_(0)
{}

// Overridden from DirectSparseSolver

const DirectSparseSolver::basis_matrix_factory_ptr_t
DirectSparseSolverMA28::basis_matrix_factory() const
{
  namespace mmp = MemMngPack;
  return Teuchos::rcp(new Teuchos::AbstractFactoryStd<BasisMatrix,BasisMatrixMA28>());
}

void DirectSparseSolverMA28::estimated_fillin_ratio(
  value_type estimated_fillin_ratio
  )
{
  estimated_fillin_ratio_ = estimated_fillin_ratio;
}

// Overridden from DirectSparseSolverImp

const Teuchos::RCP<DirectSparseSolver::FactorizationStructure>
DirectSparseSolverMA28::create_fact_struc() const
{
  return Teuchos::rcp(new FactorizationStructureMA28);
}

const Teuchos::RCP<DirectSparseSolverImp::FactorizationNonzeros>
DirectSparseSolverMA28::create_fact_nonzeros() const
{
  return Teuchos::rcp(new FactorizationNonzerosMA28);
}

void DirectSparseSolverMA28::imp_analyze_and_factor(
  const AbstractLinAlgPack::MatrixConvertToSparse   &A
  ,FactorizationStructure                         *fact_struc
  ,FactorizationNonzeros                          *fact_nonzeros
  ,DenseLinAlgPack::IVector                            *row_perm
  ,DenseLinAlgPack::IVector                            *col_perm
  ,size_type                                      *rank
  ,std::ostream                                   *out
  )
{
  using Teuchos::dyn_cast;
  typedef MatrixConvertToSparse MCTS;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  if(out)
    *out << "\nUsing MA28 to analyze and factor a new matrix ...\n";

  // Get the concrete factorization and nonzeros objects
  FactorizationStructureMA28
    &fs = dyn_cast<FactorizationStructureMA28>(*fact_struc);
  FactorizationNonzerosMA28
    &fn = dyn_cast<FactorizationNonzerosMA28>(*fact_nonzeros);

  // Set MA28 parameters
  set_ma28_parameters(&fs);
  
  // Get the dimensions of things.
  const index_type
    m = A.rows(),
    n = A.cols(),
    nz = A.num_nonzeros( MCTS::EXTRACT_FULL_MATRIX ,MCTS::ELEMENTS_ALLOW_DUPLICATES_SUM );

  // Validate input
  TEUCHOS_TEST_FOR_EXCEPTION(
    n <= 0 || m <= 0 || m > n, std::invalid_argument
    ,"DirectSparseSolverMA28::imp_analyze_and_factor(...) : Error!" );

  // Memorize the dimenstions for checks later
  fs.m_ = m; fs.n_ = n; fs.nz_ = nz;
  fs.max_n_ = my_max(fs.m_,fs.n_);

  // By default set licn and ircn equal to estimated_fillin_ratio * nz.
  if( estimated_fillin_ratio_ < 1.0 ) {
    if( out ) *out << "Warning, client set estimated_fillin_ratio = " << estimated_fillin_ratio_
             << " < 1.0.\nSetting estimated_fillin_ratio = 10.0 ...\n";
    estimated_fillin_ratio_ = 10.0;
  }
  if( fs.licn_ < fs.nz_ ) fs.licn_ = static_cast<index_type>(estimated_fillin_ratio_ * fs.nz_);
  if( fs.lirn_ < fs.nz_ ) fs.lirn_ = static_cast<index_type>(estimated_fillin_ratio_ * fs.nz_);

  // Initialize structure storage
  fs.ivect_.resize(fs.nz_); // perminatly stores nz row indexes
  fs.jvect_.resize(fs.nz_); // perminatly stores nz column indexes

  index_type iflag = 0;
  for( int num_fac = 0; num_fac < 5; ++num_fac ) {
    
    // Initialize matrix factorization storage and temporary storage
    fs.icn_.resize(fs.licn_); // first nz entries stores the column indexes
     fn.a_.resize(fs.licn_);
    fs.ikeep_.resize( fs.ma28_.lblock() ? 5*fs.max_n_ :  4*fs.max_n_ + 1 );
    Workspace<index_type>  irn_tmp_(wss,fs.lirn_), iw(wss,8*fs.max_n_);
    Workspace<value_type>  w(wss,fs.max_n_);
    
    // Fill in the matrix information
    A.coor_extract_nonzeros(
      MCTS::EXTRACT_FULL_MATRIX
      ,MCTS::ELEMENTS_ALLOW_DUPLICATES_SUM
      ,fs.nz_
      ,&fn.a_[0]
      ,fs.nz_
      ,&fs.ivect_[0]
      ,&fs.jvect_[0]
      );
    std::copy( &fs.ivect_[0], &fs.ivect_[0] + fs.nz_, &irn_tmp_[0] );
    std::copy( &fs.jvect_[0], &fs.jvect_[0] + fs.nz_, &fs.icn_[0] );
    
    // Scale the matrix
    if( fs.matrix_scaling_.get() )
      fs.matrix_scaling_->scale_matrix(
        fs.m_, fs.n_, fs.nz_, &fs.ivect_[0], &fs.jvect_[0], true
        ,&fn.a_[0]
        );
    
    // Analyze and factor the matrix
    if(out)
      *out << "\nCalling ma28ad(...) ...\n";
    fs.ma28_.ma28ad(
      fs.max_n_, fs.nz_, &fn.a_[0], fs.licn_, &irn_tmp_[0], fs.lirn_, &fs.icn_[0], fs.u_
      ,&fs.ikeep_[0], &iw[0], &w[0], &iflag
      );
    
    if(iflag != 0 && out)
      *out << "\nMA28AD returned iflag = " << iflag << " != 0!\n";

    // Print MA28 outputs
    print_ma28_outputs(true,iflag,fs,&w[0],out);

    if( iflag >= 0 ) break;

    switch( iflag ) {
      case LICN_AND_LIRN_TOO_SMALL:
      case LICN_TOO_SMALL:
      case LICN_FAR_TOO_SMALL:
      case LIRN_TOO_SMALL:
        if(out)
          *out << "\nWarning! iflag = " << iflag << ", LIRN and/or LICN are too small!\n"
             << "Increasing lirn = " << fs.lirn_ << " amd licn = " << fs.licn_
             << " by a factor of 10\n"
             << "(try increasing estimated_fillin_ratio = " << estimated_fillin_ratio_
             << " to a larger value next time)...\n";
        fs.lirn_ = 10 * fs.lirn_;
        fs.licn_ = 10 * fs.licn_;
        break;
    }
  }

  // Check for errors and throw exception if you have to.
  ThrowIFlagException(iflag);

  // Extract the basis matrix selection

  *rank = fs.ma28_.irank();

  row_perm->resize(fs.m_);
  if( *rank < fs.m_ ) {
    index_type
      *row_perm_ikeep = &fs.ikeep_[fs.max_n_],
      *row_perm_itr   = &(*row_perm)(1),
      *row_perm_last  = row_perm_itr + fs.m_;
    for(; row_perm_itr != row_perm_last;)
      *row_perm_itr++ = abs(*row_perm_ikeep++);
    // Sort partitions in assending order (required!)
    std::sort(&(*row_perm)[0]           , &(*row_perm)[0] + (*rank) );
    std::sort(&(*row_perm)[0] + (*rank)	, &(*row_perm)[0] + m       );
  }
  else {
    DenseLinAlgPack::identity_perm( row_perm );
  }

  col_perm->resize(fs.n_);
  if( *rank < fs.n_ ) {
    index_type
      *col_perm_ikeep = &fs.ikeep_[2*fs.max_n_],
      *col_perm_itr   = &(*col_perm)(1),
      *col_perm_last  = col_perm_itr + fs.n_;
    for(; col_perm_itr != col_perm_last;)
      *col_perm_itr++ = abs(*col_perm_ikeep++);
    // Sort partitions in assending order (required!)
    std::sort(&(*col_perm)[0]           , &(*col_perm)[0] + (*rank) );
    std::sort(&(*col_perm)[0] + (*rank)	, &(*col_perm)[0] + n       );
  }
  else {
    DenseLinAlgPack::identity_perm( col_perm );
  }

  // Set internal copy of basis selection
  fs.row_perm_ = *row_perm;
  fs.col_perm_ = *col_perm;
  fs.rank_     = *rank;

}

void DirectSparseSolverMA28::imp_factor(
  const AbstractLinAlgPack::MatrixConvertToSparse   &A
  ,const FactorizationStructure                   &fact_struc
  ,FactorizationNonzeros                          *fact_nonzeros
  ,std::ostream                                   *out
  )
{
  using Teuchos::dyn_cast;
  typedef MatrixConvertToSparse MCTS;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  if(out)
    *out << "\nUsing MA28 to factor a new matrix ...\n";

  // Get the concrete factorization and nonzeros objects
  const FactorizationStructureMA28
    &fs = dyn_cast<const FactorizationStructureMA28>(fact_struc);
  FactorizationNonzerosMA28
    &fn = dyn_cast<FactorizationNonzerosMA28>(*fact_nonzeros);

  // Get the dimensions of things.
  const index_type
    m = A.rows(),
    n = A.cols(),
    nz = A.num_nonzeros( MCTS::EXTRACT_FULL_MATRIX ,MCTS::ELEMENTS_ALLOW_DUPLICATES_SUM );

  // Validate input
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    m != fs.m_ || n != fs.n_ || fs.nz_ != nz, std::invalid_argument
    ,"DirectSparseSolverMA28::imp_factor(...) : Error, "
    "A is not compatible with matrix passed to imp_analyze_and_factor()!" );
#endif

  // Initialize matrix factorization storage and temporary storage
  if(fn.a_.size() < fs.licn_)  fn.a_.resize(fs.licn_);
  Workspace<index_type>   iw(wss,5*fs.max_n_);
  Workspace<value_type>   w(wss,fs.max_n_);

  // Fill in the matrix nonzeros (we already have the structure)
  A.coor_extract_nonzeros(
    MCTS::EXTRACT_FULL_MATRIX
    ,MCTS::ELEMENTS_ALLOW_DUPLICATES_SUM
    ,fs.nz_
    ,&fn.a_[0]
    ,0
    ,NULL
    ,NULL
    );

  // Scale the matrix
  if( fs.matrix_scaling_.get() )
    fs.matrix_scaling_->scale_matrix(
      fs.m_, fs.n_, fs.nz_, &cva(fs.ivect_)[0], &cva(fs.jvect_)[0], false
      ,&fn.a_[0]
      );

  // Factor the matrix
  index_type iflag = 0;
  fs.ma28_.ma28bd(
    fs.max_n_, fs.nz_, &fn.a_[0], fs.licn_, &cva(fs.ivect_)[0], &cva(fs.jvect_)[0], &cva(fs.icn_)[0]
    ,&cva(fs.ikeep_)[0], &iw[0], &w[0], &iflag
    );

  if(iflag != 0 && out)
    *out << "\nMA28BD returned iflag = " << iflag << " != 0!\n";

  // Print MA28 outputs
  print_ma28_outputs(false,iflag,fs,&w[0],out);

  // Check for errors and throw exception if you have to.
  ThrowIFlagException(iflag);

}

// private

void DirectSparseSolverMA28::set_ma28_parameters( FactorizationStructureMA28* fs )
{
  // Set common block parameters
  fs->ma28_.lblock( FortranTypes::F_FALSE ); // Do not permute to block triangular form (*** This is critical!)
  fs->u_ = u_;
  fs->ma28_.grow( grow_ ? FortranTypes::F_TRUE : FortranTypes::F_FALSE );
  fs->ma28_.tol(tol_);
  fs->ma28_.nsrch(nsrch_);
  fs->ma28_.lbig( lbig_ ? FortranTypes::F_TRUE : FortranTypes::F_FALSE );
  // Setup output file
  if( output_file_name_.length() > 0 && fs->ma28_.lp() == 0 ) {
    // Open a fortran file
    index_type iout = 25; // Unique?
    FortranTypes::f_open_file( iout, output_file_name_.c_str() );
    fs->ma28_.mp(iout);
    fs->ma28_.lp(iout);
  }
  else if( output_file_name_.length() == 0 && fs->ma28_.lp() ) {
    fs->ma28_.mp(0);
    fs->ma28_.lp(0);
  } 
}

void DirectSparseSolverMA28::print_ma28_outputs(
  bool                               ma28ad_bd
  ,index_type                        iflag
  ,const FactorizationStructureMA28  &fs
  ,const value_type                  w[]
  ,std::ostream                      *out
  )
{
  if( print_ma28_outputs_ && out ) {
    *out << "\nReturn parameters from MA28 (call number = " << ++file_output_num_ << ")\n";
    if( fs.ma28_.grow() == FortranTypes::F_TRUE )
      *out << "w(1)   = " << w[0] << std::endl;
    *out << "rmin   = " << fs.ma28_.rmin() << std::endl;
    *out << "irncp  = " << fs.ma28_.irncp() << std::endl;
    *out << "icncp  = " << fs.ma28_.icncp() << std::endl;
    *out << "minirn = " << fs.ma28_.minirn() << std::endl;
    *out << "minicn = " << fs.ma28_.minicn() << std::endl;
    *out << "irank  = " << fs.ma28_.irank() << std::endl;
    *out << "themax = " << fs.ma28_.themax() << std::endl;
    if( fs.ma28_.lbig() == FortranTypes::F_TRUE )
      *out << "big    = " << fs.ma28_.big() << std::endl;
    *out << "ndrop  = " << fs.ma28_.ndrop() << std::endl;
    if( iflag >= 0 ) {
      *out << "\nAnalysis:\n"
         << "estimated_fillin_ratio can be reduced to max(minirn,minicn)/nz = "
         << "max(" << fs.ma28_.minirn() << "," << fs.ma28_.minicn() << ")/" << fs.nz_
         << " = " << my_max( fs.ma28_.minirn(), fs.ma28_.minicn() ) / (double)fs.nz_
         << std::endl;
    }
  }
}

void DirectSparseSolverMA28::ThrowIFlagException(index_type iflag)
{
  E_IFlag e_iflag = static_cast<E_IFlag>(iflag);
  const char msg_err_head[] = "DirectSparseSolverMA28::ThrowIFlagException(iflag) : Error";
  switch(e_iflag) {
    case SLOW_ITER_CONV :
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::runtime_error
        ,msg_err_head << ", Convergence to slow" );
    case MAXIT_REACHED :
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::runtime_error
        ,msg_err_head << ", Maximum iterations exceeded");
    case MA28BD_CALLED_WITH_DROPPED :
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error
        ,msg_err_head << ", ma28bd called with elements dropped in ma28ad");
    case DUPLICATE_ELEMENTS :
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, FactorizationFailure
        ,msg_err_head << ", Duplicate elements have been detected");
    case NEW_NONZERO_ELEMENT :
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, FactorizationFailure
        ,msg_err_head << ", A new non-zero element has be passed to ma28bd that was not ot ma28ad");
    case N_OUT_OF_RANGE :
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, FactorizationFailure
        ,msg_err_head << ", 1 <=max(n,m) <= 32767 has been violated");
    case NZ_LE_ZERO :
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error
        ,msg_err_head << ", nz <= 0 has been violated");
    case LICN_LE_NZ :
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error
        ,msg_err_head << ", licn <= nz has been violated");
    case LIRN_LE_NZ :
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error
        ,msg_err_head << ", lirn <= nz has been violated");
    case ERROR_DURRING_BLOCK_TRI :
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, FactorizationFailure
        ,msg_err_head << ", An error has occured durring block triangularization");
    case LICN_AND_LIRN_TOO_SMALL :
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, FactorizationFailure
        ,msg_err_head << ", licn and lirn are to small to hold matrix factorization");
    case LICN_TOO_SMALL :
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, FactorizationFailure
        ,msg_err_head << ", licn is to small to hold matrix factorization");
    case LICN_FAR_TOO_SMALL :
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, FactorizationFailure
        ,msg_err_head << ", licn is to far small to hold matrix factorization");
    case LIRN_TOO_SMALL :
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, FactorizationFailure
        ,msg_err_head << ", lirn is to small to hold matrix factorization");
    case NUMERICALLY_SINGULAR :
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, FactorizationFailure
        ,msg_err_head << ", matrix is numerically singular, see \'abort2\'");
    case STRUCTURALLY_SINGULAR :
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, FactorizationFailure
        ,msg_err_head << ", matrix is structurally singular, see \'abort1\'");
    default:
      return; // We don't throw exceptions for other values of iflag.
  }
}

}	// end namespace AbstractLinAlgPack 


#endif // HAVE_MOOCHO_MA28
