#if 0

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

#include <ostream>
#include <iomanip>

#include "MoochoPack_PBFGS_helpers.hpp"
#include "AbstractLinAlgPack_SortByDescendingAbsValue.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_SpVectorOp.hpp"
#include "MiWorkspacePack.h"

namespace MoochoPack {

bool PBFGSPack::act_set_calmed_down( 
  const ActSetStats         &stats
  ,const value_type         act_set_frac_proj_start
  ,EJournalOutputLevel      olevel
  ,std::ostream             &out
  )
{
  typedef ActSetStats ASS;
  const size_type
    num_active_indep = stats.num_active_indep(),
    num_adds_indep   = stats.num_adds_indep(),
    num_drops_indep  = stats.num_drops_indep();
  const value_type
    frac_same
    = ( num_adds_indep == ASS::NOT_KNOWN || num_drops_indep == ASS::NOT_KNOWN || num_active_indep == 0
      ? 0.0
      : std::_MAX(((double)(num_active_indep)-num_adds_indep-num_drops_indep) / num_active_indep, 0.0 ) );
  const bool act_set_calmed_down = ( num_active_indep > 0 && frac_same >= act_set_frac_proj_start );
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out << "\nnum_active_indep = " << num_active_indep;
    if( num_active_indep ) {
      out	<< "\nmax(num_active_indep-num_adds_indep-num_drops_indep,0)/(num_active_indep) = "
        << "max("<<num_active_indep<<"-"<<num_adds_indep<<"-"<<num_drops_indep<<",0)/("<<num_active_indep<<") = "
        << frac_same;
      if( act_set_calmed_down )
        out << " >= ";
      else
        out << " < ";
      out << "act_set_frac_proj_start = " << act_set_frac_proj_start;
      if( act_set_calmed_down )
        out << "\nThe active set has calmed down enough\n";
      else
        out << "\nThe active set has not calmed down enough\n";
    }
  }
  return act_set_calmed_down;
}

void PBFGSPack::init_i_x_free_sRTsR_sRTyR(
  const SpVectorSlice        &nu_indep
  ,const DVectorSlice         &s
  ,const DVectorSlice         &y
  ,size_type                 *n_pz_R
  ,size_type                 i_x_free[]
  ,value_type                *sRTsR
  ,value_type                *sRTyR
  )
{
  const size_type
    n_pz = nu_indep.size();
  SpVectorSlice::const_iterator
    nu_indep_itr   = nu_indep.begin(),
    nu_indep_end   = nu_indep.end();
  const SpVectorSlice::difference_type
    o = nu_indep.offset();
  *n_pz_R = 0;
  *sRTsR  = 0.0;
  *sRTyR  = 0.0;
  for( size_type i = 1; i <= n_pz; ++i ) {
    if( nu_indep_itr != nu_indep_end && nu_indep_itr->indice() + o == i ) {
      ++nu_indep_itr;
    }
    else {
      ++(*n_pz_R);
      *i_x_free++ = i;
      const value_type s_i = s(i);
      (*sRTsR) += s_i*s_i;
      (*sRTyR) += s_i*y(i);
    }
  }
}

void PBFGSPack::sort_fixed_max_cond_viol(
  const SpVectorSlice        &nu_indep
  ,const DVectorSlice         &s
  ,const DVectorSlice         &y
  ,const DVectorSlice         &B_XX
  ,const value_type          sRTBRRsR
  ,const value_type          sRTyR
  ,value_type                *sXTBXXsX
  ,value_type                *sXTyX
  ,size_type                 l_x_fixed_sorted[]
  )
{
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  const size_type
    n_pz = nu_indep.size();
  
  *sXTBXXsX = 0.0;
  *sXTyX    = 0.0;

  // Initial spare vector so we can sort this stuff
  typedef SpVector::element_type ele_t;
  Workspace<ele_t> sort_array(wss,nu_indep.nz());
  {
    SpVectorSlice::const_iterator
      nu_indep_itr = nu_indep.begin();
    ele_t
      *itr         = &sort_array[0];
    for( size_type l = 1 ; l <= nu_indep.nz(); ++l, ++itr, ++nu_indep_itr ) {
      const size_type i = nu_indep_itr->indice() + nu_indep.offset();
      const value_type
        s_i          = s(i),
        y_i          = y(i),
        B_ii         = B_XX[l-1],
        s_i_B_ii_s_i = s_i*B_ii*s_i,
        s_i_y_i      = s_i*y_i;
      *sXTBXXsX += s_i_B_ii_s_i;
      *sXTyX    += s_i_y_i;
      itr->initialize( l, s_i_B_ii_s_i/sRTBRRsR + ::fabs(s_i_y_i)/sRTyR );
    }
  }
  // Sort this sparse vector in decending order
  std::sort(
    &sort_array[0], &sort_array[0] + sort_array.size()
    , AbstractLinAlgPack::SortByDescendingAbsValue()
    );
  // Extract this ordering
  {
    for( size_type l = 0; l < nu_indep.nz(); ++l )
      l_x_fixed_sorted[l] = sort_array[l].indice();
  }
}

void PBFGSPack::choose_fixed_free(
  const value_type                       project_error_tol
  ,const value_type                      super_basic_mult_drop_tol
  ,const SpVectorSlice                   &nu_indep
  ,const DVectorSlice                     &s
  ,const DVectorSlice                     &y
  ,const DVectorSlice                     &B_XX
  ,const size_type                       l_x_fixed_sorted[]
  ,EJournalOutputLevel                   olevel
  ,std::ostream                          &out
  ,value_type                            *sRTBRRsR
  ,value_type                            *sRTyR
  ,value_type                            *sXTBXXsX
  ,value_type                            *sXTyX
  ,size_type                             *n_pz_X
  ,size_type                             *n_pz_R
  ,size_type                             i_x_free[]
  ,size_type                             i_x_fixed[]
  ,ConstrainedOptPack::EBounds  bnd_fixed[]
  )
{
  using std::setw;
  using std::endl;
  using std::right;
  using AbstractLinAlgPack::norm_inf;
  namespace COP = ConstrainedOptPack;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  const size_type
    n_pz = nu_indep.size();

  // Store the status of the variables so that we can put sorted i_x_free[]
  // and i_x_fixed[] together at the end
  Workspace<long int>  i_x_status(wss,n_pz);  // free if > 0 , fixed if < 0 , error if 0
  std::fill_n( &i_x_status[0], n_pz, 0 );
  {for( size_type l = 0; l < (*n_pz_R); ++l ) {
    i_x_status[i_x_free[l]-1] = +1;
  }
  // Adjust i_x_free[] and i_x_fixed to meat the projection conditions
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out << "\nDetermining which fixed variables to put in superbasis and which to leave out (must have at least one in superbasis)...\n";
  }
  const value_type
    max_nu_indep = norm_inf(nu_indep);
  const bool
    all_fixed = n_pz == nu_indep.nz();
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ACTIVE_SET) ) {
    out << "\nmax{|nu_k(indep)|,i=r+1...n} = " << max_nu_indep            << std::endl
      << "super_basic_mult_drop_tol    = " << super_basic_mult_drop_tol << std::endl
      << "project_error_tol            = " << project_error_tol         << std::endl;
  }
  if( super_basic_mult_drop_tol > 1.0 ) {
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
      out << "super_basic_mult_drop_tol = " << super_basic_mult_drop_tol << " > 1"
        << "\nNo variables will be removed from the super basis!  (You might consider decreasing super_basic_mult_drop_tol < 1)\n";
    }
  }
  else {
    const int prec = out.precision();
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ACTIVE_SET) ) {
      out << endl
        << right << setw(10)      << "i"
        << right << setw(prec+12) << "nu_indep(i)"
        << right << setw(1)       << " "
        << right << setw(prec+12) << "s(i)*B(ii)*s(i)"
        << right << setw(prec+12) << "s_R'*B_RR*s_R"
        << right << setw(prec+12) << "s_X'*B_XX*s_X"
        << right << setw(1)       << " "
        << right << setw(prec+12) << "s(i)*y(i)"
        << right << setw(prec+12) << "s_R'*y_R"
        << right << setw(prec+12) << "s_X'*y_X"
        << right << setw(1)       << " "
        << right << setw(14)      << "status"
        << endl
        << right << setw(10)      << "--------"
        << right << setw(prec+12) << "---------------"
        << right << setw(1)       << " "
        << right << setw(prec+12) << "---------------"
        << right << setw(prec+12) << "---------------"
        << right << setw(prec+12) << "---------------"
        << right << setw(1)       << " "
        << right << setw(prec+12) << "---------------"
        << right << setw(prec+12) << "---------------"
        << right << setw(prec+12) << "---------------"
        << right << setw(1)       << " "
        << right << setw(14)      << "------------"
        << endl;
    }
    // Loop through the fixed variables in decending order of the violation.
    bool kept_one = false;
    for( size_type k = 0; k < nu_indep.nz(); ++k ) {
      const size_type
        l = l_x_fixed_sorted[k];
      const SpVectorSlice::element_type
        &nu_i = *(nu_indep.begin() + (l-1));
      const size_type
        i = nu_i.indice() + nu_indep.offset();
      const value_type
        abs_val_nu   = ::fabs(nu_i.value()),
        rel_val_nu   = abs_val_nu / max_nu_indep,
        s_i          = s(i),
        y_i          = y(i),
        B_ii         = B_XX[l-1],
        s_i_B_ii_s_i = s_i*B_ii*s_i,
        s_i_y_i      = s_i*y_i;
      const bool
        nu_cond       =  rel_val_nu < super_basic_mult_drop_tol,
        sXTBXXsX_cond = (*sXTBXXsX) / (*sRTBRRsR) > project_error_tol,
        sXTyX_cond    = ::fabs(*sXTyX) / ::fabs(*sRTyR) > project_error_tol,
        keep = ( (all_fixed && abs_val_nu == max_nu_indep && !kept_one)
             || nu_cond || sXTBXXsX_cond || nu_cond );
      if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ACTIVE_SET) ) {
        out << right << setw(10)      << i
          << right << setw(prec+12) << nu_i.value()
          << right << setw(1)       << (nu_cond ? "*" : " ")
          << right << setw(prec+12) << s_i_B_ii_s_i
          << right << setw(prec+12) << (*sRTBRRsR)
          << right << setw(prec+12) << (*sXTBXXsX)
          << right << setw(1)       << (sXTBXXsX_cond ? "*" : " ")
          << right << setw(prec+12) << s_i_y_i
          << right << setw(prec+12) << (*sRTyR)
          << right << setw(prec+12) << (*sXTyX)
          << right << setw(1)       << (sXTyX_cond ? "*" : " ")
          << right << setw(14)      << (keep ? "superbasic" : "nonbasic")
          << endl;
      }
      if(keep) {
        kept_one = true;
        *sRTBRRsR += s_i_B_ii_s_i;
        *sXTBXXsX -= s_i_B_ii_s_i;
        *sRTyR    += s_i_y_i;
        *sXTyX    -= s_i_y_i;
        i_x_status[i-1] = +1;
        ++(*n_pz_R);
      }
      else {
        i_x_status[i-1] = -1;
        ++(*n_pz_X);
      }
    }}
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
      out << "\nFinal selection: n_pz_X = " << (*n_pz_X) << " nonbasic variables and n_pz_R = " << (*n_pz_R) << " superbasic variables\n";
    }
  }
  // Set the final set of i_x_free[] and i_x_fixed[]
  SpVector::const_iterator
    nu_itr = nu_indep.begin(),
    nu_end = nu_indep.end();
  SpVector::difference_type
    nu_o = nu_indep.offset();
  size_type
    *i_x_free_itr  = i_x_free,
    *i_x_fixed_itr = i_x_fixed;
  ConstrainedOptPack::EBounds
    *bnd_fixed_itr = bnd_fixed;
  {for( size_type i = 1; i <= n_pz; ++i ) {
    long int status = i_x_status[i-1];
    TEUCHOS_TEST_FOR_EXCEPT( !( status ) ); // should not be zero!
    if( status > 0 ) {
      // A superbasic variable
      *i_x_free_itr++ = i;
    }
    else {
      // A nonbasic variable
      for( ; nu_itr->indice() + nu_o < i; ++nu_itr ); // Find the multiplier
      TEUCHOS_TEST_FOR_EXCEPT( !( nu_itr != nu_end ) );
      TEUCHOS_TEST_FOR_EXCEPT( !( nu_itr->indice() + nu_o == i  ) );
      *i_x_fixed_itr++ = i;
      *bnd_fixed_itr++ = ( nu_itr->value() > 0.0 ? COP::UPPER : COP::LOWER );
    }
  }}
  TEUCHOS_TEST_FOR_EXCEPT( !(  i_x_free_itr  - i_x_free  == *n_pz_R  ) );
  TEUCHOS_TEST_FOR_EXCEPT( !(  i_x_fixed_itr - i_x_fixed == *n_pz_X  ) );
  TEUCHOS_TEST_FOR_EXCEPT( !(  bnd_fixed_itr - bnd_fixed == *n_pz_X  ) );
}

} // end namespace MoochoPack

#endif // 0
