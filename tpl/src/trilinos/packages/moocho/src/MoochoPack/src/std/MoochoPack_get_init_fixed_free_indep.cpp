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

#include "MoochoPack_get_init_fixed_free_indep.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_SpVectorOp.hpp"

void MoochoPack::get_init_fixed_free_indep(
  const size_type                        n
  ,const size_type                       r
  ,const SpVectorSlice                   &nu_indep
  ,const value_type                      super_basic_mult_drop_tol
  ,EJournalOutputLevel                   olevel
  ,std::ostream                          &out
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

  const size_type
    n_pz = n-r;

  // Loop through and set i_x_free and i_x_fixed
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out << "\nDetermining which fixed variables to remove from rHL to form rHL_RR (can remove all but one)...\n";
  }
  const value_type
    max_nu_indep = norm_inf(nu_indep);
  const bool
    all_fixed = n_pz == nu_indep.nz();
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out << "\nmax{|nu_k(indep)|,i=r+1...n} = " << max_nu_indep << std::endl;
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
        << right << setw(prec+12) << "nu(i)"
        << right << setw(8)       << "status"
        << endl
        << right << setw(10)      << "--------"
        << right << setw(prec+12) << "--------"
        << right << setw(8)       << "------"
        << endl;
    }
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
    *n_pz_X = 0;
    *n_pz_R = 0;
    bool kept_one = false;
    {for( size_type i_indep = 1; i_indep <= n_pz; ++i_indep ) {
      if( nu_itr != nu_end && (nu_itr->indice() + nu_o) == i_indep ) {
        const value_type
          abs_val = ::fabs(nu_itr->value()),
          rel_val = abs_val / max_nu_indep;
        const bool
          keep = ( (all_fixed && abs_val == max_nu_indep && !kept_one)
               || rel_val < super_basic_mult_drop_tol );
        if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ACTIVE_SET) ) {
          out << right << setw(10)      << i_indep + r
            << right << setw(prec+12) << nu_itr->value()
            << right << setw(8)       << (keep ? "keep" : "drop")
            << endl;
        }
        if(!keep) {
          *i_x_fixed_itr++ = i_indep;
          namespace COP = ConstrainedOptPack;
          *bnd_fixed_itr++
            = ( nu_itr->value() > 0.0 ? COP::UPPER : COP::LOWER );
          // ToDo: Consider fixed variable bounds
          ++(*n_pz_X);
        }
        else {
          kept_one = true;
        }
        ++nu_itr;
        if(!keep) continue;
      }
      *i_x_free_itr++ = i_indep;
      ++(*n_pz_R);
    }}
    TEUCHOS_TEST_FOR_EXCEPT( !(  i_x_free_itr  - i_x_free  == *n_pz_R  ) );
    TEUCHOS_TEST_FOR_EXCEPT( !(  i_x_fixed_itr - i_x_fixed == *n_pz_X  ) );
    TEUCHOS_TEST_FOR_EXCEPT( !(  bnd_fixed_itr - bnd_fixed == *n_pz_X  ) );
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
      out << "\nRemoving n_pz_X = " << (*n_pz_X) << " from the superbasic set and keeping n_pz_R = " << (*n_pz_R) << std::endl;
    }
  }
}

#endif // 0
