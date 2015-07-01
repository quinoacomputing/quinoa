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

#include <assert.h>

#include <iomanip>
#include <ostream>

#include "MoochoPack_active_set_change.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_SpVectorClass.hpp"

void MoochoPack::active_set_change(
  const SpVectorSlice& nu_k, const SpVectorSlice& nu_km1, Range1D indep
  ,EJournalOutputLevel olevel, std::ostream* out
  ,size_type* num_adds, size_type* num_drops
  ,size_type* num_active_indep, size_type* num_adds_indep, size_type* num_drops_indep
  )
{
  using std::setw;
  using std::endl;
  
  const int w = 12;

  *num_adds = *num_drops = *num_adds_indep = *num_drops_indep = 0;
  *num_active_indep = nu_k(indep).nz();

  if( !nu_k.nz() && !nu_km1.nz() )
    return;

  TEUCHOS_TEST_FOR_EXCEPT( !(  nu_k.is_sorted() && nu_km1.is_sorted()  ) );

  bool dump_change = (int)olevel >= (int)PRINT_ACTIVE_SET;

  if( dump_change ) {
    *out
      << "\n*** Changes in active set\n\n"
      << setw(w) << "i"
      << setw(w) << "uplo"
      << setw(w) << "dep/indep"
      << setw(w) << "change\n"
      << setw(w) << "----------"
      << setw(w) << "----------"
      << setw(w) << "----------"
      << setw(w) << "----------\n";
  }

  SpVectorSlice::const_iterator
    nu_k_itr	= nu_k.begin(),
    nu_k_end	= nu_k.end(),
    nu_km1_itr	= nu_km1.begin(),
    nu_km1_end	= nu_km1.end();

  while( nu_k_itr != nu_k_end || nu_km1_itr != nu_km1_end ) {
    if( nu_k_itr != nu_k_end && ( nu_km1_itr == nu_km1_end || ( nu_k_itr != nu_k_end
        && nu_k_itr->indice()+nu_k.offset() < nu_km1_itr->indice()+nu_km1.offset() ) ) )
    {
      // *nu_k_itr was added to active set.
      const size_type i = nu_k_itr->indice() + nu_k.offset();
      const bool is_indep = indep.in_range(i);
      if(is_indep)
        (*num_adds_indep)++;
      (*num_adds)++;
      if(dump_change)
        *out
          << setw(w) << i
          << setw(w) << ( nu_k_itr->value() >= 0.0 ? "upper" : "lower" )
          << setw(w) << ( is_indep ? "indep" : "dep" )
          << setw(w) << "added" << endl;
      nu_k_itr++;
    }
    else if( nu_km1_itr != nu_km1_end && ( nu_k_itr == nu_k_end || ( nu_km1_itr != nu_km1_end
        && nu_k_itr->indice()+nu_k.offset() > nu_km1_itr->indice()+nu_km1.offset() ) ) )
    {
      // *nu_km1_itr was removed from the active set.
      const size_type i = nu_km1_itr->indice() + nu_km1.offset();
      const bool is_indep = indep.in_range(i);
      if(is_indep)
        (*num_drops_indep)++;
      (*num_drops)++;
      if(dump_change)
        *out
          << setw(w) << i
          << setw(w) << ( nu_km1_itr->value() >= 0.0 ? "upper" : "lower" )
          << setw(w) << ( is_indep ? "indep" : "dep" )
          << setw(w) << "dropped" << endl;
      nu_km1_itr++;
    }
    else {
      // same variable (but the bound may have changed)
      const size_type i = nu_k_itr->indice() + nu_k.offset();
      const bool is_indep = indep.in_range(i);
      if( nu_k_itr->value() * nu_km1_itr->value() < 0.0 ) {
        // Switched bounds.
        if(is_indep) {
          (*num_adds_indep)++;
          (*num_drops_indep)++;
        }
        (*num_adds)++;
        (*num_drops)++;
      if(dump_change)
        *out
          << setw(w) << i
          << setw(w) << ( nu_k_itr->value() >= 0.0 ? "upper" : "lower" )
          << setw(w) << ( is_indep ? "indep" : "dep" )
          << setw(w) << "switch bnd" << endl;
      }
      nu_k_itr++;
      nu_km1_itr++;
    }
  }

  // Output summary
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    *out
      << "\n*** Active set change summary\n"
      << "nact_old  = "			<< nu_km1.nz()			<< endl
      << "nact_new  = "			<< nu_k.nz()			<< endl
      << "num_adds  = "			<< *num_adds			<< endl
      << "num_drops = "			<< *num_drops			<< endl
      << "nact_indep_old  = "		<< nu_km1(indep).nz()	<< endl
      << "nact_indep_new  = "		<< *num_active_indep	<< endl
      << "num_indep_adds  = "		<< *num_adds_indep		<< endl
      << "num_indep_drops = "		<< *num_drops_indep		<< endl;
  }
}

#endif // 0
