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
#include <sstream>

#include "ConstrainedOptPack_DirectLineSearchArmQuad_Strategy.hpp"
#include "ConstrainedOptPack_MeritFuncCalc1D.hpp"
#include "check_nan_inf.h"

namespace ConstrainedOptPack {
inline value_type min(value_type v1, value_type v2) {
  return (v1 < v2) ? v1 : v2;
}
inline value_type max(value_type v1, value_type v2) {
  return (v1 > v2) ? v1 : v2;
}
}

ConstrainedOptPack::DirectLineSearchArmQuad_Strategy::DirectLineSearchArmQuad_Strategy(
  int           max_iter
  ,value_type   eta
  ,value_type   min_frac
  ,value_type   max_frac
  ,bool         max_out_iter
  )
  :max_iter_(max_iter)
  ,eta_(eta)
  ,min_frac_(min_frac)
  ,max_frac_(max_frac)
  ,max_out_iter_(max_out_iter)
{}

void ConstrainedOptPack::DirectLineSearchArmQuad_Strategy::set_max_iter(int max_iter)
{
  max_iter_ = max_iter;
}

int ConstrainedOptPack::DirectLineSearchArmQuad_Strategy::max_iter() const
{
  return max_iter_;
}

int ConstrainedOptPack::DirectLineSearchArmQuad_Strategy::num_iterations() const
{
  return num_iter_;
}

bool ConstrainedOptPack::DirectLineSearchArmQuad_Strategy::do_line_search(
    const MeritFuncCalc1D& phi, value_type phi_k
  , value_type* alpha_k, value_type* phi_kp1
  , std::ostream* out )
{
  using std::setw;
  using std::endl;

  if(*alpha_k < 0.0) {
    throw std::logic_error(	"DirectLineSearchArmQuad_Strategy::do_line_search(): "
      "alpha_k can't start out less than 0.0"	);
  }

  validate_parameters();
  
  int w = 20;
  int prec = 8;
  int prec_saved;
  if(out) {
    prec_saved = out->precision();
    *out	<< std::setprecision(prec)
        << "\nStarting Armijo Quadratic interpolation linesearch ...\n";
  }

  // Loop initialization (technically the first iteration)

  const value_type Dphi_k = phi.deriv();

  // output header
  if(out) {
    if(max_out_iter())
      *out << "\nmax_out_iter == true, maxing out the number of iterations!"; 
    *out	<< "\nDphi_k = "	<< Dphi_k
        << "\nphi_k = "		<< phi_k << "\n\n"
        << setw(5)			<< "itr"
        << setw(w)			<< "alpha_k"
        << setw(w)			<< "phi_kp1"
        << setw(w)			<< "phi_kp1-frac_phi\n"
        << setw(5)			<< "----"
        << setw(w)			<< "----------------"
        << setw(w)			<< "----------------"
        << setw(w)			<< "----------------\n";
  }

  // Check that this is a decent direction
  if(Dphi_k >= 0) throw DirectLineSearch_Strategy::NotDescentDirection(
    "DirectLineSearchArmQuad_Strategy::do_line_search(): "
    "The given d_k is not a descent direction for the given "
    "phi (phi.deriv() is >= 0)"			);

  // keep memory of the best value
  value_type	best_alpha = *alpha_k, best_phi = *phi_kp1;

  // Perform linesearch.
  bool success = false;
  for( num_iter_ = 0; num_iter_ < max_iter(); ++num_iter_ ) {

    // Print out this iteration.

    value_type frac_phi = phi_k + eta() * (*alpha_k) * Dphi_k;
    if(out)
      *out	<< setw(5)			<< num_iter_
          << setw(w)			<< *alpha_k
          << setw(w)			<< *phi_kp1
          << setw(w)			<< ((*phi_kp1)-frac_phi)	<< endl;
    
    // Check that this is a valid number.
    if( RTOp_is_nan_inf( *phi_kp1 ) ) {
      // Cut back the step to min_frac * alpha_k
      *alpha_k = min_frac()*(*alpha_k);
      best_alpha = 0.0;
      best_phi = phi_k;
    }
    else {		

      // Armijo condition
      if( *phi_kp1 < frac_phi ) {
        // We have found an acceptable point
        success = true;
        if( !max_out_iter() || ( max_out_iter() && num_iter_ == max_iter() - 1 ) )
          break;	// get out of the loop
      }

      // Select a new alpha to try:
      //   alpha_k = ( min_frac*alpha_k <= quadratic interpolation <= max_frac*alpha_k )

      // Quadratic interpolation of alpha_k that minimizes phi.
      // We know the values of phi at the initail point and alpha_k and
      // the derivate of phi w.r.t alpha at the initial point and
      // that's enough information for a quandratic interpolation.
      
      value_type alpha_quad =		( -0.5 * Dphi_k * (*alpha_k) * (*alpha_k) )
                    / ( (*phi_kp1) - phi_k - (*alpha_k) * Dphi_k );

      *alpha_k = min( max(min_frac()*(*alpha_k),alpha_quad), max_frac()*(*alpha_k) );

    }
    
    // Evaluate the point

    *phi_kp1 = phi(*alpha_k);

    // Save the best point found
    if(*phi_kp1 < best_phi) {
      best_phi = *phi_kp1;
      best_alpha = *alpha_k;
    }

  }

  // Be nice and reset the precision
  if(out) {
    out->precision(prec_saved);
  }

  if( success ) {
    return true;
  }

  // Line search failure.  Return the best point found and let the 
  // client decide what to do.
  *alpha_k = best_alpha;
  *phi_kp1 = phi(best_alpha);	// Make this the last call to phi(x)
  return false; 
}

void ConstrainedOptPack::DirectLineSearchArmQuad_Strategy::print_algorithm(
  std::ostream& out, const std::string& L) const
{
  out
    << L << "*** start line search using the Armijo cord test and quadratic interpolation of alpha\n"
    << L << "default: max_ls_iter = " << max_iter() << std::endl
    << L << "         eta = " << eta() << std::endl
    << L << "         min_frac = " << min_frac() << std::endl
    << L << "         max_frac = " << max_frac() << std::endl
    << L << "         max_out_iter = " << max_out_iter() << std::endl
    << L << "Dphi_k = phi.deriv()\n"
    << L << "if Dphi_k >= 0\n"
    << L << "    throw not_descent_direction\n"
    << L << "    end line search\n"
    << L << "end\n"
    << L << "best_alpha = alpha_k\n"
    << L << "best_phi = phi_kp1\n"
    << L << "for num_iter = 0... max_ls_iter\n"
    << L << "    frac_phi = phi_k + eta * alpha_k * Dphi_k\n"
    << L << "    print iteration\n"
    << L << "    if phi_kp1 is not a valid number then\n"
    << L << "        *** Cut back the step so the NLP's functions\n"
    << L << "        *** will hopefully be defined.\n"
    << L << "        alpha_k = min_frac * alpha_k\n"
    << L << "        best_alpha = 0\n"
    << L << "        best_phi = phi_k\n"
    << L << "    else\n"
    << L << "        if ( phi_kp1 < frac_phi ) then\n"
    << L << "            *** We have found an acceptable point\n"
    << L << "            if( !max_out_iter || max_out_iter && num_iter == max_ls_iter - 1 ) )\n"
    << L << "                end line search\n"
    << L << "            end\n"
    << L << "        end\n"
    << L << "        *** Use a quadratic interpoation to minimize phi(alpha)\n"
    << L << "        alpha_quad = (-0.5 * Dphi_k * alpha_k^2) / ( phi_kp1 - phi_k - alpha_k*Dphi_k )\n"
    << L << "        alpha_k = min( max( min_frac*alpha_k, alpha_quad ), max_frac*alpha_k )\n"
    << L << "    end\n"
    << L << "    phi_kp1 = phi(alpha_k)\n"
    << L << "    if phi_kp1 < best_phi\n"
    << L << "        best_phi = phi_kp1\n"
    << L << "        best_alpha = alpha_k\n"
    << L << "    end\n"
    << L << "end\n"
    << L << "*** If you get there the line search failed.\n"
    << L << "alpha_k = best_alpha\n"
    << L << "phi_kp1 = phi(alpha_k)\n"
    << L << "linesearch_failure = true\n";
}

void ConstrainedOptPack::DirectLineSearchArmQuad_Strategy::validate_parameters() const
{
  if( eta() < 0.0 || 1.0 < eta() ) {
    std::ostringstream omsg;
    omsg
      << "DirectLineSearchArmQuad_Strategy::validate_parameters() : "
      << "Error, eta = " << eta() << " is not in the range [0, 1]";
    throw std::invalid_argument( omsg.str() );
  }
  if( min_frac() < 0.0 || 1.0 < min_frac() ) {
    std::ostringstream omsg;
    omsg
      << "DirectLineSearchArmQuad_Strategy::validate_parameters() : "
      << "Error, min_frac = " << min_frac() << " is not in the range [0, 1]";
    throw std::invalid_argument( omsg.str() );
  }
  if( max_frac() < 0.0 || ( !max_out_iter() && 1.0 < max_frac() ) ) {
    std::ostringstream omsg;
    omsg
      << "DirectLineSearchArmQuad_Strategy::validate_parameters() : "
      << "Error, max_frac = " << max_frac() << " is not in the range [0, 1]";
    throw std::invalid_argument( omsg.str() );
  }
  if( max_frac() < min_frac() ) {
    std::ostringstream omsg;
    omsg
      << "DirectLineSearchArmQuad_Strategy::validate_parameters() : "
      << "Error, max_frac = " << max_frac()
      << " < min_frac = " << min_frac();;
    throw std::invalid_argument( omsg.str() );
  }
}
