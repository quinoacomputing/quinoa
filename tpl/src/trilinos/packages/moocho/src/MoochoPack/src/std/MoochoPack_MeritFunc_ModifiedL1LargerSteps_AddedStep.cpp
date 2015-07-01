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
#include <typeinfo>

#include "MoochoPack_MeritFunc_ModifiedL1LargerSteps_AddedStep.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "ConstrainedOptPack_MeritFuncPenaltyParams.hpp"
#include "ConstrainedOptPack_MeritFuncNLPDirecDeriv.hpp"
#include "ConstrainedOptPack/src/VectorWithNorms.h"
#include "DenseLinAlgPack_DVectorOp.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_DVectorOut.hpp"

namespace {

typedef MoochoPack::value_type value_type;
inline value_type max(value_type v1, value_type v2)
{	return (v1 > v2) ? v1 : v2; }

}

namespace MoochoPack {

MeritFunc_ModifiedL1LargerSteps_AddedStep::MeritFunc_ModifiedL1LargerSteps_AddedStep(
      const merit_func_ptr_t& merit_func
    , value_type	eta
    , int			after_k_iter
    , value_type	obj_increase_threshold
    , value_type	max_pos_penalty_increase
    , value_type	pos_to_neg_penalty_increase
    , value_type	incr_mult_factor )
  : merit_func_(merit_func), eta_(eta), after_k_iter_(after_k_iter)
    , obj_increase_threshold_(obj_increase_threshold)
    , max_pos_penalty_increase_(max_pos_penalty_increase)
    , pos_to_neg_penalty_increase_(pos_to_neg_penalty_increase)
    , incr_mult_factor_(incr_mult_factor)
{}

bool MeritFunc_ModifiedL1LargerSteps_AddedStep::do_step(Algorithm& _algo
  , poss_type step_poss, IterationPack::EDoStepType type
  , poss_type assoc_step_poss)
{
  using DenseLinAlgPack::norm_inf;
  using DenseLinAlgPack::dot;

  NLPAlgo	&algo	= rsqp_algo(_algo);
  NLPAlgoState	&s		= algo.rsqp_state();
  
  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
  }

  MeritFuncPenaltyParams
    *params = dynamic_cast<MeritFuncPenaltyParams*>(&merit_func());
  if( !params ) {
    std::ostringstream omsg;
    omsg
      << "MeritFunc_ModifiedL1LargerSteps_AddedStep::do_step(...), Error "
      << "The class " << typeName(&merit_func()) << " does not support the "
      << "MeritFuncPenaltyParams iterface\n";
    out << omsg.str();
    throw std::logic_error( omsg.str() );
  }

  bool consider_modifications = s.k() >= after_k_iter();
  
  if( !consider_modifications )
    return true;

  if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
    out << "\nk = " << s.k() << " >= " << " after_k_iter = " << after_k_iter()
      << "\nConsidering increasing the penalty parameters ...\n";
  }

  // /////////////////////////////////////////
  // Get references to iteration quantities

  const value_type
    &f_k		= s.f().get_k(0),
    &f_kp1		= s.f().get_k(+1);

  const DVector
    &c_k		= s.c().get_k(0).cv(),
    &c_kp1		= s.c().get_k(+1).cv();

  const DVector
    &Gf_k		= s.Gf().get_k(0).cv(),
    &d_k		= s.d().get_k(0).cv();

  // Determining if the objective increase is sufficent.

  const value_type
    very_small	= std::numeric_limits<value_type>::min(),
    obj_increase = ( f_kp1 - f_k ) / ::fabs( f_kp1 + f_k + very_small );
  bool attempt_modifications = obj_increase >= obj_increase_threshold(); 

  if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
    out << "\n(f_kp1-f_k)/|f_kp1+f_k+very_small| = " << obj_increase
      << ( attempt_modifications ? " >= " : " < " )
      << "obj_increase_threshold = " << obj_increase_threshold() << std::endl;
  }

  if( obj_increase < obj_increase_threshold() ) {
    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      out << "\nLeave the penalty parameters as they are.\n";
    }
  }
  else {
    // Compute the penalty parameters.
    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      out << "\nTry to increase the penalty parameters to allow a larger SQP step... .\n";
    }

    DVectorSlice
      mu = params->mu();

    // ///////////////////////////////////////////////////////////
    // Derivation of the modification to the penalty parameters
    //
    // Given the modified L1 merit function:
    //
    // phi(x) = f(x) + sum( mu(j) * |c(x)(j)|, j=1..m )
    // Dphi(x_k,d) = Gf_k'*d - sum( mu(j) * |c(x)(j)|, j=1..m )
    //
    // Given the armijo condition for a full step:
    //
    // phi(x_kp1) <= phi(x_k) + eta * Dphi(x_k,d)
    //
    // ->
    //
    // f_kp1 - sum(mu(j)*|c_kp1(j)|,j=1..m)
    //		<= f_k - sum(mu(j)*|c_k(j)|,j=1..m)
    //			+ eta*( Gf_k'*d - sum(mu(j)*|c_k(j)|,j=1..m) )
    //
    // -> (1)
    //
    // f_kp1 - f_k - eta*Gf_k'*d <= sum( mu(j) * del(j), j=1...m )
    //
    // where:
    //		del(j) = (1-eta) * c_k(j) - c_kp1(j)
    //
    // Define the sets:
    //
    // DelPos = { j | del(j) >  0 }		(2.a)
    // DelNeg = { j | del(j) =< 0 }		(2.b)
    //
    // Define the update expresions for the penalty parameters:
    //
    // mu(j) <- mu(j) + a * ( mu_max - mu(j) ), for j <: DelPos	    (3.a)
    //
    // mu(j) <- mu(j) + (a/b) * ( mu_max - mu(j) ), for j <: DelNeg	(3.b)
    //
    // where:
    //		a < max_pos_penalty_increase ( >= 0 )
    //		b = pos_to_neg_penalty_increase ( >= 0 )
    //		mu_max = (1.0 + incr_mult_factor) * ||mu||inf
    //		0 < a < max_pos_penalty_increase : The length to be determined
    //			so that (1) can be satsified.
    //
    // The idea there is to make (1) an equality, plug (3) into it and then
    // solve for a.
    //
    // (1), (3) -> (4)
    //
    // a = ( f_kp1 - f_k - eta*Gf_k'*d - num_term ) / ( pos_denom_term + neg_denom_term )
    //
    // where:
    //		num_term = sum( mu(j) * del(j), j=1..m )
    //		pos_denom_term = sum( (mu_max - mu(j)) * del(j), j <: DelPos )
    //		neg_denom_term = (1/b) * sum( (mu_max - mu(j)) * del(j), j <: NegPos )
    //
    // If the value of a from (4) is within 0 < a < max_pos_penalty_increase
    // then that means that can increase the penalties
    // enough and satisfy (1).  If a < 0 then we would
    // have to decrease the penalties and we are not allowed to do this.
    // If a > max_pos_penalty_increase then we are not allowed to increase
    // the penalties enough to
    // satisfy (1) but it suggests that if we increase them up to (3) for
    // a = max_pos_penalty_increase
    // that we would be able to take a larger SQP step durring our linesearch.

    // //////////////////////////
    // Compute the terms in (4)

    const value_type
      mu_max = norm_inf( mu ) * (1.0 + incr_mult_factor());

    value_type
      num_term = 0.0,
      pos_denom_term = 0.0,
      neg_denom_term = 0.0;

    typedef std::vector<bool> del_pos_t;	// Remember the sets DelPos, DelNeg
    del_pos_t
      del_pos( mu.size() );

    {
      DVectorSlice::const_iterator
        mu_itr		= const_cast<const DVectorSlice&>(mu).begin();
      DVector::const_iterator
        c_k_itr		= c_k.begin(),
        c_kp1_itr	= c_kp1.begin();

      del_pos_t::iterator
        del_pos_itr = del_pos.begin();

      for( ; c_k_itr != c_k.end(); ++mu_itr, ++c_k_itr, ++c_kp1_itr, ++del_pos_itr ) {
        TEUCHOS_TEST_FOR_EXCEPT( !(  mu_itr < const_cast<const DVectorSlice&>(mu).end()  ) );
        TEUCHOS_TEST_FOR_EXCEPT( !(  c_kp1_itr < c_kp1.end()  ) );
        TEUCHOS_TEST_FOR_EXCEPT( !(  del_pos_itr < del_pos.end()  ) );

        const value_type
          del_j = ( 1 - eta() ) * ::fabs( *c_k_itr ) - ::fabs( *c_kp1_itr );
          
        num_term += (*mu_itr) * del_j;

        if( del_j > 0 ) {
          *del_pos_itr = true;
          pos_denom_term += ( mu_max - (*mu_itr) ) * del_j;
        }
        else {
          *del_pos_itr = false;
          neg_denom_term += ( mu_max - (*mu_itr) ) * del_j;
        }
      }
      neg_denom_term /= pos_to_neg_penalty_increase();
    }

    // Compute a from (4)
    value_type
      a = ( f_kp1 - f_k - eta() * dot(Gf_k,d_k) - num_term)
        / ( pos_denom_term + neg_denom_term );

    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      out << "\nnum_term       = " << num_term
        << "\npos_denom_term = " << pos_denom_term
        << "\nneg_denom_term = " << neg_denom_term
        << "\n\na = " << a << std::endl;
    }

    if( a < 0.0 ) {
      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        out << "\na < 0 : Leave the penalty parameters alone\n";
      }
      return true;
    }
    else if( a > max_pos_penalty_increase() ) {
      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        out << "\na > max_pos_penalty_increase = " << max_pos_penalty_increase()
          << "\nSet a = max_pos_penalty_increase ...\n";
      }
      a = max_pos_penalty_increase();
    }
    else {
      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        out << "\n0 <= a <= max_pos_penalty_increase = " << max_pos_penalty_increase()
          << "\nWe should be able to take a full SQP step ...\n";
      }
    }

    // Update the penalty parameters using (3)
    {
      const value_type
        pos_step = a,
        neg_step = pos_step / pos_to_neg_penalty_increase(); 
      del_pos_t::const_iterator
        del_pos_itr = const_cast<const del_pos_t&>(del_pos).begin();
      DVectorSlice::iterator
        mu_itr		= mu.begin();
      for( ; mu_itr != mu.end(); ++del_pos_itr, ++mu_itr ) {
        TEUCHOS_TEST_FOR_EXCEPT( !(  del_pos_itr < const_cast<const del_pos_t&>(del_pos).end()  ) );
        *mu_itr = *mu_itr
              + (*del_pos_itr ?pos_step :neg_step) * (mu_max - (*mu_itr));
      }		
    }

    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      out << "\nmax(|mu(j)|) = " << (*std::max_element( mu.begin(), mu.end() ))
        << "\nmin(|mu(j)|) = " << (*std::min_element( mu.begin(), mu.end() ))
          << std::endl;
    }

    if( (int)olevel >= (int)PRINT_VECTORS ) {
      out << "\nmu = \n" << mu;
    }
  }

  return true;
}

void MeritFunc_ModifiedL1LargerSteps_AddedStep::print_step( const Algorithm& algo
  , poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  , std::ostream& out, const std::string& L ) const
{
  out
    << L << "*** Increase the penalty parameters for the modified L1 merit function\n"
    << L << "*** so as to allow for a larger SQP step.\n"
    << L << "default: eta                          = " << eta() << std::endl
    << L << "         after_k_iter                 = " << after_k_iter() << std::endl
    << L << "         obj_increase_threshold       = " << obj_increase_threshold() << std::endl
    << L << "         max_pos_penalty_increase     = " << max_pos_penalty_increase() << std::endl
    << L << "         pos_to_neg_penalty_increase  = " << pos_to_neg_penalty_increase() << std::endl
    << L << "         incr_mult_factor             = " << incr_mult_factor() << std::endl
    << L << "if k < after_k_iter then\n"
    << L << "    goto next step\n"
    << L << "end\n"
    << L << "if (f_kp1-f_k)/abs(f_kp1+f_k+very_small) >= obj_increase_threshold then\n"
    << L << "    mu = phi.mu()\n"
    << L << "    *** Try to increase to penalty parameters mu(j) to allow for a full step.\n"
    << L << "    mu_max = norm(mu,inf) * (1.0+incr_mult_factor)\n"
    << L << "    num_term = 0\n"
    << L << "    pos_denom_term = 0\n"
    << L << "    neg_denom_term = 0\n"
    << L << "    for j = 1 ... m\n"
    << L << "        del(j) = (1-eta)*abs(c_k(j)) - abs(c_kp1(k))\n"
    << L << "        num_term = num_term + mu(j) * del(j)\n"
    << L << "        if del(j) > 0 then\n"
    << L << "            del_pos(j) = true\n"
    << L << "            pos_denom_term = pos_denom_term + (mu_max - mu(j)) * del(j)\n"
    << L << "        else\n"
    << L << "            del_pos(j) = false\n"
    << L << "            neg_denom_term = neg_denom_term + (mu_max - mu(j)) * del(j)\n"
    << L << "        end\n"
    << L << "    end\n"
    << L << "    neg_denom_term = (1/pos_to_neg_penalty_increase) * neg_denom_term\n"
    << L << "    a = ( f_kp1 - f_k - eta * dot(Gf_k,d_k) - num_term)\n"
    << L << "    		/ ( pos_denom_term + neg_denom_term )\n"
    << L << "    if a < 0 then\n"
    << L << "        *** We can't take a larger SQP step by increasing mu(j)\n"
    << L << "        goto next step\n"
    << L << "    else if a > max_pos_penalty_increase then\n"
    << L << "        *** We are not allowed to increase mu(j) enough to allow a\n"
    << L << "        *** full SQP step but we will still increase mu(j) to take\n"
    << L << "        *** a hopefully larger step\n"
    << L << "        a = max_pos_penalty_increase\n"
    << L << "    else\n"
    << L << "        *** We can increase mu(j) and take a full SQP step\n"
    << L << "    end\n"
    << L << "    *** Increase the multipliers\n"
    << L << "    for j = 1...m\n"
    << L << "        if del_pos(j) == true then\n"
    << L << "            mu(j) = mu(j) + a*(mu_max - mu(j))\n"
    << L << "        else\n"
    << L << "            mu(j) = mu(j) + (a/pos_to_neg_penalty_increase)*(mu_max - mu(j))\n"
    << L << "        end\n"
    << L << "    end\n"
    << L << "end\n";
}

}	// end namespace MoochoPack

#endif // 0
