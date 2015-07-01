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

#include <math.h>

#include <ostream>
#include <fstream>
#include <typeinfo>

#include "MoochoPack_LineSearchFilter_Step.hpp"
#include "MoochoPack_Exceptions.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "NLPInterfacePack_NLPObjGrad.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_assert_print_nan_inf.hpp"
#include "AbstractLinAlgPack_VectorMutableSubView.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"

//#define FILTER_DEBUG_OUT 1


namespace {


class RemoveFilterEntryPred {
public:

  RemoveFilterEntryPred(
    const MoochoPack::value_type &f_with_boundary,
    const MoochoPack::value_type &theta_with_boundary,
    const Teuchos::RCP<std::ostream> &out
    )
    :f_with_boundary_(f_with_boundary),
     theta_with_boundary_(theta_with_boundary),
     out_(out)
    {}

  bool operator()(const MoochoPack::FilterEntry &entry)
    {
      if (
        entry.f >= f_with_boundary_
        &&
        entry.theta >= theta_with_boundary_
        )
      {
        if (!is_null(out_)) {
          *out_
            << "\nRemoving from the filter the redundant point:"
            << "\n  f_with_boundary     = " << entry.f
            << "\n  theta_with_boundary = " << entry.theta
            << "\n  iteration added     = " << entry.iter
            << std::endl;
        }
        return true;
      }
      return false;
    }

private:

  const MoochoPack::value_type f_with_boundary_;
  const MoochoPack::value_type theta_with_boundary_;
  const Teuchos::RCP<std::ostream> out_;

};


} // namespace


namespace MoochoPack {

// This must exist somewhere already, ask Ross
value_type MIN(value_type x, value_type y)
{ return (x < y) ? x : y; }

value_type MAX(value_type x, value_type y)
{ return (x > y) ? x : y; }

LineSearchFilter_Step::LineSearchFilter_Step( 
  Teuchos::RCP<NLPInterfacePack::NLP> nlp
  ,const std::string obj_iq_name
  ,const std::string grad_obj_iq_name
  ,const value_type &gamma_theta
  ,const value_type &gamma_f
  ,const value_type &f_min
  ,const value_type &gamma_alpha
  ,const value_type &delta        
  ,const value_type &s_f
  ,const value_type &s_theta
  ,const value_type &theta_small_fact
  ,const value_type &theta_max
  ,const value_type &eta_f
  ,const value_type &back_track_frac
  )
  :
  gamma_theta_(gamma_theta),
  gamma_f_(gamma_f),
  f_min_(f_min),
  gamma_alpha_(gamma_alpha),
  delta_(delta),
  s_f_(s_f),
  s_theta_(s_theta),
  theta_small_fact_(theta_small_fact),
  theta_max_(theta_max),
  eta_f_(eta_f),
  back_track_frac_(back_track_frac),
  filter_(FILTER_IQ_STRING),
  obj_f_(obj_iq_name),
  grad_obj_f_(grad_obj_iq_name),
  nlp_(nlp)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !nlp_.get(),
    std::logic_error,
    "Null nlp passed to LineSearchFilter_Step constructor"
    );

#if defined(FILTER_DEBUG_OUT)
  std::ofstream fout("filter_out.xml", std::ofstream::out | std::ofstream::trunc);
  fout << "<FilterDebugDocument>" << std::endl;
  fout.close();
#endif
}

LineSearchFilter_Step::~LineSearchFilter_Step()
{
#if defined(FILTER_DEBUG_OUT)
  std::ofstream fout("filter_out.xml", std::ofstream::out | std::ofstream::app);
  fout << "</FilterDebugDocument>" << std::endl;
  fout.close();
#endif
}
  
bool LineSearchFilter_Step::do_step(
  Algorithm& _algo, poss_type step_poss, 
  IterationPack::EDoStepType type
  ,poss_type assoc_step_poss)
{

  using Teuchos::dyn_cast;
  using IterationPack::print_algorithm_step;
  using LinAlgOpPack::Vp_StV;
  using std::setw;
    
  // Get Algorithm (cast), state, and problem
  NLPAlgo &algo = rsqp_algo(_algo);
  NLPAlgoState &s = algo.rsqp_state();

  EJournalOutputLevel olevel  = algo.algo_cntr().journal_output_level();
  std::ostream &out = algo.track().journal_out();
    
  // print step header
  if (static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS)) 
  { 
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
  }

  Teuchos::VerboseObjectTempState<NLP>
    nlpOutputTempState(
      nlp_, Teuchos::getFancyOStream(Teuchos::rcp(&out,false)),
      Teuchos::VERB_DEFAULT );
  // Above, we don't want to litter the output with any trace from the NLP.
  // However, if the user forces the verbosity level to be increased, then we
  // want to set the stream so that it knows where to print to.
    
  const size_type
    m  = nlp_->m();
    
  // Get the iteration quantity container objects
  IterQuantityAccess<value_type>
    &f_iq = obj_f_(s),
    &alpha_iq = s.alpha();
    
  IterQuantityAccess<VectorMutable>
    &x_iq   = s.x(),
    *c_iq   = m > 0 ? &s.c() : NULL,
    *h_iq   = NULL,                   // RAB: Replace latter!
    &Gf_iq  = grad_obj_f_(s);

  // check that all the pertinent information is known
  if (!s.d().updated_k(0) || !x_iq.updated_k(0))
  {
    // Dead in the water
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Error, d_k or x_k not updated." ); 		
    return false;
  }
    
  if (!alpha_iq.updated_k(0) || alpha_iq.get_k(0) > 1 || alpha_iq.get_k(0) <= 0)
  {
    // if alpha_k is not known then we would need to calculate all the new points
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::out_of_range, "Error, alpha_k not updated or out of range [0, 1)." ); 		
    return false;
  }

  // Setup some necessary parameters
  // Assuming that f_iq, Gf_iq, c_iq, h_iq are updated for k
  const value_type Gf_t_dk =
    (
      Gf_iq.updated_k(0)
      ? Gf_iq.get_k(0).inner_product( s.d().get_k(0) )
      : dyn_cast<NLPInterfacePack::NLPObjGrad>(*nlp_).calc_Gf_prod(
        x_iq.get_k(0),s.d().get_k(0)
        )
      );
  const value_type theta_k = CalculateTheta_k( c_iq, h_iq, 0);
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
    out << "\ntheta_k = ||c_k||1/c_k.dim() = " << theta_k << std::endl;
  const value_type gamma_f_used = CalculateGammaFUsed(f_iq,theta_k,olevel,out);
  const value_type theta_small = theta_small_fact_ * MAX(1.0,theta_k);
  const value_type alpha_min = CalculateAlphaMin( gamma_f_used, Gf_t_dk, theta_k, theta_small );

  value_type &alpha_k = alpha_iq.get_k(0);
  value_type theta_kp1 = 0.0;;

  // Print out some header/initial information
  int w = 15;
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
  {
    out << "\nBeginning Filter line search method.\n"
        << "\nCurrent Filter\n"
        << "-----------------------------------------------------\n"
        << "|" << setw(25) << "f_with_boundary     " 
        << "|" << setw(25) << "theta_with_boundary    "
        << "|\n"
        << "-----------------------------------------------------\n";

    IterQuantityAccess<Filter_T>& filter_iq = filter_(s);
    
    if (filter_iq.updated_k(-1)) {
      Filter_T& filter = filter_iq.get_k(-1);
      if (!filter.empty()) {
        for (
          Filter_T::iterator entry = filter.begin();
          entry != filter.end();
          ++entry
          )
        {	
          out << "|" << setw(25) << entry->f
              << " " << setw(25) << entry->theta
              << "|\n";
        }
      }
      else {
        out << "Filter is empty.\n";
      }
    }
    else {
      out << "Filter is empty.\n";
    }
  
    // dump header
    out << "\n  Iteration Status\n";
    out << "----------------------------------------------------------------------------------------------------------\n"; 

    out << "|" << setw(w) << "alpha_k    " 
        << "|" << setw(w) << "f_kp1     "
        << "|" << setw(w) << "theta_kp1   "
        << "|" << setw(w) << "pt. status   "
        << "|" << setw(40) << "comment                "
        << "|" << std::endl;
  
    out << "----------------------------------------------------------------------------------------------------------" 
        << std::endl;
  }

  // Begin the line search
  bool augment_filter = false;
  bool accepted = false;
  while (alpha_k > alpha_min && !accepted)
  {
    accepted = true;

    // Check that point is safe for calculations (no nans, infs, etc)
    if (!ValidatePoint(x_iq, f_iq, c_iq, h_iq, false))
    {
      accepted = false;
      
      // Print out some point information
      if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
      {
        int w = 15;
        // dump point
        out << "|" << setw(w) << " --- " 
            << " " << setw(w) << " --- "
            << " " << setw(w) << " --- "
            << " " << setw(w) << " failed "
            << " " << setw(40) << " nan_or_inf in calc"
            << " " << std::endl;
      }
  
      // Really, we do not need to throw an exception here, we can try and backtrack
      // alpha to get into an acceptable region
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::out_of_range, "Point Not Valid." );	
    }
      
    // Check if point satisfies filter
    if (accepted)
    {
      theta_kp1 = CalculateTheta_k(c_iq, h_iq, +1);

      // Print out some point information
      if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
      {
        // dump point
        out << "|" << setw(w) << alpha_k 
            << " " << setw(w) << f_iq.get_k(+1)
            << " " << setw(w) << theta_kp1;
      }

      accepted = CheckFilterAcceptability(f_iq.get_k(+1), theta_kp1, s);

      // Print out failure information
      if( !accepted && static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
      {
        out << " " << setw(w) << "failed"
            << " " << setw(40) << "Unacceptable to filter"
            << "|" << std::endl;
      }
    }

    // Check if point has theta less than theta_max
    if (accepted) {
      if (theta_kp1 > theta_max_) {
        accepted = false;
        // Print out failure information
        if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
          out << " " << setw(w) << "failed"
              << " " << setw(40) << "theta_kp1 > theta_max"
              << "|" << std::endl;
        }
      }
    }


    // Check if point satisfies sufficient decrease (Armijo on f if switching cond holds)
    if (accepted)
    {
      // Check for switching condition
      if (ShouldSwitchToArmijo(Gf_t_dk, alpha_k, theta_k, theta_small))
      {
        accepted = CheckArmijo(Gf_t_dk, alpha_k, f_iq);

        // Print out point information
        if(static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
        {
          if (accepted)
          { out << " " << setw(w) << "accepted"; }
          else
          { out << " " << setw(w) << "failed"; }

          out << " " << setw(40) << "Switch Cond. Holds (Armijo)" << "|" << std::endl;
        }
      }
      else
      {
        accepted = CheckFractionalReduction(f_iq, gamma_f_used, theta_kp1, theta_k);
        if (accepted)
        { augment_filter = true; }

        // Print out point information
        if(static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
        {
          if (accepted)
          { out << " " << setw(w) << "accepted"; }
          else
          { out << " " << setw(w) << "failed"; }

          out << " " << setw(40) << "Fraction Reduction (! Switch Cond )" << "|" << std::endl;
        }

      }
    }

    // if the point fails any of the tests, then backtrack
    if (!accepted)
    {
      // try a smaller alpha_k
      alpha_k = alpha_k*back_track_frac_;
      UpdatePoint(s.d().get_k(0), alpha_k, x_iq, f_iq, c_iq, h_iq, *nlp_);
    }	  

  } // end while


  if (accepted)
  {
    // Print status
    if(static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
      if (augment_filter)
        out << "\nPoint was accepted - augmenting the filter ...\n";
      else
        out << "Point was accepted - did NOT augment filter.\n";
    }

    if (augment_filter) {
      AugmentFilter( gamma_f_used, f_iq.get_k(+1), theta_kp1, s, olevel, out );
    }
    else {
      // Just update the filter from the last iteration
      UpdateFilter(s);
    }
  }
  else {
    // Could not find an acceptable alpha_k, go to restoration
    // Print status
    //if(static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) 
    //	{
    //	out << "\nCould not find acceptable alpha_k - going to restoration phase.\n"; 
    //  }
    
    //TEUCHOS_TEST_FOR_EXCEPTION( true, std::out_of_range, "Tried to go to restoration phase." );	
    
    TEUCHOS_TEST_FOR_EXCEPTION( true, LineSearchFailure
                        ,"FilterLineSearchFailure : Should go to restoration"
      );
  }

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) 
  {
    out << "\nx_kp1 =\n" << x_iq.get_k(+1);
    if (c_iq)
    { out << "\nc_kp1 =\n" << c_iq->get_k(+1); }
    if (h_iq)
    { out << "\nh_kp1 =\n" << h_iq->get_k(+1); }
  }

#if defined(FILTER_DEBUG_OUT)
  std::ofstream fout("filter_out.xml", std::ofstream::out | std::ostream::app);
  fout << "   <FilterIteration iter=\"" << s.k() << "\">" << std::endl;
  fout << "      <SelectedPoint alpha=\"" << alpha_k 
       << "\" f=\"" << f_iq.get_k(+1) 
       << "\" theta=\"" << theta_kp1
       << "\" />" << std::endl;

  // Output the filter
  fout << "      <Filter>" << std::endl;
    
  IterQuantityAccess<Filter_T>& filter_iq = filter_(s);
  if (filter_iq.updated_k(0))
  {    
    Filter_T& current_filter = filter_iq.get_k(0);
    for (
      Filter_T::iterator entry = current_filter.begin();
      entry != current_filter.end();
      ++entry
      )
    {
      fout << "         <FilterPoint iter=\"" << entry->iter 
           << "\" f=\"" << entry->f 
           << "\" theta=\"" << entry->theta << "\>" << std::endl;
    }
  }
  else
  {
    fout << "         <FilterNotUpdated/>" << std::endl;
  }

  fout << "      </Filter>" << std::endl;

    
  // Output the alpha curve
  fout << "      <AlphaCurve>" << std::endl;
  value_type alpha_tmp = 1.0;
  for (int i = 0; i < 10 || alpha_tmp > alpha_k; ++i)
  {
    UpdatePoint(s.d().get_k(0), alpha_tmp, x_iq, f_iq, c_iq, h_iq, *nlp_);
    if (ValidatePoint(x_iq, f_iq, c_iq, h_iq, false))
    {
      value_type theta = CalculateTheta_k(c_iq, h_iq, +1);
      fout << "         <AlphaPoint "
           << "alpha=\"" << alpha_tmp << "\" "
           << "f=\"" << f_iq.get_k(+1) << "\" "
           << "theta=\"" << theta << "\>" << std::endl;
    }

    alpha_tmp=alpha_tmp*back_track_frac_;
  }

  // restore alpha_k
  UpdatePoint(s.d().get_k(0), alpha_k, x_iq, f_iq, c_iq, h_iq, *nlp_);

  fout << "      </AlphaCurve>" << std::endl;

  fout << "   </FilterIteration" << std::endl;

  fout.close();

#endif
    
  return true;
}
  
void LineSearchFilter_Step::print_step(
  const Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss, std::ostream& out, const std::string& L
  ) const
{
  out
  << L << "*** Filter line search method\n"
  << L << "# Assumes initial d_k & alpha_k (0-1) is known and\n"
  << L << "# x_kp1, f_kp1, c_kp1 and h_kp1 are calculated for that alpha_k\n"
  << L << "Gf_t_dk = <Gf,dk>\n"
  << L << "theta_k = norm_1(c_k)/c_k.dim()\n"
  << L << "theta_small = theta_small_fact*max(1.0,theta_k)\n"
  << L << "if f_min != F_MIN_UNBOUNDED then\n"
  << L << "  gamma_f_used = gamma_f * (f_k-f_min)/theta_k\n"
  << L << "else\n"
  << L
  << L << "  gamma_f_used = gamma_f\n"
  << L << "end\n"
  << L << "if Gf_t_dk < 0 then\n"
  << L << "  alpha_min = min(gamma_theta, gamma_f_used*theta_k/(-Gf_t_dk))\n"
  << L << "  if theta_k <= theta_small then\n"
  << L << "    alpha_min = min(alpha_min, delta_*(theta_k^s_theta)/((-Gf_t_dk)^s_f))\n"
  << L << "  end\n"
  << L << "else\n"
  << L << "  alpha_min = gamma_theta\n"
  << L << "end\n"
  << L << "alpha_min = alpha_min*gamma_alpha\n"
  << L << "# Start the line search\n"
  << L << "accepted = false\n"
  << L << "augment = false\n"
  << L << "while alpha > alpha_min and accepted = false then\n"
  << L << "  accepted = true"
  << L << "  if any values in x_kp1, f_kp1, c_kp1, h_kp1 are nan or inf then\n"
  << L << "    accepted = false\n"
  << L << "  end\n"
  << L << "  # Check filter\n"
  << L << "  if accepted = true then\n"
  << L << "    theta_kp1 = norm_1(c_kp1)/c_kp1.dim()\n"
  << L << "    for each pt in the filter do\n"
  << L << "      if theta_kp1 >= pt.theta and f_kp1 >= pt.f then\n"
  << L << "        accepted = false\n"
  << L << "        break for\n"
  << L << "      end\n"
  << L << "    next pt\n"
  << L << "  end\n"
  << L << "  #Check Sufficient Decrease\n"
  << L << "  if accepted = true then"
  << L << "    # if switching condition is satisfied, use Armijo on f\n"
  << L << "    if theta_k < theta_small and Gf_t_dk < 0 and\n"
  << L << "        ((-Gf_t_dk)^s_f)*alpha_k < delta*theta_k^s_theta then\n"
  << L << "      if f_kp1 <= f_k + eta_f*alpha_k*Gf_t_dk then\n"
  << L << "        accepted = true\n"
  << L << "      end\n"
  << L << "    else\n"
  << L << "      # Verify factional reduction\n"
  << L << "      if theta_kp1 < (1-gamma_theta)*theta_k or f_kp1 < f_k - gamma_f_used*theta_k then\n"
  << L << "        accepted = true\n"
  << L << "        augment = true\n"
  << L << "      end\n"
  << L << "    end\n"
  << L << "  end\n"
  << L << "  if accepted = false then\n"
  << L << "    alpha_k = alpha_k*back_track_frac\n"
  << L << "    x_kp1  = x_k + alpha_k * d_k\n"
  << L << "    f_kp1  = f(x_kp1)\n"
  << L << "    c_kp1  = c(x_kp1)\n"
  << L << "    h_kp1  = h(x_kp1)\n"
  << L << "  end\n"    
  << L << "end while\n"
  << L << "if accepted = true then\n"
  << L << "  if augment = true then\n"
  << L << "    Augment the filter (use f_with_boundary and theta_with_boundary)\n"
  << L << "  end\n"
  << L << "else\n"
  << L << "  goto the restoration phase\n"
  << L << "end\n";
}
  
bool LineSearchFilter_Step::ValidatePoint( 
  const IterQuantityAccess<VectorMutable>& x,
  const IterQuantityAccess<value_type>& f,
  const IterQuantityAccess<VectorMutable>* c,
  const IterQuantityAccess<VectorMutable>* h,
  const bool throw_excpt
  ) const
{

  using AbstractLinAlgPack::assert_print_nan_inf;
  
  if (assert_print_nan_inf(x.get_k(+1), "x", throw_excpt, NULL) 
      || assert_print_nan_inf(f.get_k(+1), "f", throw_excpt, NULL)
      || (!c || assert_print_nan_inf(c->get_k(+1), "c", throw_excpt, NULL))
      || (!h || assert_print_nan_inf(h->get_k(+1), "c", throw_excpt, NULL)))
  {
    return true;
  }
  return false;
}
  
bool LineSearchFilter_Step::CheckFilterAcceptability( 
  const value_type f, 
  const value_type theta,
  const AlgorithmState& s
  ) const
{
  bool accepted = true;

  const IterQuantityAccess<Filter_T>& filter_iq = filter_(s);
    
  if (filter_iq.updated_k(-1))
  {
    const Filter_T &current_filter = filter_iq.get_k(-1);
    
    for (
      Filter_T::const_iterator entry = current_filter.begin();
      entry != current_filter.end();
      ++entry
      )
    {	
      if (f >= entry->f && theta >= entry->theta)
      {
        accepted = false;
        break;
      }
    }
  }

  return accepted;
}

bool LineSearchFilter_Step::CheckArmijo( 
  const value_type Gf_t_dk, 
  const value_type alpha_k, 
  const IterQuantityAccess<value_type>& f_iq
  ) const
{
  bool accepted = false;

  // Check Armijo on objective fn
  double f_kp1 = f_iq.get_k(+1);
  double f_k = f_iq.get_k(0);
  double lhs = f_k - f_kp1;
  double rhs = -eta_f_*alpha_k*Gf_t_dk;
  if ( lhs >= rhs )
  {
    // Accept pt, do NOT augment filter
    accepted = true;
  }

  return accepted;
}

bool LineSearchFilter_Step::CheckFractionalReduction( 
  const IterQuantityAccess<value_type>& f_iq,
  const value_type gamma_f_used,
  const value_type theta_kp1, 
  const value_type theta_k
  ) const
{
  bool accepted = false;
  if (theta_kp1 <= (1-gamma_theta_)*theta_k
      || f_iq.get_k(+1) <= f_iq.get_k(0)-gamma_f_used*theta_k )
  {
    // Accept pt and augment filter
    accepted = true;
  }
  return accepted;
}

void LineSearchFilter_Step::UpdatePoint( 
  const VectorMutable& d,
  const value_type alpha, 
  IterQuantityAccess<VectorMutable> &x,
  IterQuantityAccess<value_type>& f,
  IterQuantityAccess<VectorMutable>* c,
  IterQuantityAccess<VectorMutable>* h,
  NLP& nlp ) const
{  
  using LinAlgOpPack::Vp_StV;
  using AbstractLinAlgPack::assert_print_nan_inf;
  VectorMutable& x_kp1 = x.set_k(+1);
  x_kp1 = x.get_k(0);
  Vp_StV( &x_kp1, alpha, d);

  if (assert_print_nan_inf(x_kp1, "x", true, NULL))
  {
    // Calcuate f and c at the new point.
    nlp.unset_quantities();
    nlp.set_f( &f.set_k(+1) );
    if (c) nlp.set_c( &c->set_k(+1) );
    nlp.calc_f( x_kp1 ); 
    if (c) nlp.calc_c( x_kp1, false );
    nlp.unset_quantities();
  }
}

value_type LineSearchFilter_Step::CalculateAlphaMin( 
  const value_type gamma_f_used,
  const value_type Gf_t_dk,
  const value_type theta_k, 
  const value_type theta_small
  ) const
{
  value_type alpha_min = 0;
    
  if (Gf_t_dk < 0)
  {
    alpha_min = MIN(gamma_theta_, gamma_f_used*theta_k/(-Gf_t_dk));
    if (theta_k <= theta_small)
    {
      value_type switch_bound = delta_*pow(theta_k, s_theta_)/pow(-Gf_t_dk,s_f_);
      alpha_min = MIN(alpha_min, switch_bound);
    }
  }
  else
  {
    alpha_min = gamma_theta_;
  }
  return alpha_min * gamma_alpha_;
}

value_type LineSearchFilter_Step::CalculateTheta_k( 
  const IterQuantityAccess<VectorMutable>* c,
  const IterQuantityAccess<VectorMutable>* h,
  int k
  ) const
{
  value_type theta = 0.0;
  if (h) {
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::out_of_range, "Error, do not support inequalities yet" );
  }
  if (c) {
    const Vector &c_k = c->get_k(k);
    theta += ( c_k.norm_1() / c_k.dim() );
  }
  return theta;
}

value_type LineSearchFilter_Step::CalculateGammaFUsed(
  const IterQuantityAccess<value_type> &f,
  const value_type theta_k,
  const EJournalOutputLevel olevel,
  std::ostream &out
  ) const
{
  if( f_min() == F_MIN_UNBOUNDED) {
    if((int)olevel >= static_cast<int>(PRINT_ALGORITHM_STEPS))
      out << "\nf_min==F_MIN_UNBOUNDED:  Setting gamma_f_used = gamma_f = "<<gamma_f()<<".\n";
    return gamma_f();
  }
  const value_type f_k = f.get_k(0);
  if( f_k < f_min() ) {
    if((int)olevel >= static_cast<int>(PRINT_ALGORITHM_STEPS))
      out << "\nWarning!!! f_k = "<<f_k<<" < f_min = "<<f_min()<<"!  Setting gamma_f_used = gamma_f = "<<gamma_f()<<".\n";
    return gamma_f();
  }
  const value_type gamma_f_used = gamma_f() * ( f_k - f_min() ) / theta_k;
  if((int)olevel >= static_cast<int>(PRINT_ALGORITHM_STEPS))
    out << "\nf_min = "<<f_min()<<"!=F_MIN_UNBOUNDED: Setting gamma_f_used = gamma_f * (f_k - f_min)/theta_k = "
        <<gamma_f()<<" * ("<<f_k<<" - "<<f_min()<<")/"<<theta_k<<" = "<<gamma_f_used <<".\n";
  return gamma_f_used;
}

bool LineSearchFilter_Step::ShouldSwitchToArmijo( 
  const value_type Gf_t_dk,
  const value_type alpha_k,
  const value_type theta_k,
  const value_type theta_small
  ) const
{
  if (theta_k < theta_small && Gf_t_dk < 0) {
    if (pow(-Gf_t_dk, s_f_)*alpha_k - delta_*pow(theta_k, s_theta_) > 0) {
      return true;
    }
  }
  return false;
}


void LineSearchFilter_Step::UpdateFilter( IterationPack::AlgorithmState& s ) const
{
  
  IterQuantityAccess<Filter_T>& filter_iq = filter_(s);
    
  if (!filter_iq.updated_k(0))
  {
    if (filter_iq.updated_k(-1))
    {
      // initialize the filter from the last iteration
      filter_iq.set_k(0,-1);
    }
    else
    {
      // create an uninitialized filter
      filter_iq.set_k(0);
    }
  }

}

void LineSearchFilter_Step::AugmentFilter( 
  const value_type gamma_f_used,
  const value_type f,
  const value_type theta,
  IterationPack::AlgorithmState& s,
  const EJournalOutputLevel olevel,
  std::ostream &out
  ) const
{

  using Teuchos::as;

  const value_type
    f_with_boundary = f-gamma_f_used*theta,
    theta_with_boundary = (1.0-gamma_theta())*theta;

  if(static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    out << "\nAugmenting the filter with the point:"
        << "\n  f_with_boundary = f_kp1 - gamma_f_used*theta_kp1 = "<<f<<" - "<<gamma_f_used<<"*"<<theta<< " = " <<f_with_boundary
        << "\n  theta_with_boundary = (1-gamma_theta)*theta_kp1 = (1-"<<gamma_theta()<<")*"<<theta<<" = "<<theta_with_boundary
        << std::endl;
  }
  
  UpdateFilter(s);
  Filter_T& current_filter = filter_(s).get_k(0);

  // Remove current filter entries that are not longer needed
  current_filter.remove_if(
    RemoveFilterEntryPred(
      f_with_boundary, theta_with_boundary,
      as<int>(olevel) >= as<int>(PRINT_ALGORITHM_STEPS)
      ? Teuchos::rcpFromRef(out) : Teuchos::null
      )
    );

  // Now append the current point
  current_filter.push_front(FilterEntry(f_with_boundary, theta_with_boundary, s.k()));

}

// static

value_type LineSearchFilter_Step::F_MIN_UNBOUNDED = std::numeric_limits<value_type>::min();

} // end namespace MoochoPack
