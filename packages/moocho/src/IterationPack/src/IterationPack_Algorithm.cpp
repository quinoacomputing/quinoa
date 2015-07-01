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

#include <signal.h>

#include <iterator>
#include <numeric>

#include "IterationPack_Algorithm.hpp"
#include "StopWatchPack_stopwatch.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#ifdef HAVE_MPI
#include "mpi.h"
#endif

// Define to see MPI/interrupt deugging output
//#define ITERATION_PACK_ALGORITHM_SHOW_MPI_DEBUG_INFO

// Define of the MPI implementation receives signals on all processes
//#define ITERATION_PACK_ALGORITHM_SIGNALS_ON_ALL_PROCESSES;

extern "C" {

void sig_handler_interrupt_algorithm( int signum )
{
  IterationPack::Algorithm::interrupt();
}

} // extern "C"

namespace {

// Helper functions

template< class T >
inline
T my_max( const T& v1, const T& v2 ) { return v1 > v2 ? v1 : v2; }

// Private static data for IterationPack::Algorithm.
// I put it here so that I can modify it without affecting the
// header file and avoiding unnecessary recompilations.

enum EInterruptStatus { NOT_INTERRUPTED=0, STOP_END_STEP=1, STOP_END_ITER=2, ABORT_PROGRAM=3 };

int static_mpi_initialized = false;
int static_num_running_algorithms = 0;
int static_num_proc = 0; // Flag that no algorithm has been even allocated yet!
int static_proc_rank = 0;
bool static_interrupt_called = false;
bool static_processed_user_interrupt = false;
EInterruptStatus static_interrupt_status = NOT_INTERRUPTED;
bool static_interrupt_terminate_return = false;

} // end namespace

// ToDo: change step_itr and assoc_step_itr to just return iterators without
// asserting if the names exist.  This will be more useful.

namespace IterationPack {

// constructors / destructor

Algorithm::Algorithm()
  :running_state_(NOT_RUNNING), max_iter_(100)
  ,max_run_time_(std::numeric_limits<double>::max())
  ,next_step_name_(0), do_step_next_called_(false), reconfigured_(false)
  ,time_stats_computed_(false)
{
  // Set MPI info
  static_num_proc = 1;
  static_proc_rank = 0;
#ifdef HAVE_MPI
  // If MPI is not initialized then this must be because the code was
  // compiled with support for MPI but it not actually using it.
  // Therefore, we will initialize MPI but not bother to finialize it.
  if(!static_mpi_initialized) {
    int mpi_initialized = false;
    MPI_Initialized(&mpi_initialized);
    if(!mpi_initialized) {
      int argc = 1;
      char arg_str[] = "dummy_prg";
      char *arg_str_ptr = arg_str;
      char **argv = &arg_str_ptr;
      MPI_Init( &argc, &argv );
    }
    static_mpi_initialized = true;
  }
  // ToDo: Allow the specification of another communicator if needed!
  MPI_Comm_size( MPI_COMM_WORLD, &static_num_proc );
  MPI_Comm_rank( MPI_COMM_WORLD, &static_proc_rank );
#ifdef ITERATION_PACK_ALGORITHM_SHOW_MPI_DEBUG_INFO
  std::cerr << "\np=" << static_proc_rank << ": Algorithm::Algorithm() being called (num_proc = "<<static_num_proc<<") ... \n";
#endif
#endif // HAVE_MPI
}

Algorithm::~Algorithm()
{}

// maximum iterations

void Algorithm::max_iter(size_t max_iter)
{	max_iter_ = max_iter; }

size_t Algorithm::max_iter() const
{	return max_iter_; }

// maximum run tine

void Algorithm::max_run_time(double max_run_time)
{	max_run_time_ = max_run_time; }

double Algorithm::max_run_time() const
{	return max_run_time_; }


// step information / access

int Algorithm::num_steps() const
{	return steps_.size(); }

Algorithm::poss_type Algorithm::get_step_poss(const std::string& step_name) const
{	
  steps_t::const_iterator itr = step_itr(step_name);
  return itr == steps_.end() ? DOES_NOT_EXIST : std::distance( steps_.begin(), itr ) + 1;
}

const std::string& Algorithm::get_step_name(poss_type step_poss) const
{	return steps_[validate(step_poss) - 1].name; }

Algorithm::step_ptr_t& Algorithm::get_step(poss_type step_poss)
{	return steps_[validate(step_poss) - 1].step_ptr; }

const Algorithm::step_ptr_t& Algorithm::get_step(poss_type step_poss) const
{	return steps_[validate(step_poss) - 1].step_ptr; }

// pre/post step information / access

int Algorithm::num_assoc_steps(poss_type step_poss, EAssocStepType type) const
{	return assoc_steps_[validate(step_poss) - 1][type].size(); }

Algorithm::poss_type Algorithm::get_assoc_step_poss(poss_type step_poss, EAssocStepType type
  ,const std::string& assoc_step_name) const
{	
  // ToDo: change to return DOES_NOT_EXIST if it does not exist.
  const assoc_steps_ele_list_t &assoc_list = assoc_steps_[validate(step_poss) - 1][type];
  assoc_steps_ele_list_t::const_iterator itr = assoc_step_itr(assoc_list,assoc_step_name);
  return itr == assoc_list.end() ? DOES_NOT_EXIST : std::distance( assoc_list.begin() , itr ) + 1;
}

const std::string& Algorithm::get_assoc_step_name(poss_type step_poss, EAssocStepType type
  , poss_type assoc_step_poss) const
{
  const assoc_steps_ele_list_t &assoc_list= assoc_steps_[validate(step_poss) - 1][type];
  validate(assoc_list,assoc_step_poss);
  assoc_steps_ele_list_t::const_iterator itr = assoc_list.begin();
  std::advance( itr, assoc_step_poss - 1 );
  return (*itr).name;
}

Algorithm::step_ptr_t& Algorithm::get_assoc_step(poss_type step_poss, EAssocStepType type
  , poss_type assoc_step_poss)
{
  assoc_steps_ele_list_t &assoc_list= assoc_steps_[validate(step_poss) - 1][type];
  validate(assoc_list,assoc_step_poss);
  assoc_steps_ele_list_t::iterator itr = assoc_list.begin();
  std::advance( itr, assoc_step_poss - 1 );
  return (*itr).step_ptr;
}

const Algorithm::step_ptr_t& Algorithm::get_assoc_step(poss_type step_poss, EAssocStepType type
  , poss_type assoc_step_poss) const
{
  const assoc_steps_ele_list_t &assoc_list= assoc_steps_[validate(step_poss) - 1][type];
  validate(assoc_list,assoc_step_poss);
  assoc_steps_ele_list_t::const_iterator itr = assoc_list.begin();
  std::advance( itr, assoc_step_poss - 1 );
  return (*itr).step_ptr;
}

// step manipulation

void Algorithm::insert_step(poss_type step_poss, const std::string& step_name, const step_ptr_t& step)
{
  validate_not_in_state(RUNNING);
  TEUCHOS_TEST_FOR_EXCEPTION(
    step.get() == NULL, std::invalid_argument
    ,"Algorithm::insert_step(...) : A step with the name = \'" << step_name
    << "\' being inserted into the position  = " << step_poss
    << " has step.get() == NULL!" );
  // Make sure a step with this name does not already exist.
  steps_t::iterator itr;
  if( steps_.end() != ( itr = step_itr(step_name) ) )
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, AlreadyExists
      ,"Algorithm::insert_step(...) : A step with the name = " << step_name
      << " already exists at step_poss = " << std::distance(steps_.begin(),itr) + 1 );
  // insert the step in such a way that any container can be used for steps_
  itr = steps_.begin();
  std::advance ( itr , validate(step_poss,+1) - 1 );
  steps_.insert( itr , steps_ele_t(step,step_name) );
  // insert the assoc_step element in such a way that any container can be used for assoc_steps_
  assoc_steps_t::iterator a_itr = assoc_steps_.begin();
  std::advance ( a_itr , step_poss - 1 );
  assoc_steps_.insert( a_itr , assoc_steps_ele_t() );
}

void Algorithm::change_step_name(poss_type step_poss, const std::string& new_name)
{
  validate_not_in_state(RUNNING);
  if(running_state() == RUNNING_BEING_CONFIGURED) {
    validate_not_curr_step(validate(step_poss));
    validate_not_next_step(steps_[step_poss - 1].name);
  }
  steps_[step_poss - 1].name = new_name;
}

void Algorithm::replace_step(poss_type step_poss, const step_ptr_t& step)
{
  validate_not_in_state(RUNNING);
  if(running_state() == RUNNING_BEING_CONFIGURED)	validate_not_curr_step(validate(step_poss));
  steps_[step_poss - 1].step_ptr = step;
}

void Algorithm::remove_step(poss_type step_poss)
{
  validate_not_in_state(RUNNING);
  if(running_state() == RUNNING_BEING_CONFIGURED) {
    validate_not_curr_step(validate(step_poss));
    validate_not_next_step(steps_[step_poss - 1].name);
  }
  // remove the step in such a way that any container can be used for steps_
  steps_t::iterator itr = steps_.begin();
  std::advance ( itr , validate(step_poss) - 1 );
  steps_.erase( itr );
  // remove the assoc_step element in such a way that any container can be used for assoc_steps_
  assoc_steps_t::iterator a_itr = assoc_steps_.begin();
  std::advance ( a_itr , step_poss - 1 );
  assoc_steps_.erase( a_itr );
}

// pre/post step manipulation

void Algorithm::insert_assoc_step(poss_type step_poss, EAssocStepType type, poss_type assoc_step_poss
  , const std::string& assoc_step_name, const step_ptr_t& assoc_step)
{
  validate_not_in_state(RUNNING);
  TEUCHOS_TEST_FOR_EXCEPTION(
    assoc_step.get() == NULL, std::invalid_argument
    ,"Algorithm::insert_assoc_step(...) : A step with the name = \'" << assoc_step_name
    << "\' being inserted into the position  = " << step_poss
    << "." << ( type == PRE_STEP
          ? (int)assoc_step_poss - num_assoc_steps(step_poss,type) - 1
          : assoc_step_poss )
    << " has assoc_step.get() == NULL!" );
  if(running_state() == RUNNING_BEING_CONFIGURED) validate_not_curr_step(validate(step_poss));
  // Make sure an associated step with this name does not already exist.
  assoc_steps_ele_list_t &assoc_list = assoc_steps_[step_poss - 1][type];
  validate(assoc_list,assoc_step_poss,+1);
  assoc_steps_ele_list_t::iterator itr = assoc_list.begin();
  char assoc_type_name[2][10] = { "PRE_STEP" , "POST_STEP" };
  if( assoc_list.end() != ( itr = assoc_step_itr(assoc_list,assoc_step_name) ) )
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, AlreadyExists
      ,"Algorithm::insert_assoc_step(...) : An associated step of type = "
      <<	assoc_type_name[type]
      << " with the name = " << assoc_step_name
      << " already exists at step_poss = " << step_poss
      << " and assoc_step_poss = " <<  std::distance(assoc_list.begin(),itr) + 1 );
  // insert an associated step in such a way that any container could be used.
  itr = assoc_list.begin();
  std::advance( itr, assoc_step_poss - 1 );
  assoc_list.insert( itr , assoc_steps_ele_list_ele_t(assoc_step,assoc_step_name) );
}

void Algorithm::remove_assoc_step(poss_type step_poss, EAssocStepType type, poss_type assoc_step_poss)
{
  validate_not_in_state(RUNNING);
  if(running_state() == RUNNING_BEING_CONFIGURED) validate_not_curr_step(validate(step_poss));
  validate(step_poss);
  assoc_steps_ele_list_t &assos_list = assoc_steps_[step_poss - 1][type];
  validate(assos_list,assoc_step_poss);
  assoc_steps_ele_list_t::iterator itr = assos_list.begin();
  std::advance( itr, assoc_step_poss - 1 );
  assos_list.erase( itr );
}

//  runtime configuration updating control

void Algorithm::begin_config_update()
{
  validate_in_state(RUNNING);
  saved_next_step_name_ = *next_step_name_;
  saved_curr_step_name_ = steps_[curr_step_poss_ - 1].name;
  change_running_state(RUNNING_BEING_CONFIGURED);
}

void Algorithm::end_config_update()
{
  validate_in_state(RUNNING_BEING_CONFIGURED);

  // update next_step_poss_ and next_step_name_.
  steps_t::iterator itr = step_itr(saved_next_step_name_);
  TEUCHOS_TEST_FOR_EXCEPT( !(  itr != steps_.end()  ) );	// the step with this name should not have been deleted or changed.
  next_step_poss_ = std::distance( steps_.begin() , itr ) + 1;
  next_step_name_ = &(*itr).name;

  // update curr_step_poss_
  itr = step_itr(saved_curr_step_name_);
  TEUCHOS_TEST_FOR_EXCEPT( !(  itr != steps_.end()  ) );	// the step with this name should not have been deleted or changed.
  curr_step_poss_ = std::distance( steps_.begin() , itr ) + 1;

  // inform the step objects that *this has changes.
  imp_inform_steps( &AlgorithmStep::inform_updated );

  change_running_state(RUNNING);
  reconfigured_ = true;
}

// algorithmic control

void Algorithm::do_step_next(const std::string& step_name)
{
  validate_in_state(RUNNING);
  steps_t::iterator itr = step_itr_and_assert(step_name);
  next_step_poss_ = std::distance( steps_.begin() , itr ) + 1;
  next_step_name_ = &(*itr).name;
  do_step_next_called_ = true;
}

void Algorithm::do_step_next(poss_type step_poss)
{
  validate_in_state(RUNNING);
  const steps_ele_t &ele = steps_[validate(step_poss) - 1];
  next_step_poss_ = step_poss;
  next_step_name_ = &ele.name;
  do_step_next_called_ = true;
}

const std::string& Algorithm::what_is_next_step_name() const
{
  validate_in_state(RUNNING);
  return *next_step_name_;
}

Algorithm::poss_type Algorithm::what_is_next_step_poss() const
{	
  validate_in_state(RUNNING);
  return next_step_poss_;
}

bool Algorithm::do_step(const std::string& step_name)
{
  validate_in_state(RUNNING);
  return imp_do_step( std::distance( steps_.begin() , step_itr_and_assert(step_name) ) + 1 );
}

bool Algorithm::do_step(poss_type step_poss)
{
  validate_in_state(RUNNING);
  return imp_do_step(step_poss);
}

void Algorithm::terminate(bool success)
{	
  validate_in_state(RUNNING);
  terminate_status_ = success ? STATUS_TERMINATE_TRUE : STATUS_TERMINATE_FALSE;
}

// start iterations

EAlgoReturn Algorithm::do_algorithm(poss_type step_poss)
{
  using StopWatchPack::stopwatch;

  validate_in_state(NOT_RUNNING);

  track().initialize();

  try{
  
  terminate_status_ = STATUS_KEEP_RUNNING;
  change_running_state(RUNNING);

  first_k_ = state().k();
  next_step_poss_ = validate(step_poss);
  next_step_name_ = &steps_[step_poss - 1].name;
  
  // Prepair for timing algorithm
  step_times_.resize( algo_timing_ ? (num_steps()+1) * (max_iter()+1+NUM_STEP_TIME_STATS) : 0 );
  if( algo_timing_ ) {
//		step_times_[ max_iter() ] = 0.0;	// flag for statistics not calc. yet.
//		// set iteration totals to zero
//		if( step_times_[(max_iter() + 1 + 5) * num_steps()] != 0.0 )
//			std::fill_n( step_times_.begin() + (max_iter() + 1 + 5) * num_steps(), max_iter(), 0.0 );
    std::fill_n( step_times_.begin(), step_times_.size(), 0.0 );	// Try setting everything to zero?
    time_stats_computed_ = false;
  }
  stopwatch step_timer;
  stopwatch overall_timer;

  imp_inform_steps( &AlgorithmStep::initialize_step );

  overall_timer.start();
  for(;;) {

    curr_step_poss_ = next_step_poss_;
    // Note that curr_step_poss_ may change if there is a runtime
    // change in the configuration of the steps.

    bool keep_on = true;

    // Execute the steps for this step
    
    if( algo_timing_ ) {
      step_timer.reset();
      step_timer.start();
    }

    keep_on = imp_do_step(curr_step_poss_);

    if( algo_timing_ ) {
      const double time = my_max(step_timer.stop(),-1e-50);	// negative somehow (g++ -O2 ?)
      // time for step k for the iteration
      step_times_[state().k()-first_k_+(curr_step_poss_-1)*(max_iter()+1+NUM_STEP_TIME_STATS)] = time;
      // Add to time for the full iteration
      step_times_[state().k()-first_k_+(num_steps())*(max_iter()+1+NUM_STEP_TIME_STATS)] += time;
    }

    // See if a step object called terminate(...)
    if(terminate_status_ != STATUS_KEEP_RUNNING) {
      EAlgoReturn algo_return;
      if( static_interrupt_status == STOP_END_STEP ) {
        algo_return = ( terminate_status_ == STATUS_TERMINATE_TRUE
                        ? INTERRUPTED_TERMINATE_TRUE
                        : INTERRUPTED_TERMINATE_FALSE );
        static_interrupt_status = NOT_INTERRUPTED;
      }
      else {
        algo_return = ( terminate_status_ == STATUS_TERMINATE_TRUE
                        ? TERMINATE_TRUE
                        : TERMINATE_FALSE );
      }
      return finalize_algorithm(algo_return);
    }

    if(keep_on) {

      // All the step objects returned true so increment the step and loop around

      if( curr_step_poss_ == static_cast<poss_type>(num_steps()) ) {

        //
        // This is the last step in the algorithm
        //
  
        // Output this iteration
        track().output_iteration(*this);

        // Check if the maximum number of iterations has been exceeded.
        if( state().k() - first_k_ >= max_iter() ) {
          return finalize_algorithm(MAX_ITER_EXCEEDED);
        }

        // Check if the maximum runtime has been exceeded.
        if( ( overall_timer.read() / 60 ) >= max_run_time() ) {
          return finalize_algorithm(MAX_RUN_TIME_EXCEEDED);
        }

        // Set if the algorithm was interrupted
        if( static_interrupt_status == STOP_END_ITER ) {
          static_interrupt_status = NOT_INTERRUPTED;
          const EAlgoReturn algo_return = ( static_interrupt_terminate_return
                                            ? INTERRUPTED_TERMINATE_TRUE
                                            : INTERRUPTED_TERMINATE_FALSE );
          return finalize_algorithm(algo_return);
        }

        // Transition the iteration quantities to k = k + 1
        state().next_iteration();

        // Setup to start the major loop over again
        next_step_poss_ = 1;
        next_step_name_ = &steps_[0].name;

      }
      else {

        // else just increment the step
        ++next_step_poss_;
        next_step_name_ = &steps_[next_step_poss_ - 1].name;

      }

      continue;	// loop around

    }
    else {
      // some step object returned false from its do_step(..) operation so it
      // should have called do_step_next(...) to request a jump to
      // a specific operation.
      if(!do_step_next_called_)
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, InvalidControlProtocal
          ,"EAlgoReturn Algorithm::do_algorithm(...) :"
          " A step object returned false from its do_step(...) operation"
          " without calling do_step_next(...) to request jump to a specific"
          " step." );
      do_step_next_called_ = false;
      // just loop around and do the step that the step object requested
      // by changing next_step_poss_ by its call to do_step_next(...).
    }
  }	// end for(;;)

  }	// end try
  catch(...) {
    try {
      finalize_algorithm(TERMINATE_FALSE);
    }
    catch(...) {
      // We tried to finalize gracefully but we failed!
    }
    throw;
  }
}

// algorithm information output

void Algorithm::print_steps(std::ostream& out) const
{
  out << "\n*** Algorithm Steps ***\n\n";
  imp_print_algorithm(out,false);
  out << std::endl;
}

void Algorithm::print_algorithm(std::ostream& out) const
{
  out << "\n*** Iteration Quantities ***\n\n";
  state().dump_iter_quant(out);
  out << std::endl;
  out << "\n*** Algorithm Description ***\n\n";
  imp_print_algorithm(out,true);
  out << std::endl;
}

// Algorithm Timing.

void Algorithm::set_algo_timing( bool algo_timing ) {
   validate_not_in_state(RUNNING);
  algo_timing_ = algo_timing;
}

bool Algorithm::algo_timing() const {
  return algo_timing_;
}

void Algorithm::print_algorithm_times( std::ostream& out ) const
{
  using std::setw;
  using std::endl;

  validate_not_in_state(RUNNING);

  if( step_times_.size() == 0 ) {
    out << "No step timing was performed\n";
    return;
  }

  const int w = 10;
  const int prec = 4;
  const int n = num_steps();					// Total steps
  const int m = state().k() - first_k_ + 1;	// Total number of iterations performed
  const int mm = max_iter()+1;				// Total number of possible iterations
  const int mmm = mm + NUM_STEP_TIME_STATS;	// total entries in a step_i row
   
  // Print the header.
  out	<< "\n\n**************************************\n"
    << "*** Algorithm step CPU times (sec) ***\n";

  // Print the step names.
  out	<< "\nStep names"
    << "\n----------\n";
  {for( int i = 1; i <= n; ++i ) {
    out	<< i << ") \"" << get_step_name(i) << "\"\n";	
  }}
  out	<< n+1 << ") Iteration total\n";	
  out << endl;

  out << std::right << std::setprecision(prec);

  // Print table header
  out << setw(w) << "" << "  steps 1..." << n+1 << " ->\n\n";
  
  // print step numbers
  out	<< setw(w)	<< " iter k";
  {for( int i = 1; i <= n+1; ++i ) {
    out	<< setw(w) << i;	
  }}
  out << endl;
  out	<< setw(w)	<< "--------";
  {for( int i = 1; i <= n+1; ++i ) {
    out	<< setw(w) << "--------";	
  }}
  out << endl;
  // Print the step times.
  {for( int k = 0; k < m; ++k ) {
    out	<< setw(w)	<< first_k_ + k;
    {for( int i = 0; i < n+1; ++i ) {
      out	<< setw(w) << step_times_[k+i*mmm];	
    }}
    out << endl;
  }}

  // Compute the (1) totals for each step, the (2) average, (3) min and (4) max times
  // per iteration for each step and the (5) precentages for each step.

  compute_final_time_stats();

  // Ouput time statistics.
  
  out	<< setw(w)	<< "--------";
  {for( int i = 1; i <= n+1; ++i ) {
    out	<< setw(w) << "--------";	
  }}

  // Output the total times for each step.
  out << endl;
  out	<< setw(w)	<< "total(sec)";
  {for( int i = 0; i < n+1; ++i ) {
    const double *step_i_times = &step_times_[i*mmm];
    out	<< setw(w) << step_i_times[ mm + TIME_STAT_TOTALS_OFFSET ];	
  }}
  out << endl;

  // Output the average times per iteration
  out	<< setw(w)	<< "av(sec)/k";
  {for( int i = 0; i < n+1; ++i ) {
    const double *step_i_times = &step_times_[i*mmm];
    out	<< setw(w) << step_i_times[ mm + TIME_STAT_AV_OFFSET ];	
  }}
  out << endl;

  // Output the min times per iteration
  out	<< setw(w)	<< "min(sec)";
  {for( int i = 0; i < n+1; ++i ) {
    const double *step_i_times = &step_times_[i*mmm];
    out	<< setw(w) << step_i_times[ mm + TIME_STAT_MIN_OFFSET ];	
  }}
  out << endl;

  // Output the max times per iteration
  out	<< setw(w)	<< "max(sec)";
  {for( int i = 0; i < n+1; ++i ) {
    const double *step_i_times = &step_times_[i*mmm];
    out	<< setw(w) << step_i_times[ mm + TIME_STAT_MAX_OFFSET ];	
  }}
  out << endl;

  // Output the precentage times for each step.
  out	<< setw(w)	<< "% total";
  {for( int i = 0; i < n+1; ++i ) {
    const double *step_i_times = &step_times_[i*mmm];
    out	<< setw(w) << step_i_times[ mm + TIME_STAT_PERCENT_OFFSET ] * 100.0;	
  }}
  out << endl;


  // Print total time for entire algorithm.
  out << "------------------------------" << endl
    << "total CPU time = " << total_time_ << " sec\n";;
}


void Algorithm::get_step_times_k( int offset, double step_times[] ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    step_times_.size() == 0, std::logic_error
    ,"Algorithm::get_step_times_k(...) : times requested, but no times calculated!"
    );
  TEUCHOS_TEST_FOR_EXCEPTION(
    offset > 0, std::invalid_argument
    ,"Algorithm::get_step_times_k(...) : Can\'t get times for an iteratin that has not occured yet!."
    );

  const int n = num_steps();					// Total steps
  //const int m = state().k() - first_k_ + 1;	// Total number of iterations performed
  const int mm = max_iter()+1;				// Total number of possible iterations
  const int mmm = mm + NUM_STEP_TIME_STATS;	// total entries in a step_i row
  
  const int k = state().k() + offset;
  {for (int step = 0; step < n+1; ++step) {
    step_times[step] = step_times_[step*mmm + k];		
  }}

}

void Algorithm::get_final_step_stats( size_t step, double* total, double* average, double* min, double* max, double* percent) const
{
  // Compute the (1) totals for each step, the (2) average, (3) min and (4) max times
  // per iteration for each step and the (5) precentages for each step.
  compute_final_time_stats();

  //const int n = num_steps();					// Total steps
  //const int m = state().k() - first_k_ + 1;	// Total number of iterations performed
  const int mm = max_iter()+1;				// Total number of possible iterations
  const int mmm = mm + NUM_STEP_TIME_STATS;	// total entries in a step_i row

  double* step_i_times = &const_cast<step_times_t&>(step_times_)[step*mmm];
  if (total) {
    *total   = step_i_times[mm + TIME_STAT_TOTALS_OFFSET];
  }
  if (average) {
    *average = step_i_times[mm + TIME_STAT_AV_OFFSET];
  }
  if (min) {
    *min     = step_i_times[mm + TIME_STAT_MIN_OFFSET];
  }
  if (max) {
    *max     = step_i_times[mm + TIME_STAT_MAX_OFFSET];
  }
  if (percent) {
    *percent = step_i_times[mm + TIME_STAT_PERCENT_OFFSET];
  }
}

EAlgoReturn Algorithm::finalize_algorithm( EAlgoReturn algo_return )
{
  change_running_state(NOT_RUNNING);
  imp_inform_steps( &AlgorithmStep::finalize_step );
  track().output_final(*this,algo_return);
  return algo_return;
}

void Algorithm::compute_final_time_stats() const
{
  if (!time_stats_computed_) {
    time_stats_computed_ = true;
    
    const int n = num_steps();					// Total steps
    const int m = state().k() - first_k_ + 1;	// Total number of iterations performed
    const int mm = max_iter()+1;				// Total number of possible iterations
    const int mmm = mm + NUM_STEP_TIME_STATS;	// total entries in a step_i row
    
    // compute totals for each step (1...n) and the full iteration (n+1)
    double &_total_time = const_cast<double&>(total_time_);
    _total_time = 0.0;
    
    {for( int i = 0; i < n+1; ++i ) {
      double *step_i_times = &const_cast<step_times_t&>(step_times_)[i*mmm];
      // compute total step times (and total algorithm time)
      const double
        step_time = std::accumulate( step_i_times, step_i_times + m, (double)0.0 );
      if(i < n)
        _total_time += step_time;
      step_i_times[ mm + TIME_STAT_TOTALS_OFFSET ] = step_time;
      // compute average per step.
      step_i_times[ mm + TIME_STAT_AV_OFFSET ] = step_time / m;
      // compute min per step
      step_i_times[ mm + TIME_STAT_MIN_OFFSET ]= *std::min_element( step_i_times, step_i_times + m );
      // compute max per step
      step_i_times[ mm + TIME_STAT_MAX_OFFSET ]= *std::max_element( step_i_times, step_i_times + m );
    }}
    
    {for( int i = 0; i < n+1; ++i ) {
      double *step_i_times = &const_cast<step_times_t&>(step_times_)[i*mmm];
      // compute fractions for each step.
      step_i_times[ mm + TIME_STAT_PERCENT_OFFSET ]
        = step_i_times[ mm + TIME_STAT_TOTALS_OFFSET ] / total_time_;
    }}
  }
}

// private

void Algorithm::change_running_state(ERunningState _running_state)
{
  if( running_state() != RUNNING && _running_state == RUNNING ) {
    if( static_num_running_algorithms == 0 ) {
      // Register the signal handler for the SIGINT
      signal( SIGINT, &sig_handler_interrupt_algorithm );
      static_interrupt_called = false;
      static_processed_user_interrupt = false;
    }
    ++static_num_running_algorithms;
  }
  else if( running_state() != NOT_RUNNING && _running_state == NOT_RUNNING ) {
    --static_num_running_algorithms;
    if( static_num_running_algorithms == 0 ) {
      // Put back the default signal handler
      signal( SIGINT, SIG_DFL );
      static_interrupt_called = false;
      static_processed_user_interrupt = false;
    }
  }
  running_state_ = _running_state;
}

void Algorithm::validate_in_state(ERunningState _running_state) const {
  const char running_state_name[3][25] = { "NOT_RUNNING" , "RUNNING", "RUNNING_BEING_CONFIGURED" };
  if(running_state() != _running_state)
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, InvalidRunningState
      ,"Algorithm::validate_in_state(...) : The condition running_state() == "
      << running_state_name[running_state()] << " has been violated with "
      << " running_state = " << running_state_name[_running_state] );
}

void Algorithm::validate_not_in_state(ERunningState _running_state) const {
  const char running_state_name[3][25] = { "NOT_RUNNING" , "RUNNING", "RUNNING_BEING_CONFIGURED" };
  if(running_state() == _running_state)
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, InvalidRunningState
      ,"Algorithm::validate_not_in_state(...) : The condition running_state() != "
      << running_state_name[running_state()] << " has been violated" );
}

void Algorithm::validate_not_curr_step(poss_type step_poss) const {
  if(step_poss == curr_step_poss_)
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, InvalidConfigChange
      ,"Algorithm::validate_not_curr_step(step_poss="<<step_poss<<") : "
      "Error, You can not modify the step being currently executed" );
}

void Algorithm::validate_not_next_step(const std::string& step_name) const {
  if( step_name == saved_next_step_name_ )
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, InvalidConfigChange,
      "Algorithm::validate_not_next_step(step_name): "
      "Error, You can not modify name or remove the step given by "
      "step_name = what_is_next_name() = " << step_name );
}

Algorithm::steps_t::iterator Algorithm::step_itr_and_assert(const std::string& step_name)
{
  steps_t::iterator itr = step_itr(step_name);
  if(itr == steps_.end())
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, DoesNotExist
      ,"Algorithm::step_itr(...) : A step with the name "
      << step_name << " does not exist." );
  return itr;	
}

Algorithm::steps_t::const_iterator Algorithm::step_itr_and_assert(const std::string& step_name) const
{
  steps_t::const_iterator itr = step_itr(step_name);
  if(itr == steps_.end())
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, DoesNotExist
      ,"Algorithm::step_itr(...) : A step with the name "
      << step_name << " does not exist." );
  return itr;	
}

bool Algorithm::imp_do_step(poss_type step_poss) {
  curr_step_poss_ = step_poss;
  // do the pre steps in order
  if( !imp_do_assoc_steps(PRE_STEP) ) return false;
  // do the main step
  if( !steps_[curr_step_poss_-1].step_ptr->do_step(*this, curr_step_poss_, DO_MAIN_STEP, 0) ) return false;
  // do the post steps in order
  if( !imp_do_assoc_steps(POST_STEP) ) return false;
  // if you get here all the pre steps, step, and post steps returned true.
  if( static_interrupt_status == NOT_INTERRUPTED )
    look_for_interrupt();
  if( static_interrupt_status == STOP_END_STEP ) {
    terminate( static_interrupt_terminate_return );
    return false;
  }
  return true;
}

bool Algorithm::imp_do_assoc_steps(EAssocStepType type) {
  assoc_steps_ele_list_t				*assoc_list	= &assoc_steps_[curr_step_poss_ - 1][type];
  assoc_steps_ele_list_t::iterator	itr			= assoc_list->begin();
  int									n			= assoc_list->size();
  for(int i = 1; i <= n; ++itr, ++i) {
    if(reconfigured_) {
      // The associated step just has reconfigured *this
      // so we must update our pointers and iterators.
      // Since it is not allowed for this step or its associated steps
      // to have been changed, the next associated step to
      // execute will not change.
      assoc_list	= &assoc_steps_[curr_step_poss_ - 1][type];
      itr			= assoc_list->begin();
      std::advance( itr, i - 1 );
      reconfigured_ = false;	// This works as long as no one else needs to know
                  // if *this has been reconfigured.
    }
    if( !(*(*itr).step_ptr).do_step(*this, curr_step_poss_, do_step_type(type), i) ) return false;
  }
  return true;	// All the associated steps returned true.
}

void Algorithm::imp_inform_steps(inform_func_ptr_t inform_func_ptr)
{
  steps_t::const_iterator         s_itr = steps_.begin();
  assoc_steps_t::const_iterator   a_itr = assoc_steps_.begin();
  poss_type step_i = 1;
  for(; step_i <= static_cast<poss_type>(num_steps()); ++step_i, ++s_itr, ++a_itr) {
    // pre_steps (e.q. 2.-3, 2.-2, 2.-1)
    const assoc_steps_ele_list_t &pre_steps = (*a_itr)[PRE_STEP];
    assoc_steps_ele_list_t::const_iterator pre_step_itr = pre_steps.begin();
    for(int pre_step_i = - pre_steps.size(); pre_step_i < 0; ++pre_step_i, ++pre_step_itr) {
      ((&*(*pre_step_itr).step_ptr)->*inform_func_ptr)(
        *this, step_i, DO_PRE_STEP, pre_steps.size()+pre_step_i+1
        );
    }
    // The main step.
    ((&*(*s_itr).step_ptr)->*inform_func_ptr)( *this, step_i, DO_MAIN_STEP, 0 );
    // post_steps (e.q. 2.1, 2.2, 2.3)
    const assoc_steps_ele_list_t &post_steps = (*a_itr)[POST_STEP];
    assoc_steps_ele_list_t::const_iterator post_step_itr = post_steps.begin();
    for(int post_step_i = 1; post_step_i <= static_cast<int>(post_steps.size()); ++post_step_i, ++post_step_itr) {
      ((&*(*post_step_itr).step_ptr)->*inform_func_ptr)(
        *this, step_i, DO_POST_STEP, post_step_i
        );
    }
  }
}

void Algorithm::imp_print_algorithm(std::ostream& out, bool print_steps) const
{
  using Teuchos::typeName;
  const std::string leading_str = "    ";
  
  steps_t::const_iterator          s_itr = steps_.begin();
  assoc_steps_t::const_iterator    a_itr = assoc_steps_.begin();
  poss_type step_i = 1;
  for(; step_i <= static_cast<poss_type>(num_steps()); ++step_i, ++s_itr, ++a_itr) {
    // list pre_steps (e.q. 2.-3, 2.-2, 2.-1)
    const assoc_steps_ele_list_t &pre_steps = (*a_itr)[PRE_STEP];
    assoc_steps_ele_list_t::const_iterator pre_step_itr = pre_steps.begin();
    for(int pre_step_i = - pre_steps.size(); pre_step_i < 0; ++pre_step_i, ++pre_step_itr) {
      out		<< step_i << "." << pre_step_i << ". \""
          << (*pre_step_itr).name << "\"\n"
          << leading_str << "(" << typeName(*(*pre_step_itr).step_ptr) << ")\n";
      if(print_steps) {
        (*(*pre_step_itr).step_ptr).print_step( *this, step_i, DO_PRE_STEP
          , pre_steps.size()+pre_step_i+1, out, leading_str );
        out << std::endl;
      }
    }
    // The main step.
    out		<< step_i << ". \"" << (*s_itr).name
        << "\"\n"
        << leading_str << "(" << typeName(*(*s_itr).step_ptr) << ")\n";
    if(print_steps) {
      (*(*s_itr).step_ptr).print_step( *this, step_i, DO_MAIN_STEP, 0, out, leading_str );
      out << std::endl;
    }
    // list post_steps (e.q. 2.1, 2.2, 2.3)
    const assoc_steps_ele_list_t &post_steps = (*a_itr)[POST_STEP];
    assoc_steps_ele_list_t::const_iterator post_step_itr = post_steps.begin();
    for(int post_step_i = 1; post_step_i <= static_cast<poss_type>(post_steps.size()); ++post_step_i, ++post_step_itr) {
      out		<< step_i << "." << post_step_i << ". \""
          << (*post_step_itr).name << "\"\n"
          << leading_str << "(" << typeName(*(*post_step_itr).step_ptr) << ")\n";
      if(print_steps) {
        (*(*post_step_itr).step_ptr).print_step( *this, step_i, DO_POST_STEP, post_step_i
          , out, leading_str );
        out << std::endl;
      }
    }
  }
  if(print_steps) {
    out
      << step_i << ". \"Major Loop\" :\n"
      << "    if k >= max_iter then\n"
      << "        terminate the algorithm\n"
      << "    elseif run_time() >= max_run_time then\n"
      << "        terminate the algorithm\n"
      << "    else\n"
      << "        k = k + 1\n"
      << "        goto 1\n"
      << "    end\n";
  }
}

// validate poss

Algorithm::poss_type Algorithm::validate(poss_type step_poss, int past_end) const
{
    
  TEUCHOS_TEST_FOR_EXCEPTION(
    step_poss < 1 || steps_.size() + past_end < step_poss, DoesNotExist
    ,"Algorithm::validate(step_poss) : The step_poss = " << step_poss
    << " is not in range of 1 to " << steps_.size() + past_end );
  return step_poss;
}	

Algorithm::poss_type Algorithm::validate(const assoc_steps_ele_list_t& assoc_list
  , poss_type assoc_step_poss, int past_end) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    assoc_step_poss < 1 || assoc_list.size() + past_end < assoc_step_poss, DoesNotExist
    ,"Algorithm::validate(assoc_list,assoc_step_poss) : The assoc_step_poss = "
    << assoc_step_poss << " is not in range of 1 to " << assoc_list.size() + past_end );
  return assoc_step_poss;
}

void Algorithm::look_for_interrupt()
{
  //
  // Get the mode of aborting from the user!
  //
  if( static_interrupt_called && !static_processed_user_interrupt && static_proc_rank == 0 ) {
    // Allow for another interrupt possibly
    static_interrupt_called = false;
    //
    // Get the response from the user
    //
    enum EResponse { R_ABORT_NOW, R_CONTINUE, R_STOP_END_STEP, R_STOP_END_ITER };
    EResponse response = R_ABORT_NOW;
    const int max_tries = 3;
    bool valid_response = false;
    for( int tries = 0; !valid_response && tries < max_tries; ++tries ) {
      std::cerr
        << "\nIterationPack::Algorithm: Received signal SIGINT."
        << "\nJust completed current step curr_step_name = \""
        << get_step_name(curr_step_poss_) << "\",  curr_step_poss = "
        << curr_step_poss_ << " of steps [1..." << num_steps() << "]."
        << "\nDo you want to:\n"
        << "  (a) Abort the program immediately?\n"
        << "  (c) Continue with the algorithm?\n"
        << "  (s) Gracefully terminate the algorithm at the end of this step?\n"
        << "  (i) Gracefully terminate the algorithm at the end of this iteration?\n"
        << "Answer a, c, s or i ? ";
      char abort_mode = 'a';
      std::cin >> abort_mode;
      if( abort_mode == 'a' ) {
        response = R_ABORT_NOW;
        valid_response = true;
      }
      else if( abort_mode == 'c' ) {
        response = R_CONTINUE;
        valid_response = true;
      }
      else if( abort_mode == 's' || abort_mode == 'i' ) {
        if( abort_mode == 's')
          response = R_STOP_END_STEP;
        else
          response = R_STOP_END_ITER;
        std::cerr
          << "\nTerminate the algorithm with true (t) or false (f) ? ";
        std::cin >> abort_mode;
        if( abort_mode == 't' ) {
          static_interrupt_terminate_return = true;
          valid_response = true;
        }
        else if( abort_mode == 'f' ) {
          static_interrupt_terminate_return = false;
          valid_response = true;
        }
        else {
          std::cerr	<< "Invalid response! Expecting \'t\' or \'f\'\n";
        }
      }
      else {
        std::cerr	<< "\nInvalid response! Expecting \'a\', \'c\', \'s\' or \'i\'\n";
      }
      std::cerr << std::endl;
    }
    if(!valid_response) {
      std::cerr << "Three strikes, you are out!\n";
    }
    //
    // Interpret the response
    //
    switch(response) {
      case R_ABORT_NOW: {
        static_interrupt_status = ABORT_PROGRAM;
        break;
      }
      case R_CONTINUE: {
        static_interrupt_status = NOT_INTERRUPTED;
        break;
      }
      case R_STOP_END_STEP: {
        static_interrupt_status = STOP_END_STEP;
        break;
      }
      case R_STOP_END_ITER: {
        static_interrupt_status = STOP_END_ITER;
        break;
      }
      default: {
        TEUCHOS_TEST_FOR_EXCEPT(true);
      }
    }
    static_processed_user_interrupt = true;
  }
  else if( interrupt_file_name().length() && !static_processed_user_interrupt && static_proc_rank == 0 ) {
    //
    // If there was not an interactive interrupt then look for an
    // interrupt file if we have not already done this
    // (static_processed_user_interrupt).
    //
    std::ifstream interrupt_file(interrupt_file_name().c_str());
    if(interrupt_file) {
      std::cerr
        << "\nIterationPack::Algorithm: Found the interrupt file \""<<interrupt_file_name()<<"\"!"
        << "\nJust completed current step curr_step_name = \""
        << get_step_name(curr_step_poss_) << "\",  curr_step_poss = "
        << curr_step_poss_ << " of steps [1..." << num_steps() << "].\n";
      char abort_mode = 0;
      interrupt_file >> abort_mode;
      std::cerr	<< "Read a value of abort_mode = \'"<<abort_mode<<"\': ";
      if( abort_mode == 'a' ) {
        std::cerr	<< "Will abort the program immediatly!\n";
        static_interrupt_status = ABORT_PROGRAM;
      }
      else if( abort_mode == 's' || abort_mode == 'i' ) {
        if( abort_mode == 's') {
          std::cerr	<< "Will abort the program gracefully at the end of this step!\n";
          static_interrupt_status = STOP_END_STEP;
        }
        else {
          std::cerr	<< "Will abort the program gracefully at the end of this iteration!\n";
          static_interrupt_status = STOP_END_ITER;
        }
        TEUCHOS_TEST_FOR_EXCEPTION(
          interrupt_file.eof(), std::logic_error,
          "IterationPack::Algorithm: Error, expected input for terminate_bool option from the "
          "file \""<<interrupt_file_name()<<"\"!"
          );
        char terminate_bool = 0;
        interrupt_file >> terminate_bool;
        std::cerr	<< "Read a value of terminate_bool = \'"<<terminate_bool<<"\': ";
        if( terminate_bool == 't' ) {
          std::cerr	<< "Will return a success flag!\n";
          static_interrupt_terminate_return = true;
        }
        else if( terminate_bool == 'f' ) {
          std::cerr	<< "Will return a failure flag!\n";
          static_interrupt_terminate_return = false;
        }
        else {
          TEUCHOS_TEST_FOR_EXCEPTION(
            true, std::logic_error
            ,"Error, the value of terminate_bool = \'"<<terminate_bool<<"\' is not "
            "valid!  Valid values include only \'t\' or \'f\'\n"
            );
        }
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error
          ,"Error, the value of abort_mode = \'"<<abort_mode<<"\' is not "
          "valid!  Valid values include only \'a\', \'s\' or \'i\'\n"
          );
      }
      std::cerr << std::endl;
      static_processed_user_interrupt = true;
    }
  }
  //
  // Make sure that all of the processes get the same
  // response
  //
#ifdef HAVE_MPI
  const bool query_for_interrupt = true; // ToDo: Make this an external option!
  if( static_num_proc > 1 && query_for_interrupt ) {
    //
    // Here we will do a global reduction to see of a processor has
    // recieved an interrupt.  Here we will do a sum operation since only the
    // root process should be getting these options.
    //
    int sendbuf[2] = { 0, 0 };
    int recvbuf[2] = { 0, 0 };
    if(static_proc_rank == 0) {
      sendbuf[0] = (int)static_interrupt_status;
      sendbuf[1] = static_interrupt_terminate_return ? 1 : 0;
    }
    // Note: this global reduction will synchronize all of the processors!
#ifdef ITERATION_PACK_ALGORITHM_SHOW_MPI_DEBUG_INFO
    std::cerr	<< "\np="<<static_proc_rank<<": IterationPack::Algorithm::interrupt(): Calling MPI_Allreduce(...) ...\n";
#endif
    MPI_Allreduce(
      sendbuf                  // sendbuf
      ,recvbuf                 // recvbuf
      ,2                       // count
      ,MPI_INT                 // datatype
      ,MPI_SUM                 // op
      ,MPI_COMM_WORLD          // comm (ToDo: Make more general?)
      );
#ifdef ITERATION_PACK_ALGORITHM_SHOW_MPI_DEBUG_INFO
    std::cerr
      << "\np="<<static_proc_rank<<": IterationPack::Algorithm::interrupt(): After MPI_Allreduce(...)"
      << "\np="<<static_proc_rank<<": recvbuf[0] = " << recvbuf[0] << ", recvbuf[1] = " << recvbuf[1] << std::endl;
#endif
    // Set static_interrupt_status
    switch( (EInterruptStatus)recvbuf[0] ) {
      case NOT_INTERRUPTED:
        static_interrupt_status = NOT_INTERRUPTED;
        break;
      case STOP_END_STEP:
        static_interrupt_status = STOP_END_STEP;
        break;
      case STOP_END_ITER:
        static_interrupt_status = STOP_END_ITER;
        break;
      case ABORT_PROGRAM:
        static_interrupt_status = ABORT_PROGRAM;
        break;
      default:
        std::cerr
          << "p=" << static_proc_rank << ": Algorithm::look_for_interrupt(): Error, the globally reduced value of "
          "recvbuf[0] = " << recvbuf[0] << " is not valid!";
        std::abort();
    }
    // Set static_interrupt_terminate_return
    static_interrupt_terminate_return = ( recvbuf[1] == 0 ? false : true );
  }
  //
  // Abort the program now if the user did not already press Ctrl-C again!
  //
  if( static_interrupt_status == ABORT_PROGRAM ) {
    if( static_proc_rank == 0 ) {
      std::cerr << "\nAborting the program now!\n";
    }
    std::abort();
  }
#endif
}

// static

void Algorithm::interrupt()
{
  //
  // This function assumes that every process will recieve the same
  // signal which I found to be the case with MPICH.  I am not clear
  // what the MPI standard says about interrupts so I can not
  // guarantee that this is 100% portable.  If other behavior is
  // needed, this will have to be compiled in differently.
  //
  // Note: I have found that on MPICH that you can not guarantee that
  // only a single signal will be sent to a slave process so this
  // function will ignore interupts for slave processes.
  //
  // Note that you have to be very careful what you do inside of a
  // signal handler and in general you should only be setting flags or
  // aborting.
  //
  static_processed_user_interrupt = false;
#ifdef ITERATION_PACK_ALGORITHM_SHOW_MPI_DEBUG_INFO
  std::cerr	<< "\np="<<static_proc_rank<<": IterationPack::Algorithm::interrupt() called!\n";
#endif
  //
  // See if an algorithm is possibly even running yet!
  //
  if( static_num_proc == 0 ) {
    if( static_proc_rank == 0 )
      std::cerr
        << "\nIterationPack::Algorithm::interrupt(): Received signal SIGINT but an Algorithm "
        << "object has not been allocated yet and no algorithm is running.\n"
        << "\nAborting the program now!\n";
    std::abort();
    return;  // Should not be called!
  }
  //
  // See if we are going to query for an interrupt when running in MPI mode
  //
  const bool query_for_interrupt = true; // ToDo: Make this an external option!
  if( !query_for_interrupt && static_num_proc > 1 ) {
    if( static_proc_rank == 0 ) {
      std::cerr
        << "\nIterationPack::Algorithm::interrupt(): Received signal SIGINT but num_proc = "
        << static_num_proc << " > 1 and query_for_interrupt = false so:\n"
        << "\nAborting the program now!\n";
    }
    std::abort();
    return;  // Should not be called!
  }
  //
  // Remember that this interrupt has been called!
  //
  if( static_proc_rank == 0 ) {
    std::cerr
      << "\nIterationPack::Algorithm::interrupt(): Received signal SIGINT.  "
      << "Wait for the end of the current step and respond to an interactive query,  "
      << "kill the process by sending another signal (i.e. SIGKILL).\n";
  }
  static_interrupt_called = true;
}

} // end namespace IterationPack
