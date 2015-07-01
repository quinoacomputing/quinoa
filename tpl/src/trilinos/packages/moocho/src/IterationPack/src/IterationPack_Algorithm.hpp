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

#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <assert.h>

#include <string>
#include <deque>
#include <list>
#include <vector>
#include <algorithm>

#include "IterationPack_AlgorithmState.hpp"
#include "IterationPack_AlgorithmTracker.hpp"
#include "IterationPack_AlgorithmStep.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace IterationPack {

// ToDo: 7/31/98: Finish documentation.

/** \brief Acts as the central hub for an iterative algorithm.
  *
  * This class is the center for a framework for iterative algorithms.
  * These iterative algorithms are of the form:
  *
  \verbatim

  Step1------>Step2------>Step3------>Step4------>Step5
   /|\         /|\         /|\          |           |
    |           |           |_Minor L 1_|           |
    |           |                       |           |
    |           |_____Minor Loop 2______|           |
    |                                               |
    |_______________Major Loop (k = k+1)____________|
  \endverbatim
  *
  * For the typical iteration the steps are executed sequantially from Step1 to Step2
  * and then control loops around the Major Loop to Step1 again.
  * Durring some iterations however
  * Minor Loop 1 may be executed several times before control is continued alone the
  * major loop.  The same may also apply to Minor Loop 2.
  *
  * To allow for greater algorithmic control any step object can take over the role of
  * <tt>Algorithm</tt> and take over complete control the algorithm.  For examle, Step 4 may
  * need to Execute Step1->Step3->Step2->Step5 before returning algorithmic control to
  * <tt>Algorithm</tt>.
  *
  * <tt>Algorithm</tt> executes the steps of the algorithm through step objects of the
  * base type <tt>AlgorithmStep</tt>.  In addition to major step objects as shown above
  * there are also PreStep and PostStep objects.
  * These are steps that are intimatly associated with a major step object and
  * will always (well almost always) be exectuted alone with a major step.
  *
  * ToDo: Finish documentation.
  *
  * These functions provide information as to the number of major steps
  * , their possitions given their names and their names given their
  * possitions.
  *
  * In addition, access is given to the step objects themselves through
  * the RCP<...> objects that are used to manage thier memory.
  * Using this type of direct access allows clients to take over memory
  * management if needed and to call the step objects in any order
  * and thereby taking over control of the algorithm.
  *
  * These functions can be invoked in any state of the algorithm.
  *
  * Note that the currently running algorithm can always be interupted
  * by invoking the SIGINT signal (i.e. Ctrl-C) from the console.
  */
class Algorithm {
public:

  /** @name Public types */
  //@{

  /** \brief . */
  typedef Teuchos::RCP<AlgorithmState>      state_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<AlgorithmTracker>    track_ptr_t;
  /** \brief . */
  typedef Teuchos::RCP<AlgorithmStep>       step_ptr_t;
  /** \brief . */
  typedef size_t                                    poss_type;
  /** \brief . */
  enum { DOES_NOT_EXIST = 1000 };  // never be that many steps
  /** \brief . */
  enum ERunningState { NOT_RUNNING = 0, RUNNING = 1, RUNNING_BEING_CONFIGURED = 2 };

  /// Thrown if name or id does not exist
  class DoesNotExist : public std::logic_error
  {public: DoesNotExist(const std::string& what_arg) : std::logic_error(what_arg) {}};

  /// Thrown if name already exists
  class AlreadyExists : public std::logic_error
  {public: AlreadyExists(const std::string& what_arg) : std::logic_error(what_arg) {}};

  /// Thrown if an invalid control protocal is used.
  class InvalidControlProtocal : public std::logic_error
  {public: InvalidControlProtocal(const std::string& what_arg) : std::logic_error(what_arg) {}};

  /// Thrown if a member function is called while <tt>this</tt> is in an invalid running state..
  class InvalidRunningState : public std::logic_error
  {public: InvalidRunningState(const std::string& what_arg) : std::logic_error(what_arg) {}};

  /// Thrown if a member function is called while <tt>this</tt> is in an invalid running state..
  class InvalidConfigChange : public std::logic_error
  {public: InvalidConfigChange(const std::string& what_arg) : std::logic_error(what_arg) {}};

  /// Thrown if Algorithm was interrupted by the user
  class AlgorithmInterrupted : public std::runtime_error
  {public: AlgorithmInterrupted(const std::string& what_arg) : std::runtime_error(what_arg) {}};

  //@}

  /** @name Constructors & destructors */
  //@{

  /** \brief Constructs an algorithm with no steps and a default of max_iter() == 100.
    *
    */
  Algorithm();

  /** \brief . */
  virtual ~Algorithm();

  //@}

  /** \brief Name of an file that will cause the algorithm to terminate.
   *
   * If <tt>interrupt_file_name()!=""</tt> then the file <tt>interrupt_file_name()</tt>
   * will be looked for in the current directory and if it exists, it will be read
   * to determine how the algorithm should be terminated.  The format of this file
   * is as follows:
   *
   \verbatim
   
   abort_mode  terminate_bool

   \endverbatim
   *
   * Above, <tt>abort_mode</tt> can be:
   * <tt>a</tt> for "abort the program immediately",
   * <tt>s</tt> for "Gracefully terminate the algorithm at the end of this step",
   * <tt>i</tt> for "Gracefully terminate the algorithm at the end of this iteration"
   *
   * If the option <tt>abort_mode</tt> is set to <tt>s</tt> or <tt>i</tt> then the the
   * value of <tt>terminate_bool</tt> must be set to <tt>t</tt> for "true" or
   * <tt>f</tt> for "false".  If this <tt>abort_mode</tt> is set ot <tt>a</tt> then the
   * value of <tt>terminate_bool</tt> is not read.
   *
   * Note that the option values <tt>abort_mode</tt> and
   * <tt>terminate_bool</tt> are simple <tt>char</tt> data objects and
   * the only requirement is that they be seperated by whitespace.
   *
   * If the format of this file is not as described above, then an
   * exception will be thrown (which will have the affect of aborting
   * the process most likely).
   *
   * The default is <tt>interrupt_file_name()==""</tt> which means
   * that this file will not be looked for.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( std::string, interrupt_file_name );

  /** @name «std comp» members for state */
  //@{

  /** \brief . */
  void set_state(const state_ptr_t& state);
  /** \brief . */
  state_ptr_t& get_state();
  /** \brief . */
  const state_ptr_t& get_state() const;
  /** \brief . */
  AlgorithmState& state();
  /** \brief . */
  const AlgorithmState& state() const;

  //@}

  /** @name «std comp» members for track */
  //@{

  /** \brief . */
  void set_track(const track_ptr_t& track);
  /** \brief . */
  track_ptr_t& get_track();
  /** \brief . */
  const track_ptr_t& get_track() const;
  /** \brief . */
  AlgorithmTracker& track();
  /** \brief . */
  const AlgorithmTracker& track() const;

  //@}

  /** @name Maximum iterations */
  //@{
  
  /** \brief . */
  virtual void max_iter(size_t max_iter);
  /** \brief . */
  virtual size_t max_iter() const;

  //@}

  /** @name Maximum runtime (in minutes) */
  //@{
  
  /** \brief Set the maximum runtime (in minues)
    * The runtime is checked at the end of each iteration and if it exceeds
    * this value then the algorithm is terminated.
    */
  virtual void max_run_time(double max_iter);
  /** \brief . */
  virtual double max_run_time() const;

  //@}

  /** @name Step information & access */
  //@{

  /// Return the number of main steps
  virtual int num_steps() const;

  /** \brief Return the possition in the major loop of a named step.
    *
    * If a step with this name does not exist then the value
    * DOES_NOT_EXIST will be returned.
    */
  virtual poss_type get_step_poss(const std::string& step_name) const;

  /** \brief Return the name of a step given its possition.
    *
    * Preconditions:<ul>
    * <li> <tt>1 <= step_poss && step_poss <= num_steps()</tt> (throw <tt>DoesNotExist</tt>)
    * </ul>
    */
  virtual const std::string& get_step_name(poss_type step_poss) const;

  /** \brief Return the RCP<...> object for the step object at step_poss.
    *
    * Preconditions:<ul>
    * <li> <tt>1 <= step_poss && step_poss <= num_steps()</tt> (throw <tt>DoesNotExist</tt>)
    * </ul>
    */
  virtual step_ptr_t& get_step(poss_type step_poss);

  /** \brief . */
  virtual const step_ptr_t& get_step(poss_type step_poss) const;

  //@}

  /** @name Pre/post step information & access */
  //@{

  /** \brief Return the number of pre or post steps for the main step step_poss.
    *
    * Preconditions:<ul>
    * <li> <tt>1 <= step_poss && step_poss <= num_steps()</tt> (throw <tt>DoesNotExist</tt>)
    * </ul>
    */
  virtual int num_assoc_steps(poss_type step_poss, EAssocStepType type) const;

  /** \brief Return the possition of the pre or post step for the main step_poss.
    *
    * If a pre or post step does not exist with the name <tt>assoc_step_name</tt>
    * then a value of DOES_NOT_EXIST will be retruned.
    *
    * Preconditions:<ul>
    * <li> <tt>1 <= step_poss && step_poss <= num_steps()</tt> (throw <tt>DoesNotExist</tt>)
    * </ul>
    */
  virtual poss_type get_assoc_step_poss(poss_type step_poss, EAssocStepType type
    ,const std::string& assoc_step_name) const;

  /** \brief Return the name of the pre or post step at step_poss and at assoc_step_poss.
    *
    * Preconditions:<ul>
    * <li> <tt>1 <= step_poss && step_poss <= num_steps()</tt> (throw <tt>DoesNotExist</tt>)
    * <li> <tt>1 <= assoc_step_poss && assoc_step_poss <= num_assoc_steps(step_poss,type)</tt>
    *    (throw <tt>DoesNotExist</tt>)
    * </ul>
    */
  virtual const std::string& get_assoc_step_name(poss_type step_poss, EAssocStepType type
    , poss_type assoc_step_poss) const;

  /** \brief Return the RCP<...> object for the associated step object at step_poss
    * and assoc_step_poss.
    *
    * Preconditions:<ul>
    * <li> <tt>1 <= step_poss && step_poss <= num_steps()</tt> (throw <tt>DoesNotExist</tt>)
    * <li> <tt>1 <= assoc_step_poss && assoc_step_poss <= num_assoc_steps(step_poss,type)</tt>
    *    (throw <tt>DoesNotExist</tt>)
    * </ul>
    */
  virtual step_ptr_t& get_assoc_step(poss_type step_poss, EAssocStepType type
    , poss_type assoc_step_poss);

  /** \brief . */
  virtual const step_ptr_t& get_assoc_step(poss_type step_poss, EAssocStepType type
    , poss_type assoc_step_poss) const;

  //@}

  /** @name Step manipulation */
  //@{

  /** \brief Insert a step object with the name <tt>step_name</tt> into the possition <tt>step_poss</tt>.
    *
    * All the steps at and after <tt>step_poss</tt> are pushed back one possition unless
    * <tt>step_poss == num_steps() + 1</tt> in which case the new step is appended to the end.
    * Initiaily this step will have no pre or post steps associated with it.
    * 
    * Preconditions:<ul>
    * <li> <tt>running_state() != RUNNING]</tt> (throw <tt>InvalidRunningState</tt>)
    * <li> <tt>1 <= ste_poss && step_poss <= num_steps() + 1</tt> (throw <tt>DoesNotExist</tt>)
    * <li> <tt>step.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
    * </ul>
    */
  virtual void insert_step(poss_type step_poss, const std::string& step_name, const step_ptr_t& step);

  /** \brief Change the name of an existing step.
    *
    * None of the pre or post steps for the existing step are changes.
    *
    * Preconditions:<ul>
    * <li> <tt>running_state() != RUNNING]</tt> (throw <tt>InvalidRunningState</tt>)
    * <li> <tt>1 <= poss && poss <= num_steps()</tt> (throw <tt>DoesNotExist</tt>)
    * </ul>
    */
  virtual void change_step_name(poss_type step_poss, const std::string& new_name);

  /** \brief Replace the step object of an existing step.
    *
    * None of the pre or post steps for the existing step are changes.
    *
    * Preconditions:<ul>
    * <li> <tt>running_state() != RUNNING]</tt> (throw <tt>InvalidRunningState</tt>)
    * <li> <tt>1 <= poss && poss <= num_steps()</tt> (throw <tt>DoesNotExist</tt>)
    * </ul>
    */
  virtual void replace_step(poss_type step_poss, const step_ptr_t& step);

  /** \brief Remove an existing step object and all of its pre and post steps.
    *
    * All of the steps after <tt>step_poss</tt> will have thier possitions
    * decreased by one.
    *
    * Preconditions:<ul>
    * <li> <tt>running_state() != RUNNING]</tt> (throw <tt>InvalidRunningState</tt>)
    * <li> <tt>1 <= poss && poss <= num_steps()</tt> (throw <tt>DoesNotExist</tt>)
    * </ul>
    */
  virtual void remove_step(poss_type step_poss);

  //@}

  /** @name Pre/post step manipulation */
  //@{

  /** \brief Insert an pre or post step into for the main step step_poss into the possition
    * assoc_step_poss.
    *
    * All of the pre or post steps at and after <tt>assoc_step_poss</tt> will be pushed back
    * one possition.
    *
    * Preconditions:<ul>
    * <li> <tt>running_state() != RUNNING]</tt> (throw <tt>InvalidRunningState</tt>)
    * <li> <tt>1 <= step_poss && step_poss <= num_steps() + 1</tt> (throw <tt>DoesNotExist</tt>)
    * <li> <tt>1 <= assoc_step_poss && assoc_step_poss <= num_assoc_steps(step_poss,type) + 1</tt>
    *    (throw <tt>DoesNotExist</tt>)
    * <li> <tt>assoc_step.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
    * </ul>
    */
  virtual void insert_assoc_step(poss_type step_poss, EAssocStepType type, poss_type assoc_step_poss
    , const std::string& assoc_step_name, const step_ptr_t& assoc_step);

  /** \brief Remove an pre or post step for the main step step_poss in the possition
    * assoc_step_poss.
    *
    * All of the pre or post steps after <tt>assoc_step_poss</tt> will be pushed forward
    * one possition to fill in the hole.
    *
    * Preconditions:<ul>
    * <li> <tt>running_state() != RUNNING]</tt> (throw <tt>InvalidRunningState</tt>)
    * <li> <tt>1 <= step_poss && step_poss <= num_steps()</tt> (throw <tt>DoesNotExist</tt>)
    * <li> <tt>1 <= assoc_step_poss && assoc_step_poss <= num_assoc_steps(step_poss,type)</tt>
    *    (throw <tt>DoesNotExist</tt>)
    * </ul>
    */
  virtual void remove_assoc_step(poss_type step_poss, EAssocStepType type, poss_type assoc_step_poss);

  //@}

  /** @name Runtime configuration updating control */
  //@{

  /// Return the current running state of \c this algorithm object.
  ERunningState running_state() const;

  /** \brief Changes from running_state() == RUNNING to running_state() == RUNNING_BEING_CONFIGURED.
    *
    * Must be called before the algorithm's configuration can be changed while it is running.
    *
    * Preconditions:<ul>
    * <li> <tt>running_state() == RUNNING]</tt> (throw <tt>InvalidRunningState</tt>)
    * </ul>
    *
    * Postconditions:<ul>
    * <li> <tt>running_state() == RUNNING_BEING_CONFIGURED]</tt>
    * </ul>
    */
  virtual void begin_config_update();

  /** \brief Changes from running_state() == RUNNING_BEING_CONFIGURED to running_state() == RUNNING.
    *
    * Must be called after the algorithm's configuration can be changed while it is running.
    *
    * Preconditions:<ul>
    * <li> <tt>running_state() == RUNNING_BEING_CONFIGURED]</tt> (throw <tt>InvalidRunningState</tt>)
    * </ul>
    *
    * Postconditions:<ul>
    * <li> <tt>running_state() == RUNNING]</tt>
    * </ul>
    */
  virtual void end_config_update();

  //@}

  /** @name Algorithmic control */
  //@{

  /** \brief Called by step objects to set the step (given its name) that <tt>this</tt> will envoke the next time
    * <tt>this</tt> calls a step.
    *
    * Preconditions:<ul>
    * <li> <tt>running_state() == RUNNING<tt> (throw <tt>InvalidRunningState</tt>)
    * <li> <tt>get_step_poss(step_name) != DOES_NOT_EXIST</tt> (throw <tt>DoesNotExist</tt>)
    * </ul>
    */
  virtual void do_step_next(const std::string& step_name);

  /** \brief Called by step objects to set the step (given its possition) that <tt>this</tt> will envoke the next time
    * <tt>this</tt> calls a step.
    *
    * Preconditions:<ul>
    * <li> <tt>running_state() == RUNNING</tt> (throw <tt>InvalidRunningState</tt>)
    * <li> <tt>1 <= step_poss && step_poss <= num_steps()</tt> (throw <tt>DoesNotExist</tt>)
    * </ul>
    */
  virtual void do_step_next(poss_type step_poss);

  /** \brief Returns the name of the next step <tt>this</tt> will call the next time it calls a step.
    *
    * Preconditions:<ul>
    * <li> running_state() != NOT_RUNNING</tt> (throw <tt>InvalidRunningState</tt>)
    * </ul>
    */
  virtual const std::string& what_is_next_step_name() const;

  /** \brief Returns the possition of the next step <tt>this</tt> will call the next time it calls a step.
    *
    * Preconditions:<ul>
    * <li> running_state() != NOT_RUNNING</tt> (throw <tt>InvalidRunningState</tt>)
    * </ul>
    */
  virtual poss_type what_is_next_step_poss() const;

  /** \brief Calls <tt>do_step()</tt> on all of the pre step objects the step object and the post step objects
    * in order for the step named <tt>step_name</tt>.
    *
    * This operation is called by step objects that need to take over control of the algorithm
    * at some point.
    *
    * If any of the of the pre or post objects or the step object returns false, then this
    * operation immediatly returns false.  It is assumed that if any step object returns
    * false from its <tt>do_step()</tt> that it has either also called <tt>terminate()</tt>
    * or <tt>do_step_next()</tt>.
    *
    * Preconditions:<ul>
    * <li> <tt>running_state() == RUNNING</tt> (throw <tt>InvalidRunningState</tt>)
    * <li> <tt>get_step_poss(step_name) != DOES_NOT_EXIST</tt> (throw <tt>DoesNotExist</tt>)
    * </ul>
    */
  virtual bool do_step(const std::string& step_name);

  /** \brief Call <tt>do_step()</tt> on all of the pre step objects the step object and the post step objects
    * in order for the step in the possition <tt>step_poss</tt>.
    *
    * This operation is called by step objects that need to take over control of the algorithm
    * at some point.
    *
    * If any of the of the pre or post objects or the step object returns false, then this
    * operation immediatly returns false.  It is assumed that if any step object returns
    * false from do step that it has either also called <tt>terminate()</tt> or <tt>do_step_next()</tt>. 
    *
    * Preconditions:<ul>
    * <li> <tt>running_state() == RUNNING</tt> (throw <tt>InvalidRunningState</tt>)
    * <li> <tt>1 <= step_poss && step_poss <= num_steps()</tt> (throw <tt>DoesNotExist</tt>)
    * </ul>
    */
  virtual bool do_step(poss_type step_poss);

  /** \brief Called by step objects to terminate the algorithm.
    *
    * Calling with <tt>success == true</tt> cause <tt>do_algorithm()</tt> to completely return <tt>TERMINATE_TRUE</tt>
    * and with <tt>success == false</tt> return  <tt>TERMINATE_FALSE</tt>.
    *
    * Preconditions:<ul>
    * <li> <tt>running_state() == RUNNING</tt> (throw <tt>InvalidRunningState</tt>)
    * </ul>
    */
  virtual void terminate(bool success);

  //@}

  /** @name Start iterations */
  //@{

  /** \brief Called by clients to begin an algorithm.
    *
    * Preconditions:<ul>
    * <li> <tt>this->get_track() != NULL</tt> (throw <tt>???</tt>)
    * <li> <tt>running_state() == NOT_RUNNING</tt> (throw <tt>InvalidRunningState</tt>)
    * <li> <tt>1 <= step_poss && step_poss <= num_steps()</tt> (throw <tt>DoesNotExist</tt>)
    * </ul>
    *
    * This operation acts as the central hub for the algorithm.  It calls the <tt>do_step(i)</tt>
    * each <tt>i</tt> = 1,...,<tt>num_steps()</tt> and then loops around again for the major loop.  If
    * <tt>do_step(i)</tt> returns false then it goes executes the step specified by the 
    * <tt>do_step_next()</tt> operation which the step object supposivly called.  If a step
    * object returns false but does not call <tt>do_step_next()</tt> to specify a step to
    * jump to, then <tt>this</tt> will throw an <tt>InvalidControlProtocal</tt> exception.
    *
    * Before the algorithm is started, <tt>this</tt> calls <tt>track().initialize()</tt>.
    * At the end of each iteration <tt>this</tt> calls <tt>track().output_iteration(*this)</tt> and
    * <tt>state().next_iteration()</tt>.  It then checks if <tt>state.k() - k_start</tt> >= <tt>max_iter()</tt>.
    * If it is then the <tt>do_algorithm()</tt> immediatly terminates with a value of
    * <tt>MAX_ITER_EXCEEDED</tt> after it passes <tt>MAX_ITER_EXCEEDED</tt> to <tt>track().output_final()</tt>.
    * If the maxinum runtime <tt>this->max_run_time()</tt> is exceeded then <tt>MAX_RUN_TIME_EXCEEDED</tt>
    * will be passed to <tt>track().output_final()</tt> and this function will return
    * <tt>track().output_final()</tt>.  If the algorithm throws any exception then
    * <tt>track().output_final()</tt> will be called with the value <tt>TERMINATE_FALSE</tt>
    * and this exception will be rethrown clean out of here.
    *
    * Any step object can cause the algorithm to terminate by calling <tt>terminate(success)</tt>.
    * This operation will then immediatly return <tt>TERMINATE_TRUE</tt> if <tt>success == true</tt>
    * and <tt>TERMINATE_FALSE</tt>  if <tt>success == false</tt>.
    *
    * The algorithm starts on the step specified with <tt>step_poss</tt>.
    *
    */
  virtual EAlgoReturn do_algorithm(poss_type step_poss = 1);

  //@}

  /** @name Algorithm information output */
  //@{

  /** \brief Print out just a listing of the steps, their positions in the algorithm and the subclasses.
    */
  virtual void print_steps(std::ostream& out) const;

  /** \brief Print out the entire algorithm by calling <tt>print_step(...)</tt> on the step objects.
    */
  virtual void print_algorithm(std::ostream& out) const;

  //@}

  /** @name Algorithm Timing
    */
  //@{

  /** \brief Causes algorithm to be timed.
    *
    * Call with <tt>algo_timing == true</tt> before <tt>do_algorithm()</tt> to have the algorithm timed.
    *
    * Do not call when algorithm is running.
    */
  virtual void set_algo_timing( bool algo_timing );

  /** \brief . */
  virtual bool algo_timing() const;

  /** \brief Outputs table of times for each step, cummulative times and other
    * statistics.
    *
    * Call after <tt>do_algorithm()</tt> has executed to get a table
    * of times.
    *
    * Do not call when algorithm is running.
    */
  virtual void print_algorithm_times( std::ostream& out ) const;

  /** \brief Returns the step_times for iteration offset.
   *
   * @param  offset       [in] The interation offset to retrieve times for.
   * @param  step_times   [out] Array (size <tt>this->num_steps() + 1</tt>) with the
   *                      output step times (in seconds) for iteration <tt>k+offset</tt>.
   *                      The last entry <tt>step_times[this->num_steps()]</tt> gives
   *                      the total time for the entire iteration.
   *
   * Preconditions:<ul>
   * <li> <tt>offset <= 0</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Note that <tt>offset</tt> must be a nonpositive number since we can only retreve
   * timings from the current or previous iterations.
   */  
  void get_step_times_k( int offset, double step_times[] ) const;

  /** \brief Returns the final statistics for a given step
    *  Do not call when algorithm is running
    */
  void get_final_step_stats( size_t step, double* total, double* average, double* min, double* max, double* percent) const;

  //@}

private:

  // /////////////////////////////////////////////////////
  // Private types

  /** \brief . */
  template<class T_Step_ptr>
  struct step_ptr_and_name {
    /** \brief . */
    step_ptr_and_name(const T_Step_ptr& _step_ptr
        , const std::string& _name )
      : step_ptr(_step_ptr), name(_name)
    {}
    /** \brief . */
    T_Step_ptr step_ptr;
    //
    std::string name;
  };  // end struct step_ptr_and_name

  /** \brief . */
  typedef step_ptr_and_name<step_ptr_t>           steps_ele_t;
  /** \brief . */
  typedef std::deque<steps_ele_t>                 steps_t;

  /** \brief . */
  typedef step_ptr_and_name<step_ptr_t>           assoc_steps_ele_list_ele_t;
  /** \brief . */
  typedef std::list<assoc_steps_ele_list_ele_t>   assoc_steps_ele_list_t;
  /** \brief . */
  struct assoc_steps_ele_t {
    /** \brief . */
    assoc_steps_ele_list_t& operator[](int i)
    { return assoc_steps_lists_[i]; }
    /** \brief . */
    const assoc_steps_ele_list_t& operator[](int i) const
    { return assoc_steps_lists_[i]; }
  private:
    assoc_steps_ele_list_t assoc_steps_lists_[2];
  };

  //typedef assoc_steps_ele_list_t[2]        assoc_steps_ele_t; // PRE_STEP, POST_STEP
  /** \brief . */
  typedef std::deque<assoc_steps_ele_t>      assoc_steps_t;
  
  /** \brief . */
  enum ETerminateStatus { STATUS_KEEP_RUNNING, STATUS_TERMINATE_TRUE, STATUS_TERMINATE_FALSE };

  /** \brief . */
  template<class T_ele>
  class name_comp {
  public:
    /** \brief . */
    name_comp(const std::string& name) : name_(name) {}
    /** \brief . */
    bool operator()(const T_ele& ele) { return ele.name == name_; }
  private:
    const std::string& name_;
  };  // end class name_comp

  typedef std::vector<double> step_times_t;
  
  /** \brief . */
  static const int
    TIME_STAT_TOTALS_OFFSET  = 0,
    TIME_STAT_AV_OFFSET    = 1,
    TIME_STAT_MIN_OFFSET    = 2,
    TIME_STAT_MAX_OFFSET    = 3,
    TIME_STAT_PERCENT_OFFSET  = 4;
  /** \brief . */
  enum { NUM_STEP_TIME_STATS = 5 };

  /** \brief . */
  typedef void (AlgorithmStep::*inform_func_ptr_t)(
    Algorithm&             algo
    ,poss_type             step_poss
    ,EDoStepType           type
    ,poss_type             assoc_step_poss
    );

  // /////////////////////////////////////////////////////
  // Private data members

  // aggregate members

#ifdef DOXYGEN_COMPILE
  AlgorithmState       *state;
  AlgorithmTracker     *track;
  AlgorithmStep        *steps;
#else
  state_ptr_t        state_;
  // RCP<...> object for the aggragate AlgorithmState object.

  track_ptr_t        track_;
  // RCP<...> object for the aggragate AlgorithmTracker object.
#endif

  // algorithm control etc.
  
  ERunningState      running_state_;
  // The state this Algorithm object is in:
  //
  // NOT_RUNNING    do_algorithm() is not active.
  // RUNNING        do_algorithm() has been called and is active.
  // RUNNING_BEING_CONFIGURED
  //                do_algorithm() is active and begin_config_update() has been called
  //                but end_config_update() has not.
  //
  // Note: Only change this variable through the private function change_running_state(...)
  // unless you are 100% sure that you know what you are doing!

  size_t          first_k_;
  // The first iteration from state().k().
  
  size_t          max_iter_;
  // The maximum number of iterations that <tt>this</tt> will execute.

  double          max_run_time_;
  // The maximum amount of time the algorithm is allowed to execute.

  ETerminateStatus    terminate_status_;
  // Flag for if it is time to terminate do_algorithm().

  poss_type        next_step_poss_;
  // The next step possition that <tt>this</tt> will execute when control is returned to do_algorithm().

  const std::string*    next_step_name_;
  // The name of the next step that <tt>this</tt> will execute when control is returned to do_algorithm().
  
  bool          do_step_next_called_;
  // Flag for if do_step_next() was called so that <tt>do_algorithm()</tt> can validate
  // that if a step object returned <tt>false</tt> from its <tt>do_step()</tt> operation that it
  // must also call <tt>do_step_next()</tt> to specify a step to jump to.

  poss_type        curr_step_poss_;
  // The current step being executed in do_algorithm().
  // If the current step being executed is changed during the imp_do_step() operation, then
  // imp_do_step() must adjust to this step.

  std::string        saved_curr_step_name_;
  // The name of the current step that is saved when begin_config_update() is called
  // so that curr_step_poss_ can be reset when end_config_update() is called.

  std::string        saved_next_step_name_;
  // The name of the next step to call so that when begin_config_update() is called
  // so that next_step_poss_ and next_step_name_ can be reset when end_config_update()
  // is called.

  bool          reconfigured_;
  // A flag that is set to true when a runtime configuration has been preformed.  It
  // is used to flag this event for imp_do_assoc_steps().

  // step and associated step object data structures

  steps_t          steps_;
  // Array of std::pair<RCP<step_ptr_t>,std::string> objects.
  //
  // *steps_[step_poss].first returns the step object for step_poss = 1...steps_.size().
  // steps_[step_poss].second returns the name of the step for step_poss = 1...steps_.size().

  assoc_steps_t      assoc_steps_;
  // Array of two lists of std::pair<step_ptr_t,std::string> objects
  //
  // *(assoc_steps_[step_poss][PRE_STEP].begin() + pre_step_poss).first gives the pre step object.
  // (assoc_steps_[step_poss][PRE_STEP].begin() + pre_step_poss).second gives the name of the pre step
  // *(assoc_steps_[step_poss][POST_STEP].begin() + post_step_poss).first gives the post step object.
  // (assoc_steps_[step_poss][POST_STEP].begin() + post_step_poss).second gives the name of the post step

  bool algo_timing_;
  // If true each step will be timed.

  mutable step_times_t step_times_;
  // Array of step times ( size (max_iter() + 1 + NUM_STEP_TIME_STATS) * (num_steps() + 1) ).
  //  The time in sec. for step step_i (one based)
  // for iteration iter_k (zero based) is:
  //   step_times_[ iter_k + (step_i - 1) * (max_iter() + 1 + NUM_STEP_TIME_STATS) ].
  // Therefore the times for each step are stored by column (consecutive elements)
  // so that statistics will be easy to compute at the end.
  // The last five elements after max_iter() for each step are reserved for:
  // * total time for the step
  // * average time for the step
  // * min step time
  // * max step time
  // * percentage for each step to the total.
  // The last column is for the total times for each iteration with the last five
  // elements being for the statistics for each iteration.   

  mutable bool time_stats_computed_;
  // A flag for if the timing statistics have already been computed or not.
  
  mutable double total_time_;
  // Records the total computed time for the algorithm.

  // /////////////////////////////////////////////////////
  // Private member functions

  /// Validate a step_poss and throw a DoesNotExist exception if it does not.
  poss_type validate(poss_type step_poss, int past_end = 0) const;

  /// Validate an assoc_step_poss and throw a DoesNotExist exception if it does not.
  poss_type validate(const assoc_steps_ele_list_t& assoc_list, poss_type assoc_step_poss, int past_end = 0) const;

  /// Change the running state
  void change_running_state(ERunningState running_state);

  /// Validate that <tt>this</tt> is in a specific running state.
  void validate_in_state(ERunningState running_state) const;

  /// Validate that <tt>this</tt> is not in a specific running state.
  void validate_not_in_state(ERunningState running_state) const;

  /// Validate that the step_poss in not the current step.
  void validate_not_curr_step(poss_type step_poss) const;

  /// Validate that the step_name in not the next step.
  void validate_not_next_step(const std::string& step_name) const;

  /** \brief Find a step given its name and throw a DoesNotExist exception if not found.
    */
  steps_t::iterator step_itr(const std::string& step_name);

  /** \brief . */
  steps_t::const_iterator step_itr(const std::string& step_name) const;

  /** \brief Find a step given its name and throw a DoesNotExist exception if not found.
    */
  steps_t::iterator step_itr_and_assert(const std::string& step_name);

  /** \brief . */
  steps_t::const_iterator step_itr_and_assert(const std::string& step_name) const;

  /** \brief Find a an associated step given its name and throw a DoesNotExist exception if not found.
    */
  static assoc_steps_ele_list_t::iterator assoc_step_itr(assoc_steps_ele_list_t& assoc_list
    , const std::string& assoc_step_name);

  /** \brief . */
  static assoc_steps_ele_list_t::const_iterator assoc_step_itr(const assoc_steps_ele_list_t& assoc_list
    , const std::string& assoc_step_name);

  /** \brief . */
  bool imp_do_step(poss_type step_poss);

  /** \brief . */
  bool imp_do_assoc_steps(EAssocStepType type);

  /** \brief . */
  void imp_inform_steps(inform_func_ptr_t inform_func_ptr);

  /** \brief . */
  void imp_print_algorithm(std::ostream& out, bool print_steps) const;

  /// EAssocStepType -> EDoStepType
  EDoStepType do_step_type(EAssocStepType assoc_step_type);

  /** \brief . */
  EAlgoReturn finalize_algorithm( EAlgoReturn algo_return );

  /** \brief . */
  void compute_final_time_stats() const;

  // Look for interrup
  void look_for_interrupt();

public:

  // This is put here out of need.  Not for any user to call!
  static void interrupt();

};  // end class Algorithm

// //////////////////////////////////////////////////////////////////////////////////////////////////
// Inline member function definitions for Algorithm

// «std comp» members for state 

inline
void Algorithm::set_state(const state_ptr_t& state)
{  state_ = state; }

inline
Algorithm::state_ptr_t& Algorithm::get_state()
{  return state_; }

inline
const Algorithm::state_ptr_t& Algorithm::get_state() const
{  return state_; }

inline
AlgorithmState& Algorithm::state()
{  TEUCHOS_TEST_FOR_EXCEPT( !( state_.get() ) ); return *state_; }

inline
const AlgorithmState& Algorithm::state() const
{  TEUCHOS_TEST_FOR_EXCEPT( !( state_.get() ) ); return *state_; }

// «std comp» members for track 

inline
void Algorithm::set_track(const track_ptr_t& track)
{  track_ = track; }

inline
Algorithm::track_ptr_t& Algorithm::get_track()
{  return track_; }

inline
const Algorithm::track_ptr_t& Algorithm::get_track() const
{  return track_; }

inline
AlgorithmTracker& Algorithm::track()
{  TEUCHOS_TEST_FOR_EXCEPT( !( track_.get() ) ); return *track_; }

inline
const AlgorithmTracker& Algorithm::track() const
{  TEUCHOS_TEST_FOR_EXCEPT( !( track_.get() ) ); return *track_; }

// running state

inline
Algorithm::ERunningState Algorithm::running_state() const
{  return running_state_; }

// lookup iterator given name

inline
Algorithm::steps_t::iterator Algorithm::step_itr(const std::string& step_name)
{
  return std::find_if( steps_.begin() , steps_.end()
    , name_comp<steps_ele_t>(step_name) );
}

inline
Algorithm::steps_t::const_iterator Algorithm::step_itr(const std::string& step_name) const
{
  return std::find_if( steps_.begin() , steps_.end()
    , name_comp<steps_ele_t>(step_name) );
}

inline
Algorithm::assoc_steps_ele_list_t::iterator Algorithm::assoc_step_itr(
  assoc_steps_ele_list_t& assoc_list, const std::string& assoc_step_name)
{
  return std::find_if( assoc_list.begin() , assoc_list.end()
    , name_comp<assoc_steps_ele_list_ele_t>(assoc_step_name) );
}

inline
Algorithm::assoc_steps_ele_list_t::const_iterator Algorithm::assoc_step_itr(
  const assoc_steps_ele_list_t& assoc_list, const std::string& assoc_step_name)
{
  return std::find_if( assoc_list.begin() , assoc_list.end()
    , name_comp<assoc_steps_ele_list_ele_t>(assoc_step_name) );
}

inline
EDoStepType Algorithm::do_step_type(EAssocStepType assoc_step_type) {
  switch(assoc_step_type) {
    case PRE_STEP  : return DO_PRE_STEP;
    case POST_STEP  : return DO_POST_STEP;
  }
  TEUCHOS_TEST_FOR_EXCEPT( !( true ) );
  return DO_PRE_STEP;  // will never execute.
}

}  // end namespace IterationPack 

#endif // ALGORITHM_H
