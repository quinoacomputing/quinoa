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

#include "MoochoPack_MoochoSolver.hpp"

#include "MoochoPack_NLPAlgoConfigMamaJama.hpp"
#include "MoochoPack_NLPAlgoConfigIP.hpp"

#include "MoochoPack_NLPSolverClientInterfaceSetOptions.hpp"
#include "MoochoPack_NLPAlgoClientInterface.hpp"
#include "MoochoPack_NLPAlgoContainer.hpp"
#include "MoochoPack_NLPAlgoState.hpp"
#include "MoochoPack_MoochoTrackerSummaryStd.hpp"
#include "MoochoPack_MoochoTrackerConsoleStd.hpp"
#include "MoochoPack_MoochoTrackerStatsStd.hpp"
#include "IterationPack_AlgorithmTrackerComposite.hpp"
#include "NLPInterfacePack_NLPFirstOrder.hpp"
#include "NLPInterfacePack_NLPDirect.hpp"
#include "NLPInterfacePack_test_nlp_first_order.hpp"
#include "NLPInterfacePack_test_nlp_direct.hpp"
#include "OptionsFromStreamPack_OptionsFromStream.hpp"
#include "StopWatchPack_stopwatch.hpp"
#include "OptionsFromStreamPack_StringToIntMap.hpp"
#include "OptionsFromStreamPack_StringToBool.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_Utils.hpp"

namespace MoochoPack {

// Initialization and algorithm configuration

MoochoSolver::MoochoSolver(
  const std::string &options_file_name
  ,const std::string &extra_options_str
  )
  :reconfig_solver_(true)
  ,workspace_MB_(-1.0)
  ,obj_scale_(1.0)
  ,test_nlp_(true)
  ,print_algo_(true)
  ,algo_timing_(true)
  ,generate_stats_file_(false)
  ,print_opt_grp_not_accessed_(true)
  ,throw_exceptions_(false)
  ,output_file_tag_("")
  ,do_console_outputting_(true)
  ,do_summary_outputting_(true)
  ,do_journal_outputting_(true)
  ,do_algo_outputting_(true)
  ,configuration_(MAMA_JAMA)
  ,output_to_black_hole_(OUTPUT_TO_BLACK_HOLE_DEFAULT)
  ,file_context_postfix_("")
  ,file_proc_postfix_("")
{
  
  Teuchos::RCP<Teuchos::FancyOStream>
    defaultOut = Teuchos::VerboseObjectBase::getDefaultOStream();
  error_out_used_ = defaultOut;

  set_output_context("");

  commandLineOptionsFromStreamProcessor_.options_file_name_opt_name(
    "moocho-options-file");
  commandLineOptionsFromStreamProcessor_.options_file_name_opt_doc(
    "Set the name of the MOOCHO options file in the OptionsFromStream format."
    "  File is ignored if it does not exist!");
  commandLineOptionsFromStreamProcessor_.options_file_name(options_file_name);

  commandLineOptionsFromStreamProcessor_.extra_options_str_opt_name(
    "moocho-extra-options");
  commandLineOptionsFromStreamProcessor_.extra_options_str_opt_doc(
    "Extra MOOCHO options specified in the format"
    " \"OptGrp1{name1=val1,...,namen=valn}:OptGr2{name1=val1,...,namen=valn}:...\"");
  commandLineOptionsFromStreamProcessor_.extra_options_str(extra_options_str);

}

OptionsFromStreamPack::CommandLineOptionsFromStreamProcessor&
MoochoSolver::commandLineOptionsFromStreamProcessor()
{
  return commandLineOptionsFromStreamProcessor_;
}

const OptionsFromStreamPack::CommandLineOptionsFromStreamProcessor&
MoochoSolver::commandLineOptionsFromStreamProcessor() const
{
  return commandLineOptionsFromStreamProcessor_;
}

void MoochoSolver::setup_commandline_processor(
  Teuchos::CommandLineProcessor *clp
  )
{
  commandLineOptionsFromStreamProcessor_.setup_commandline_processor(clp);
}

void MoochoSolver::set_output_context(
  const std::string    &file_context_postfix
  ,EOutputToBlackHole  output_to_black_hole
  ,const int           procRank_in
  ,const int           numProcs_in
  )
{

  file_context_postfix_ = file_context_postfix;

  output_to_black_hole_ = output_to_black_hole;
  
  Teuchos::RCP<Teuchos::FancyOStream>
    defaultOut = Teuchos::VerboseObjectBase::getDefaultOStream();
  
  const int procRank
    = ( procRank_in >= 0 ? procRank_in : Teuchos::GlobalMPISession::getRank() );
  const int numProcs
    = ( numProcs_in > 0 ? numProcs_in : Teuchos::GlobalMPISession::getNProc() );

  if(output_to_black_hole_==OUTPUT_TO_BLACK_HOLE_DEFAULT) {
    if( numProcs > 1 && defaultOut->getOutputToRootOnly() >= 0 ) {
      output_to_black_hole_
        = ( defaultOut->getOutputToRootOnly()==procRank
            ? OUTPUT_TO_BLACK_HOLE_FALSE
            : OUTPUT_TO_BLACK_HOLE_TRUE
          );
    }
    else {
      output_to_black_hole_ = OUTPUT_TO_BLACK_HOLE_FALSE;
    }
  }

  if( numProcs > 1 ) {
    file_proc_postfix_ = Teuchos::Utils::getParallelExtension(procRank,numProcs);
  }
  else {
    file_proc_postfix_ = "";
  }
  
  summary_out_ = Teuchos::null;
  journal_out_ = Teuchos::null;
  algo_out_ = Teuchos::null;
  
  summary_out_used_ = Teuchos::null;
  journal_out_used_ = Teuchos::null;
  algo_out_used_ = Teuchos::null;
  stats_out_used_ = Teuchos::null;

}

void MoochoSolver::set_nlp(const nlp_ptr_t& nlp)
{
  nlp_ = nlp;
  reconfig_solver_ = true;
}
  
const MoochoSolver::nlp_ptr_t&
MoochoSolver::get_nlp() const
{
  return nlp_;
}

void MoochoSolver::set_track(const track_ptr_t& track)
{
  track_ = track;
  solver_.set_track(Teuchos::null); // Force the track objects to be rebuilt and added!
}
  
const MoochoSolver::track_ptr_t&
MoochoSolver::get_track() const
{
  return track_;
}
  
void MoochoSolver::set_config( const config_ptr_t& config )
{
  config_ = config;
  solver_.set_config(Teuchos::null); // Must unset the config object.
  reconfig_solver_ = true;
}

const MoochoSolver::config_ptr_t&
MoochoSolver::get_config() const
{
  return config_;
}

void MoochoSolver::set_options( const options_ptr_t& options )
{
  options_ = options;                    // Must totally free all of the references we
  const config_ptr_t                     // have to the current options.  That includes
    &config = solver_.get_config();      // removing the options object for the configuration
  if(config.get())                       // object.
    config->set_options(Teuchos::null);  // ...
  options_used_ = options;
  commandLineOptionsFromStreamProcessor_.set_options(options);
  reconfig_solver_ = true;
}

const MoochoSolver::options_ptr_t&
MoochoSolver::get_options() const
{
  return options_;
}

void MoochoSolver::set_error_handling(
  bool                    throw_exceptions
  ,const ostream_ptr_t&   error_out
  )

{
  if( error_out_.get() != NULL ) {
    if( error_out.get() == NULL )
      error_out_used_ = Teuchos::VerboseObjectBase::getDefaultOStream();
    else 
      error_out_used_ = error_out;
  }
  else if( error_out.get() != NULL ) {
    error_out_used_ = error_out;
  }
  throw_exceptions_ = throw_exceptions;
  error_out_       = error_out;
}

bool MoochoSolver::throw_exceptions() const
{
  return throw_exceptions_;
}

const MoochoSolver::ostream_ptr_t&
MoochoSolver::error_out() const
{
  return error_out_;
}

void MoochoSolver::set_console_out( const ostream_ptr_t& console_out )
{
  console_out_      = console_out;
  console_out_used_ = Teuchos::null;  // Remove every reference to this ostream object!
  solver_.set_track(Teuchos::null);
}

const MoochoSolver::ostream_ptr_t&
MoochoSolver::get_console_out() const
{
  return console_out_;
}

void MoochoSolver::set_summary_out( const ostream_ptr_t& summary_out )
{
  summary_out_      = summary_out;
  summary_out_used_ = Teuchos::null;
  solver_.set_track(Teuchos::null);     // Remove every reference to this ostream object!
}
  
const MoochoSolver::ostream_ptr_t&
MoochoSolver::get_summary_out() const
{
  return summary_out_;
}

void MoochoSolver::set_journal_out( const ostream_ptr_t& journal_out )
{
  journal_out_      = journal_out;
  journal_out_used_ = Teuchos::null;
  solver_.set_track(Teuchos::null);     // Remove every reference to this ostream object!
}
  
const MoochoSolver::ostream_ptr_t&
MoochoSolver::get_journal_out() const
{
  return journal_out_;
}

void MoochoSolver::set_algo_out( const ostream_ptr_t& algo_out )
{
  algo_out_      = algo_out;
  algo_out_used_ = Teuchos::null;
}
  
const MoochoSolver::ostream_ptr_t&
MoochoSolver::get_algo_out() const
{
  return algo_out_;
}

Teuchos::RCP<std::ostream>
MoochoSolver::generate_output_file(const std::string &fileNameBase) const
{
  if( output_to_black_hole_ == OUTPUT_TO_BLACK_HOLE_TRUE )
    return Teuchos::rcp(new Teuchos::oblackholestream());
  std::string fileName = fileNameBase;
  if(file_context_postfix_.length())
    fileName += "." + file_context_postfix_;
  if(file_proc_postfix_.length())
    fileName += "." + file_proc_postfix_;
  if(output_file_tag_.length())
    fileName += ( "." + output_file_tag_ );
  fileName += ".out";
  return Teuchos::rcp(new std::ofstream(fileName.c_str()));
}

// Solve the NLP

void MoochoSolver::update_solver() const
{

  using std::endl;
  using std::setw;
  using StopWatchPack::stopwatch;
  using Teuchos::RCP;
  namespace ofsp = OptionsFromStreamPack;
  using ofsp::OptionsFromStream;
  using ofsp::StringToIntMap;
  using ofsp::StringToBool;

  /** \brief . */
  // Validate the input
  //
  
  TEUCHOS_TEST_FOR_EXCEPTION(
    nlp_.get() == NULL, std::logic_error
    ,"MoochoSolver::update_solver() : Error, this->get_nlp().get() can not be NULL!" );

  
  //
  // Get the options (or lack of)
  //

  bool MoochoSolver_opt_grp_existed = true;
    
  if(reconfig_solver_) {
    
    if( options_used_.get() == NULL ) {
      if( options_.get() == NULL ) {
        options_used_ = commandLineOptionsFromStreamProcessor_.process_and_get_options();
      }
      else
        options_used_ = options_;
    }
    
    //
    // Read in some options for "MoochoSolver" if needed
    //
    
    if( options_used_.get() ) {
      
      options_used_->reset_unaccessed_options_groups();
      
      const std::string optgrp_name = "MoochoSolver";
      OptionsFromStream::options_group_t optgrp = options_used_->options_group( optgrp_name );
      if( OptionsFromStream::options_group_exists( optgrp ) ) {
        
        const int num_opt = 12;
        enum EOptions {
          WORKSPACE_MB
          ,OBJ_SCALE
          ,TEST_NLP
          ,CONSOLE_OUTPUTTING
          ,SUMMARY_OUTPUTTING
          ,JOURNAL_OUTPUTTING
          ,ALGO_OUTPUTTING
          ,PRINT_ALGO
          ,ALGO_TIMING
          ,GENERATE_STATS_FILE
          ,PRINT_OPT_GRP_NOT_ACCESSED
          ,CONFIGURATION
        };
        const char* SOptions[num_opt] = {
          "workspace_MB"
          ,"obj_scale"
          ,"test_nlp"
          ,"console_outputting"
          ,"summary_outputting"
          ,"journal_outputting"
          ,"algo_outputting"
          ,"print_algo"
          ,"algo_timing"
          ,"generate_stats_file"
          ,"print_opt_grp_not_accessed"
          ,"configuration"
        };
        
        const int num_config_opt = 2;
        
        const char* SConfigOptions[num_config_opt] = {
          "mama_jama"
          ,"interior_point"
        };
        
        StringToIntMap	config_map( optgrp_name, num_config_opt, SConfigOptions );
        
        StringToIntMap	opt_map( optgrp_name, num_opt, SOptions );
        
        OptionsFromStream::options_group_t::const_iterator itr = optgrp.begin();
        for( ; itr != optgrp.end(); ++itr ) {
          switch( (EOptions)opt_map( ofsp::option_name(itr) ) ) {
            case WORKSPACE_MB:
              workspace_MB_ = std::atof( ofsp::option_value(itr).c_str() );
              break;
            case OBJ_SCALE:
              obj_scale_ = std::atof( ofsp::option_value(itr).c_str() );
              break;
            case TEST_NLP:
              test_nlp_ = StringToBool( "test_nlp", ofsp::option_value(itr).c_str() );
              break;
            case CONSOLE_OUTPUTTING:
              do_console_outputting_ = StringToBool( "console_outputting", ofsp::option_value(itr).c_str() );
              break;
            case SUMMARY_OUTPUTTING:
              do_summary_outputting_ = StringToBool( "summary_outputting", ofsp::option_value(itr).c_str() );
              break;
            case JOURNAL_OUTPUTTING:
              do_journal_outputting_ = StringToBool( "journal_outputting", ofsp::option_value(itr).c_str() );
              break;
            case ALGO_OUTPUTTING:
              do_algo_outputting_ = StringToBool( "algo_outputting", ofsp::option_value(itr).c_str() );
              break;
            case PRINT_ALGO:
              print_algo_ = StringToBool( "print_algo", ofsp::option_value(itr).c_str() );
              break;
            case ALGO_TIMING:
              algo_timing_ = StringToBool( "algo_timing", ofsp::option_value(itr).c_str() );
              break;
            case GENERATE_STATS_FILE:
              generate_stats_file_ = StringToBool( "generate_stats_file", ofsp::option_value(itr).c_str() );
              break;
            case PRINT_OPT_GRP_NOT_ACCESSED:
              print_opt_grp_not_accessed_ = StringToBool( "algo_timing", ofsp::option_value(itr).c_str() );
              break;
            case CONFIGURATION:
              configuration_ = config_map( ofsp::option_value(itr).c_str() );
              break;
            default:
              TEUCHOS_TEST_FOR_EXCEPT(true);	// this would be a local programming error only.
          }
        }
      }
      else {
        MoochoSolver_opt_grp_existed = false;
      }

    }
    
  }
  
  //
  // Get the output streams if needed
  //

  generate_output_streams();  
  
  if( do_algo_outputting() && !MoochoSolver_opt_grp_existed )
    *algo_out_used_
      << "\nWarning!  The options group \'MoochoSolver\' was not found.\n"
      "Using a default set of options ...\n";

  //
  // Configure the algorithm
  //
    
  if(reconfig_solver_) {
    
    //
    // Print the headers for the output files
    //
      
    int prec = 8;
    if(do_summary_outputting())
      summary_out_used_->precision(prec);
    if(do_journal_outputting())
      *journal_out_used_ << std::setprecision(prec) << std::scientific;

    if(do_algo_outputting())
      *algo_out_used_
        << "\n********************************************************************"
        << "\n*** Algorithm information output                                 ***"
        << "\n***                                                              ***"
        << "\n*** Below, information about how the the MOOCHO algorithm is     ***"
        << "\n*** setup is given and is followed by detailed printouts of the  ***"
        << "\n*** contents of the algorithm state object (i.e. iteration       ***"
        << "\n*** quantities) and the algorithm description printout           ***"
        << "\n*** (if the option MoochoSolver::print_algo = true is set).      ***"
        << "\n********************************************************************\n";
    if(do_summary_outputting())
      *summary_out_used_
        << "\n********************************************************************"
        << "\n*** Algorithm iteration summary output                           ***"
        << "\n***                                                              ***"
        << "\n*** Below, a summary table of the SQP iterations is given as     ***"
        << "\n*** well as a table of the CPU times for each step (if the       ***"
        << "\n*** option MoochoSolver::algo_timing = true is set).             ***"
        << "\n********************************************************************\n";
    if(do_journal_outputting())
      *journal_out_used_
        << "\n********************************************************************"
        << "\n*** Algorithm iteration detailed journal output                  ***"
        << "\n***                                                              ***"
        << "\n*** Below, detailed information about the SQP algorithm is given ***"
        << "\n*** while it is running.  The amount of information that is      ***"
        << "\n*** produced can be specified using the option                   ***"
        << "\n*** NLPSolverClientInterface::journal_output_level (the default ***"
        << "\n*** is PRINT_NOTHING and produces no output                      ***"
        << "\n********************************************************************\n";
    
    // Echo options.
    if(do_summary_outputting()) {
      *summary_out_used_ << "\n*** Echoing input options ...\n";
      if(options_used_.get())
        options_used_->print_options( *summary_out_used_ );
      summary_out_used_->flush();
    }
    if(do_algo_outputting()) {
      *algo_out_used_ << "\n*** Echoing input options ...\n";
      if(options_used_.get())
        options_used_->print_options( *algo_out_used_ );
      algo_out_used_->flush();
    }
    if(do_journal_outputting()) {
      *journal_out_used_ << "\n*** Echoing input options ...\n";
      if(options_used_.get())
        options_used_->print_options( *journal_out_used_ );
      journal_out_used_->flush();
    }
      
    //
    // Allocate the workspace
    //
      
    nlp_->set_options(options_used_);
    nlp_->initialize();
    const int default_ws_scale = 10;
    if( workspace_MB_ < 0.0 ) {
      workspace_MB_ = nlp_->n() * default_ws_scale * 1e-6 * sizeof(value_type);
      if(do_algo_outputting())
        *algo_out_used_
          << "\nworkspace_MB < 0.0:\n"
          << "Setting workspace_MB = n * default_ws_scale * 1e-6 * sizeof(value_type) = "
          << nlp_->n() << " * " << default_ws_scale << " * 1e-6 * " << sizeof(value_type)
          << " = " << workspace_MB_ << " MB\n";
    }
    if(do_summary_outputting())
      *summary_out_used_
        << "\nAllocating workspace_MB = " << workspace_MB_ << " megabytes of temporary "
        "workspace for automatic arrays only ...\n";
    Teuchos::set_default_workspace_store(
      Teuchos::rcp(new Teuchos::WorkspaceStoreInitializeable(static_cast<size_t>(1e+6*workspace_MB_)))
      );
    
    //
    // Reconfigure the algorithm
    //
      
    // Get and set up the configuration object
    config_ptr_t _config;
    if(config_.get() == NULL) {
      if (configuration_ == (EConfigOptions) INTERIOR_POINT) {
          _config = Teuchos::rcp(new NLPAlgoConfigIP());
            }
      else {
          _config = Teuchos::rcp(new NLPAlgoConfigMamaJama());
      }
        }
    else
      _config = config_;
    _config->set_options( options_used_ );
      
    if( do_summary_outputting() || do_journal_outputting() || do_algo_outputting() ) {
      std::ostringstream msg;
      msg << "\n*** Setting up to run MOOCHO on the NLP using a "
        << "configuration object of type \'" << typeName(*_config) << "\' ...\n";
      if(do_summary_outputting())
        *summary_out_used_ << msg.str();
      if(do_journal_outputting())
        *journal_out_used_ << msg.str();
      if(do_algo_outputting())
        *algo_out_used_ << msg.str();
    }
      
    // Set up the solver
      
    solver_.set_nlp(nlp_);           // Set the NLP object
    solver_.set_config(_config);     // Set the configuration object
      
    // Set the client interface options
    if(options_used_.get()) {
      NLPSolverClientInterfaceSetOptions
        solver_options_setter( &solver_ );
      solver_options_setter.set_options( *options_used_ );
    }
    
    reconfig_solver_ = false;

    //
    // Set up the track objects
    //
    
    RCP<AlgorithmTrackerComposite>
      composite_track = Teuchos::rcp(new AlgorithmTrackerComposite(journal_out_used_));
    if(do_console_outputting())
      composite_track->tracks().push_back(
        Teuchos::rcp(new MoochoTrackerConsoleStd(console_out_used_,journal_out_used_)) );
    if(do_summary_outputting())
      composite_track->tracks().push_back(
        Teuchos::rcp(new MoochoTrackerSummaryStd(summary_out_used_,journal_out_used_)) );
    if(stats_out_used_.get()) {
      composite_track->tracks().push_back(
        Teuchos::rcp(new MoochoTrackerStatsStd(stats_out_used_,stats_out_used_)) );
    }
    if( track_.get() ) {
      track_->set_journal_out(journal_out_used_);
      composite_track->tracks().push_back( track_ );
    }
    solver_.set_track( composite_track );
    
    //
    // Configure and print the algorithm if needed
    //
    
    solver_.configure_algorithm(algo_out_used_.get());
    if( do_algo_outputting() && print_algo_ )
      solver_.print_algorithm(*algo_out_used_);

  }
  
}

MoochoSolver::ESolutionStatus MoochoSolver::solve_nlp() const
{
  using std::endl;
  using std::setw;
  using std::flush;
  using StopWatchPack::stopwatch;
  using Teuchos::RCP;

  stopwatch timer;
  bool threw_exception = false;
  ESolutionStatus solve_return = SOLVE_RETURN_EXCEPTION;
  NLPSolverClientInterface::EFindMinReturn
    r_find_min = NLPSolverClientInterface::SOLUTION_FOUND;
  
  try {
    
    update_solver();

    {
      std::ostringstream os;
      os
        << "\n*****************************"
        << "\n*** MoochoSolver::solve() ***"
        << "\n*****************************\n";
      *this->get_console_out() << os.str();
      *this->get_summary_out() << os.str();
      *this->get_journal_out() << os.str();
    }
    //    
    // Direct any output from the NLP to the journal output file
    //

    EJournalOutputLevel olevel = solver_.journal_output_level();
    Teuchos::VerboseObjectTempState<NLP>
      nlpOutputTempState(nlp_,Teuchos::getFancyOStream(journal_out_used_),convertToVerbLevel(olevel));

    //
    // Scale the NLP objective function
    //
  
    nlp_->scale_f(obj_scale_);

    //
    // Test the nlp if needed
    //
    
    if(test_nlp_) {
      
      const char msg1[] = "\ntest_nlp = true: Testing the NLP! ...\n";
      if(do_console_outputting())
        *console_out_used_ << msg1;
      if(do_summary_outputting())
        *summary_out_used_ << msg1;
      if(do_journal_outputting())
        *journal_out_used_ << msg1;
      if(NLPFirstOrder* nlp_foi = dynamic_cast<NLPFirstOrder*>(nlp_.get())) {
        const char msg[] = "\nTesting the supported NLPFirstOrder interface ...\n";
        if(do_console_outputting())
          *console_out_used_ << msg;
        if(do_summary_outputting())
          *summary_out_used_ << msg;
        if(do_journal_outputting())
          *journal_out_used_ << msg;
        const bool
          result = NLPInterfacePack::test_nlp_first_order(
            nlp_foi,options_used_.get()
            ,do_journal_outputting() ? journal_out_used_.get() : NULL
            );
        if(!result) {
          const char msg[] = "\nNLPFirstOrder test failed (see journal file)!  exiting!\n";
          if(do_console_outputting())
            *console_out_used_ << msg;
          if(do_summary_outputting())
            *summary_out_used_ << msg;
          if(do_journal_outputting())
            *journal_out_used_ << msg;
          solve_return = SOLVE_RETURN_NLP_TEST_FAILED;
          return solve_return;
        }
      }
      else if(NLPDirect* nlp_fod = dynamic_cast<NLPDirect*>(nlp_.get())) {
        const char msg[] = "\nTesting the supported NLPDirect interface ...\n";
        if(do_console_outputting())
          *console_out_used_ << msg;
        if(do_summary_outputting())
          *summary_out_used_ << msg;
        if(do_journal_outputting())
          *journal_out_used_ << msg;
        const bool
          result = NLPInterfacePack::test_nlp_direct(
            nlp_fod,options_used_.get()
            ,do_journal_outputting() ? journal_out_used_.get() : NULL
            );
        if(!result) {
          const char msg[] = "\nNLPDirect test failed (see journal file)!  exiting!\n";
          if(do_console_outputting())
            *console_out_used_ << msg;
          if(do_summary_outputting())
            *summary_out_used_ << msg;
          if(do_journal_outputting())
            *journal_out_used_ << msg;
          solve_return = SOLVE_RETURN_NLP_TEST_FAILED;
          return solve_return;
        }
      }
      const char msg2[] = "\nSuccessful end of testing of the nlp\n";
      if(do_console_outputting())
        *console_out_used_ << msg2;
      if(do_summary_outputting())
        *summary_out_used_ << msg2;
      if(do_journal_outputting())
        *journal_out_used_ << msg2;

    }
    
    //
    // Solve the NLP
    //
    
    if(do_journal_outputting())
      *journal_out_used_
        << "\n************************************"
        << "\n*** MoochoSolver::solve_nlp()    ***"
        << "\n************************************\n"	
        << "\n*** Starting iterations ...\n\n";
    
    solver_.set_algo_timing(algo_timing_);
    timer.start();
    r_find_min = solver_.find_min();
    
    }
  catch(const std::exception& excpt) {
    std::ostringstream msg;
    msg << "\nMoochoSolver: Caught an std::exception of type "
        << typeName(excpt) << " described as : " << excpt.what() << endl;
    *error_out_used_  << msg.str() << flush;
    if(do_summary_outputting())
      *summary_out_used_ << msg.str() << flush;
    if(do_journal_outputting())
      *journal_out_used_ << msg.str() << flush;
    if(throw_exceptions_)
      throw;
    threw_exception = true;
  }
  catch(...) {
    std::ostringstream msg;
    msg << "\nMoochoSolver: Caught an unknown exception (i.e. ...)\n";
    *error_out_used_   << msg.str() << flush;
    if(do_summary_outputting())
      *summary_out_used_ << msg.str() << flush;
    if(do_journal_outputting())
      *journal_out_used_ << msg.str() << flush;
    if(throw_exceptions_)
      throw;
    threw_exception = true;
  }
  
  timer.stop();
  
  if(threw_exception) {
    if(do_summary_outputting())
      *summary_out_used_	<< "\n\n****************************\n"
                << "**** Threw an exception ****\n";
    solve_return = SOLVE_RETURN_EXCEPTION;
  }
  else {
    switch( r_find_min ) {
        case NLPSolverClientInterface::SOLUTION_FOUND: {
        if(do_summary_outputting())
          *summary_out_used_	<< "\n\n************************\n"
                    << "**** Solution Found ****\n";
        *error_out_used_    << "Solution Found!\n";
        solve_return = SOLVE_RETURN_SOLVED;
        break;
      }
        case NLPSolverClientInterface::MAX_ITER_EXCEEDED: {
        if(do_summary_outputting())
          *summary_out_used_	<< "\n\n**********************************************\n"
                    << "**** Maximun number of iteration exceeded ****\n";
        *error_out_used_    << "Maximun number of iteration exceeded!\n";
        solve_return = SOLVE_RETURN_MAX_ITER;
        break;
      }
        case NLPSolverClientInterface::MAX_RUN_TIME_EXCEEDED: {
        if(do_summary_outputting())
          *summary_out_used_	<< "\n\n**********************************\n"
                    << "**** Maximun runtime exceeded ****\n";
        *error_out_used_    << "Maximun runtime exceeded!\n";
        solve_return = SOLVE_RETURN_MAX_RUN_TIME;
        break;
      }
        case NLPSolverClientInterface::ALGORITHMIC_ERROR: {
        if(do_summary_outputting())
          *summary_out_used_	<< "\n\n*********************************************\n"
                    << "**** Some error occurred in the algorithm ****\n";
        *error_out_used_    << "Some algorithmic error occurred!\n";
        solve_return = SOLVE_RETURN_EXCEPTION;
        break;
      }
    }
  }
  
  if(do_summary_outputting()) {
    *summary_out_used_	<< "\n  total time = " << timer.read() << " sec.\n";
    if( solver_.algo_timing() ) {
      Teuchos::TimeMonitor::format().setRowsBetweenLines(100);
      solver_.print_algorithm_times( *summary_out_used_ );
      *summary_out_used_ << "\n\n";
      Teuchos::TimeMonitor::summarize(*summary_out_used_);
    }
  }
  
  // Print workspace usage statistics
  if(do_summary_outputting())
    *summary_out_used_
      << "\n*** Statistics for automatic array workspace:"
      << "\nNumber of megabytes of preallocated workspace                = "
      << workspace_MB_
      << "\nNumber of allocations using preallocated workspace           = "
      << Teuchos::get_default_workspace_store()->num_static_allocations()
      << "\nNumber of dynamic allocations beyond preallocated workspace  = "
      << Teuchos::get_default_workspace_store()->num_dyn_allocations();
  
  // Print which options groups were not read
  if( do_algo_outputting() && print_opt_grp_not_accessed_ ) {
    *algo_out_used_
      <<	"\n***************************************************************\n"
      "Warning, the following options groups where not accessed.\n"
      "An options group may not be accessed if it is not looked for\n"
      "or if an \"optional\" options group was looked from and the user\n"
      "spelled it incorrectly:\n\n";
    if(options_used_.get())
      options_used_->print_unaccessed_options_groups(*algo_out_used_);
  }

  if(do_console_outputting())
    console_out_used_->flush();
  if(do_summary_outputting())
    summary_out_used_->flush();
  if(do_journal_outputting())
    journal_out_used_->flush();
  if(do_algo_outputting())
    algo_out_used_->flush();
  
  return solve_return;
  
}

// Get the underlying solver object

NLPSolverClientInterface& MoochoSolver::get_solver()
{
  update_solver();
  return solver_;
}

const NLPSolverClientInterface& MoochoSolver::get_solver() const
{
  update_solver();
  return solver_;
}

// private

void MoochoSolver::generate_output_streams() const
{
  if( do_console_outputting() && console_out_used_.get() == NULL ) {
    if( console_out_.get() == NULL ) {
      console_out_used_ = Teuchos::VerboseObjectBase::getDefaultOStream();
      console_out_ = console_out_used_;
    }
    else {
      console_out_used_ = console_out_;
    }
  }
  if( do_summary_outputting() && summary_out_used_.get()==NULL ) {
    if( summary_out_.get() == NULL ) {
      summary_out_used_ = generate_output_file("MoochoSummary");
      summary_out_ = summary_out_used_;
    }
    else {
      summary_out_used_ = summary_out_;
    }
  }
  if( do_journal_outputting() && journal_out_used_.get() == NULL ) {
    if( journal_out_.get() == NULL ) {
      journal_out_used_ = generate_output_file("MoochoJournal");
      journal_out_ = journal_out_used_;
    }
    else {
      journal_out_used_ = journal_out_;
    }
  }
  else {
    journal_out_used_ = Teuchos::rcp(new Teuchos::oblackholestream());
  }
  if( do_algo_outputting() && algo_out_used_.get() == NULL ) {
    if( algo_out_.get() == NULL ) {
      algo_out_used_ = generate_output_file("MoochoAlgo");
      algo_out_ = algo_out_used_;
    }
    else {
      algo_out_used_ = algo_out_;
    }
  }
  if( generate_stats_file_ && stats_out_used_.get() == NULL ) {
    stats_out_used_ = generate_output_file("MoochoStats");
  }
}

} // end namespace MoochoPack
