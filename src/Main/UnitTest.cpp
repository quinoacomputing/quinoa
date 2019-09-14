// *****************************************************************************
/*!
  \file      src/Main/UnitTest.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     UnitTest's Charm++ main chare and main().
  \details   UnitTest's Charm++ main chare and main(). This file contains
    the definition of the Charm++ main chare, equivalent to main() in Charm++-
    land, running the serial and Charm++ unit tests as well as the ordinary
    main() function, running the MPI unit test suite.
*/
// *****************************************************************************

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <utility>
#include <cstddef>

#include "NoWarning/tut_runner.hpp"

#include "NoWarning/tutsuite.decl.h"
#include "NoWarning/unittest.decl.h"

#include "Print.hpp"
#include "Timer.hpp"
#include "Tags.hpp"
#include "Exception.hpp"
#include "Init.hpp"
#include "QuinoaConfig.hpp"
#include "HelpFactory.hpp"
#include "Assessment.hpp"
#include "ProcessException.hpp"
#include "UnitTest/CmdLine/CmdLine.hpp"
#include "UnitTestPrint.hpp"
#include "UnitTestDriver.hpp"
#include "UnitTest/CmdLine/Parser.hpp"
#include "TUTConfig.hpp"
#include "ChareStateCollector.hpp"
#include "QuietCerr.hpp"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wmissing-variable-declarations"
#endif

//! \brief Charm handle to the main proxy, facilitates call-back to finalize,
//!    etc., must be in global scope, unique per executable
CProxy_Main mainProxy;

//! Chare state collector Charm++ chare group proxy
tk::CProxy_ChareStateCollector stateProxy;

//! If true, call and stack traces are to be output with exceptions
//! \note This is true by default so that the trace is always output between
//!   program start and the Main ctor in which the user-input from command line
//!   setting for this overrides this true setting.
bool g_trace = true;

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

//! UnitTest declarations and definitions
namespace unittest {

//! Global-scope data. Initialized by the main chare and distibuted to all PEs
//! by the Charm++ runtime system. Though semantically not const, all these
//! global data should be considered read-only. See also
//! http://charm.cs.illinois.edu/manuals/html/charm++/manual.html. The data
//! below is global-scope because they must be available to all PEs which could
//! be on different machines.

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wmissing-variable-declarations"
#endif

//! Template Unit Test test runner
tut::test_runner_singleton g_runner;

//! Test suite Charm++ proxy facilitating call-back to unit test suite by
//! individual unit tests spawning Charm++ chares
CProxy_TUTSuite g_suiteProxy;

//! UnitTest executable name. So that FileParser's unit tests can access a file
//! for opening.
std::string g_executable;

//! Max number of tests in every group
int g_maxTestsInGroup = tut::MAX_TESTS_IN_GROUP;

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

//! Pack/Unpack test runner. This Pack/Unpack method (re-)creates the
//! test runner singleton on all processing elements. Therefore we circumvent
//! Charm's usual pack/unpack for this type, and thus sizing does not make
//! sense: sizing is a no-op. We could initialize the stack in UnitTestDriver's
//! constructor and let this function re-create the runner only when unpacking,
//! but that leads to repeating the same code twice: once in UnitTestDriver's
//! constructor, once here. Another option is to use this pack/unpack routine to
//! both initially create (when packing) and to re-create (when unpacking) the
//! runner, which eliminates the need for pre-creating the object in
//! UnitTestDriver's constructor and therefore eliminates the repeated code.
//! This explains the guard for sizing: the code below is called for packing
//! only (in serial) and packing and unpacking (in parallel).
inline void operator|( PUP::er& p, tut::test_runner_singleton& runner )
{ if (!p.isSizing()) runner = tut::test_runner_singleton(); }

} // unittest::

//! \brief Charm++ main chare for the unit test suite executable, unittest.
//! \details Note that this object should not be in a namespace.
// cppcheck-suppress noConstructor
class Main : public CBase_Main {

  public:
    //! \brief Constructor
    //! \details UnitTest's main chare constructor is the entry point of the
    //!   program, called by the Charm++ runtime system. The constructor does
    //!   basic initialization steps, e.g., parser the command-line, prints out
    //!   some useful information to screen (in verbose mode), and instantiates
    //!   a driver. Since Charm++ is fully asynchronous, the constructor
    //!   usually spawns asynchronous objects and immediately exits. Thus in the
    //!   body of the main chare constructor we fire up an 'execute' chare,
    //!   which then calls back to Main::execute(). Finishing the main chare
    //!   constructor the Charm++ runtime system then starts the
    //!   network-migration of all global-scope data (if any). The execute chare
    //!   calling back to Main::execute() signals the end of the migration of
    //!   the global-scope data. Then we are ready to execute the driver. Since
    //!   the unit test suite is parallel and asynchronous, its driver fires up
    //!   additional Charm++ chare objects which then call back to
    //!   Main::finalize() at some point in the future when all work has been
    //!   finished. finalize() then exits by calling Charm++'s CkExit(),
    //!   shutting down the runtime system.
    //! \see http://charm.cs.illinois.edu/manuals/html/charm++/manual.html
    Main( CkArgMsg* msg )
    try :
      m_signal( tk::setSignalHandlers() ),
      m_helped( false ),
      m_cmdline(),
      // Parse command line into m_cmdline using default simple pretty printer
      m_cmdParser( msg->argc, msg->argv, tk::Print(), m_cmdline, m_helped ),
      // Create pretty printer initializing output streams based on command line
      m_print( m_cmdline.get< tag::verbose >() ? std::cout : std::clog ),
      // Create UnitTest driver
      m_driver( tk::Main< unittest::UnitTestDriver >
                        ( msg->argc, msg->argv,
                          m_cmdline,
                          tk::HeaderType::UNITTEST,
                          tk::unittest_executable(),
                          m_print ) ),
      m_timer(1), // Start new timer measuring the serial+Charm++ runtime
      m_timestamp()
    {
      g_trace = m_cmdline.get< tag::trace >();
      // Immediately exit if any help was requested; help is printed in main()
      if (m_helped) CkExit();
      // Save executable name to global-scope string so FileParser can access it
      unittest::g_executable = msg->argv[0];
      delete msg;
      // Call generic mainchare contructor
      tk::MainCtor( mainProxy, thisProxy, m_timer, m_cmdline,
                    CkCallback( CkIndex_Main::quiescence(), thisProxy ) );
      // If quiescence detection is on or user requested it, create chare state
      // collector Charm++ chare group
      if ( m_cmdline.get< tag::chare >() || m_cmdline.get< tag::quiescence >() )
        stateProxy = tk::CProxy_ChareStateCollector::ckNew();
      // Fire up an asynchronous execute object, which when created at some
      // future point in time will call back to this->execute(). This is
      // necessary so that this->execute() can access already migrated
      // global-scope data.
      CProxy_execute::ckNew();
      // Quiet std::cerr
      tk::CProxy_QuietCerr::ckNew();
    } catch (...) { tk::processExceptionCharm(); }

    void execute() {
      try {
        m_timestamp.emplace_back("Migrate global-scope data", m_timer[1].hms());
        m_driver.execute();       // fires up async chares
      } catch (...) { tk::processExceptionCharm(); }
    }

    //! Towards normal exit but collect chare state first (if any)
    void finalize( bool pass ) {
      tk::finalize( m_cmdline, m_timer, m_print, stateProxy, m_timestamp,
                    CkCallback( CkIndex_Main::dumpstate(nullptr), thisProxy ),
                    pass );
    }

    //! Entry method triggered when quiescence is detected
    void quiescence() {
      try {
        stateProxy.collect( /* error= */ true,
          CkCallback( CkIndex_Main::dumpstate(nullptr), thisProxy ) );
      } catch (...) { tk::processExceptionCharm(); }
    }

    //! Dump chare state
    void dumpstate( CkReductionMsg* msg ) {
      tk::dumpstate( m_cmdline, m_print, msg );
    }

  private:
    int m_signal;                               //!< Used to set signal handlers
    bool m_helped;      //!< Indicates if help was requested on the command line
    unittest::ctr::CmdLine m_cmdline;                   //!< Command line
    unittest::CmdLineParser m_cmdParser;                //!< Command line parser
    unittest::UnitTestPrint m_print;                    //!< Pretty printer
    unittest::UnitTestDriver m_driver;                  //!< Driver
    std::vector< tk::Timer > m_timer;                   //!< Timers

    //! Time stamps in h:m:s with labels
    std::vector< std::pair< std::string, tk::Timer::Watch > > m_timestamp;
};

//! \brief Charm++ chare execute
//! \details By the time this object is constructed, the Charm++ runtime system
//!    has finished migrating all global-scoped read-only objects which happens
//!    after the main chare constructor has finished.
class execute : public CBase_execute {
  public: execute() { mainProxy.execute(); }
};

#include "NoWarning/unittest.def.h"
