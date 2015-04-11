//******************************************************************************
/*!
  \file      src/Main/RegTest.C
  \author    J. Bakosi
  \date      Sat 11 Apr 2015 07:50:33 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     RegTest's Charm++ main chare and main().
  \details   RegTest's Charm++ main chare and main(). This file contains
    the definition of the Charm++ main chare, equivalent to main() in Charm++-
    land, running the regression test suite.
*/
//******************************************************************************
#include <Config.h>
#include <RegTestPrint.h>
#include <RegTestDriver.h>
#include <RegTest/CmdLine/Parser.h>
#include <ProcessException.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <regtest.decl.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#include <Init.h>

//! \brief Charm handle to the main proxy, facilitates call-back to finalize,
//!    etc., must be in global scope, unique per executable
CProxy_Main mainProxy;

//! RegTest declarations and definitions
namespace regtest {

//! Global-scope data. Initialized by the main chare and distibuted to all PEs
//! by the Charm++ runtime system. Though semantically not const, all these
//! global data should be considered read-only. See also
//! http://charm.cs.illinois.edu/manuals/html/charm++/manual.html. The data
//! below is global-scope because they must be available to all PEs which could
//! be on different machines.

// ...

} // regtest::

//! \brief Charm++ main chare for the regression test suite executable, regtest.
//! \details Note that this object should not be in a namespace.
class Main : public CBase_Main {

  public:
    //! \brief Constructor
    //! \details The main chare constructor is the main entry point of the
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
    //!   the regression test suite is parallel and asynchronous, its driver
    //!   fires up additional Charm++ chare objects which then call back to
    //!   Main::finalize() at some point in the future when all work has been
    //!   finished. finalize() then exits by calling Charm++'s CkExit(),
    //!   shutting down the runtime system.
    //! \see http://charm.cs.illinois.edu/manuals/html/charm++/manual.html
    Main( CkArgMsg* msg )
    try :
      // Parse command line into m_cmdline using default simple pretty printer
      m_cmdParser( msg->argc, msg->argv, tk::Print(), m_cmdline ),
      // Create pretty printer initializing output streams based on command line
      m_print( m_cmdline.get< tag::verbose >() ? std::cout : std::clog ),
      // Create RegTest driver
      m_driver( tk::Main< regtest::RegTestDriver >
                        ( msg->argc, msg->argv,
                          m_cmdline,
                          tk::HeaderType::REGTEST,
                          REGTEST_EXECUTABLE,
                          m_print ) ),
      m_timer(1)  // Start new timer measuring the serial+Charm++ runtime
    {
      delete msg;
      mainProxy = thisProxy;
      // Fire up an asynchronous execute object, which when created at some
      // future point in time will call back to this->execute(). This is
      // necessary so that this->execute() can access already migrated
      // global-scope data.
      CProxy_execute::ckNew();
      // Start new timer measuring the migration of global-scope data
      m_timer.emplace_back();
    } catch (...) { tk::processExceptionCharm(); }

    void execute() {
      try {
        m_timestamp.emplace_back(
          "Migration of global-scope data + fire up all tests",
           m_timer[1].hms() );
        m_driver.execute();       // fires up async chares
      } catch (...) { tk::processExceptionCharm(); }
    }

    void finalize( bool worked ) {
      try {
        if (worked && !m_timer.empty()) {
          m_timestamp.emplace_back( "Total runtime", m_timer[0].hms() );
          m_print.time( "Timers (h:m:s)", m_timestamp );
          m_print.endpart();
        }
      } catch (...) { tk::processExceptionCharm(); }
      // Tell the Charm++ runtime system to exit
      CkExit();
    }

  private:
    regtest::ctr::CmdLine m_cmdline;                    //!< Command line
    regtest::CmdLineParser m_cmdParser;                 //!< Command line parser
    regtest::RegTestPrint m_print;                      //!< Pretty printer
    regtest::RegTestDriver m_driver;                    //!< Driver
    std::vector< tk::Timer > m_timer;                   //!< Timers

    //! Time stamps in h:m:s with labels
    std::vector< std::pair< std::string, tk::Timer::Watch > > m_timestamp;
};

//! \brief Charm++ chare execute
//! \details By the time this object is constructed, the Charm++ runtime system
//!    has finished migrating all global-scoped read-only objects which happens
//!    after the main chare constructor has finished.
//! \author J. Bakosi
struct execute : CBase_execute { execute() { mainProxy.execute(); } };

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <regtest.def.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif
