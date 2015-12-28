//******************************************************************************
/*!
  \file      src/Main/Inciter.C
  \author    J. Bakosi
  \date      Fri 11 Dec 2015 12:40:25 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Inciter, computational shock hydrodynamics tool, Charm++ main
    chare.
  \details   Inciter, computational shock hydrodynamics tool, Charm++ main
    chare. This file contains the definition of the Charm++ main chare,
    equivalent to main() in Charm++-land.
*/
//******************************************************************************

#include <unordered_map>
#include <vector>
#include <iostream>

#include "Types.h"
#include "Init.h"
#include "Config.h"
#include "Timer.h"
#include "Exception.h"
#include "ProcessException.h"
#include "InciterPrint.h"
#include "InciterDriver.h"
#include "Inciter/CmdLine/Parser.h"
#include "Inciter/CmdLine/CmdLine.h"
#include "Inciter/InputDeck/InputDeck.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "inciter.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

//! \brief Charm handle to the main proxy, facilitates call-back to finalize,
//!    etc., must be in global scope, unique per executable
CProxy_Main mainProxy;

//! Inciter declarations and definitions
namespace inciter {

//! Global-scope data. Initialized by the main chare and distibuted to all PEs
//! by the Charm++ runtime system. Though semantically not const, all these
//! global data should be considered read-only. See also
//! http://charm.cs.illinois.edu/manuals/html/charm++/manual.html. The data
//! below is global-scope because they must be available to all PEs which could
//! be on different machines.

//! Defaults of input deck, facilitates detection what is set by user
ctr::InputDeck g_inputdeck_defaults;
//! Input deck filled by parser, containing all input data
ctr::InputDeck g_inputdeck;

} // inciter::

//! \brief Charm++ main chare for the shock hydroddynamics executable, inciter.
//! \details In inciter the Charm++ runtime system is initialized only after the
//!   mesh has been read in, partitioned, and the necessary data structures,
//!   e.g., communication maps, have been generated. This delayed initialization
//!   of the Charm++ runtime system is required since the mesh partitioning is
//!   done by Zoltan, an MPI library. Note that this Charm++ main chare object
//!   should not be in a namespace.
//! \author J. Bakosi
class Main : public CBase_Main {

  public:
    //! \brief Constructor
    //! \details The main chare constructor is the main entry point of the
    //!   Charm++ portion of inciter, called by the Charm++ runtime system. The
    //!   constructor does basic initialization steps, prints out some useful
    //!   information to screen (in verbose mode), and instantiates a driver.
    //!   Since Charm++ is fully asynchronous, the constructor usually spawns
    //!   asynchronous objects and immediately exits. Thus in the body of the
    //!   main chare constructor we fire up an 'execute' chare, which then calls
    //!   back to Main::execute(). Finishing the main chare constructor the
    //!   Charm++ runtime system then starts the network-migration of all
    //!   global-scope data (if any). The execute chare calling back to
    //!   Main::execute() signals the end of the migration of the global-scope
    //!   data. Then we are ready to execute the driver. Since inciter is
    //!   parallel and asynchronous, its driver fires up additional Charm++
    //!   chare objects which then call back to Main::finalize() at some point
    //!   in the future when all work has been finished. finalize() then exits
    //!   by calling Charm++'s CkExit(), shutting down the runtime system.
    //! \see http://charm.cs.illinois.edu/manuals/html/charm++/manual.html
    Main( CkArgMsg* msg )
    try :
      // Parse command line into m_cmdline using default simple pretty printer
      m_cmdParser( msg->argc, msg->argv, tk::Print(), m_cmdline ),
      // Create pretty printer initializing output streams based on command line
      m_print( m_cmdline.get< tag::verbose >() ? std::cout : std::clog ),
      // Create Inciter driver
      m_driver( tk::Main< inciter::InciterDriver >
                        ( msg->argc, msg->argv,
                          m_cmdline,
                          tk::HeaderType::INCITER,
                          INCITER_EXECUTABLE,
                          m_print ) ),
      // Start new timer measuring the total runtime
      m_timer(1)
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

    //! Execute driver created and initialized by constructor
    void execute() {
      try {
        m_timestamp.emplace_back("Migrate global-scope data", m_timer[1].hms());
        m_driver.execute();
      } catch (...) { tk::processExceptionCharm(); }
    }

    //! Normal exit point
    void finalize() {
      try {
        if (!m_timer.empty()) {
          m_timestamp.emplace_back( "Total runtime", m_timer[0].hms() );
          m_print.time( "Timers (h:m:s)", m_timestamp );
          m_print.perf( "Performance statistics", m_perf );
          m_print.endpart();
        }
      } catch (...) { tk::processExceptionCharm(); }
      // Tell the Charm++ runtime system to exit
      CkExit();
    }

    //! Add a time stamp contributing to final timers output
    void timestamp( std::string label, tk::real stamp ) {
      try {
        m_timestamp.emplace_back( label, tk::hms( stamp ) );
      } catch (...) { tk::processExceptionCharm(); }
    }
    //! Add multiple time stamps contributing to final timers output
    void timestamp( const std::vector< std::pair< std::string, tk::real > >& s )
    { for (const auto& t : s) timestamp( t.first, t.second ); }

    //! Add a performance statistic contributing to final perfstat output
    void perfstat( std::string label, tk::real value ) {
      try {
        m_perf.emplace_back( label, value );
      } catch (...) { tk::processExceptionCharm(); }
    }
    //! Add multiple performance statistics contributing to final perf output
    void perfstat( const std::vector< std::pair< std::string, tk::real > >& p )
    { for (const auto& s : p) perfstat( s.first, s.second ); }

  private:
    inciter::ctr::CmdLine m_cmdline;            //!< Command line
    inciter::CmdLineParser m_cmdParser;         //!< Command line parser
    inciter::InciterPrint m_print;              //!< Pretty printer
    inciter::InciterDriver m_driver;            //!< Driver
    std::vector< tk::Timer > m_timer;           //!< Timers

    //! Time stamps in h:m:s with labels
    std::vector< std::pair< std::string, tk::Timer::Watch > > m_timestamp;

    //! Performance statistics
    std::vector< std::pair< std::string, tk::real > > m_perf;
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

#include "inciter.def.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif
