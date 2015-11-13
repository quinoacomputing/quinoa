//******************************************************************************
/*!
  \file      src/Main/Inciter.C
  \author    J. Bakosi
  \date      Tue 10 Nov 2015 07:47:01 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Inciter, computational shock hydrodynamics tool, Charm++ main
    chare.
  \details   Inciter, computational shock hydrodynamics tool, Charm++ main
    chare. This file contains the definition of the Charm++ main chare,
    equivalent to main() in Charm++-land.
*/
//******************************************************************************

#include <map>
#include <vector>
#include <iostream>
#include <utility>
#include <type_traits>
#include <cstddef>

#include <boost/format.hpp>

#include "Types.h"
#include "Writer.h"
#include "Timer.h"
#include "Exception.h"
#include "ProcessException.h"
#include "HelpFactory.h"
#include "InciterPrint.h"
#include "InciterDriver.h"
#include "InciterSetup.h"
#include "Inciter/CmdLine/CmdLine.h"
#include "Inciter/InputDeck/InputDeck.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <mpi.h>
#include <mpi-interoperate.h>
#include <pup_stl.h>
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

//! Command line filled by command line parser, with command line user input
ctr::CmdLine g_cmdline;
//! Defaults of input deck, facilitates detection what is set by user
ctr::InputDeck g_inputdeck_defaults;
//! Input deck filled by parser, containing all input data
ctr::InputDeck g_inputdeck;

//! \brief Total number of points in computational mesh
//! \details While this data is declared in global scope (so that Charm++
//!   chares can access it), it is intentionally NOT declared in the Charm++
//!   main module interface file for Inciter in Main/inciter.ci, so that the
//!   Charm++ runtime system does not migrate it across all PEs. This is okay,
//!   since there is no need for any of the other Charm++ chares, other than the
//!   main chare below, to access it in the future.
std::size_t g_npoin;

//! \brief Total number of chares
//! \details While this data is declared in global scope (so that Charm++
//!   chares can access it), it is intentionally NOT declared in the Charm++
//!   main module interface file for Inciter in Main/inciter.ci, so that the
//!   Charm++ runtime system does not migrate it across all PEs. This is okay,
//!   since there is no need for any of the other Charm++ chares, other than the
//!   main chare below, to access it in the future.
uint64_t g_nchare;

//! \brief Global mesh element ids owned by each chare (associated to chare IDs)
//! \details This data holds different element IDs on different MPI ranks. It
//!   holds the global mesh element IDs assigned to each Charm++ chare. While
//!   this data is declared in global scope (so that Charm++ chares can access
//!   it), it is intentionally NOT declared in the Charm++ main module interface
//!   file for Inciter in Main/inciter.ci, so that the Charm++ runtime system
//!   does not migrate it across all PEs. This is intentional, since it will be
//!   used to initialize elements of a Charm++ group of which a single one is
//!   created on each PE. This data is generated in the initial MPI portion and
//!   transfer of this data from the initial MPI portion to the Charm++ portion
//!   is facilitated by it being in global-scope. While this data could be
//!   passed down to the point where the consumer Charm++ chare group is fired
//!   up, that would be the wrong thing to do, because that way a copy of the
//!   data generated on PE 0 would be sent to each PE by Charm++. Keeping this
//!   untouched in global scope, i.e., not listing in Main/inciter.ci as
//!   readonly data, the Charm++ chare group elements, created on different PEs,
//!   will simply access this data, correctly a different one in their own
//!   global scope, as intended.
std::map< int, std::vector< std::size_t > > g_element;

//! \brief Time stamps in h:m:s for the initial MPI portion
//! \details Time stamps collected here are those collected by the initial MPI
//!   portion and are displayed by the Charm++ main chare at the end. While this
//!   map of timers is declared in global scope (so that Charm++ chares can
//!   access it), it is intentionally NOT declared in the Charm++ main module
//!   interface file for Inciter in Main/inciter.ci, so that the Charm++ runtime
//!   system does not migrate it across all PEs. This is okay, since there is no
//!   need for any of the other Charm++ chares to access it in the future. In
//!   fact, the main chare grabs it and swallows it right away during its
//!   constructor for output at the end of the run.
std::vector< std::pair< std::string, tk::Timer::Watch > > g_timestamp;

//! Conductor Charm++ proxy facilitating call-back to Conductor by the
//! individual performers
CProxy_Conductor g_ConductorProxy;

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
      // Create pretty printer initializing output streams based on command line
      m_print( inciter::g_cmdline.get< tag::verbose >()
                 ? std::cout : std::clog ),
      // Create Inciter driver
      m_driver( inciter::InciterDriver( m_print ) ),
      // Start new timer measuring the total runtime
      m_timer(1),
      // Import, i.e., swallow, timers from the initial MPI portion
      m_timestamp( std::move(inciter::g_timestamp) )
    {
      const auto helpcmd = inciter::g_cmdline.get< tag::help >();
      const auto helpctr = inciter::g_cmdline.get< tag::helpctr >();
      const auto helpkw = inciter::g_cmdline.get< tag::helpkw >();
      // Exit if help was requested or exectuable was called without argument
      if (msg->argc == 1 || helpcmd || helpctr || !helpkw.keyword.empty()) {
        CkExit();
      } else {  // business as usual
        mainProxy = thisProxy;
        // Fire up an asynchronous execute object, which when created at some
        // future point in time will call back to this->execute(). This is
        // necessary so that this->execute() can access already migrated
        // global-scope data.
        CProxy_execute::ckNew();
        // Start new timer measuring the migration of global-scope data
        m_timer.emplace_back();
      }
      delete msg;
      m_print.diagstart( "Migrating read-only global-scope data ..." );
    } catch (...) { tk::processExceptionCharm(); }

    //! Execute driver created and initialized by constructor
    void execute() {
      try {
        m_print.diagend( "done" );
        m_timestamp.emplace_back( "Migrate Charm++ read-only global-scope data",
                                  m_timer[1].hms());
        m_driver.execute( inciter::g_npoin, inciter::g_nchare );
      } catch (...) { tk::processExceptionCharm(); }
    }

    //! Normal exit point
    void finalize() {
      try {
        if (!m_timer.empty()) {
          m_timestamp.emplace_back(
           "Total Charm++ runtime . . . . . . . . . . . . . . . . . . . . . . "
           ". . . . .",
           m_timer[0].hms() );
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
    inciter::InciterPrint m_print;                    //!< Pretty printer
    inciter::InciterDriver m_driver;                  //!< Driver
    std::vector< tk::Timer > m_timer;                 //!< Timers

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

//! \brief Inciter main()
//! \details Inciter does have a main() function so that Zoltan, an MPI library,
//!   can partition the mesh. Then we initialize the Charm++ runtime system and
//!   do a calculation. This is necessary, since MPI_Init() is a bit adamant
//!   about capturing resources it wants and hence it has to be called before
//!   Charm is initialized.
//! \author J. Bakosi
int main( int argc, char **argv ) {

  using namespace inciter;

  tk::Timer mpi;        // start timing the MPI portion

  // Initialize MPI
  MPI_Init( &argc, &argv );

  try {
    // Parse command line
    parseCmdLine( argc, argv, g_cmdline );
   
    // Instantiate inciter's pretty printer
    InciterPrint
      iprint( g_cmdline.get< tag::verbose >() ? std::cout : std::clog );

    // Parse input deck, echo info
    init( g_cmdline, iprint, g_inputdeck, argc, argv );

    // Prepare computational mesh and fill some global-scope data (all given at
    // the top of this file) so that Charm++ chares can access them
    prepareMesh( g_cmdline, iprint, g_inputdeck,  // <- const
                 g_timestamp, g_nchare, g_npoin, g_element );

  } catch (...) { tk::processExceptionMPI(); }

  g_timestamp.emplace_back(
    "Total MPI runtime . . . . . . . . . . . . . . . . . . . . . . . . . . "
    ". . .",
    mpi.hms() );

  // Run Charm++ main chare using the partitioned graph
  CharmLibInit( MPI_COMM_WORLD, argc, argv );
  //MPI_Barrier( MPI_COMM_WORLD );
  // ...
  //MPI_Barrier( MPI_COMM_WORLD );
  CharmLibExit();

  // Finalize MPI
  MPI_Barrier( MPI_COMM_WORLD );
  MPI_Finalize();

  return tk::ErrCode::SUCCESS;
}

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "inciter.def.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif
