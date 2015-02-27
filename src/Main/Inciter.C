//******************************************************************************
/*!
  \file      src/Main/Inciter.C
  \author    J. Bakosi
  \date      Fri 27 Feb 2015 04:01:37 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Inciter, computational shock hydrodynamics tool, Charm++ main
    chare.
  \details   Inciter, computational shock hydrodynamics tool, Charm++ main
    chare. This file contains the definition of the Charm++ main chare,
    equivalent to main() in Charm++-land.
*/
//******************************************************************************

#include <mpi.h>

#include <mpi-interoperate.h>   // For interoperation of MPI and Charm++
#include <pup_stl.h>

#include <Config.h>
#include <RNG.h>
#include <RNGStack.h>
#include <InciterPrint.h>
#include <InciterDriver.h>
#include <Inciter/CmdLine/Parser.h>
#include <ZoltanInterOp.h>
#include <inciter.decl.h>
#include <Init.h>

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
//! \brief Mesh is global so that it is accessible to Zoltan
//! \details Zoltan is an MPI library interoperating with Charm++. However,
//!    g_mesh is NOT declared in Inciter's Charm++ main chare in Main/inciter.ci
//!    so that Charm does not migrate the entire mesh across all PEs. The mesh
//!    is read in on MPI rank 0 and is partitioned by Zoltan, then Charm++
//!    chares will be created containing pieces of the partitioned mesh.
tk::UnsMesh g_mesh;

//! Distributor Charm++ proxy facilitating call-back to Distributor by the
//! individual integrators
//CProxy_Distributor g_DistributorProxy;

} // inciter::

//! \brief Charm++ main chare for the computational fluid dynamics executable,
//!    inciter.
//! \details Note that this object should not be in a namespace.
class Main : public CBase_Main {

  public:
    //! \brief Constructor
    //! \details The main chare constructor is the main entry point of the
    //!   program, called by the Charm++ runtime system. The constructor does
    //!   basic initialization steps, e.g., parser the command-line, prints out
    //!   some useful information to screen (in verbose mode), and instantiates
    //!   a driver. Since Charm++ is fully asynchronous, the constructure
    //!   usually spawns asynchronous objects and immediately exits. Thus in the
    //!   body of the main chare constructor we fire up an 'execute' chare,
    //!   which then calls back to Main::execute(). Finishing the main chare
    //!   constructor the Charm++ runtime system then starts the
    //!   network-migration of all global-scope data (if any). The execute chare
    //!   calling back to Main::execute() signals the end of the migration of
    //!   the global-scope data. Then we are ready to execute the driver. Since
    //!   the fluid dynamics tool is parallel and asynchronous, its driver fires
    //!   up additional Charm++ chare objects which then call back to
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
      // Create Inciter driver
      m_driver( tk::Main< inciter::InciterDriver >
                        ( msg->argc, msg->argv,
                          m_cmdline,
                          tk::HeaderType::INCITER,
                          INCITER_EXECUTABLE,
                          m_print ) ),
      m_timer(1)        // Start new timer measuring the total runtime
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
    } catch (...) { processException(); }

    //! Execute driver created and initialized by constructor
    void execute() {
      m_timestamp.emplace("Migration of global-scope data", m_timer[1].hms());
      m_driver.execute();       // fires up async chares
    }

    //! Normal exit point
    void finalize() {
      if (!m_timer.empty()) {
        m_timestamp.emplace( "Total runtime", m_timer[0].hms() );
        m_print.time( "Timers (h:m:s)", m_timestamp );
        m_print.endpart();
      }
      CkExit();
    }

    //! Add time stamp contributing to final timers output
    void timestamp( std::string label, tk::real stamp )
    { m_timestamp.emplace( label, tk::hms( stamp ) ); }

    //! Process an exception
    void processException() {
      try {
        throw;      // rethrow exception to deal with it here
      }
        // Catch Quina::Exceptions
        catch ( tk::Exception& qe ) {
          qe.handleException();
        }
        // Catch std::exception and transform it into tk::Exception without
        // file:line:func information
        catch ( std::exception& se ) {
          tk::Exception qe( se.what() );
          qe.handleException();
        }
        // Catch uncaught exception
        catch (...) {
          tk::Exception qe( "Non-standard exception" );
          qe.handleException();
        }

      // Tell the runtime system to exit
      finalize();
    }

  private:
    inciter::ctr::CmdLine m_cmdline;                  //!< Command line
    inciter::CmdLineParser m_cmdParser;               //!< Command line parser
    inciter::InciterPrint m_print;                    //!< Pretty printer
    inciter::InciterDriver m_driver;                  //!< Drive
    std::vector< tk::Timer > m_timer;                 //!< Timers

    //! Time stamps in h:m:s with labels
    std::map< std::string, tk::Timer::Watch > m_timestamp;
};

//! \brief Charm++ chare execute
//! \details By the time this object is constructed, the Charm++ runtime system
//!    has finished migrating all global-scoped read-only objects which happens
//!    after the main chare constructor has finished.
//! \author J. Bakosi
struct execute : CBase_execute { execute() { mainProxy.execute(); } };

//! \brief Inciter main()
//! \details Inciter does have a main() function so that Charm++ can
//!   interoperate with Zoltan, which is an MPI library. This is necessary,
//!   since MPI_Init() is a bit adamant about capturing resources it wants and
//!   hence it has to be called before Charm is initialized.
//! \author J. Bakosi
int main( int argc, char **argv ) {
  int peid, numpes;

  // Initialize MPI
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &peid );
  MPI_Comm_size( MPI_COMM_WORLD, &numpes );

  // Initialize Charm++
  CharmLibInit( MPI_COMM_WORLD, argc, argv );
  MPI_Barrier( MPI_COMM_WORLD );

  // Partition mesh using Zoltan
  tk::zoltan::partitionMesh( inciter::g_mesh );

  // Finalize Charm++
  CharmLibExit();
  MPI_Barrier( MPI_COMM_WORLD );

  // Finalize MPI
  MPI_Finalize();

  return 0;  
}

#include <inciter.def.h>
