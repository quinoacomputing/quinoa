//******************************************************************************
/*!
  \file      src/Main/UnitTest.C
  \author    J. Bakosi
  \date      Thu 12 Mar 2015 10:20:24 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     UnitTest's unit test suite Charm++ main chare.
  \details   UnitTest's unit test suite Charm++ main chare. This file contains
    the definition of the Charm++ main chare, equivalent to main() in Charm++-
    land.
*/
//******************************************************************************
#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <mpi.h>
#include <mpi-interoperate.h>   // for interoperation of MPI and Charm++

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#include <pup_stl.h>

#include <Config.h>
#include <UnitTestPrint.h>
#include <UnitTestDriver.h>
#include <UnitTest/CmdLine/Parser.h>
#include <ProcessException.h>
#include <Assessment.h>
#include <unittest.decl.h>
#include <Init.h>

// Unit test groups to be tested. Each file defines a different test group.
#include <tests/Base/flip_map.h>
#include <tests/Base/make_list.h>
#include <tests/Base/Timer.h>
#include <tests/Base/CharmUtil.h>
#include <tests/Base/Has.h>
#include <tests/Base/ParticleProperties.h>
#include <tests/Base/Factory.h>
#include <tests/Base/Print.h>
#include <tests/Base/TaggedTuple.h>
#include <tests/Base/Exception.h>
#include <tests/Base/ExceptionMPI.h>
#include <tests/Base/PUPUtil.h>
#include <tests/Base/Reader.h>
#include <tests/Base/StrConvUtil.h>
#include <tests/Base/Writer.h>
#include <tests/Base/LoadDistributor.h>

#include <tests/Control/Components.h>
#include <tests/Control/Control.h>
#include <tests/Control/FileParser.h>
#include <tests/Control/StringParser.h>
#include <tests/Control/Toggle.h>
#ifdef HAS_MKL
  #include <tests/Control/Options/MKLGaussianMethod.h>
  #include <tests/Control/Options/MKLUniformMethod.h>
#endif
#include <tests/Control/Options/RNG.h>

//! \brief Charm handle to the main proxy, facilitates call-back to finalize,
//!    etc., must be in global scope, unique per executable
CProxy_Main mainProxy;

//! UnitTest declarations and definitions
namespace unittest {

//! Global-scope data. Initialized by the main chare and distibuted to all PEs
//! by the Charm++ runtime system. Though semantically not const, all these
//! global data should be considered read-only. See also
//! http://charm.cs.illinois.edu/manuals/html/charm++/manual.html. The data
//! below is global-scope because they must be available to all PEs which could
//! be on different machines.

//! Template Unit Test test runner. This Pack/Unpack method (re-)creates the
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
tut::test_runner_singleton g_runner;

//! Test suite Charm++ proxy facilitating call-back to unit test suite by
//! individual unit tests spawning Charm++ chares
CProxy_TUTSuite g_suiteProxy;

//! UnitTest executable name. So that FileParser's unit tests can access a file
//! for opening.
std::string g_executable;

//! \brief Max number tests in every test group to attempt to run
//! \details If any of the unit test groups have larger than this number than
//!   this should be increased.
int g_maxTestsInGroup = 50;

//! Pack/Unpack test runner
inline void operator|( PUP::er& p, tut::test_runner_singleton& runner )
{ if (!p.isSizing()) runner = tut::test_runner_singleton(); }

} // unittest::

//! \brief Charm++ main chare for the unit test suite executable, unittest.
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
    //!   the unit test suite is parallel and asynchronous, its driver fires up
    //!   additional Charm++ chare objects which then call back to
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
      // Create UnitTest driver
      m_driver( tk::Main< unittest::UnitTestDriver >
                        ( msg->argc, msg->argv,
                          m_cmdline,
                          tk::HeaderType::UNITTEST,
                          UNITTEST_EXECUTABLE,
                          m_print ) ),
      m_timer(1)  // Start new timer measuring the serial+Charm++ runtime
    {
      // Save executable name to global-scope string so FileParser can access it
      unittest::g_executable = msg->argv[0];
      delete msg;
      mainProxy = thisProxy;
      // Fire up an asynchronous execute object, which when created at some
      // future point in time will call back to this->execute(). This is
      // necessary so that this->execute() can access already migrated
      // global-scope data.
      CProxy_execute::ckNew();
      // Start new timer measuring the migration of global-scope data
      m_timer.emplace_back();
    } catch (...) { tk::processException(); }

    void execute() {
      try {
        m_timestamp.emplace(
          "Migration of global-scope data + fire up all tests",
           m_timer[1].hms() );
        m_driver.execute();       // fires up async chares
      } catch (...) { tk::processException(); }
    }

    void finalize() {
      try {
        if (!m_timer.empty()) {
          m_timestamp.emplace( "Serial and Charm++ tests runtime",
                               m_timer[0].hms() );
          m_print.time( "Serial and Charm++ test suite timers (h:m:s)",
                        m_timestamp );
          m_print.endpart();
        }
      } catch (...) { tk::processException(); }
      // Tell the Charm++ runtime system to exit
      CkExit();
    }

  private:
    unittest::ctr::CmdLine m_cmdline;                   //!< Command line
    unittest::CmdLineParser m_cmdParser;                //!< Command line parser
    unittest::UnitTestPrint m_print;                    //!< Pretty printer
    unittest::UnitTestDriver m_driver;                  //!< Driver
    std::vector< tk::Timer > m_timer;                   //!< Timers

    //! Time stamps in h:m:s with labels
    std::map< std::string, tk::Timer::Watch > m_timestamp;
};

//! \brief Charm++ chare execute
//! \details By the time this object is constructed, the Charm++ runtime system
//!    has finished migrating all global-scoped read-only objects which happens
//!    after the main chare constructor has finished.
//! \author J. Bakosi
struct execute : CBase_execute { execute() { mainProxy.execute(); } };

//! \brief UnitTest main()
//! \details UnitTest does have a main() function so that we can have tests
//!   calling MPI functions. Thus we are using Charm++'s MPI-interoperation
//!   capability as would have to be done with interoperation with an MPI
//!   library. This is necessary, since MPI_Init() is a bit adamant about
//!   capturing resources it wants and hence it has to be called before Charm is
//!   initialized.
//! \author J. Bakosi
int main( int argc, char **argv ) {

  int peid, numpes;

  // Initialize MPI
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &peid );
  MPI_Comm_size( MPI_COMM_WORLD, &numpes );

  // Run serial and Charm++ unit test suite
  CharmLibInit( MPI_COMM_WORLD, argc, argv );
  CharmLibExit();

  // Run MPI test suite
  try {

    tk::Print print;    // quiet output by default using print, see ctr
    unittest::ctr::CmdLine cmdline;
    unittest::CmdLineParser cmdParser( argc, argv, print, cmdline );
    unittest::UnitTestPrint
      uprint( cmdline.get< tag::verbose >() ? std::cout : std::clog );

    if (peid == 0) {
      uprint.endpart();
      uprint.part( "MPI unit test suite" );
      uprint.unithead( "Unit tests computed" );
    }

    tk::Timer timer;  // start new timer measuring the MPI-suite runtime

    // Fire up all tests in all test groups exercising MPI on rank 0
    std::size_t nrun=0, ncomplete=0, nwarn=0, nskip=0, nexcp=0, nfail=0;
    for (const auto& g : unittest::g_runner.get().list_groups())
      if (g.find("MPI") != std::string::npos) // only start MPI test groups
        for (int t=1; t<=unittest::g_maxTestsInGroup; ++t) {
          tut::test_result tr;
          unittest::g_runner.get().run_test( g, t, tr );
          if (peid == 0) {
            ++nrun;
            std::vector< std::string > status
              { tr.group, tr.name, std::to_string(tr.result), tr.message,
                tr.exception_typeid };
            unittest::evaluate( status, ncomplete, nwarn, nskip, nexcp, nfail );
            uprint.test( ncomplete, nfail, status );
          }
        }

    if (peid == 0) {
      unittest::assess( uprint, "MPI", nfail, nwarn, nskip, nexcp, ncomplete );
      std::map< std::string, tk::Timer::Watch > timestamp;
      timestamp.emplace( "MPI tests runtime", timer.hms() );
      uprint.time( "MPI test suite timers (h:m:s)", timestamp );
    }

  } catch (...) { tk::processException(); }

  // Finalize MPI
  MPI_Finalize();

  return tk::ErrCode::SUCCESS;
}

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <charmchild.def.h>
#include <migrated.def.h>
#include <unittest.def.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif
