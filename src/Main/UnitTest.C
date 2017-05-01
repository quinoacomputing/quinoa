// *****************************************************************************
/*!
  \file      src/Main/UnitTest.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
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

#include "NoWarning/tut_result.h"
#include "NoWarning/tut_runner.h"

#include "NoWarning/charm.h"
#include "NoWarning/mpi.h"
#include "NoWarning/mpi-interoperate.h"

#include "NoWarning/tutsuite.decl.h"
#include "NoWarning/unittest.decl.h"

#include "Print.h"
#include "Timer.h"
#include "Tags.h"
#include "Exception.h"
#include "Init.h"
#include "QuinoaConfig.h"
#include "HelpFactory.h"
#include "Assessment.h"
#include "ProcessException.h"
#include "UnitTest/CmdLine/CmdLine.h"
#include "TestArray.h"
#include "UnitTestPrint.h"
#include "UnitTestDriver.h"
#include "UnitTest/CmdLine/Parser.h"

namespace tut {

//! \brief Maximum number of tests in every test group to attempt to run
//! \details If any of the unit test groups have more tests than this number,
//!   this should be increased. All test groups included below will use this
//!   value to override the default template argument for tut::test_group<>.
const int MAX_TESTS_IN_GROUP = 80;

} // tut::

// // Unit test groups to be tested. Each file defines a different test group.
#include "tests/Base/TestCharmUtil.h"
#include "tests/Base/TestFactory.h"
#include "tests/Base/TestTimer.h"
#include "tests/Base/TestPUPUtil.h"

#include "tests/Base/TestFlip_map.h"
#include "tests/Base/TestMake_list.h"
#include "tests/Base/TestHas.h"
#include "tests/Base/TestData.h"
#include "tests/Base/TestPrint.h"
#include "tests/Base/TestTaggedTuple.h"
#include "tests/Base/TestException.h"
#include "tests/Base/TestExceptionMPI.h"
#include "tests/Base/TestReader.h"
#include "tests/Base/TestStrConvUtil.h"
#include "tests/Base/TestWriter.h"
#include "tests/Base/TestProcessControl.h"
#include "tests/Base/TestVector.h"
#include "tests/Base/TestContainerUtil.h"

#include "tests/Control/TestSystemComponents.h"
#include "tests/Control/TestControl.h"
#include "tests/Control/TestFileParser.h"
#include "tests/Control/TestStringParser.h"
#include "tests/Control/TestToggle.h"
#ifdef HAS_MKL
  #include "tests/Control/Options/TestMKLUniformMethod.h"
  #include "tests/Control/Options/TestMKLGaussianMethod.h"
  #include "tests/Control/Options/TestMKLBetaMethod.h"
#endif
#include "tests/Control/Options/TestRNG.h"

#include "tests/IO/TestMesh.h"
#include "tests/IO/TestExodusIIMeshReader.h"

#include "tests/Mesh/TestDerivedData.h"

#include "tests/RNG/TestRNG.h"
#ifdef HAS_MKL
  #include "tests/RNG/TestMKLRNG.h"
#endif
#ifdef HAS_RNGSSE2
#include "tests/RNG/TestRNGSSE.h"
#endif
#include "tests/RNG/TestRandom123.h"

#include "tests/LoadBalance/TestLoadDistributor.h"
// Disabled due to API change.
// See https://lists.cs.illinois.edu/lists/arc/charm/2017-01/msg00018.html.
#include "tests/LoadBalance/TestLinearMap.h"
#include "tests/LoadBalance/TestUnsMeshMap.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wmissing-variable-declarations"
#endif

//! \brief Charm handle to the main proxy, facilitates call-back to finalize,
//!    etc., must be in global scope, unique per executable
CProxy_Main mainProxy;

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

//! Bool indicating whether all Charm++ and serial tests have passed
bool g_charmpass = true;

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
                          UNITTEST_EXECUTABLE,
                          m_print ) ),
      m_timer(1), // Start new timer measuring the serial+Charm++ runtime
      m_timestamp()
    {
      // Immediately exit if any help was requested; help is printed in main()
      if (m_helped) CkExit();
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
    } catch (...) { tk::processExceptionCharm(); }

    void execute() {
      try {
        m_timestamp.emplace_back("Migrate global-scope data", m_timer[1].hms());
        m_driver.execute();       // fires up async chares
      } catch (...) { tk::processExceptionCharm(); }
    }

    void finalize( bool worked, bool pass ) {
      try {
        if (worked && !m_timer.empty()) {
          m_timestamp.emplace_back( "Serial and Charm++ tests runtime",
                                    m_timer[0].hms() );
          m_print.time( "Serial and Charm++ test suite timers (h:m:s)",
                        m_timestamp );
          m_print.endpart();
        }
      } catch (...) { tk::processExceptionCharm(); }
      // Set global bool indicating whether all tests have passed
      unittest::g_charmpass = pass;
      // Tell the Charm++ runtime system to exit
      CkExit();
    }

  private:
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
//! \author J. Bakosi
class execute : public CBase_execute {
 public: execute() { mainProxy.execute(); }
};

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
  
  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wold-style-cast"
  #endif

  // Initialize MPI
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &peid );
  MPI_Comm_size( MPI_COMM_WORLD, &numpes );

  // Run serial and Charm++ unit test suite
  CharmLibInit( MPI_COMM_WORLD, argc, argv );
  CharmLibExit();

  bool mpipass = true;

  // Lambda to compute exit code based on test failures and exit. This is the
  // single exit point and we must exit from the program.
  auto stop = [numpes](int pass) {
    // Combine pass-status from Charm++/serial and MPI suites
    int mypass = unittest::g_charmpass && pass ? 1 : 0;
    // Add up every PE's pass status
    int g = 0;
    MPI_Allreduce( &mypass, &g, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
    // Exit code: g==numpes: pass, g<numpes: fail
    g /= numpes;
    if (g == 1)
      MPI_Finalize();   // will pass exit code 0 to shell
    else
      MPI_Abort( MPI_COMM_WORLD, tk::ErrCode::FAILURE ); // nonzero exit code
    // Since this is an MPI program, the exit code passed to exit() does not
    // matter, however, calling exit() here is important, because we must exit.
    exit( tk::ErrCode::SUCCESS );
  };

  #if defined(__clang__)
    #pragma clang diagnostic pop
  #endif

  // Run MPI test suite
  try {

    tk::Print print;    // quiet output by default using print, see ctor
    unittest::ctr::CmdLine cmdline;
    bool helped;
    unittest::CmdLineParser cmdParser( argc, argv, print, cmdline, helped );

    // Print out help on all command-line arguments if help was requested
    const auto helpcmd = cmdline.get< tag::help >();
    if (peid == 0 && helpcmd)
      print.help< tk::QUIET >( UNITTEST_EXECUTABLE,
                               cmdline.get< tag::cmdinfo >(),
                               "Command-line Parameters:", "-" );

    // Print out verbose help for a single keyword if requested
    const auto helpkw = cmdline.get< tag::helpkw >();
    if (peid == 0 && !helpkw.keyword.empty())
      print.helpkw< tk::QUIET >( UNITTEST_EXECUTABLE, helpkw );

    // Immediately exit if any help was output
    if (helpcmd || !helpkw.keyword.empty()) stop( mpipass );

    unittest::UnitTestPrint
      uprint( cmdline.get< tag::verbose >() ? std::cout : std::clog );

    const auto& groups = unittest::g_runner.get().list_groups();

    // Get group name string passed in by -g
    const auto grp = cmdline.get< tag::group >();

    // If only select groups to be run, see if there is any that will run
    bool work = false;
    if (grp.empty())
      work = true;
    else
      for (const auto& g : groups)
        if ( g.find("MPI") != std::string::npos &&  // only consider MPI groups
             g.find(grp) != std::string::npos )
          work = true;

    // Quit if there is no work to be done
    if (!work) {
      if (peid == 0)
        uprint.note( "\nNo MPI test groups to be executed because no test "
                     "group names match '" + grp + "'.\n" );
      stop( mpipass );
    }

    if (peid == 0) {
      uprint.endpart();
      uprint.part( "MPI unit test suite" );
      uprint.unithead( "Unit tests computed", cmdline.get< tag::group >() );
    }

    std::size_t nrun=0, ncomplete=0, nwarn=0, nskip=0, nexcp=0, nfail=0;
    tk::Timer timer;  // start new timer measuring the MPI-suite runtime

    // Lambda to fire up all tests in a test group
    auto spawngrp = [&]( const std::string& g ) {
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
    };

    // Fire up all tests in all test groups exercising MPI on rank 0
    for (const auto& g : groups)
      if (g.find("MPI") != std::string::npos) { // only start MPI test groups
        if (grp.empty()) {                        // consider all test groups
          spawngrp( g );
        } else if (g.find(grp) != std::string::npos) {
          // spawn only the groups that match the string specified via -g string
          spawngrp( g );
        }
      }

    if (peid == 0) {
      mpipass =
       unittest::assess( uprint, "MPI", nfail, nwarn, nskip, nexcp, ncomplete );
      std::vector< std::pair< std::string, tk::Timer::Watch > > timestamp;
      timestamp.emplace_back( "MPI tests runtime", timer.hms() );
      uprint.time( "MPI test suite timers (h:m:s)", timestamp );
    }

  } catch (...) { tk::processExceptionMPI(); }

  stop( mpipass );
}

#include "NoWarning/charmchild.def.h"
#include "NoWarning/charmtimer.def.h"
#include "NoWarning/migrated.def.h"
#include "NoWarning/testarray.def.h"
#include "NoWarning/unittest.def.h"
