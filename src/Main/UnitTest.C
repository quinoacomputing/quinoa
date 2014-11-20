//******************************************************************************
/*!
  \file      src/Main/UnitTest.C
  \author    J. Bakosi
  \date      Wed 19 Nov 2014 04:55:17 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     UnitTest: Quinoa's unit test suite
  \details   UnitTest: Quinoa's unit test suite
*/
//******************************************************************************

#include <pup_stl.h>

#include <Config.h>
#include <UnitTestPrint.h>
#include <UnitTestDriver.h>
#include <UnitTest/CmdLine/Parser.h>
#include <unittest.decl.h>
#include <Init.h>

// Unit test groups to be tested
#include <tests/Base/flip_map.h>
#include <tests/Base/make_list.h>
#include <tests/Base/StrConvUtil.h>
#include <tests/Base/Timer.h>
#include <tests/Base/CharmUtil.h>
#include <tests/Base/Factory.h>
#include <tests/Base/Print.h>
#include <tests/Base/TaggedTuple.h>
#include <tests/Base/Exception.h>
#include <tests/Base/PUPUtil.h>
#include <tests/Control/FileParser.h>
#include <tests/Control/StringParser.h>
#include <tests/Control/Toggle.h>
#include <tests/Control/Options/MKLGaussianMethod.h>
#include <tests/Control/Options/MKLUniformMethod.h>
#include <tests/Control/Options/RNG.h>
#include <tests/IO/Reader.h>

//! Charm handle to the main proxy, facilitates call-back to finalize, etc.,
//! must be in global scope, unique per executable
CProxy_Main mainProxy;

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

//! Pack/Unpack test runner
inline void operator|( PUP::er& p, tut::test_runner_singleton& runner )
{ if (!p.isSizing()) runner = tut::test_runner_singleton(); }

} // unittest::

//! Charm++ main chare
class Main : public CBase_Main {

  public:
    Main( CkArgMsg* msg )
    try :
      // Parse command line into m_cmdline using default simple pretty printer
      m_cmdParser( msg->argc, msg->argv, tk::Print(), m_cmdline ),
      // Create pretty printer initializing output streams based on command line
      m_print( m_cmdline.get< tk::tag::verbose >() ? std::cout : std::clog ),
      // Create UnitTest driver
      m_driver( tk::Main< unittest::UnitTestDriver >
                        ( msg->argc, msg->argv,
                          m_cmdline,
                          tk::HeaderType::UNITTEST,
                          UNITTEST_EXECUTABLE,
                          m_print ) ),
      m_timer(1)        // Start new timer measuring the total runtime
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
    } catch (...) { processException(); }

    void execute() {
      m_timestamp.emplace( "Migration of global-scope data + fire up all tests",
                           m_timer[1].hms() );
      m_driver.execute();       // fires up async chares
    }

    void finalize() {
      if (!m_timer.empty()) {
        m_timestamp.emplace( "Total runtime", m_timer[0].hms() );
        m_print.time( "Timers (h:m:s)", m_timestamp );
        m_print.endpart();
      }
      CkExit();
    }

    //! Process an exception
    void processException() {
      try {
        throw;      // rethrow exception to deal with it here
      }
        // Catch Quina::Exceptions
        catch ( tk::Exception& qe ) {
          qe.handleException();
        }
        // Catch std::exception and transform it into Quinoa::Exception without
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
    unittest::ctr::CmdLine m_cmdline;                   //!< Command line
    unittest::CmdLineParser m_cmdParser;                //!< Command line parser
    unittest::UnitTestPrint m_print;                    //!< Pretty printer
    unittest::UnitTestDriver m_driver;                  //!< Driver
    std::vector< tk::Timer > m_timer;                   //!< Timers

    //! Time stamps in h:m:s with labels
    std::map< std::string, tk::Timer::Watch > m_timestamp;
};

//! Charm++ chare execute: by the time this object is constructed, the Charm++
//! runtime system has finished migrating all global-scoped read-only objects
//! which happens after the main chare constructor has finished.
struct execute : CBase_execute { execute() { mainProxy.execute(); } };

#include <charmchild.def.h>
#include <migrated.def.h>
#include <unittest.def.h>
