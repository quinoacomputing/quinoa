//******************************************************************************
/*!
  \file      src/Main/UnitTest.C
  \author    J. Bakosi
  \date      Fri 25 Jul 2014 11:54:47 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     UnitTest: Quinoa's unit test suite
  \details   UnitTest: Quinoa's unit test suite
*/
//******************************************************************************

#include <Config.h>
#include <Paradigm.h>
#include <UnitTestPrint.h>
#include <UnitTestDriver.h>
#include <UnitTest/CmdLine/Parser.h>
#include <TPLInfo/Boost.h>
#include <TPLInfo/OpenMP.h>
#include <unittest.decl.h>
#include <Init.h>

// Unit test groups to be tested
#include <tests/Base/flip_map.h>
#include <tests/Base/make_list.h>
#include <tests/Base/StrConvUtil.h>

//! Charm handle to the main proxy, facilitates call-back to finalize, etc.,
//! must be in global scope, unique per executable
CProxy_Main mainProxy;

namespace unittest {

void echoTPL(const tk::Print& print)
//******************************************************************************
//  Echo TPL version informaion for libs specific to UnitTest
//! \author  J. Bakosi
//******************************************************************************
{
  echoOpenMP( print, "OpenMP runtime" );
  echoBoost( print, "Boost C++ Libraries" );
  print.endpart();
}

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

//! Pack/Unpack test runner
inline void operator|( PUP::er& p, tut::test_runner_singleton& runner )
{ if (!p.isSizing()) runner = tut::test_runner_singleton(); }

} // unittest::

//! Charm++ main chare
class Main : public CBase_Main {

  public:
    Main( CkArgMsg* msg ) :
      // Parse command line into m_cmdline using default simple pretty printer
      m_cmdParser( msg->argc, msg->argv, tk::Print(), m_cmdline ),
      // Create pretty printer initializing output streams based on command line
      m_print( m_cmdline.get< tk::tag::verbose >() ? std::cout : tk::null ),
      // Create UnitTest driver
      m_driver( tk::Main< unittest::UnitTestDriver >
                        ( msg->argc, msg->argv,
                          m_cmdline,
                          tk::HeaderType::UNITTEST,
                          UNITTEST_EXECUTABLE,
                          m_print,
                          unittest::echoTPL ) ),
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
    }

    void execute() {
      m_timestamp.emplace( "Migration of global-scope data", m_timer[1].hms() );
      m_driver.execute();       // fires up async chares
    }

    void finalize() {
      m_timestamp.emplace( "Total runtime", m_timer[0].hms() );
      m_print.time( "Timers (h:m:s)", m_timestamp );
      m_print.endpart();
      CkExit();
    }

  private:
    unittest::ctr::CmdLine m_cmdline;                   //!< Command line
    unittest::CmdLineParser m_cmdParser;                //!< Command line parser
    unittest::UnitTestPrint m_print;                    //!< Pretty printer
    unittest::UnitTestDriver m_driver;                  //!< Driver
    std::vector< tk::Timer > m_timer;                   //!< Timers

    //! Time stamps in h:m:s with labels
    std::map< std::string, tk::Watch > m_timestamp;
};

//! Charm++ chare execute: by the time this object is constructed, the Charm++
//! runtime system has finished migrating all global-scoped read-only objects
//! which happens after the main chare constructor has finished.
struct execute : CBase_execute { execute() { mainProxy.execute(); } };

#include <unittest.def.h>
