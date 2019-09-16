// *****************************************************************************
/*!
  \file      src/Main/FileConv.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     File converter Charm++ main chare
  \details   File converter Charm++ main chare. This file contains the
    definition of the Charm++ main chare, equivalent to main() in Charm++-land.
*/
// *****************************************************************************

#include <vector>
#include <utility>
#include <iostream>

#include "Print.hpp"
#include "Timer.hpp"
#include "Types.hpp"
#include "QuinoaConfig.hpp"
#include "Init.hpp"
#include "Tags.hpp"
#include "FileConvDriver.hpp"
#include "FileConv/CmdLine/CmdLine.hpp"
#include "FileConv/CmdLine/Parser.hpp"
#include "ProcessException.hpp"
#include "ChareStateCollector.hpp"

#include "NoWarning/charm.hpp"
#include "NoWarning/fileconv.decl.h"

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

//! \brief Charm++ main chare for the file converter executable, fileconv.
//! \details Note that this object should not be in a namespace.
// cppcheck-suppress noConstructor
class Main : public CBase_Main {

  public:
    //! \brief Constructor
    //! \details FileConv's main chare constructor is the entry point of the
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
    //!   the global-scope data. Then we are ready to execute the driver which
    //!   calls back to Main::finalize() when it finished. Then finalize() exits
    //!   by calling Charm++'s CkExit(), shutting down the runtime system.
    //! \see http://charm.cs.illinois.edu/manuals/html/charm++/manual.html
    Main( CkArgMsg* msg )
    try :
      m_signal( tk::setSignalHandlers() ),
      m_cmdline(),
      // Parse command line into m_cmdline using default simple pretty printer
      m_cmdParser( msg->argc, msg->argv, tk::Print(), m_cmdline ),
      // Create pretty printer initializing output streams based on command line
      m_print( m_cmdline.get< tag::verbose >() ? std::cout : std::clog ),
      // Create FileConv driver
      m_driver( tk::Main< fileconv::FileConvDriver >
                        ( msg->argc, msg->argv,
                          m_cmdline,
			  tk::HeaderType::FILECONV,
			  tk::fileconv_executable(),
			  m_print ) ),
      m_timer(1),       // Start new timer measuring the total runtime
      m_timestamp()
    {
      delete msg;
      g_trace = m_cmdline.get< tag::trace >();
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
    } catch (...) { tk::processExceptionCharm(); }

    void execute() {
      try {
        m_timestamp.emplace_back("Migrate global-scope data", m_timer[1].hms());
        m_driver.execute();
      } catch (...) { tk::processExceptionCharm(); }
    }

    //! Towards normal exit but collect chare state first (if any)
    void finalize() {
      tk::finalize( m_cmdline, m_timer, m_print, stateProxy, m_timestamp,
                    CkCallback( CkIndex_Main::dumpstate(nullptr), thisProxy ) );
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

    //! Add a time stamp contributing to final timers output
    void timestamp( std::string label, tk::real stamp ) {
      try {
        m_timestamp.emplace_back( label, tk::hms( stamp ) );
      } catch (...) { tk::processExceptionCharm(); }
    }
    //! Add multiple time stamps contributing to final timers output
    void timestamp( const std::vector< std::pair< std::string, tk::real > >& s )
    { for (const auto& t : s) timestamp( t.first, t.second ); }

  private:
    int m_signal;                               //!< Used to set signal handlers
    fileconv::ctr::CmdLine m_cmdline;           //!< Command line
    fileconv::CmdLineParser m_cmdParser;        //!< Command line parser
    tk::Print m_print;                          //!< Pretty printer
    fileconv::FileConvDriver m_driver;          //!< Driver
    std::vector< tk::Timer > m_timer;           //!< Timers

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

#include "NoWarning/fileconv.def.h"
