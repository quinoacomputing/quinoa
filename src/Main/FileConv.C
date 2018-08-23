// *****************************************************************************
/*!
  \file      src/Main/FileConv.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     File converter Charm++ main chare
  \details   File converter Charm++ main chare. This file contains the
    definition of the Charm++ main chare, equivalent to main() in Charm++-land.
*/
// *****************************************************************************

#include <vector>
#include <utility>
#include <iostream>

#include "Print.h"
#include "Timer.h"
#include "Types.h"
#include "QuinoaConfig.h"
#include "Init.h"
#include "Tags.h"
#include "FileConvDriver.h"
#include "FileConv/CmdLine/CmdLine.h"
#include "FileConv/CmdLine/Parser.h"
#include "ProcessException.h"
#include "ChareStateCollector.h"

#include "NoWarning/charm.h"
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

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

//! \brief Charm++ main chare for the file converter executable, fileconv.
//! \details Note that this object should not be in a namespace.
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
      mainProxy = thisProxy;
      // If quiescence detection is on or user requested it, create chare state
      // collector Charm++ chare group
      if (m_cmdline.get< tag::chare >() || m_cmdline.get< tag::quiescence >())
        stateProxy = tk::CProxy_ChareStateCollector::ckNew();
      // Optionally enable quiscence detection
      if (m_cmdline.get< tag::quiescence >())
        CkStartQD( CkCallback( CkIndex_Main::quiescence(), thisProxy ) );
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
        m_driver.execute();
      } catch (...) { tk::processExceptionCharm(); }
    }

    //! Towards normal exit but collect chare state first (if any)
    void finalize() {
      try {
        if (!m_timer.empty()) {
          m_timestamp.emplace_back( "Total runtime", m_timer[0].hms() );
          m_print.time( "Timers (h:m:s)", m_timestamp );
          m_print.endpart();
         // If quiescence detection is on or user requested it, collect chare
         // state
         if (m_cmdline.get< tag::chare >()||m_cmdline.get< tag::quiescence >())
           stateProxy.collect( /* error = */ false,
             CkCallback( CkIndex_Main::dumpstate(nullptr), thisProxy ) );
         else
           CkExit(); // tell the Charm++ runtime system to exit
        }
      } catch (...) { tk::processExceptionCharm(); }
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
      try {
        std::unordered_map< int, std::vector< tk::ChareState > > state;
        PUP::fromMem creator( msg->getData() );
        creator | state;
        delete msg;
        // find out if chare state collection was triggered due to an error
        auto it = state.find( -1 );
        bool error = it != end(state);
        if (error) state.erase( it );
        // pretty-print collected chare state (only if user requested it or
        // quiescence was detected which is and indication of a logic error)
        if (m_cmdline.get< tag::chare >() || error)
          m_print.charestate( state );
        // exit differently depending on how we were called
        if (error)
          Throw( "Quiescence detected" );
        else
          CkExit(); // tell the Charm++ runtime system to exit
      } catch (...) { tk::processExceptionCharm(); }
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
