// *****************************************************************************
/*!
  \file      src/Main/Walker.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Random walker Charm++ main chare
  \details   Random walker Charm++ main chare. This file contains the definition
    of the Charm++ main chare, equivalent to main() in Charm++-land.
*/
// *****************************************************************************

#include <map>
#include <iostream>
#include <utility>
#include <vector>

#include "NoWarning/format.hpp"

#include "NoWarning/charm.hpp"
#include "NoWarning/pup.hpp"

#include "Print.hpp"
#include "Timer.hpp"
#include "Types.hpp"
#include "Init.hpp"
#include "QuinoaConfig.hpp"
#include "Tags.hpp"
#include "ProcessException.hpp"
#include "RNG.hpp"
#include "RNGStack.hpp"
#include "DiffEq.hpp"
#include "DiffEqStack.hpp"
#include "Options/RNG.hpp"
#include "WalkerPrint.hpp"
#include "WalkerDriver.hpp"
#include "Walker/CmdLine/Parser.hpp"
#include "Walker/CmdLine/CmdLine.hpp"
#include "Walker/InputDeck/InputDeck.hpp"
#include "ChareStateCollector.hpp"

#include "NoWarning/walker.decl.h"

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

//! Walker declarations and definitions
namespace walker {

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

//! Defaults of input deck, facilitates detection what is set by user
//! \details This object is in global scope, it contains the default of all
//!   possible user input, and thus it is made available to all PEs for
//!   convenience reasons. The runtime system distributes it to all PEs during
//!   initialization. Once distributed, the object does not change.
ctr::InputDeck g_inputdeck_defaults;
//! Input deck filled by parser, containing all input data
//! \details This object is in global scope, it contains all of user input, and
//!   thus it is made available to all PEs for convenience reasons. The runtime
//!   system distributes it to all PEs during initialization. Once distributed,
//!   the object does not change.
ctr::InputDeck g_inputdeck;
//! Random number generators selected by user
//! \details This map is in global scope, because it holds polymorphic
//!   objects, and thus must be distributed to all PEs during initialization.
//!   Once distributed by the runtime system, the objects do not change.
std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;
//! Differential equations selected by user
//! \details This vector is in global scope, because it holds polymorphic
//!   objects, and thus must be distributed to all PEs during initialization.
//!   Once distributed by the runtime system, the objects do not change.
std::vector< DiffEq > g_diffeqs;

//! Distributor Charm++ proxy facilitating call-back to Distributor by the
//! individual integrators
CProxy_Distributor g_DistributorProxy;

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

//! Pack/Unpack selected RNGs. This Pack/Unpack method (re-)creates the full RNG
//! stack since it needs to (re-)bind function pointers on different processing
//! elements. Therefore we circumvent Charm's usual pack/unpack for this type,
//! and thus sizing does not make sense: sizing is a no-op. We could initialize
//! the stack in RNGTestDriver's constructor and let this function re-create the
//! stack only when unpacking, but that leads to repeating the same code twice:
//! once in RNGTestDriver's constructor, once here. Another option is to use
//! this pack/unpack routine to both initially create (when packing) and to
//! re-create (when unpacking) the stack, which eliminates the need for
//! pre-creating the object in RNGTestDriver's constructor and therefore
//! eliminates the repeated code. This explains the guard for sizing: the code
//! below is called for packing only (in serial) and packing and unpacking (in
//! parallel).
inline
void operator|( PUP::er& p, std::map< tk::ctr::RawRNGType, tk::RNG >& rng ) {
  try {
    if (!p.isSizing()) {
      tk::RNGStack stack(
        #ifdef HAS_MKL
        g_inputdeck.get< tag::param, tag::rngmkl >(),
        #endif
        #ifdef HAS_RNGSSE2
        g_inputdeck.get< tag::param, tag::rngsse >(),
        #endif
        g_inputdeck.get< tag::param, tag::rng123 >() );
      rng = stack.selected( g_inputdeck.get< tag::selected, tag::rng >() );
    }
  } catch (...) { tk::processExceptionCharm(); }
}

//! Pack/Unpack selected differential equations. This Pack/Unpack method
//! (re-)creates the DiffEq factory since it needs to (re-)bind function
//! pointers on different processing elements. Therefore we circumvent Charm's
//! usual pack/unpack for this type, and thus sizing does not make sense: sizing
//! is a no-op. We could initialize the factory in WalkerDriver's constructor
//! and let this function re-create the stack only when unpacking, but that
//! leads to repeating the same code twice: once in WalkerDriver's constructor,
//! once here. Another option is to use this pack/unpack routine to both
//! initially create (when packing) and to re-create (when unpacking) the
//! factory, which eliminates the need for pre-creating the object in
//! WalkerDriver's constructor and therefore eliminates the repeated code. This
//! explains the guard for sizing: the code below is called for packing only (in
//! serial) and packing and unpacking (in parallel).
inline
void operator|( PUP::er& p, std::vector< DiffEq >& eqs ) {
  try {
    if (!p.isSizing()) eqs = DiffEqStack().selected();
  } catch (...) { tk::processExceptionCharm(); }
}

} // walker::

//! \brief Charm++ main chare for the random walker executable, walker.
//! \details Note that this object should not be in a namespace.
// cppcheck-suppress noConstructor
class Main : public CBase_Main {

  public:
    //! \brief Constructor
    //! \details Walker's main chare constructor is the entry point of the
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
    //!   the random walker is parallel and asynchronous, its driver fires up
    //!   additional Charm++ chare objects which then call back to
    //!   Main::finalize() at some point in the future when all work has been
    //!   finished. finalize() then exits by calling Charm++'s CkExit(),
    //!   shutting down the runtime system.
    //! \see http://charm.cs.illinois.edu/manuals/html/charm++/manual.html
    Main( CkArgMsg* msg )
    try :
      m_signal( tk::setSignalHandlers() ),
      m_cmdline(),
      // Parse command line into m_cmdline using default simple pretty printer
      m_cmdParser( msg->argc, msg->argv, tk::Print(), m_cmdline ),
      // Create pretty printer initializing output streams based on command line
      m_print( m_cmdline.get< tag::verbose >() ? std::cout : std::clog ),
      // Create Walker driver
      m_driver( tk::Main< walker::WalkerDriver >
                        ( msg->argc, msg->argv,
                          m_cmdline,
                          tk::HeaderType::WALKER,
                          tk::walker_executable(),
                          m_print ) ),
      m_timer(1),       // start new timer measuring the total runtime
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

    //! Execute driver created and initialized by constructor
    void execute() {
      try {
        m_timestamp.emplace_back("Migrate global-scope data", m_timer[1].hms());
        m_driver.execute();       // fires up async chares
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

    //! Add time stamp contributing to final timers output
    void timestamp( std::string label, tk::real stamp ) {
      try{
        m_timestamp.emplace_back( label, tk::hms( stamp ) );
      } catch (...) { tk::processExceptionCharm(); }
    }

  private:
    int m_signal;                               //!< Used to set signal handlers
    walker::ctr::CmdLine m_cmdline;             //!< Command line
    walker::CmdLineParser m_cmdParser;          //!< Command line parser
    walker::WalkerPrint m_print;                //!< Pretty printer
    walker::WalkerDriver m_driver;              //!< Driver
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

#include "NoWarning/walker.def.h"
