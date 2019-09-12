// *****************************************************************************
/*!
  \file      src/Main/Inciter.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter, computational shock hydrodynamics tool, Charm++ main
    chare.
  \details   Inciter, computational shock hydrodynamics tool, Charm++ main
    chare. This file contains the definition of the Charm++ main chare,
    equivalent to main() in Charm++-land.
*/
// *****************************************************************************

#include <unordered_map>
#include <vector>
#include <iostream>

#include "Types.hpp"
#include "Init.hpp"
#include "QuinoaConfig.hpp"
#include "Timer.hpp"
#include "Exception.hpp"
#include "CGPDE.hpp"
#include "DGPDE.hpp"
#include "PDEStack.hpp"
#include "ProcessException.hpp"
#include "InciterPrint.hpp"
#include "InciterDriver.hpp"
#include "Inciter/CmdLine/Parser.hpp"
#include "Inciter/CmdLine/CmdLine.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "ChareStateCollector.hpp"
#include "LBSwitch.hpp"

#include "NoWarning/inciter.decl.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wmissing-variable-declarations"
#endif

//! \brief Charm handle to the main proxy, facilitates call-back to finalize,
//!    etc., must be in global scope, unique per executable
CProxy_Main mainProxy;

//! Chare state collector Charm++ chare group proxy
tk::CProxy_ChareStateCollector stateProxy;

//! Load balancer switch group proxy
tk::CProxy_LBSwitch LBSwitchProxy;

//! If true, call and stack traces are to be output with exceptions
//! \note This is true by default so that the trace is always output between
//!   program start and the Main ctor in which the user-input from command line
//!   setting for this overrides this true setting.
bool g_trace = true;

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

//! Inciter declarations and definitions
namespace inciter {

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
//! Partial differential equations using continuous Galerkin selected by user
//! \details This vector is in global scope, because it holds polymorphic
//!   objects, and thus must be distributed to all PEs during initialization.
//!   Once distributed by the runtime system, the objects do not change.
std::vector< CGPDE > g_cgpde;
//! Partial differential equations using discontinuous Galerkin selected by user
//! \details This vector is in global scope, because it holds polymorphic
//!   objects, and thus must be distributed to all PEs during initialization.
//!   Once distributed by the runtime system, the objects do not change.
std::vector< DGPDE > g_dgpde;

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

//! \brief Pack/Unpack selected partial differential equations using continuous
//!   Galerkin discretization.
//! \details This Pack/Unpack method (re-)creates the PDE factory since it needs
//!   to (re-)bind function pointers on different processing elements. Therefore
//!   we circumvent Charm's usual pack/unpack for this type, and thus sizing
//!   does not make sense: sizing is a no-op. We could initialize the factory in
//!   InciterDriver's constructor and let this function re-create the stack only
//!   when unpacking, but that leads to repeating the same code twice: once in
//!   InciterDriver's constructor, once here. Another option is to use this
//!   pack/unpack routine to both initially create (when packing) and to
//!   re-create (when unpacking) the factory, which eliminates the need for
//!   pre-creating the object in InciterDriver's constructor and therefore
//!   eliminates the repeated code. This explains the guard for sizing: the code
//!   below is called for packing only (in serial) and packing and unpacking (in
//!   parallel).
inline
void operator|( PUP::er& p, std::vector< CGPDE >& eqs ) {
  try {
    if (!p.isSizing()) eqs = PDEStack().selectedCG();
  } catch (...) { tk::processExceptionCharm(); }
}

//! \brief Pack/Unpack selected partial differential equations using
//!   discontinuous Galerkin discretization.
//! \details This Pack/Unpack method (re-)creates the PDE factory since it needs
//!   to (re-)bind function pointers on different processing elements. Therefore
//!   we circumvent Charm's usual pack/unpack for this type, and thus sizing
//!   does not make sense: sizing is a no-op. We could initialize the factory in
//!   InciterDriver's constructor and let this function re-create the stack only
//!   when unpacking, but that leads to repeating the same code twice: once in
//!   InciterDriver's constructor, once here. Another option is to use this
//!   pack/unpack routine to both initially create (when packing) and to
//!   re-create (when unpacking) the factory, which eliminates the need for
//!   pre-creating the object in InciterDriver's constructor and therefore
//!   eliminates the repeated code. This explains the guard for sizing: the code
//!   below is called for packing only (in serial) and packing and unpacking (in
//!   parallel).
inline
void operator|( PUP::er& p, std::vector< DGPDE >& eqs ) {
  try {
    if (!p.isSizing()) eqs = PDEStack().selectedDG();
  } catch (...) { tk::processExceptionCharm(); }
}

} // inciter::

//! \brief Charm++ main chare for the shock hydrodynamics executable, inciter.
//! \details In inciter the Charm++ runtime system is initialized only after the
//!   mesh has been read in, partitioned, and the necessary data structures,
//!   e.g., communication maps, have been generated. This delayed initialization
//!   of the Charm++ runtime system is required since the mesh partitioning is
//!   done by Zoltan, an MPI library. Note that this Charm++ main chare object
//!   should not be in a namespace.
// cppcheck-suppress noConstructor
class Main : public CBase_Main {

  public:
    //! \brief Constructor
    //! \details Inciter's main chare constructor is the entry point of the
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
      m_signal( tk::setSignalHandlers() ),
      m_cmdline(),
      // Parse command line into m_cmdline using default simple pretty printer
      m_cmdParser( msg->argc, msg->argv, tk::Print(), m_cmdline ),
      // Create pretty printer initializing output streams based on command line
      m_print( m_cmdline.get< tag::verbose >() ? std::cout : std::clog ),
      // Create Inciter driver
      m_driver( tk::Main< inciter::InciterDriver >
                        ( msg->argc, msg->argv,
                          m_cmdline,
                          tk::HeaderType::INCITER,
                          tk::inciter_executable(),
                          m_print ) ),
      // Start new timer measuring the total runtime
      m_timer(1),
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

    //! Migrate constructor: returning from a checkpoint
    explicit Main( CkMigrateMessage* msg ) : CBase_Main( msg ),
      m_signal( tk::setSignalHandlers() ),
      m_cmdline(),
      m_cmdParser( reinterpret_cast<CkArgMsg*>(msg)->argc,
                   reinterpret_cast<CkArgMsg*>(msg)->argv,
                   tk::Print(),
                   m_cmdline ),
      m_print( m_cmdline.get< tag::verbose >() ? std::cout : std::clog ),
      m_driver( tk::Main< inciter::InciterDriver >
                        ( reinterpret_cast<CkArgMsg*>(msg)->argc,
                          reinterpret_cast<CkArgMsg*>(msg)->argv,
                          m_cmdline,
                          tk::HeaderType::INCITER,
                          tk::inciter_executable(),
                          m_print ) ),
      m_timer(1),
      m_timestamp()
    {
      g_trace = m_cmdline.get< tag::trace >();
      tk::MainCtor( mainProxy, thisProxy, m_timer, m_cmdline,
                    CkCallback( CkIndex_Main::quiescence(), thisProxy ) );
    }

    //! Execute driver created and initialized by constructor
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

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \note This is a Charm++ mainchare, pup() is thus only for
    //!    checkpoint/restart.
    void pup( PUP::er &p ) override {
      p | m_timer;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] m Mainchare object reference
    friend void operator|( PUP::er& p, Main& m ) { m.pup(p); }
    //@}

  private:
    int m_signal;                               //!< Used to set signal handlers
    inciter::ctr::CmdLine m_cmdline;            //!< Command line
    inciter::CmdLineParser m_cmdParser;         //!< Command line parser
    inciter::InciterPrint m_print;              //!< Pretty printer
    inciter::InciterDriver m_driver;            //!< Driver
    std::vector< tk::Timer > m_timer;           //!< Timers
    //! Time stamps in h:m:s with labels
    std::vector< std::pair< std::string, tk::Timer::Watch > > m_timestamp;
};

//! \brief Charm++ chare execute
//! \details By the time this object is constructed, the Charm++ runtime system
//!    has finished migrating all global-scoped read-only objects which happens
//!    after the main chare constructor has finished.
class execute : public CBase_execute {
  public:
    //! Constructor
    execute() { mainProxy.execute(); }
    //! Migrate constructor
    explicit execute( CkMigrateMessage* m ) : CBase_execute( m ) {}
};

#include "NoWarning/inciter.def.h"
