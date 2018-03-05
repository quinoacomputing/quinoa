// *****************************************************************************
/*!
  \file      src/Main/RNGTest.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     RNGTest's random number generator test suite's Charm++ main chare.
  \details   RNGTest's random number generator test suite's Charm++ main chare.
    This file contains the definition of the Charm++ main chare, equivalent to
    main() in Charm++-land.
*/
// *****************************************************************************

#include <map>
#include <vector>
#include <iostream>
#include <utility>

#include "NoWarning/format.h"

#include "QuinoaConfig.h"
#include "Print.h"
#include "Init.h"
#include "Timer.h"
#include "Tags.h"
#include "ProcessException.h"
#include "RNG.h"
#include "RNGStack.h"
#include "TestStack.h"
#include "Options/RNG.h"
#include "RNGTestPrint.h"
#include "RNGTestDriver.h"
#include "RNGTest/CmdLine/CmdLine.h"
#include "RNGTest/CmdLine/Parser.h"
#include "RNGTest/InputDeck/InputDeck.h"

#include "NoWarning/charm.h"
#include "NoWarning/pup.h"
#include "NoWarning/rngtest.decl.h"

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

//! RNGTest declarations and definitions
namespace rngtest {

//! Global-scope data. Initialized by the main chare and distibuted to all PEs by
//! the Charm++ runtime system. Though semantically not const, all these global
//! data should be considered read-only. See also http://charm.cs.illinois.edu/
//! manuals/html/charm++/manual.html. The data below is global-scope because they
//! must be available to all PEs which could be on different machines. In a
//! previous non-Charm++ design, most of this data was held at class-level, but
//! since the generators in g_rng must be possible to be called from
//! global-scope, as external generators to TestU01, it is easier to make g_rng
//! global-scope, as well the additional data required to initialize it,
//! contained in g_inputdeck (storing all parsed user input).  This is required
//! for serializing during migration of g_rng across the network.
//!
//! Note that the container (std::map) holding tk::RNG objects uses value
//! semantics which is safer and generally less error-prone than reference
//! semantics. At the same time tk::RNG is used in a polymorphic fashion with
//! various classes that adhere to the concepts required by Concept defined
//! inside tk::RNG. tk::RNG does not define a default, i.e., non-templated
//! constructor, since then the "derived" class object could not be initialized
//! rendering the class tk::RNG empty-constructed, which invites abuse and
//! ill-defined behavior. As such, the "derived" class type comes through the
//! constructors and thus would not be available for a pack/unpack migrator
//! required by Charm++ from within. Templating the class tk::RNG is not an
//! option since then we could not hold tk::RNG objects in a simple std::vector.
//! As a result of the above requirements, the tk::RNG objects in g_rng are
//! migrated (here in global-scope) by reinstantiating RNGStack, which
//! reinstatiates the RNG factory, from which the RNGs selected by the user are
//! instantiated.
//!
//! Note also that RNGFactory associates tk::ctr::RNG ids (enum class values) to
//! function pointers (std::function objects pointing to tk::RNG constructors
//! bound with their arguments). Since function pointers cannot simply be
//! serialized and migrated via the network, they must also be recreated on
//! remote machines. This initial migration of global-scope data is done by the
//! Charm++ runtime once the main chare constructor is finished -- see the
//! RNGTestDriver constructor, which initializes the data required for the
//! migration).

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
//!   Once distributed by the runtime system, the objects do not change and
//!   available to all PEs.
std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;
//! Statistical test wrappers, grouped by test libraries
//! \details This object is in global scope, because it holds function pointers
//!   to test runners, and thus must be distributed to all PEs during
//!   initialization. Once distributed by the runtime system, the objects do not
//!   change.
TestStack g_testStack;

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
}

//! Pack/Unpack test stack. This Pack/Unpack method (re-)creates the full test
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
inline void operator|( PUP::er& p, TestStack& stack )
{ if (!p.isSizing()) stack = TestStack(); }

} // rngtest::

//! \brief Charm++ main chare for the random number generator test suite
//!    executable, rngtest.
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
    //!   the random number generator test suite is parallel and asynchronous,
    //!   its driver fires up additional Charm++ chare objects which then call
    //!   back to Main::finalize() at some point in the future when all work has
    //!   been finished. finalize() then exits by calling Charm++'s CkExit(),
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
      // Create RNGTest driver
      m_driver( tk::Main< rngtest::RNGTestDriver >
                        ( msg->argc, msg->argv,
                          m_cmdline,
                          tk::HeaderType::RNGTEST,
                          RNGTEST_EXECUTABLE,
                          m_print ) ),
      m_timer(1),       // Start new timer measuring the total runtime
      m_timestamp()
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
    } catch (...) { tk::processExceptionCharm(); }

    void execute() {
      try {
        m_timestamp.emplace_back("Migrate global-scope data", m_timer[1].hms());
        m_driver.execute();       // fires up async chares
      } catch (...) { tk::processExceptionCharm(); }   
    }

    void finalize() {
      try {
        if (!m_timer.empty()) {
          m_timestamp.emplace_back( "Total runtime", m_timer[0].hms() );
          m_print.time( "Timers (h:m:s)", m_timestamp );
          m_print.endpart();
        }
      } catch (...) { tk::processExceptionCharm(); }
      // Tell the Charm++ runtime system to exit
      CkExit();
    }

  private:
    int m_signal;                               //!< Used to set signal handlers
    rngtest::ctr::CmdLine m_cmdline;            //!< Command line
    rngtest::CmdLineParser m_cmdParser;         //!< Command line parser
    rngtest::RNGTestPrint m_print;              //!< Pretty printer
    rngtest::RNGTestDriver m_driver;            //!< Driver
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

#include "NoWarning/rngtest.def.h"
