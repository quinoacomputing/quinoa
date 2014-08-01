//******************************************************************************
/*!
  \file      src/Main/RNGTest.C
  \author    J. Bakosi
  \date      Fri 01 Aug 2014 11:41:04 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     RNGTest: Quinoa's random number generator test suite
  \details   RNGTest: Quinoa's random number generator test suite
*/
//******************************************************************************

#include <Config.h>
#include <Paradigm.h>
#include <RNG.h>
#include <RNGStack.h>
#include <RNGTestPrint.h>
#include <RNGTestDriver.h>
#include <RNGTest/CmdLine/Parser.h>
#include <RNGTest/InputDeck/InputDeck.h>
#include <TestStack.h>
#include <PUPUtil.h>
#include <TPLInfo/MKL.h>
#include <TPLInfo/Boost.h>
#include <TPLInfo/OpenMP.h>
#include <rngtest.decl.h>
#include <Init.h>

//! Charm handle to the main proxy, facilitates call-back to finalize, etc.,
//! must be in global scope, unique per executable
CProxy_Main mainProxy;

namespace rngtest {

void echoTPL( const tk::Print& print )
//******************************************************************************
//  Echo TPL version informaion
//! \author  J. Bakosi
//******************************************************************************
{
  echoOpenMP( print, "OpenMP runtime" );
#ifdef HAS_MKL
  echoMKL( print, "Intel Math Kernel Library" );
#else
  print.item( "Intel Math Kernel Library", "n/a" );
#endif
  echoBoost( print, "Boost C++ Libraries" );
  print.endpart();
}

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
//! remote machines.  This initial migration of global-scope data is done by the
//! Charm++ runtime once the main chare constructor is finished -- see the
//! RNGTestDriver constructor, which initializes the data required for the
//! migration).

//! Defaults of input deck, facilitates detection what is set by user
ctr::InputDeck g_inputdeck_defaults;
//! Input deck filled by parser, containing all input data
ctr::InputDeck g_inputdeck;
//! Random number generators selected by user
std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;
//! Statistical test wrappers, grouped by test libraries
TestStack g_testStack;

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
    tk::RNGFactory factory;
    tk::RNGStack stack;
    stack.initFactory( factory, CkNumPes(),
                       #ifdef HAS_MKL
                       g_inputdeck.get< tag::param, tk::tag::rngmkl >(),
                       #endif
                       g_inputdeck.get< tag::param, tk::tag::rngsse >() );
    rng = stack.createSelected(
            factory, g_inputdeck.get< tag::selected, tk::tag::rng >() );
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

//! Charm++ main chare
class Main : public CBase_Main {

  public:
    Main( CkArgMsg* msg ) :
      // Parse command line into m_cmdline using default simple pretty printer
      m_cmdParser( msg->argc, msg->argv, tk::Print(), m_cmdline ),
      // Create pretty printer initializing output streams based on command line
      m_print( m_cmdline.get< tk::tag::verbose >() ? std::cout : std::clog ),
      // Create RNGTest driver
      m_driver( tk::Main< rngtest::RNGTestDriver >
                        ( msg->argc, msg->argv,
                          m_cmdline,
                          tk::HeaderType::RNGTEST,
                          RNGTEST_EXECUTABLE,
                          m_print,
                          rngtest::echoTPL ) ),
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
    rngtest::ctr::CmdLine m_cmdline;                    //!< Command line
    rngtest::CmdLineParser m_cmdParser;                 //!< Command line parser
    rngtest::RNGTestPrint m_print;                      //!< Pretty printer
    rngtest::RNGTestDriver m_driver;                    //!< Driver
    std::vector< tk::Timer > m_timer;                   //!< Timers

    //! Time stamps in h:m:s with labels
    std::map< std::string, tk::Timer::Watch > m_timestamp;
};

//! Charm++ chare execute: by the time this object is constructed, the Charm++
//! runtime system has finished migrating all global-scoped read-only objects
//! which happens after the main chare constructor has finished.
struct execute : CBase_execute { execute() { mainProxy.execute(); } };

#include <rngtest.def.h>
