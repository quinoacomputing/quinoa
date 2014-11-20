//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Wed 19 Nov 2014 04:56:42 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Quinoa main
  \details   Quinoa main
*/
//******************************************************************************

#include <pup_stl.h>

#include <Config.h>
#include <RNG.h>
#include <RNGStack.h>
#include <DiffEqStack.h>
#include <QuinoaPrint.h>
#include <QuinoaDriver.h>
#include <Quinoa/CmdLine/Parser.h>
#include <quinoa.decl.h>
#include <Init.h>

//! Charm handle to the main proxy, facilitates call-back to finalize, etc.,
//! must be in global scope, unique per executable
CProxy_Main mainProxy;

namespace quinoa {

//! Global-scope data. Initialized by the main chare and distibuted to all PEs
//! by the Charm++ runtime system. Though semantically not const, all these
//! global data should be considered read-only. See also
//! http://charm.cs.illinois.edu/manuals/html/charm++/manual.html. The data
//! below is global-scope because they must be available to all PEs which could
//! be on different machines.

//! Defaults of input deck, facilitates detection what is set by user
ctr::InputDeck g_inputdeck_defaults;
//! Input deck filled by parser, containing all input data
ctr::InputDeck g_inputdeck;
//! Random number generators selected by user
std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;
//! Differential equations selected by user
std::vector< DiffEq > g_diffeqs;

//! Distributor Charm++ proxy facilitating call-back to Distributor by the
//! individual integrators
CProxy_Distributor g_DistributorProxy;

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
      g_inputdeck.get< tag::param, tk::tag::rngmkl >(),
      #endif
      g_inputdeck.get< tag::param, tk::tag::rngsse >() );
    rng = stack.selected( g_inputdeck.get< tag::selected, tk::tag::rng >() );
  }
}

//! Pack/Unpack selected differential equations. This Pack/Unpack method
//! (re-)creates the DiffEq factory since it needs to (re-)bind function
//! pointers on different processing elements. Therefore we circumvent Charm's
//! usual pack/unpack for this type, and thus sizing does not make sense: sizing
//! is a no-op. We could initialize the factory in QuinoaDriver's constructor
//! and let this function re-create the stack only when unpacking, but that
//! leads to repeating the same code twice: once in QuinoaDriver's constructor,
//! once here. Another option is to use this pack/unpack routine to both
//! initially create (when packing) and to re-create (when unpacking) the
//! factory, which eliminates the need for pre-creating the object in
//! QuinoaDriver's constructor and therefore eliminates the repeated code. This
//! explains the guard for sizing: the code below is called for packing only (in
//! serial) and packing and unpacking (in parallel).
inline
void operator|( PUP::er& p, std::vector< DiffEq >& eqs ) {
  if (!p.isSizing()) eqs = DiffEqStack().selected();
}

} // quinoa::

//! Charm++ main chare
class Main : public CBase_Main {

  public:
    Main( CkArgMsg* msg )
    try :
      // Parse command line into m_cmdline using default simple pretty printer
      m_cmdParser( msg->argc, msg->argv, tk::Print(), m_cmdline ),
      // Create pretty printer initializing output streams based on command line
      m_print( m_cmdline.get< tk::tag::verbose >() ? std::cout : std::clog ),
      // Create Quinoa driver
      m_driver( tk::Main< quinoa::QuinoaDriver >
                        ( msg->argc, msg->argv,
                          m_cmdline,
                          tk::HeaderType::QUINOA,
                          QUINOA_EXECUTABLE,
                          m_print ) ),
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
    } catch (...) { processException(); }

    //! Execute driver created and initialized by constructor
    void execute() {
      m_timestamp.emplace("Migration of global-scope data", m_timer[1].hms());
      m_driver.execute();       // fires up async chares
    }

    //! Normal exit point
    void finalize() {
      if (!m_timer.empty()) {
        m_timestamp.emplace( "Total runtime", m_timer[0].hms() );
        m_print.time( "Timers (h:m:s)", m_timestamp );
        m_print.endpart();
      }
      CkExit();
    }

    //! Add time stamp contributing to final timers output
    void timestamp( std::string label, tk::real stamp )
    { m_timestamp.emplace( label, tk::hms( stamp ) ); }

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
    quinoa::ctr::CmdLine m_cmdline;                   //!< Command line
    quinoa::CmdLineParser m_cmdParser;                //!< Command line parser
    quinoa::QuinoaPrint m_print;                      //!< Pretty printer
    quinoa::QuinoaDriver m_driver;                    //!< Driver
    std::vector< tk::Timer > m_timer;                 //!< Timers

    //! Time stamps in h:m:s with labels
    std::map< std::string, tk::Timer::Watch > m_timestamp;
};

//! Charm++ chare execute: by the time this object is constructed, the Charm++
//! runtime system has finished migrating all global-scoped read-only objects
//! which happens after the main chare constructor has finished.
struct execute : CBase_execute { execute() { mainProxy.execute(); } };

#include <quinoa.def.h>
