//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 09:42:03 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Quinoa main
  \details   Quinoa main
*/
//******************************************************************************

#include <pup_stl.h>

#include <Config.h>
#include <QuinoaPrint.h>
#include <QuinoaDriver.h>
#include <Quinoa/CmdLine/Parser.h>
#include <TPLInfo/Silo.h>
#include <TPLInfo/HDF5.h>
#include <TPLInfo/Zlib.h>
#include <TPLInfo/MKL.h>
#include <TPLInfo/Boost.h>
#include <quinoa.decl.h>
#include <Init.h>

//! Charm handle to the main proxy, facilitates call-back to finalize, etc.,
//! must be in global scope, unique per executable
CProxy_Main mainProxy;

namespace quinoa {

void echoTPL( const tk::Print& print )
//******************************************************************************
//  Echo TPL version informaion for libs specific to Quinoa
//! \author  J. Bakosi
//******************************************************************************
{
  #ifdef HAS_MKL
  echoMKL( print, "Intel Math Kernel Library" );
  #else
  print.item( "Intel Math Kernel Library", "n/a" );
  #endif
  echoBoost( print, "Boost C++ Libraries" );
  tk::echoSilo(print, "Silo library");
  tk::echoHDF5(print, "HDF5 library");
  tk::echoZlib(print, "Zlib compression library");
  print.endpart();
}

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

} // quinoa::

//! Charm++ main chare
class Main : public CBase_Main {

  public:
    Main( CkArgMsg* msg ) :
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
                          m_print,
                          quinoa::echoTPL ) ),
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
      m_driver.execute();       // does not (yet) fire up async chares
      finalize();               // so call finalize right away
    }

    void finalize() {
      m_timestamp.emplace( "Total runtime", m_timer[0].hms() );
      m_print.time( "Timers (h:m:s)", m_timestamp );
      m_print.endpart();
      CkExit();
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
