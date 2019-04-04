// *****************************************************************************
/*!
  \file      src/Main/Init.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Common initialization routines for main() functions for multiple
     exectuables
  \details   Common initialization routines for main() functions for multiple
     exectuables. The functions in this file are used by multiple execitables
     to ensure code-reuse and a uniform screen-output.
*/
// *****************************************************************************
#ifndef Init_h
#define Init_h

#include <string>
#include <unordered_map>

#include "NoWarning/charm++.h"

#include "QuinoaConfig.h"
#include "Exception.h"
#include "Print.h"
#include "ChareStateCollector.h"
#include "ProcessException.h"

namespace tk {

//! Executable types for which an ascii logo is available in tk::Print
enum class HeaderType : uint8_t { INCITER=0,
                                  RNGTEST,
                                  UNITTEST,
                                  MESHCONV,
                                  FILECONV,
                                  WALKER };

//! Wrapper for the standard C library's gettimeofday() from
std::string curtime();

//! Echo program header
void echoHeader( const Print& print, HeaderType header );

//! Echo build environment
void echoBuildEnv( const Print& print, const std::string& executable );

//! Echo runtime environment
void echoRunEnv( const Print& print, int argc, char** argv,
                 bool verbose, bool quiescence, bool charestate, bool trace );

//! \brief Generic Main() used for all executables for code-reuse and a uniform
//!    output
//! \details The template arguments configure this Main class that is
//!   practically used instead of the usual main(). This allows code-reuse and a
//!   unfirom screen-output. The template arguments are:
//!   - Driver, specializaing the driver type to be created, see tk::Driver
//!   - Printer, specializaing the pretty printer type to use, see tk::Print
//!   - CmdLine, specializing the command line object storing data parsed from
//!     the command line
//! \param[in] argc Number of command-line arguments to executable
//! \param[in] argv C-style string array to command-line arguments to executable
//! \param[in] cmdline Command line object storing data parsed from the command
//!   line arguments
//! \param[in] header Header type enum indicating which executable header to
//!   print
//! \param[in] executable Name of the executable
//! \param[in] print Pretty printer to use
//! \return Instantiated driver object which can then be used to execute()
//!   whatever it is intended to drive
template< class Driver, class Printer, class CmdLine >
Driver Main( int argc, char* argv[],
             const CmdLine& cmdline,
             HeaderType header,
             const std::string& executable,
             const Printer& print )
{
  // Echo program header
  echoHeader( print, header );

  // Echo environment
  print.part( "Environment" );
  // Build environment
  echoBuildEnv( print, executable );
  // Runtime environment
  echoRunEnv( print, argc, argv, cmdline.template get< tag::verbose >(),
              cmdline.template get< tag::quiescence >(),
              cmdline.template get< tag::chare >(),
              cmdline.template get< tag::trace >() );

  // Create and return driver
  return Driver( print, cmdline );
}

//! Generic Main Charm++ module constructor for all executables
//! \tparam ExecuteProxy Charm++ proxy type for the 'excecute' chare, see
//!    src/Main/\<executable\>.C
//! \tparam MainProxy Main Charm++ chare proxy for the executable
//! \tparam CmdLine Executable-specific tagged tuple storing the rusult of the
//!    command line parser
//! \param[in] msg Charm++ CkArgMsg pointer passed (by Charm++) to the main
//!   chare proxy
//! \param[in,out] mp MainProxy to set for the main chare
//! \param[in] thisProxy 'thisProxy' to set as MainProxy
//! \param[in,out] state Chare state collector proxy
//! \param[in,out] timer Vector of timers, held by the main chare, in which to
//!   start the first timer, measuring the migration of global-scope data
//! \param[in] cmdline Command line grammar stack for the executable (assumed
//!   already parsed)
//! \param[in] quiescenceTarget Pre-created Charm++ callback to use as the
//!   target function to call if quiescence is detected
template< class ExecuteProxy, class MainProxy, class CmdLine >
void MainCtor( CkArgMsg* msg,
               MainProxy& mp,
               const MainProxy& thisProxy,
               tk::CProxy_ChareStateCollector& state,
               std::vector< tk::Timer >& timer,
               const CmdLine& cmdline,
               const CkCallback& quiescenceTarget )
{
  delete msg;

  // Set Charm++ main proxy
  mp = thisProxy;

  // If quiescence detection is on or user requested it, create chare state
  // collector Charm++ chare group
  if ( cmdline.template get< tag::chare >() ||
       cmdline.template get< tag::quiescence >() )
    state = tk::CProxy_ChareStateCollector::ckNew();

  // Optionally enable quiscence detection
  if (cmdline.template get< tag::quiescence >()) CkStartQD( quiescenceTarget );

  // Fire up an asynchronous execute object, which when created at some
  // future point in time will call back to this->execute(). This is
  // necessary so that this->execute() can access already migrated
  // global-scope data.
  ExecuteProxy::ckNew();

  // Start new timer measuring the migration of global-scope data
  timer.emplace_back();
}

//! Generic function to dump the Charm++ chare state (if collected)
//! \tparam CmdLine Executable-specific tagged tuple storing the rusult of the
//!    command line parser
//! \param[in] cmdline Command line grammar stack for the executable
//! \param[in] print Pretty printer
//! \param[in] msg Charm++ reduction message containing the chare state
//!   aggregated from all PEs
template< class CmdLine >
void dumpstate( const CmdLine& cmdline,
                const tk::Print& print,
                CkReductionMsg* msg )
{
  try {

    // unpack chare state
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
    if (cmdline.template get< tag::chare >() || error)
      print.charestate( state );

    // exit differently depending on how we were called
    if (error)
      Throw( "Quiescence detected" );
    else
      CkExit(); // tell the Charm++ runtime system to exit with zero exit code

  } catch (...) { tk::processExceptionCharm(); }
}

//! Generic finalize function for different executables
//! \param[in] cmdline Command line grammar stack for the executable
//! \param[in] timer Vector of timers, held by the main chare
//! \param[in] print Pretty printer
//! \param[in,out] state Chare state collector proxy
//! \param[in,out] timestamp Vector of time stamps in h:m:s with labels
//! \param[in] dumpstateTarget Pre-created Charm++ callback to use as the
//!   target function for dumping chare state
//! \param[in] clean True if we should exit with a zero exit code, false to
//!   exit with a nonzero exit code
template< class CmdLine >
void finalize( const CmdLine& cmdline,
               const std::vector< tk::Timer >& timer,
               const tk::Print& print,
               tk::CProxy_ChareStateCollector& state,
               std::vector< std::pair< std::string,
                                       tk::Timer::Watch > >& timestamp,
               const CkCallback& dumpstateTarget,
               bool clean = true )
{
  try {

   if (!timer.empty()) {
     timestamp.emplace_back( "Total runtime", timer[0].hms() );
     print.time( "Timers (h:m:s)", timestamp );
     print.endpart();
     // if quiescence detection is on or user requested it, collect chare state
     if ( cmdline.template get< tag::chare >() ||
          cmdline.template get< tag::quiescence >() ) {
       state.collect( /* error = */ false, dumpstateTarget );
     }
     // tell the Charm++ runtime system to exit with zero exit code
     if (clean) CkExit(); else CkAbort("Failed");
   }

  } catch (...) { tk::processExceptionCharm(); }
}

} // tk::

#endif // Init_h
