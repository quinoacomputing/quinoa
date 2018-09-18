// *****************************************************************************
/*!
  \file      src/Main/Init.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
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

#include "Print.h"

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
                 bool verbose, bool quiescence, bool charestate );

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
              cmdline.template get< tag::chare >() );

  // Create and return driver
  return Driver( print, cmdline );
}

} // tk::

#endif // Init_h
