// *****************************************************************************
/*!
  \file      src/Main/FileConvDriver.C
  \author    A. Pakki
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     File converter driver
  \details   File converter driver.
*/
// *****************************************************************************

#include <utility>
#include <iostream>

#include "Types.h"
#include "Tags.h"
#include "FileConvDriver.h"
#include "FileConvWriter.h"

#include "NoWarning/fileconv.decl.h"

using fileconv::FileConvDriver;

extern CProxy_Main mainProxy;

FileConvDriver::FileConvDriver( const tk::Print& print,
                                const ctr::CmdLine& cmdline )
  : m_print( print ),
    m_input(),
    m_output()
// *****************************************************************************
//  Constructor
//! \param[in] print Pretty printer
//! \param[in] cmdline Command line object storing data parsed from the command
//!   line arguments
//! \author A. Pakki
// *****************************************************************************
{
  // Save input file name
  m_input = cmdline.get< tag::io, tag::input >();
  // Save output file name
  m_output = cmdline.get< tag::io, tag::output >();
}

void
FileConvDriver::execute() const
// *****************************************************************************
//  Execute: Convert the file layout
//! \author A. Pakki
// *****************************************************************************
{

  std::vector< std::pair< std::string, tk::real > > times( 1 );

  tk::FileConvWriter *fcw = new tk::FileConvWriter( m_input, m_output );

  fcw->convertFiles();

  mainProxy.timestamp( times );

  delete fcw;

  mainProxy.finalize();
}
