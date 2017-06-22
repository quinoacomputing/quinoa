// *****************************************************************************
/*!
  \file      src/Main/FileDiffDriver.C
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
#include "FileDiffDriver.h"
#include "DiffWriterFiles.h"

#include "NoWarning/filediff.decl.h"

using filediff::FileDiffDriver;

extern CProxy_Main mainProxy;

FileDiffDriver::FileDiffDriver( const tk::Print& print,
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
FileDiffDriver::execute() const
// *****************************************************************************
//  Execute: Compute Differences
//! \author A. Pakki
// *****************************************************************************
{

  std::vector< std::pair< std::string, tk::real > > times( 1 );

  tk::DiffWriterFiles *dwf = new tk::DiffWriterFiles( m_input, m_output );

  dwf->convertFiles();

  mainProxy.timestamp( times );

  delete dwf;

  mainProxy.finalize();
}
