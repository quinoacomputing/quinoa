// *****************************************************************************
/*!
  \file      src/Main/FileConvDriver.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     File converter driver
  \details   File converter driver.
*/
// *****************************************************************************

#include <utility>
#include <memory>

#include "Types.hpp"
#include "Tags.hpp"
#include "FileConvDriver.hpp"
#include "FileConvWriter.hpp"

#include "NoWarning/fileconv.decl.h"

using fileconv::FileConvDriver;

extern CProxy_Main mainProxy;

FileConvDriver::FileConvDriver( const tk::Print&,
                                const ctr::CmdLine& cmdline )
  : m_input(), m_output()
// *****************************************************************************
//  Constructor
//! \param[in] cmdline Command line object storing data parsed from the command
//!   line arguments
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
// *****************************************************************************
{

  std::vector< std::pair< std::string, tk::real > > times( 1 );

  std::unique_ptr< tk::FileConvWriter > fcw(new tk::FileConvWriter
					  ( m_input, m_output ) );
  fcw->convertFiles();

  mainProxy.timestamp( times );
  mainProxy.finalize();

}
