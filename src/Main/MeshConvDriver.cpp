// *****************************************************************************
/*!
  \file      src/Main/MeshConvDriver.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Mesh converter driver
  \details   Mesh converter driver.
*/
// *****************************************************************************

#include "Types.hpp"
#include "Tags.hpp"
#include "MeshConvDriver.hpp"
#include "MeshFactory.hpp"
#include "TaggedTupleDeepPrint.hpp"
#include "Writer.hpp"

#include "NoWarning/meshconv.decl.h"

using meshconv::MeshConvDriver;

extern CProxy_Main mainProxy;

MeshConvDriver::MeshConvDriver( const ctr::CmdLine& cmdline ) :
  m_print( tk::meshconv_executable() + "_screen.log",
           cmdline.get< tag::verbose >() ? std::cout : std::clog,
           std::ios_base::app ),
  m_reorder( cmdline.get< tag::reorder >() ),
  m_input(),
  m_output()
// *****************************************************************************
//  Constructor
//! \param[in] print Pretty printer
//! \param[in] cmdline Command line object storing data parsed from the command
//!   line arguments
// *****************************************************************************
{
  // Save input file name
  m_input = cmdline.get< tag::io, tag::input >();
  // Save output file name
  m_output = cmdline.get< tag::io, tag::output >();


  // Output command line object to file
  auto logfilename = tk::meshconv_executable() + "_input.log";
  tk::Writer log( logfilename );
  tk::print( log.stream(), "cmdline", cmdline );
}

void
MeshConvDriver::execute() const
// *****************************************************************************
//  Execute: Convert mesh file
// *****************************************************************************
{
  m_print.endsubsection();

  std::vector< std::pair< std::string, tk::real > > times( 1 );

  auto mesh = tk::readUnsMesh( m_print, m_input, times[0] );
  auto wtimes = tk::writeUnsMesh( m_print, m_output, mesh, m_reorder );

  times.insert( end(times), begin(wtimes), end(wtimes) );
  mainProxy.timestamp( times );

  mainProxy.finalize();
}
