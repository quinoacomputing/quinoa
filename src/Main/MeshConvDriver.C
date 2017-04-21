// *****************************************************************************
/*!
  \file      src/Main/MeshConvDriver.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Mesh converter driver
  \details   Mesh converter driver.
*/
// *****************************************************************************

#include <utility>

#include "Types.h"
#include "Tags.h"
#include "MeshConvDriver.h"
#include "MeshFactory.h"

#include "NoWarning/meshconv.decl.h"

using meshconv::MeshConvDriver;

extern CProxy_Main mainProxy;

MeshConvDriver::MeshConvDriver( const tk::Print& print,
                                const ctr::CmdLine& cmdline )
  : m_print( print ),
    m_reorder( cmdline.get< tag::reorder >() ),
    m_input(),
    m_output()
// *****************************************************************************
//  Constructor
//! \param[in] print Pretty printer
//! \param[in] cmdline Command line object storing data parsed from the command
//!   line arguments
//! \author J. Bakosi
// *****************************************************************************
{
  // Save input file name
  m_input = cmdline.get< tag::io, tag::input >();
  // Save output file name
  m_output = cmdline.get< tag::io, tag::output >();
}

void
MeshConvDriver::execute() const
// *****************************************************************************
//  Execute: Convert mesh file
//! \author J. Bakosi
// *****************************************************************************
{
  m_print.endsubsection();

  std::vector< std::pair< std::string, tk::real > > times( 1 );

  auto mesh = tk::readUnsMesh( m_print, m_input, times[0] );
  auto wtimes = tk::writeUnsMesh( m_print,
                                  m_output,
                                  mesh,
                                  m_reorder );

  times.insert( end(times), begin(wtimes), end(wtimes) );
  mainProxy.timestamp( times );

  mainProxy.finalize();
}
