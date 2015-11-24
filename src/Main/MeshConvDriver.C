//******************************************************************************
/*!
  \file      src/Main/MeshConvDriver.C
  \author    J. Bakosi
  \date      Tue 24 Nov 2015 11:33:32 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Mesh converter driver
  \details   Mesh converter driver.
*/
//******************************************************************************

#include <utility>

#include "Types.h"
#include "Tags.h"
#include "MeshConvDriver.h"
#include "MeshFactory.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "meshconv.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

using meshconv::MeshConvDriver;

extern CProxy_Main mainProxy;

MeshConvDriver::MeshConvDriver( const tk::Print& print,
                                const ctr::CmdLine& cmdline )
  : m_print( print ),
    m_reorder( cmdline.get< tag::reorder >() )
//******************************************************************************
//  Constructor
//! \param[in] print Pretty printer
//! \param[in] cmdline Command line object storing data parsed from the command
//!   line arguments
//! \author J. Bakosi
//******************************************************************************
{
  // Save input file name
  m_input = cmdline.get< tag::io, tag::input >();
  // Save output file name
  m_output = cmdline.get< tag::io, tag::output >();
}

void
MeshConvDriver::execute() const
//******************************************************************************
//  Execute: Convert mesh file
//! \author J. Bakosi
//******************************************************************************
{
  m_print.endsubsection();

  std::vector< std::pair< std::string, tk::real > > times( 1 );

  auto wtimes = tk::writeUnsMesh( m_print,
                                  m_output,
                                  tk::readUnsMesh( m_print, m_input, times[0] ),
                                  m_reorder );

  times.insert( end(times), begin(wtimes), end(wtimes) );
  mainProxy.timestamp( times );

  mainProxy.finalize();
}
