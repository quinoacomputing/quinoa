//******************************************************************************
/*!
  \file      src/Main/MeshConvDriver.C
  \author    J. Bakosi
  \date      Sun 31 May 2015 06:26:20 AM MDT
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
  std::pair< std::string, tk::real > rtime, wtime;
  tk::writeUnsMesh( m_output, tk::readUnsMesh(m_input,rtime), wtime );
  mainProxy.timestamp( { rtime, wtime } );
  mainProxy.finalize();
}
