//******************************************************************************
/*!
  \file      src/Main/MeshConvDriver.C
  \author    J. Bakosi
  \date      Sat 04 Apr 2015 07:42:08 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Mesh converter driver
  \details   Mesh converter driver.
*/
//******************************************************************************

#include <Print.h>
#include <Exception.h>
#include <Factory.h>
#include <MeshConvDriver.h>
#include <MeshConv/CmdLine/Parser.h>
#include <MeshFactory.h>

using meshconv::MeshConvDriver;

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
  tk::writeUnsMesh( m_output, tk::readUnsMesh( m_input ) );
}
