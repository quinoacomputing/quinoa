//******************************************************************************
/*!
  \file      src/Geometry/DiscreteGeometry.C
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 05:10:09 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Discrete geometry definition
  \details   Discrete geometry definition
*/
//******************************************************************************

#include <DiscreteGeometry.h>
#include <STLTxtMeshReader.h>
#include <SiloWriter.h>

using quinoa::DiscreteGeometry;

DiscreteGeometry::DiscreteGeometry()
//******************************************************************************
//  Constructor
//! \details Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
//   // Instantiate ASCII STL mesh reader
//   STLTxtMeshReader
//     reader( base.control.get< tag::cmd, tag::io, tag::input >(), m_mesh );
// 
//   // Read in STL mesh
//   reader.read();
// 
//   // Instantiate Silo writer
//   SiloWriter writer( base.control.get< tag::cmd, tag::io, tag::output >(),
//                      m_mesh,
//                      DB_ALL_AND_DRVR);
// 
//   // Write out STL geometry to Silo file
//   writer.write();
}
