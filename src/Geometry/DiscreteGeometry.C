//******************************************************************************
/*!
  \file      src/Geometry/DiscreteGeometry.C
  \author    J. Bakosi
  \date      Thu 16 Jan 2014 10:02:40 PM MST
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Discrete geometry definition
  \details   Discrete geometry definition
*/
//******************************************************************************

#include <DiscreteGeometry.h>
#include <STLTxtMeshReader.h>
#include <SiloWriter.h>

using quinoa::DiscreteGeometry;

DiscreteGeometry::DiscreteGeometry(const Base& base) : Geometry(base)
//******************************************************************************
//  Constructor
//! \param[in] base      Essentials
//! \details Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  // Instantiate ASCII STL mesh reader
  STLTxtMeshReader
    reader( base.control.get< tag::cmd, tag::io, tag::input >(), m_mesh );

  // Read in STL mesh
  reader.read();

  // Instantiate Silo writer
  SiloWriter writer( base.control.get< tag::cmd, tag::io, tag::output >(),
                     m_mesh,
                     DB_ALL_AND_DRVR);

  // Write out STL geometry to Silo file
  writer.write();
}
