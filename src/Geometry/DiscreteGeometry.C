//******************************************************************************
/*!
  \file      src/Geometry/DiscreteGeometry.C
  \author    J. Bakosi
  \date      Mon 11 Nov 2013 08:59:15 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
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
    reader( base.control.get< ctr::cmd, ctr::io, ctr::input >(), m_mesh );

  // Read in STL mesh
  reader.read();

  // Instantiate Silo writer
  SiloWriter writer( base.control.get< ctr::cmd, ctr::io, ctr::output >(),
                     m_mesh,
                     DB_ALL_AND_DRVR);

  // Write out STL geometry to Silo file
  writer.write();
}
