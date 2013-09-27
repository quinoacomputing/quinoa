//******************************************************************************
/*!
  \file      src/Geometry/DiscreteGeometry.C
  \author    J. Bakosi
  \date      Fri Sep 27 15:28:23 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Discrete geometry definition
  \details   Discrete geometry definition
*/
//******************************************************************************

#include <DiscreteGeometry.h>
#include <Exception.h>
#include <QuinoaControl.h>
#include <STLTxtMeshReader.h>
#include <SiloWriter.h>

using namespace quinoa;

DiscreteGeometry::DiscreteGeometry(const Base& base) :
  Geometry(base)
//******************************************************************************
//  Constructor
//! \param[in] base      Essentials
//! \details Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  // Instantiate ASCII STL mesh reader object
  STLTxtMeshReader reader(base.control.get<ctr::io,ctr::input>(), m_mesh);

  // Read in STL mesh
  reader.read();

  // Instantiate Silo writer object
  SiloWriter writer(base.control.get<ctr::io,ctr::geomoutput>(),
                    m_mesh,
                    DB_ALL_AND_DRVR);

  // Write out STL geometry to Silo file
  writer.write();
}
