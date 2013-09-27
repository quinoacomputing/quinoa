//******************************************************************************
/*!
  \file      src/Geometry/DiscreteGeometry.C
  \author    J. Bakosi
  \date      Fri Sep 27 11:07:53 2013
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
  Geometry(base),
  m_mesh(new STLMesh)
//******************************************************************************
//  Constructor
//! \param[in] base      Essentials
//! \details Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  using namespace ctr;

  // Instantiate ASCII STL mesh reader object
  STLTxtMeshReader reader(base.control.get<io,input>(), m_mesh.get());

  // Read in STL mesh
  reader.read();

  // Instantiate Silo writer object
  SiloWriter writer(base.control.get<io,output>(),
                    m_mesh.get(),
                    DB_ALL_AND_DRVR);

  // Write out STL geometry to Silo file
  writer.write();
}
