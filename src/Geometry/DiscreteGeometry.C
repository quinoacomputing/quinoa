//******************************************************************************
/*!
  \file      src/Geometry/DiscreteGeometry.C
  \author    J. Bakosi
  \date      Tue 17 Sep 2013 10:45:12 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Discrete geometry definition
  \details   Discrete geometry definition
*/
//******************************************************************************

#include <DiscreteGeometry.h>
#include <Exception.h>
#include <QuinoaControl.h>
#include <STLMesh.h>
#include <STLTxtMeshReader.h>
#include <SiloWriter.h>

using namespace quinoa;

DiscreteGeometry::DiscreteGeometry(const Base& base)
//******************************************************************************
//  Constructor
//! \param[in] memory    Memory oject pointer
//! \param[in] paradigm  Parallel programming paradigm object pointer
//! \param[in] control   Control object
//! \param[in] timer     Timer object
//! \details Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
try :
  Geometry(base),
  m_mesh(nullptr)
{
  using namespace control;

  //! Echo information on discrete geometry
  echo();

  // Instantiate mesh object
  m_mesh = new(std::nothrow) STLMesh(m_base.memory);
  ErrChk(m_mesh != nullptr, ExceptType::FATAL,
         "Cannot allocate memory for STL mesh object");

  // Instantiate ASCII STL mesh reader object
  STLTxtMeshReader reader(base.control.get<io,input>(), m_mesh);

  // Read in STL mesh
  reader.read();

  // Instantiate Silo writer object
  SiloWriter writer(base.control.get<io,output>(), m_mesh, DB_ALL_AND_DRVR);

  // Write out STL geometry to Silo file
  writer.write();

} // Roll back changes and rethrow on error
  catch (std::exception&) {
    finalize();
    throw;
  }
  // Catch uncaught exceptions
  catch (...) {
    finalize();
    Throw(ExceptType::UNCAUGHT, "Non-standard exception");
  }

DiscreteGeometry::~DiscreteGeometry() noexcept
//******************************************************************************
//  Destructor
//! \details Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  finalize();
}

void
DiscreteGeometry::finalize() noexcept
//******************************************************************************
//  Finalize
//! \details Cleanup either at the end of business as usual or due to an
//!          exception. No-throw guarantee: this member function never throws
//!          exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  if (m_mesh) { delete m_mesh;  m_mesh  = nullptr; }
}

void
DiscreteGeometry::echo()
//******************************************************************************
//  Echo information on discrete geometry
//! \author J. Bakosi
//******************************************************************************
{
  //! Echo information on geometry in general
  Geometry::echo();

  //! Echo information on discrete geometry
}

void
DiscreteGeometry::init()
//******************************************************************************
//  Initialize discrete geometry
//! \author J. Bakosi
//******************************************************************************
{
}
