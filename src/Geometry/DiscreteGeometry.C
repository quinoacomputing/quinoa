//******************************************************************************
/*!
  \file      src/Geometry/DiscreteGeometry.C
  \author    J. Bakosi
  \date      Fri 26 Jul 2013 08:52:17 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Discrete geometry definition
  \details   Discrete geometry definition
*/
//******************************************************************************

#include <DiscreteGeometry.h>
#include <Exception.h>
#include <Control.h>
#include <STLMesh.h>
#include <STLTxtMeshReader.h>
#include <SiloWriter.h>

using namespace Quinoa;

DiscreteGeometry::DiscreteGeometry(Memory* const memory,
                                   Paradigm* const paradigm,
                                   Control* const control,
                                   Timer* const timer)
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
  Geometry(memory, paradigm, control, timer),
  m_mesh(nullptr)
{

  // Instantiate mesh object
  m_mesh = new(std::nothrow) STLMesh(memory);
  ErrChk(m_mesh != nullptr, ExceptType::FATAL,
         "Cannot allocate memory for STL mesh object");

  // Instantiate ASCII STL mesh reader object
  STLTxtMeshReader reader(control->get<control::INPUT>(), m_mesh);

  // Read in STL mesh
  reader.read();

  // Instantiate Silo writer object
  SiloWriter writer(control->get<control::OUTPUT>(), m_mesh, DB_ALL_AND_DRVR);

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
DiscreteGeometry::init()
//******************************************************************************
//  Initialize discrete geometry
//! \author J. Bakosi
//******************************************************************************
{
}
