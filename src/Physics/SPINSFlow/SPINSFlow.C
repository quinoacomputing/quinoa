//******************************************************************************
/*!
  \file      src/Physics/SPINSFlow/SPINSFlow.C
  \author    J. Bakosi
  \date      Tue May  7 10:52:20 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Standalone-Particle Incompressible Navier-Stokes Flow
  \details   Standalone-Particle Incompressible Navier-Stokes Flow
*/
//******************************************************************************

#include <iostream>

#include <sys/time.h>
#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP

#include <Memory.h>
#include <Control.h>
#include <MKLRandom.h>
#include <MKLRndStream.h>
#include <UnsMesh.h>
#include <GmshTxtMeshReader.h>
#include <SimplifiedLangevin.h>
#include <GeneralizedLangevin.h>
#include <SPINSFlow.h>

using namespace Quinoa;

SPINSFlow::SPINSFlow(Memory* const memory,
                     Paradigm* const paradigm,
                     Control* const control,
                     Timer* const timer,
                     const string& filename)
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  control  Control object pointer
//! \param[in]  timer    Timer object pointer
//! \param[in]  filename Mesh filename
//! \author  J. Bakosi
//******************************************************************************
try:
  Physics(memory, paradigm, control, timer),
  m_npar(control->get<control::NPAR>()),
  m_str(nullptr),
  m_filename(filename)
{

  // Instantiate random number generator
  m_random = new (nothrow) MKLRandom(m_memory, m_paradigm);
  ErrChk(m_random != nullptr, FATAL, "Cannot allocate memory");

  // Create random number leapfrog stream
  m_rndStr = m_random->addStream(VSL_BRNG_MCG59, 0);
  // Get array of MKL VSL stream state pointers right away
  m_str = m_random->getStr(m_rndStr);

  // Instantiate hydrodynamics model
  switch (control->get<control::HYDRO>()) {

    case control::HydroType::SLM :
      m_hydro = new (nothrow) SimplifiedLangevin(memory, paradigm, control);
      ErrChk(m_hydro != nullptr, FATAL, "Cannot allocate memory");
      break;

    case control::HydroType::GLM :
      m_hydro = new (nothrow) GeneralizedLangevin(memory, paradigm, control);
      ErrChk(m_hydro != nullptr, FATAL, "Cannot allocate memory");
      break;

    default:
      Throw(FATAL, "No such hydrodynamics model");
  }

  // Instantiate Eulerian mesh
  m_mesh = new (nothrow) UnsMesh(m_memory);
  ErrChk(m_mesh != nullptr, FATAL, "Cannot allocate memory");

} // Roll back changes and rethrow on error
  catch (exception&) {
    finalize();
    throw;
  }
  // Catch uncaught exceptions
  catch (...) {
    finalize();
    Throw(UNCAUGHT, "Non-standard exception");
  }

SPINSFlow::~SPINSFlow() noexcept
//******************************************************************************
//  Destructor
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//! \author  J. Bakosi
//******************************************************************************
{
  finalize();
}

void
SPINSFlow::finalize() noexcept
//******************************************************************************
//  Finalize
//! \details Single exit point, called implicitly from destructor or explicitly
//!          from anywhere else. Exception safety: no-throw guarantee: never
//!          throws exceptions.
//! \author  J. Bakosi
//******************************************************************************
{
  //m_memory->freeEntry(m_MEscalar);

  if (m_mesh) { delete m_mesh; m_mesh = nullptr; }
  if (m_hydro) { delete m_hydro; m_hydro = nullptr; }
  if (m_random) { delete m_random; m_random = nullptr; }
}

void
SPINSFlow::solve()
//******************************************************************************
//  Solve
//! \author  J. Bakosi
//******************************************************************************
{
}

void
SPINSFlow::echo() const
//******************************************************************************
//  Echo information on the physics
//! \author  J. Bakosi
//******************************************************************************
{
}

void
SPINSFlow::init()
//******************************************************************************
//  Initialize the physics
//! \author  J. Bakosi
//******************************************************************************
{
  // Read in mesh
  GmshTxtMeshReader inMesh(m_filename, m_mesh, m_memory);
  inMesh.read();
}
