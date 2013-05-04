//******************************************************************************
/*!
  \file      src/Physics/SPINSFlow/SPINSFlow.C
  \author    J. Bakosi
  \date      Fri 03 May 2013 06:38:32 AM MDT
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
                     const string& filename) :
  Physics(memory, paradigm, control, timer),
  m_npar(control->get<control::NPAR>()),
  m_filename(filename)
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  control  Control object pointer
//! \param[in]  timer    Timer object pointer
//! \param[in]  filename Mesh filename
//! \author  J. Bakosi
//******************************************************************************
{
  // Instantiate random number generator
  m_random = new (nothrow) MKLRandom(m_memory, m_paradigm);
  if (m_random == nullptr) throw Exception(FATAL, "Cannot allocate memory");
  // Create random number leapfrog stream
  m_rndStr = m_random->addStream(VSL_BRNG_MCG59, 0);
  // Get array of MKL VSL stream state pointers right away
  m_str = m_random->getStr(m_rndStr);

  // Instantiate hydrodynamics model
  switch (control->get<control::HYDRO>()) {

    case control::HydroType::SLM :
      m_hydro = new (nothrow) SimplifiedLangevin(memory, paradigm, control);
      if (m_hydro == nullptr) throw Exception(FATAL, "Cannot allocate memory");
      break;

    case control::HydroType::GLM :
      m_hydro = new (nothrow) GeneralizedLangevin(memory, paradigm, control);
      if (m_hydro == nullptr) throw Exception(FATAL, "Cannot allocate memory");
      break;

    default:
      throw Exception(FATAL, "No such hydrodynamics model");
  }

  // Instantiate Eulerian mesh
  m_mesh = new (nothrow) UnsMesh(m_memory);
  if (m_mesh == nullptr) throw Exception(FATAL, "Cannot allocate memory");
}

SPINSFlow::~SPINSFlow() noexcept
//******************************************************************************
//  Destructor
//! \author  J. Bakosi
//******************************************************************************
{
  // Free memory entries held
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
