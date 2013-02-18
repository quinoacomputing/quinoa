//******************************************************************************
/*!
  \file      src/Physics/SPINSFlow/SPINSFlow.C
  \author    J. Bakosi
  \date      Mon 18 Feb 2013 01:54:26 PM MST
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
#include <MemoryException.h>
#include <Control.h>
#include <MKLRandom.h>
#include <MKLRndStream.h>
#include <UnsMesh.h>
#include <GmshTxtMeshReader.h>
#include <SimplifiedLangevin.h>
#include <GeneralizedLangevin.h>
#include <HydroException.h>
#include <SPINSFlow.h>

using namespace Quinoa;
using namespace control;

SPINSFlow::SPINSFlow(Memory* const memory,
                     Paradigm* const paradigm,
                     Control* const control,
                     const string& filename) :
  Physics(memory, paradigm, control),
  m_npar(control->get<NPAR>()),
  m_filename(filename)
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  control  Control object pointer
//! \param[in]  filename Mesh filename
//! \author  J. Bakosi
//******************************************************************************
{
  // Instantiate random number generator
  m_random = new (nothrow) MKLRandom(m_memory, m_paradigm);
  Assert(m_random != nullptr, MemoryException,FATAL,BAD_ALLOC);
  // Create random number leapfrog stream
  m_rndStr = m_random->addStream(VSL_BRNG_MCG59, 0);
  // Get array of MKL VSL stream state pointers right away
  m_str = m_random->getStr(m_rndStr);

  // Instantiate hydrodynamics model
  switch (control->get<HYDRO>()) {

    case HydroType::SLM :
      m_hydro = new (nothrow) SimplifiedLangevin(memory, paradigm, control);
      Assert(m_hydro != nullptr, MemoryException,FATAL,BAD_ALLOC);
      break;

    case HydroType::GLM :
      m_hydro = new (nothrow) GeneralizedLangevin(memory, paradigm, control);
      Assert(m_hydro != nullptr, MemoryException,FATAL,BAD_ALLOC);
      break;

    default:
      Throw(HydroException,FATAL,NO_SUCH_HYDRO);
  }

  // Instantiate Eulerian mesh
  m_mesh = new (nothrow) UnsMesh(m_memory);
  Assert(m_mesh != nullptr, MemoryException,FATAL,BAD_ALLOC);
}

SPINSFlow::~SPINSFlow()
//******************************************************************************
//  Destructor
//! \author  J. Bakosi
//******************************************************************************
{
  // Free memory entries held
#ifndef NDEBUG  // Error checking and exceptions only in debug mode
  try {
#endif // NDEBUG
    //m_memory->freeEntry(m_MEscalar);
#ifndef NDEBUG
  } catch (...)
    { cout << "WARNING: Exception in SPINSFlow destructor" << endl; }
#endif // NDEBUG

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
SPINSFlow::echo()
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
