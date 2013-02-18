//******************************************************************************
/*!
  \file      src/Physics/HomMix/HomMix.C
  \author    J. Bakosi
  \date      Mon 18 Feb 2013 12:34:51 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous material mixing
  \details   Homogeneous material mixing
*/
//******************************************************************************

#include <sys/time.h>

#include <Memory.h>
#include <MemoryException.h>
#include <HomMix.h>
#include <MixException.h>
#include <PDFWriter.h>
#include <ControlTypes.h>
#include <Dirichlet.h>

using namespace Quinoa;
using namespace control;

HomMix::HomMix(Memory* memory,
               Paradigm* paradigm,
               const string& name,
               const MixType mix,
               const int& nscalar,
               const int& npar,
               const real time,
               const int echo,
               const int nstep) :
  Physics(memory, paradigm, name, time, echo, nstep)
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  name     Physics name
//! \param[in]  mix      Mix model
//! \param[in]  nscalar  Number of mixing scalars
//! \param[in]  npar     Number of particles
//! \param[in]  time     Maximum time to simulate
//! \param[in]  echo     One-line info in every few time step
//! \param[in]  nstep    Maximum number of time steps to take
//! \author  J. Bakosi
//******************************************************************************
{
  // Instantiate selected mix model
  switch (mix) {

    case MixType::NO_MIX :
      Throw(MixException,FATAL,NO_MIX);
      break;

    case MixType::DIRICHLET :
      m_mix = new (nothrow) Dirichlet(memory, paradigm, nscalar, npar);
      Assert(m_mix != nullptr, MemoryException,FATAL,BAD_ALLOC);
      break;

    default :
      Throw(MixException,FATAL,MIX_UNIMPLEMENTED);
  }
}

HomMix::~HomMix()
//******************************************************************************
//  Destructor
//! \author  J. Bakosi
//******************************************************************************
{
  if (m_mix) { delete m_mix; m_mix = nullptr; }
}

void
HomMix::solve()
//******************************************************************************
//  Solve
//! \author  J. Bakosi
//******************************************************************************
{
  int it=0;
  real t=0.0;
  long int hrs2end=0, mins2end=0, secs2end=0, hrs2beg=0, mins2beg=0, secs2beg=0;

  // Get start time
  gettimeofday(&m_startTime, static_cast<struct timezone*>(0));

  // Set initial time step size
  real dt = 0.05;

  // Time stepping loop
  while (fabs(t-m_time) > numeric_limits<real>::epsilon() && it < m_nstep) {

    // Advance particles
    m_mix->advance(dt);

    // Echo one-liner info
    if (!(it % m_echo)) {
      report(it, t, dt,
             hrs2beg, mins2beg, secs2beg, hrs2end, mins2end, secs2end);
    }

    // Increase timestep and iteration counter
    t += dt;
    ++it;
    if (t > m_time) t = m_time;
  }

  outJPDF();
}

void
HomMix::outJPDF()
//******************************************************************************
//  Output joint scalar PDF
//! \author  J. Bakosi
//******************************************************************************
{
  JPDF jpdf(2, 0.01);  

  m_mix->jpdf(jpdf);

  PDFWriter jpwt(m_memory,"jpdf.txt");
  jpwt.writeTxt(&jpdf);
  PDFWriter jpwg(m_memory,"jpdf.msh");
  jpwg.writeGmsh(&jpdf);
}

void
HomMix::echo()
//******************************************************************************
//  Echo informaion on homogeneous material mix
//! \author  J. Bakosi
//******************************************************************************
{
}

void
HomMix::init()
//******************************************************************************
//  Initialize homogeneous material mix
//! \author  J. Bakosi
//******************************************************************************
{
  m_mix->init();
}
