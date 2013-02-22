//******************************************************************************
/*!
  \file      src/Physics/HomMix/HomMix.C
  \author    J. Bakosi
  \date      Fri Feb 22 16:30:10 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous material mixing
  \details   Homogeneous material mixing
*/
//******************************************************************************

#include <sstream>

//#include <sys/time.h>

#include <Memory.h>
#include <MemoryException.h>
#include <ControlTypes.h>
#include <Control.h>
#include <HomMix.h>
#include <MixException.h>
#include <PDFWriter.h>
#include <Dirichlet.h>
#include <Timer.h>

using namespace Quinoa;
using namespace control;

HomMix::HomMix(Memory* const memory,
               Paradigm* const paradigm,
               Control* const control,
               Timer* const timer) :
  Physics(memory, paradigm, control, timer),
  m_term(control->get<TERM>()),
  m_jpdf_filename_base(control->get_jpdf_filename_base()),
  m_totalTime(timer->create("Total solution"))
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  control  Control object pointer
//! \author  J. Bakosi
//******************************************************************************
{
  // Instantiate selected mix model
  switch (control->get<MIX>()) {

    case MixType::NO_MIX :
      Throw(MixException,FATAL,MixExceptType::NO_MIX);
      break;

    case MixType::DIRICHLET :
      m_mix = new (nothrow) Dirichlet(memory, paradigm, control);
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

  const auto nstep = m_control->get<NSTEP>();
  const auto ttyi  = m_control->get<TTYI>();
  const auto pdfi  = m_control->get<PDFI>();
  const auto dt    = m_control->get<DT>();

  m_timer->start(m_totalTime);

  // Time stepping loop
  while (fabs(t-m_term) > numeric_limits<real>::epsilon() && it < nstep) {

    // Advance particles
    m_mix->advance(dt);

    // Echo one-liner info
    if (!(it % ttyi)) {
      report(it, nstep, t, dt);
    }

    // Output pdf at selected times
    if (!(it % pdfi)) {
      outJPDF(t);
    }

    // Increase timestep and iteration counter
    t += dt;
    ++it;
    if (t > m_term) t = m_term;
  }
}

void
HomMix::report(const int it, const int nstep, const real t, const real dt)
//******************************************************************************
//  One-liner report
//! \param[in]  it        Iteration counter
//! \param[in]  t         Time counter
//! \param[in]  dt        Time step size
//! \author  J. Bakosi
//******************************************************************************
{
  HMS elapsed, estimated;
  m_timer->query(m_totalTime, elapsed);
  m_timer->eta(m_totalTime, m_term, t, dt, nstep, it, estimated);

  cout << "it = " << it << ", t = " << t << "\t dt = " << dt << ", "
       << elapsed.h.count() << ":"
       << elapsed.m.count() << ":"
       << elapsed.s.count() << " "
       << estimated.h.count() << ":"
       << estimated.m.count() << ":"
       << estimated.s.count() << endl;
}

void
HomMix::outJPDF(const real t)
//******************************************************************************
//  Output joint scalar PDF
//! \param[in]  t    Time stamp
//! \author  J. Bakosi
//******************************************************************************
{
  // Echo
  cout << "JPDF at t = " << t << endl;

  // Contruct filename
  stringstream ss;
  ss << m_jpdf_filename_base << "." << t << ".msh";
  string filename = ss.str();

  // Create joint PDF
  JPDF jpdf(2, 0.02);

  // Estimate joint PDF
  m_mix->jpdf(jpdf);

  // Output joint PDF
  PDFWriter jpdfFile(m_memory, filename);
  jpdfFile.writeGmsh(&jpdf);
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
