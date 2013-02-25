//******************************************************************************
/*!
  \file      src/Physics/HomMix/HomMix.C
  \author    J. Bakosi
  \date      Sun 24 Feb 2013 07:14:37 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous material mixing
  \details   Homogeneous material mixing
*/
//******************************************************************************

#include <sstream>
#include <iomanip>

#include <Memory.h>
#include <MemoryException.h>
#include <ControlTypes.h>
#include <Control.h>
#include <HomMix.h>
#include <MixException.h>
#include <PDFWriter.h>
#include <Dirichlet.h>
#include <GeneralizedDirichlet.h>
#include <Timer.h>

using namespace Quinoa;
using namespace control;

HomMix::HomMix(Memory* const memory,
               Paradigm* const paradigm,
               Control* const control,
               Timer* const timer) :
  Physics(memory, paradigm, control, timer),
  m_nscalar(control->get<NSCALAR>()),
  m_term(control->get<TERM>()),
  m_jpdf_filename_base(control->get_jpdf_filename_base()),
  m_totalTime(timer->create("Total solution"))
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  control  Control object pointer
//! \param[in]  timer    Timer object pointer
//! \author  J. Bakosi
//******************************************************************************
{
  // Instantiate selected mix model
  switch (control->get<MIX>()) {

    case MixType::NO_MIX :
      Throw(MixException,FATAL,MixExceptType::NO_MIX);
      break;

    case MixType::GENERALIZED_DIRICHLET :
      m_mix = new (nothrow) GeneralizedDirichlet(memory, paradigm, control);
      Assert(m_mix != nullptr, MemoryException,FATAL,BAD_ALLOC);
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
  bool wroteJPDF=false;

  const auto nstep = m_control->get<NSTEP>();
  const auto ttyi  = m_control->get<TTYI>();
  const auto pdfi  = m_control->get<PDFI>();
  const auto dt    = m_control->get<DT>();

  m_timer->start(m_totalTime);

  if (nstep) reportHeader();

  // Time stepping loop
  while (fabs(t-m_term) > numeric_limits<real>::epsilon() && it < nstep) {

    // Advance particles
    m_mix->advance(dt);

    // Output pdf at selected times
    if (!(it % pdfi)) {
      outJPDF(t);
      wroteJPDF = true;
    }

    // Echo one-liner info
    if (!(it % ttyi)) {
      report(it, nstep, t, dt, wroteJPDF);
      wroteJPDF = false;
    }

    // Increase timestep and iteration counter
    t += dt;
    ++it;
    if (t > m_term) t = m_term;
  }
}

void
HomMix::reportHeader()
//******************************************************************************
//  Echo report header
//! \author  J. Bakosi
//******************************************************************************
{
  cout << "      it             t            dt"
          "        ETE        ETA   out\n"
       << "------------------------------------"
          "----------------------------" << endl;
}

void
HomMix::report(const int it,
               const int nstep,
               const real t,
               const real dt,
               const bool wroteJPDF)
//******************************************************************************
//  One-liner report
//! \param[in]  it         Iteration counter
//! \param[in]  nstep      Terminate time
//! \param[in]  t          Time
//! \param[in]  dt         Time step size
//! \param[in]  wroteJPDF  True if joint PDF was written
//! \author  J. Bakosi
//******************************************************************************
{
  Watch ete, eta;
  m_timer->eta(m_totalTime, m_term, t, nstep, it, ete, eta);

  cout << setfill(' ') << setw(8) << it << "  "
       << scientific << setprecision(6) << setw(12) << t << "  " << dt << "  "
       << setfill('0') << setw(3) << ete.h.count() << ":"
                       << setw(2) << ete.m.count() << ":"
                       << setw(2) << ete.s.count() << "  "
                       << setw(3) << eta.h.count() << ":"
                       << setw(2) << eta.m.count() << ":"
                       << setw(2) << eta.s.count() << "  ";

  if (wroteJPDF) cout << "J";

  cout << endl;
}

void
HomMix::outJPDF(const real t)
//******************************************************************************
//  Output joint scalar PDF
//! \param[in]  t    Time stamp
//! \author  J. Bakosi
//******************************************************************************
{
  // Contruct filename
  stringstream ss;
  ss << m_jpdf_filename_base << "." << t << ".msh";
  string filename = ss.str();

  // Create joint PDF
  JPDF jpdf(m_nscalar, 0.02);

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
