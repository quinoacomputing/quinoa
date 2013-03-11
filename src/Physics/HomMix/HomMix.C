//******************************************************************************
/*!
  \file      src/Physics/HomMix/HomMix.C
  \author    J. Bakosi
  \date      Sun 10 Mar 2013 07:23:49 PM MDT
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
#include <GlobWriter.h>
#include <Dirichlet.h>
#include <GeneralizedDirichlet.h>
#include <Timer.h>
#include <Statistics.h>

using namespace Quinoa;

HomMix::HomMix(Memory* const memory,
               Paradigm* const paradigm,
               Control* const control,
               Timer* const timer) :
  Physics(memory, paradigm, control, timer),
  m_nscalar(control->get<control::NSCALAR>()),
  m_term(control->get<control::TERM>()),
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
  switch (control->get<control::MIX>()) {

    case control::MixType::NO_MIX :
      Throw(MixException,FATAL,MixExceptType::NO_MIX);
      break;

    case control::MixType::DIRICHLET :
      m_mix = new (nothrow) Dirichlet(memory, paradigm, control);
      Assert(m_mix != nullptr, MemoryException,FATAL,BAD_ALLOC);
      break;

    case control::MixType::GENERALIZED_DIRICHLET :
      m_mix = new (nothrow) GeneralizedDirichlet(memory, paradigm, control);
      Assert(m_mix != nullptr, MemoryException,FATAL,BAD_ALLOC);
      break;

    default :
      Throw(MixException,FATAL,MIX_UNIMPLEMENTED);
  }

  // Instantiate statistics estimator
  m_statistics = new (nothrow) Statistics(memory, paradigm, control, m_mix);
  Assert(m_statistics != nullptr, MemoryException,FATAL,BAD_ALLOC);

  // Instantiate glob file writer
  m_glob = new (nothrow) GlobWriter(m_control->get<control::GLOBNAME>());
  Assert(m_glob != nullptr, MemoryException,FATAL,BAD_ALLOC);
}

HomMix::~HomMix()
//******************************************************************************
//  Destructor
//! \author  J. Bakosi
//******************************************************************************
{
  if (m_glob) { delete m_glob; m_glob = nullptr; }
  if (m_statistics) { delete m_statistics; m_statistics = nullptr; }
  if (m_mix) { delete m_mix; m_mix = nullptr; }
}

void
HomMix::solve()
//******************************************************************************
//  Solve
//! \author  J. Bakosi
//******************************************************************************
{
  int it = 0;
  real t = 0.0;
  bool wroteJPDF = false;
  bool wroteGlob = false;

  const auto nstep = m_control->get<control::NSTEP>();
  const auto ttyi  = m_control->get<control::TTYI>();
  const auto pdfi  = m_control->get<control::PDFI>();
  const auto glob  = m_control->get<control::GLOB>();
  const auto dt    = m_control->get<control::DT>();

  m_timer->start(m_totalTime);

  // Echo header
  if (nstep) reportHeader();

  // Time stepping loop
  while (fabs(t-m_term) > numeric_limits<real>::epsilon() && it < nstep) {

    // Advance particles
    m_mix->advance(dt);

    // Accumulate statistics
    m_statistics->accumulate();

    // Output pdf at selected times
    if (!(it % pdfi)) {
      outJPDF(t);
      wroteJPDF = true;
    }

    // Output domain-average statistics at selected times
    if (!(it % glob)) {
      outGlob(it,t);
      wroteGlob = true;
    }

    // Echo one-liner info
    if (!(it % ttyi)) {
      report(it, nstep, t, dt, wroteJPDF, wroteGlob);
      wroteJPDF = wroteGlob = false;
    }

    // Increase timestep and iteration counter
    t += dt;
    ++it;
    if (t > m_term) t = m_term;
  } // Time stepping loop
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
               const bool wroteJPDF,
               const bool wroteGlob)
//******************************************************************************
//  One-liner report
//! \param[in]  it         Iteration counter
//! \param[in]  nstep      Terminate time
//! \param[in]  t          Time
//! \param[in]  dt         Time step size
//! \param[in]  wroteJPDF  True if joint PDF was output
//! \param[in]  wroteGlob  True if domain-average statistics were output
//! \author  J. Bakosi
//******************************************************************************
{
  Watch ete, eta;       // estimated time elapsed and to accomplishment
  m_timer->eta(m_totalTime, m_term, t, nstep, it, ete, eta);

  cout << setfill(' ') << setw(8) << it << "  "
       << scientific << setprecision(6) << setw(12) << t << "  " << dt << "  "
       << setfill('0') << setw(3) << ete.h.count() << ":"
                       << setw(2) << ete.m.count() << ":"
                       << setw(2) << ete.s.count() << "  "
                       << setw(3) << eta.h.count() << ":"
                       << setw(2) << eta.m.count() << ":"
                       << setw(2) << eta.s.count() << "  ";

  if (wroteGlob) cout << "G";
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
  ss << m_control->get<control::JPDFNAME>() << "." << t << ".msh";
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
HomMix::outGlob(const int it, const real t)
//******************************************************************************
//  Output domain-average statistics
//! \param[in]  t    Time stamp
//! \author  J. Bakosi
//******************************************************************************
{
  m_glob->write(it, t, m_statistics->nord(), m_statistics->ordinary());
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
