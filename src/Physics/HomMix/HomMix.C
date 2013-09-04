//******************************************************************************
/*!
  \file      src/Physics/HomMix/HomMix.C
  \author    J. Bakosi
  \date      Wed Sep  4 08:09:01 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous material mixing
  \details   Homogeneous material mixing
*/
//******************************************************************************

#include <cmath>
#include <iomanip>
#include <sstream>

#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP

#include <Macro.h>
#include <Memory.h>
#include <QuinoaControl.h>
#include <Mix.h>
#include <HomMix.h>
#include <PDFWriter.h>
#include <GlobWriter.h>
#include <TxtStatWriter.h>
#include <Statistics.h>
#include <Dirichlet.h>
#include <GenDirichlet.h>

using namespace quinoa;

HomMix::HomMix(Memory* const memory,
               Paradigm* const paradigm,
               const QuinoaControl& control,
               Timer* const timer) :
  Physics(memory, paradigm, control, timer),
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
}

void
HomMix::solve()
//******************************************************************************
//  Solve
//! \author  J. Bakosi
//******************************************************************************
{
  uint64_t it = 0;
  real t = 0.0;
  bool wroteJpdf = false;
  bool wroteGlob = false;
  bool wroteStat = false;

  const auto nstep = control().get<control::incpar>().get<control::nstep>();
  const auto dt    = control().get<control::incpar>().get<control::dt>();
  const auto ttyi  = control().get<control::interval>().get<control::tty>();
  const auto pdfi  = control().get<control::interval>().get<control::pdf>();
  const auto glbi  = control().get<control::interval>().get<control::glob>();
  const auto stai  = control().get<control::interval>().get<control::plot>();

  timer()->start(m_totalTime);

  // Echo headers
  if (nstep) {
    reportHeader();
    statWriter()->header();
  }

  // Time stepping loop
  while (fabs(t-m_term) > std::numeric_limits<real>::epsilon() && it < nstep) {

    // Advance particles
    advance(dt);

    // Accumulate statistics
    statistics()->accumulate();

    // Output pdf at selected times
    if (!(it % pdfi)) { outJpdf(t); wroteJpdf = true; }

    // Append glob file at selected times
    if (!(it % glbi)) { globWriter()->write(it,t); wroteGlob = true; }

    // Append statistics file at selected times
    if (!(it % stai)) { statWriter()->write(it,t); wroteStat = true; }

    // Echo one-liner info
    if (!(it % ttyi)) {
      report(it, nstep, t, dt, wroteJpdf, wroteGlob, wroteStat);
      wroteJpdf = wroteGlob = wroteStat = false;
    }

    // Increase timestep and iteration counter
    t += dt;
    ++it;
    if (t > m_term) t = m_term;
  } // Time stepping loop
}

void
HomMix::advance(real dt)
//******************************************************************************
//  Advance particles
//! \author  J. Bakosi
//******************************************************************************
{
  uint64_t p;
  int tid;

  #ifdef _OPENMP
  #pragma omp parallel private(tid, p)
  #endif
  {
    #ifdef _OPENMP
    tid = omp_get_thread_num();
    #else
    tid = 0;
    #endif

    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (p=0; p<m_npar; ++p) {

      mix()->advance(p, tid, dt);

    } // m_npar
  } // omp parallel
}

void
HomMix::reportHeader() const
//******************************************************************************
//  Echo report header
//! \author  J. Bakosi
//******************************************************************************
{
  std::cout << "      it             t            dt"
               "        ETE        ETA   out\n"
            << "------------------------------------"
               "----------------------------" << std::endl;
}

void
HomMix::report(const uint64_t it,
               const uint64_t nstep,
               const real t,
               const real dt,
               const bool wroteJpdf,
               const bool wroteGlob,
               const bool wroteStat)
//******************************************************************************
//  One-liner report
//! \param[in]  it         Iteration counter
//! \param[in]  nstep      Terminate time
//! \param[in]  t          Time
//! \param[in]  dt         Time step size
//! \param[in]  wroteJpdf  True if joint PDF was output
//! \param[in]  wroteGlob  True if glob was output
//! \param[in]  wroteStat  True if statistics was output
//! \author  J. Bakosi
//******************************************************************************
{
  Watch ete, eta;       // estimated time elapsed and to accomplishment
  timer()->eta(m_totalTime, m_term, t, nstep, it, ete, eta);

  std::cout << std::setfill(' ') << std::setw(8) << it << "  "
            << std::scientific << std::setprecision(6) << std::setw(12) << t
            << "  " << dt << "  " << std::setfill('0')
            << std::setw(3) << ete.h.count() << ":"
            << std::setw(2) << ete.m.count() << ":"
            << std::setw(2) << ete.s.count() << "  "
            << std::setw(3) << eta.h.count() << ":"
            << std::setw(2) << eta.m.count() << ":"
            << std::setw(2) << eta.s.count() << "  ";

  if (wroteGlob) std::cout << "G";
  if (wroteJpdf) std::cout << "J";
  if (wroteStat) std::cout << "P";

  std::cout << std::endl;
}

void
HomMix::outJpdf(const real t)
//******************************************************************************
//  Output joint scalar PDF
//! \param[in]  t    Time stamp
//! \author  J. Bakosi
//******************************************************************************
{
IGNORE(t);
//   // Contruct filename
//   stringstream ss;
//   ss << control()->get<control::PDFNAME>() << "." << t << ".msh";
//   string filename = ss.str();
// 
//   // Create joint PDF
//   JPDF jpdf(m_nscalar, 0.02);
// 
//   // Estimate joint PDF
//   m_mix->jpdf(jpdf);
// 
//   // Output joint PDF
//   PDFWriter jpdfFile(filename);
//   jpdfFile.writeGmsh(&jpdf);
}

void
HomMix::init()
//******************************************************************************
//  Initialize homogeneous material mix
//! \author  J. Bakosi
//******************************************************************************
{
  uint64_t p;
  int tid;

  #ifdef _OPENMP
  #pragma omp parallel private(tid, p)
  #endif
  {
    #ifdef _OPENMP
    tid = omp_get_thread_num();
    #else
    tid = 0;
    #endif

    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (p=0; p<m_npar; ++p) {

      mix()->init(p, tid);

    } // m_npar
  } // omp parallel

}
