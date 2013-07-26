//******************************************************************************
/*!
  \file      src/Physics/HomRT/HomRT.C
  \author    J. Bakosi
  \date      Fri Jul 26 15:34:57 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous material mixing
  \details   Homogeneous material mixing
*/
//******************************************************************************

#include <sstream>
#include <iomanip>

#include <Memory.h>
#include <ControlTypes.h>
#include <Control.h>
#include <HomRT.h>
#include <PDFWriter.h>
#include <GlobWriter.h>
#include <TxtStatWriter.h>
#include <Statistics.h>
#include <Beta.h>
#include <SimplifiedLangevin.h>

using namespace Quinoa;

HomRT::HomRT(Memory* const memory,
             Paradigm* const paradigm,
             Control* const control,
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
HomRT::solve()
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

  const auto nstep = control()->get<control::NSTEP>();
  const auto ttyi  = control()->get<control::TTYI>();
  const auto pdfi  = control()->get<control::PDFI>();
  const auto glbi  = control()->get<control::GLBI>();
  const auto stai  = control()->get<control::STAI>();
  const auto dt    = control()->get<control::DT>();

  timer()->start(m_totalTime);

  // Echo headers
  if (nstep) {
    reportHeader();
    statWriter()->header();
  }

  // Time stepping loop
  while (fabs(t-m_term) > std::numeric_limits<real>::epsilon() && it < nstep) {

    // Advance particles
    mass()->advance(dt);
    hydro()->advance(dt);

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
HomRT::reportHeader() const
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
HomRT::report(const uint64_t it,
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
HomRT::outJpdf(const real t)
//******************************************************************************
//  Output joint scalar PDF
//! \param[in]  t    Time stamp
//! \author  J. Bakosi
//******************************************************************************
{
  // Contruct filename
  std::stringstream ss;
  ss << control()->get<control::PDFNAME>() << "." << t << ".msh";
  std::string filename = ss.str();

  // Create joint PDF
  JPDF jpdf(m_nscalar, 0.02);

  // Estimate joint PDF
  //mix()->jpdf(jpdf);

  // Output joint PDF
  PDFWriter jpdfFile(filename);
  jpdfFile.writeGmsh(&jpdf);
}

void
HomRT::init()
//******************************************************************************
//  Initialize homogeneous Rayleigh-Taylor
//! \author  J. Bakosi
//******************************************************************************
{
  mass()->init();
  hydro()->init();
}
