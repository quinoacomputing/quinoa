//******************************************************************************
/*!
  \file      src/Physics/HomMix/HomMix.C
  \author    J. Bakosi
  \date      Sun 12 May 2013 08:14:25 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous material mixing
  \details   Homogeneous material mixing
*/
//******************************************************************************

#include <cmath>
#include <iomanip>
#include <sstream>

#include <Memory.h>
#include <Control.h>
#include <Mix.h>
#include <HomMix.h>
#include <PDFWriter.h>
#include <GlobWriter.h>
#include <TxtPlotWriter.h>
#include <Statistics.h>
#include <Dirichlet.h>
#include <GeneralizedDirichlet.h>

using namespace Quinoa;

HomMix::HomMix(Memory* const memory,
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
HomMix::solve()
//******************************************************************************
//  Solve
//! \author  J. Bakosi
//******************************************************************************
{
  int it = 0;
  real t = 0.0;
  bool wroteJpdf = false;
  bool wroteGlob = false;
  bool wrotePlot = false;

  const auto nstep = control()->get<control::NSTEP>();
  const auto ttyi  = control()->get<control::TTYI>();
  const auto pdfi  = control()->get<control::PDFI>();
  const auto glob  = control()->get<control::GLOB>();
  const auto plti  = control()->get<control::PLTI>();
  const auto dt    = control()->get<control::DT>();

  timer()->start(m_totalTime);

  // Echo headers
  if (nstep) {
    reportHeader();
    plotWriter()->header();
  }

  // Time stepping loop
  while (fabs(t-m_term) > numeric_limits<real>::epsilon() && it < nstep) {

    // Advance particles
    mix()->advance(dt);

    // Accumulate statistics
    statistics()->accumulate();

    // Output pdf at selected times
    if (!(it % pdfi)) { outJpdf(t); wroteJpdf = true; }

    // Append glob file at selected times
    if (!(it % glob)) { globWriter()->write(it,t); wroteGlob = true; }

    // Append plot file at selected times
    if (!(it % plti)) { plotWriter()->write(it,t); wrotePlot = true; }

    // Echo one-liner info
    if (!(it % ttyi)) {
      report(it, nstep, t, dt, wroteJpdf, wroteGlob, wrotePlot);
      wroteJpdf = wroteGlob = wrotePlot = false;
    }

    // Increase timestep and iteration counter
    t += dt;
    ++it;
    if (t > m_term) t = m_term;
  } // Time stepping loop
}

void
HomMix::reportHeader() const
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
               const bool wroteJpdf,
               const bool wroteGlob,
               const bool wrotePlot)
//******************************************************************************
//  One-liner report
//! \param[in]  it         Iteration counter
//! \param[in]  nstep      Terminate time
//! \param[in]  t          Time
//! \param[in]  dt         Time step size
//! \param[in]  wroteJpdf  True if joint PDF was output
//! \param[in]  wroteGlob  True if glob was output
//! \param[in]  wrotePlot  True if plot was output
//! \author  J. Bakosi
//******************************************************************************
{
  Watch ete, eta;       // estimated time elapsed and to accomplishment
  timer()->eta(m_totalTime, m_term, t, nstep, it, ete, eta);

  cout << setfill(' ') << setw(8) << it << "  "
       << scientific << setprecision(6) << setw(12) << t << "  " << dt << "  "
       << setfill('0') << setw(3) << ete.h.count() << ":"
                       << setw(2) << ete.m.count() << ":"
                       << setw(2) << ete.s.count() << "  "
                       << setw(3) << eta.h.count() << ":"
                       << setw(2) << eta.m.count() << ":"
                       << setw(2) << eta.s.count() << "  ";

  if (wroteGlob) cout << "G";
  if (wroteJpdf) cout << "J";
  if (wrotePlot) cout << "P";

  cout << endl;
}

void
HomMix::outJpdf(const real t)
//******************************************************************************
//  Output joint scalar PDF
//! \param[in]  t    Time stamp
//! \author  J. Bakosi
//******************************************************************************
{
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
  mix()->init();
}
