//******************************************************************************
/*!
  \file      src/MonteCarlo/TestSDE.C
  \author    J. Bakosi
  \date      Sat 01 Feb 2014 11:03:54 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SDE testbed
  \details   SDE testbed
*/
//******************************************************************************

#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP

#include <TestSDE.h>
#include <GlobWriter.h>
#include <TxtStatWriter.h>
#include <Statistics.h>
#include <Macro.h>

using quinoa::TestSDE;

void
TestSDE::run()
//******************************************************************************
//  Run
//! \author  J. Bakosi
//******************************************************************************
{
  uint64_t it = 0;
  tk::real t = 0.0;
  bool wroteJpdf = false;
  bool wroteGlob = false;
  bool wroteStat = false;

  const auto nstep = control().get<tag::incpar, tag::nstep>();
  const auto dt    = control().get<tag::incpar, tag::dt>();
  const auto ttyi  = control().get<tag::interval, tag::tty>();
  const auto pdfi  = control().get<tag::interval, tag::pdf>();
  const auto glbi  = control().get<tag::interval, tag::glob>();
  const auto stai  = control().get<tag::interval, tag::plot>();

  timer().start( m_totalTime );

  // Echo headers
  if (nstep) {
    header();
    statWriter().header();
  }

  // Time stepping loop
  tk::real eps = std::numeric_limits< tk::real >::epsilon();
  while (fabs(t - m_term) > eps && it < nstep) {

    // Advance particles
    advance(dt);

    // Accumulate statistics
    statistics().accumulate();

    // Output pdf at selected times
    if (!(it % pdfi)) { outJpdf(t); wroteJpdf = true; }

    // Append glob file at selected times
    if (!(it % glbi)) { globWriter().write(it,t); wroteGlob = true; }

    // Append statistics file at selected times
    if (!(it % stai)) { statWriter().write(it,t); wroteStat = true; }

    // Echo one-liner info
    if (!(it % ttyi)) {
      report(it, nstep, t, dt, wroteJpdf, wroteGlob, wroteStat);
      wroteJpdf = wroteGlob = wroteStat = false;
    }

    // Increase timestep and iteration counter
    t += dt;
    ++it;
    if (t > m_term) t = m_term;
  }
}

void
TestSDE::advance(tk::real dt)
//******************************************************************************
//  Advance particles
//! \author  J. Bakosi
//******************************************************************************
{
  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
    #ifdef _OPENMP
    int tid = omp_get_thread_num();
    #else
    int tid = 0;
    #endif

    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (uint64_t p=0; p<m_npar; ++p) {
      //sde()->advance( p, tid, dt );
    }
  }
}

void
TestSDE::outJpdf(tk::real t)
//******************************************************************************
//  Output joint scalar PDF
//! \param[in]  t    Time stamp
//! \author  J. Bakosi
//******************************************************************************
{
IGNORE(t);
//   // Contruct filename
//   stringstream ss;
//   ss << control()->get<ctr::PDFNAME>() << "." << t << ".msh";
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
