//******************************************************************************
/*!
  \file      src/MonteCarlo/HomRT.C
  \author    J. Bakosi
  \date      Wed Apr 23 11:21:53 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous material mixing
  \details   Homogeneous material mixing
*/
//******************************************************************************

#include <iomanip>

#include <Control.h>
#include <HomRT.h>
#include <PDFWriter.h>
#include <GlobWriter.h>
#include <TxtStatWriter.h>
#include <Statistics.h>

using quinoa::HomRT;

void
HomRT::run()
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

  //timer().start( m_totalTime );

  // Echo headers
  if (nstep) {
    header();
    statWriter().header();
  }

  // Time stepping loop
  while (fabs(t-m_term) > std::numeric_limits<tk::real>::epsilon() && it<nstep) {

    // Advance particles
    //advance(dt);

    // Accumulate statistics
    statistics().accumulate();

    // Output pdf at selected times
    if (!(it % pdfi)) { outJpdf(t); wroteJpdf = true; }

    // Append glob file at selected times
    if (!(it % glbi)) { globWriter().writeGlob(it,t); wroteGlob = true; }

    // Append statistics file at selected times
    if (!(it % stai)) { statWriter().writeStat(it,t); wroteStat = true; }

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
HomRT::advance(tk::real dt)
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
      mass()->advance(p, tid, dt);
      hydro()->advance(p, tid, dt);
    }
  }
}

void
HomRT::outJpdf(const tk::real t)
//******************************************************************************
//  Output joint scalar PDF
//! \param[in]  t    Time stamp
//! \author  J. Bakosi
//******************************************************************************
{
  // Contruct filename
  std::stringstream ss;
  ss << control().get< tag::cmd, tag::io, tag::pdf >() << "." << t << ".msh";
  std::string filename = ss.str();

  // Create joint PDF
  tk::JPDF jpdf(m_nscalar, 0.02);

  // Estimate joint PDF
  //mix()->jpdf(jpdf);

  // Output joint PDF
  PDFWriter jpdfFile(filename);
  jpdfFile.writeGmsh(jpdf);
}
