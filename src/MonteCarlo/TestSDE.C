//******************************************************************************
/*!
  \file      src/MonteCarlo/TestSDE.C
  \author    J. Bakosi
  \date      Sat 22 Feb 2014 06:26:08 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SDE testbed
  \details   SDE testbed
*/
//******************************************************************************
#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP

#include <boost/mpl/cartesian_product.hpp>

#include <TestSDE.h>
#include <GlobWriter.h>
#include <TxtStatWriter.h>
#include <Statistics.h>
#include <Macro.h>
#include <Factory.h>
#include <Dirichlet.h>
#include <GenDirichlet.h>

using quinoa::TestSDE;

TestSDE::TestSDE( const Base& base ) : MonteCarlo( base )
//******************************************************************************
//  Constructor
//! \param[in]  base     Essentials
//! \author  J. Bakosi
//******************************************************************************
{
  //! Initialize factories
  initFactories( print() );

  //! Instantiate SDEs
  instantiateSDE< tag::dirichlet >( ctr::SDEType::DIRICHLET );
  instantiateSDE< tag::gendir >( ctr::SDEType::GENDIR );

  //! Echo information on test bed
  echo();
}

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
    //int tid = omp_get_thread_num();
    #else
    //int tid = 0;
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

void
TestSDE::initFactories(const QuinoaPrint& print)
//******************************************************************************
//  Initialize factories
//! \author  J. Bakosi
//******************************************************************************
{
  // Register SDEs
  namespace mpl = boost::mpl;

  // Dirichlet SDE
  // Construct vector of vectors for all possible policies
  using DirPolicies = mpl::vector< InitPolicies, DirCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< DirPolicies >(
    registerSDE< TestSDE, Dirichlet, tag::dirichlet, ctr::SDEType >
               ( this, ctr::SDEType::DIRICHLET ) );

  // Lochner's generalized Dirichlet SDE
  // Construct vector of vectors for all possible policies
  using GenDirPolicies = mpl::vector< InitPolicies, GenDirCoeffPolicies >;
  // Register SDE for all combinations of policies
  mpl::cartesian_product< GenDirPolicies >(
    registerSDE< TestSDE, GenDirichlet, tag::gendir, ctr::SDEType >
               ( this, ctr::SDEType::GENDIR ) );

  print.optionlist( "Registered SDEs and their policies", m_SDEFactory );
}

void
TestSDE::echo()
//******************************************************************************
//  Echo information on test bed
//! \author  J. Bakosi
//******************************************************************************
{
  print().endpart();
  print().part( "Problem" );

  print().section( "Title", control().get< tag::title >() );

  echoRNGs();

  print().Section< ctr::MonteCarlo, tag::selected, tag::montecarlo >();

  echoIO();

  for (std::size_t i = 0, end = m_sde.size(); i<end; ++i) {
    print().Model< ctr::SDE, tag::selected, tag::sde >( *m_sde[i], i );
  }

  echoIncpar();

  echoIntervals();

  echoStatistics();

  print().endpart();
}
