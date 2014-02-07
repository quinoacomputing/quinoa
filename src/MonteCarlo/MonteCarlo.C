//******************************************************************************
/*!
  \file      src/MonteCarlo/MonteCarlo.C
  \author    J. Bakosi
  \date      Sat 01 Feb 2014 11:00:00 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MonteCarlo base
  \details   MonteCarlo base
*/
//******************************************************************************

#include <MonteCarlo.h>

using quinoa::MonteCarlo;

void
MonteCarlo::header() const
//******************************************************************************
//  Echo standard header
//! \author  J. Bakosi
//******************************************************************************
{
  tk::Option< ctr::MonteCarlo > mc;
  auto& name = mc.name( control().get< tag::selected, tag::montecarlo >() );

  print().raw( "Start solving " + name + "\n\n" );
  print().raw( "      it             t            dt"
               "        ETE        ETA   out\n"
               "------------------------------------"
               "----------------------------\n" );
}

void
MonteCarlo::report( uint64_t it,
                    uint64_t nstep,
                    tk::real t,
                    tk::real dt,
                    bool wroteJpdf,
                    bool wroteGlob,
                    bool wroteStat )
//******************************************************************************
//  Echo standard one-liner report
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
  tk::Watch ete, eta;       // estimated time elapsed and to accomplishment
  timer().eta( m_totalTime, m_term, t, nstep, it, ete, eta );

  print().stream() << std::setfill(' ') << std::setw(8) << it << "  "
                   << std::scientific << std::setprecision(6)
                   << std::setw(12) << t << "  "
                   << dt << "  "
                   << std::setfill('0')
                   << std::setw(3) << ete.h.count() << ":"
                   << std::setw(2) << ete.m.count() << ":"
                   << std::setw(2) << ete.s.count() << "  "
                   << std::setw(3) << eta.h.count() << ":"
                   << std::setw(2) << eta.m.count() << ":"
                   << std::setw(2) << eta.s.count() << "  ";

  if (wroteGlob) print().raw( 'G' );
  if (wroteJpdf) print().raw( 'J' );
  if (wroteStat) print().raw( 'P' );

  print().raw( '\n' );
}
