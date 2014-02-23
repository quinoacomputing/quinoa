//******************************************************************************
/*!
  \file      src/MonteCarlo/MonteCarlo.C
  \author    J. Bakosi
  \date      Sat 22 Feb 2014 06:23:42 PM MST
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

void
MonteCarlo::echoRNGs() const
//******************************************************************************
//  Standard echo RNGs
//! \author  J. Bakosi
//******************************************************************************
{
  print().section("Random number generators");
  print().MKLParams( control().get< tag::selected, tk::tag::rng >(),
                     control().get< tag::param, tk::tag::mklrng >() );
  print().RNGSSEParams( control().get< tag::selected, tk::tag::rng >(),
                        control().get< tag::param, tk::tag::rngsse >() );
}

void
MonteCarlo::echoIO() const
//******************************************************************************
//  Standard echo IO filenames
//! \author  J. Bakosi
//******************************************************************************
{
  print().subsection( "Output filenames" );
  print().item( "Input", control().get< tag::cmd, tag::io, tag::input >() );
  print().item( "Output", control().get< tag::cmd, tag::io, tag::output >() );
  print().item( "Glob", control().get< tag::cmd, tag::io, tag::glob >() );
  print().item( "Statistics", control().get< tag::cmd, tag::io, tag::stat >() );
  print().item( "PDF", control().get< tag::cmd, tag::io, tag::pdf >() );
  print().endsubsection();
}

void
MonteCarlo::echoIncpar() const
//******************************************************************************
//  Standard echo increment parameters
//! \author  J. Bakosi
//******************************************************************************
{
  print().subsection( "Increment parameters" );
  print().item( "Number of particles",
                control().get< tag::incpar, tag::npar >() );
  print().item( "Number of time steps",
                control().get< tag::incpar, tag::nstep >() );
  print().item( "Terminate time",
                control().get< tag::incpar, tag::term >() );
  print().item( "Initial time step size",
                control().get< tag::incpar, tag::dt >() );
  print().endsubsection();
}

void
MonteCarlo::echoIntervals() const
//******************************************************************************
//  Standard echo intervals
//! \author  J. Bakosi
//******************************************************************************
{
  print().subsection( "Output intervals" );
  print().item( "TTY", control().get< tag::interval, tag::tty>() );
  print().item( "Dump", control().get< tag::interval, tag::dump>() );
  print().item( "Glob", control().get< tag::interval, tag::glob >() );
  print().item( "Statistics", control().get< tag::interval, tag::plot >() );
  print().item( "PDF", control().get< tag::interval, tag::pdf >() );
  print().endsubsection();
}

void
MonteCarlo::echoStatistics() const
//******************************************************************************
//  Standard echo statistics
//! \author  J. Bakosi
//******************************************************************************
{
  print().subsection( "Statistics" );
  print().RequestedStats( "Requested" );
  print().EstimatedStats( "Estimated" );
}
