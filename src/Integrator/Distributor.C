//******************************************************************************
/*!
  \file      src/Integrator/Distributor.C
  \author    J. Bakosi
  \date      Fri 22 Aug 2014 10:41:17 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Distributor drives the time integration of differential equations
  \details   Distributor drives the time integration of differential equations
*/
//******************************************************************************

#include <Distributor.h>
#include <Integrator.h>
#include <DiffEqStack.h>
#include <quinoa.decl.h>

extern CProxy_Main mainProxy;

using quinoa::Distributor;

Distributor::Distributor( const ctr::CmdLine& cmdline ) :
  m_print( cmdline.get< tk::tag::verbose >() ? std::cout : std::clog ),
  m_ninit( 0 ),
  m_nstep( 0 ),
  m_nfinish( 0 ),
  m_it( 0 ),
  m_t( 0.0 )
//******************************************************************************
// Constructor
//! \author  J. Bakosi
//******************************************************************************
{
  // Print out info data layout
  m_print.list( "Particle properties data layout policy (CMake: LAYOUT)",
                std::list< std::string >{ ParProps().major() } );

  // Re-create differential equations stack for output
  DiffEqStack stack;

  // Print out information on factory
  m_print.eqlist( "Registered differential equations", stack.factory(),
                  stack.ntypes() );
  m_print.endpart();

  // Compute load distribution given total work and specified virtualization
  uint64_t chunksize, remainder;
  computeLoadDistribution( chunksize, remainder );

  // Print out information on problem
  m_print.part( "Problem" );

  // Print out info on problem title
  if ( !g_inputdeck.get< tag::title >().empty() )
    m_print.section( "Title", g_inputdeck.get< tag::title >() );

  // Print out info on settings of selected differential equations
  m_print.diffeqs( "Differential equations integrated", stack.info() );

  // Print out info on load distirubtion
  m_print.section( "Load distribution" );
  m_print.item( "Virtualization [0.0...1.0]",
                g_inputdeck.get< tag::cmd, tag::virtualization >() );
  m_print.item( "Load (number of particles)",
                g_inputdeck.get< tag::discr, tag::npar >() );
  m_print.item( "Number of processing elements", CkNumPes() );
  m_print.item( "Number of work units",
                std::to_string( m_numchares ) + " (" +
                std::to_string( m_numchares-1 ) + "x" +
                std::to_string( chunksize ) + "+" +
                std::to_string( chunksize+remainder ) + ")" );

  // Print out time integration header
  if (g_inputdeck.get< tag::discr, tag::nstep >()) {
    header();
    //statWriter().header();
  }

  // Start timer measuring total integration time
  m_timer.emplace_back();

  // Fire up integrators
  for (auto i=decltype(m_numchares){1}; i<m_numchares; ++i)
    m_proxy.push_back( CProxyInt::ckNew( thisProxy, chunksize ) );
  m_proxy.push_back( CProxyInt::ckNew( thisProxy, chunksize+remainder ) );
}

void
Distributor::computeLoadDistribution( uint64_t& chunksize, uint64_t& remainder )
//******************************************************************************
//  Compute load distribution for given total work and virtualization
//! \author J. Bakosi
//******************************************************************************
{
  // Compute load distibution (number of chares and chunksize) based on total
  // work (total number of particles) and virtualization

  // The virtualization parameter, specified by the user, is a real number
  // between 0.0 and 1.0, inclusive, which controls the degree of virtualization
  // or over-decomposition. Independent of the value of virtualization the work
  // is approximately evenly distributed among the available processing
  // elements. For zero virtualization (no over-decomposition), the work is
  // simply decomposed into total_work/numPEs, which yields the smallest number
  // of Charm++ chares and the largest chunks of work units. The other extreme
  // is unity virtualization, which decomposes the total work into the smallest
  // size work units possible, yielding the largest number of Charm++ chares.
  // Obviously, the optimum will be between 0.0 and 1.0, depending on the
  // problem.
  //
  // The formula below uses the simplest (linear) relationship between the
  // virtualization parameter and the number of work units with the extremes
  // described above. The formula is given by
  //
  // chunksize = (1 - n) * v + n;
  //
  // where
  //         n = npar/npes
  //      npar = number of particles, representing the total work
  //      npes = number of hardware processing elements

  // Get virtualization parameter
  const auto v = g_inputdeck.get< tag::cmd, tag::virtualization >();
  Assert( v > -std::numeric_limits< tk::real >::epsilon() &&
          v < 1.0+std::numeric_limits< tk::real >::epsilon(),
          "Virtualization parameter must be between [0.0...1.0]" );

  // Get total number of particles (represents total work)
  const auto npar = g_inputdeck.get< tag::discr, tag::npar >();

  // Query number of processing elements
  const auto npe = CkNumPes();

  // Compute minimum number of work units
  const auto n = npar/npe;

  // Compute work unit size based on the linear formula above
  chunksize = (1.0-n)*v + n;

  // Compute number of work units with size computed ignoring remainder
  m_numchares = npar/chunksize;

  // Compute remainder of work if the above number of units were to be created
  remainder = npar - m_numchares*chunksize;

  // Redistribute remainder among the work units for a more equal distribution
  chunksize += remainder/m_numchares;

  // Compute new remainder (after redistribution of the previous remainder)
  remainder = npar - m_numchares*chunksize;
}

void
Distributor::init()
//******************************************************************************
// Wait for all integrators to finish initialization
//! \author  J. Bakosi
//******************************************************************************
{
  // Increase number of integrators completing initialization
  ++m_ninit;

  // Wait for all integrators completing initialization
  if ( m_ninit == m_numchares ) {
    mainProxy.timestamp( "Initial conditions", m_timer[0].dsec() );
  }
}

void
Distributor::step()
//******************************************************************************
// Wait for all integrators to finish a step in time
//! \author  J. Bakosi
//******************************************************************************
{
  // Increase number integrators completing the time step
  ++m_nstep;

  // Wait for all integrators completing the time step and continue next one
  if ( m_nstep == m_numchares ) {
    // Echo one-liner info on time step
    report();
    // Continue with next time step on all integrators
    m_nstep = 0;
    for (auto& p : m_proxy) p.advance();
  }
}

void
Distributor::header() const
//******************************************************************************
// Print out time integration header
//! \author  J. Bakosi
//******************************************************************************
{
  m_print.inthead( "Time integration",
    ctr::MonteCarlo().name(g_inputdeck.get< tag::selected, tag::montecarlo >()),
    "Legend: it - iteration count\n"
    "         t - time\n"
    "        dt - time step size\n"
    "       ETE - estimated time elapsed (h:m:s)\n"
    "       ETA - estimated time for accomplishment (h:m:s)\n"
    "       out - output saved\n",
    "\n      it             t            dt        ETE        ETA   out\n"
      " ---------------------------------------------------------------\n" );
}

void
Distributor::report()
//******************************************************************************
// Print out one-liner report on time step
//! \author  J. Bakosi
//******************************************************************************
{
  //const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  const auto dt = g_inputdeck.get< tag::discr, tag::dt >();
  ++m_it;
  m_t += dt;

  if (!(m_it % g_inputdeck.get< tag::interval, tag::tty >())) {
    tk::Timer::Watch ete, eta;  // estimated time elapsed and for accomplishment
    m_timer[0].eta( g_inputdeck.get< tag::discr, tag::term >(), m_t,
                    g_inputdeck.get< tag::discr, tag::nstep >(), m_it,
                    ete, eta );

    m_print << std::setfill(' ') << std::setw(8) << m_it << "  "
            << std::scientific << std::setprecision(6)
            << std::setw(12) << m_t << "  "
            << g_inputdeck.get< tag::discr, tag::dt >() << "  "
            << std::setfill('0')
            << std::setw(3) << ete.hrs.count() << ":"
            << std::setw(2) << ete.min.count() << ":"
            << std::setw(2) << ete.sec.count() << "  "
            << std::setw(3) << eta.hrs.count() << ":"
            << std::setw(2) << eta.min.count() << ":"
            << std::setw(2) << eta.sec.count() << "  ";

//     if (wroteGlob) print().raw( 'G' );
//     if (wroteJpdf) print().raw( 'J' );
//     if (wroteStat) print().raw( 'P' );

    m_print << '\n';
  }
}

void
Distributor::finish()
//******************************************************************************
// Wait for all integrators to finish time stepping
//! \author  J. Bakosi
//******************************************************************************
{
  // Increase number integrators completing the time step
  ++m_nfinish;

  // Wait for all integrators completing time stepping
  if ( m_nfinish == m_numchares ) {
    // Echo one-liner info on last time step
    report();
    // Normal finish, print out reason
    const auto term = g_inputdeck.get< tag::discr, tag::term >();
    const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
    m_print.endsubsection();
    if (std::fabs(m_t - term) > std::numeric_limits<tk::real>::epsilon())
       m_print.note( "Normal finish, maximum number of iterations reached: " +
                     std::to_string( nstep ) );
     else 
       m_print.note( "Normal finish, maximum time reached: " +
                     std::to_string( term ) );
    // Quit
    mainProxy.finalize();
  }
}

#include <distributor.def.h>
