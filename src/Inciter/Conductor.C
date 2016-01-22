//******************************************************************************
/*!
  \file      src/Inciter/Conductor.C
  \author    J. Bakosi
  \date      Fri 22 Jan 2016 09:23:07 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Conductor drives the time integration of a PDE
  \details   Conductor drives the time integration of a PDE
    The implementation uses the Charm++ runtime system and is fully asynchronous,
    overlapping computation, communication as well as I/O. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/Inciter/conductor.ci.
*/
//******************************************************************************

#include <string>
#include <iostream>
#include <cstddef>
#include <unordered_set>

#include "Conductor.h"
#include "ContainerUtil.h"
#include "LoadDistributor.h"
#include "ExodusIIMeshReader.h"
#include "Inciter/InputDeck/InputDeck.h"

#include "inciter.decl.h"

extern CProxy_Main mainProxy;

using inciter::Conductor;

Conductor::Conductor() :
  m_print( g_inputdeck.get<tag::cmd,tag::verbose>() ? std::cout : std::clog ),
  m_it( 0 ),
  m_t( g_inputdeck.get< tag::discr, tag::t0 >() ),
  m_dt( computedt() ),
  m_stage( 0 )
//******************************************************************************
//  Constructor
//! \author J. Bakosi
//******************************************************************************
{
  // Print out information on problem
  m_print.part( "Problem" );

  // Print out info on problem title
  if ( !g_inputdeck.get< tag::title >().empty() )
    m_print.title( g_inputdeck.get< tag::title >() );

  // Print discretization parameters
  m_print.section( "Discretization parameters" );
  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  const auto t0 = g_inputdeck.get< tag::discr, tag::t0 >();
  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto dt = g_inputdeck.get< tag::discr, tag::dt >();
  m_print.item( "Number of time steps", nstep );
  m_print.item( "Start time", t0 );
  m_print.item( "Terminate time", term );
  m_print.item( "Initial time step size", dt );

  // If the desired max number of time steps is larger than zero, and the
  // termination time is larger than the initial time, and the initial time step
  // is smaller than the duration of the time to be simulated, we have work to
  // do, otherwise, finish right away
  if ( nstep != 0 && term > t0 && dt < term-t0 ) {

    // Print I/O filenames
    m_print.section( "Output filenames" );
    m_print.item( "Field", g_inputdeck.get< tag::cmd, tag::io, tag::output >()
                           + "_<PEid>" );

    // Print output intervals
    m_print.section( "Output intervals" );
    m_print.item( "TTY", g_inputdeck.get< tag::interval, tag::tty>() );
    m_print.item( "Field", g_inputdeck.get< tag::interval, tag::field >() );
    m_print.endsubsection();

    // Create (empty) worker array
    m_performer = PerformerProxy::ckNew();
    // Create linear system merger chare group
    m_linsysmerger = LinSysMergerProxy::ckNew( thisProxy, m_performer );
    // Create mesh partitioner Charm++ chare group and start partitioning mesh
    m_print.diagstart( "Reading mesh graph ..." );
    m_partitioner =
      PartitionerProxy::ckNew( thisProxy, m_performer, m_linsysmerger );

  } else finish();      // stop if no time stepping requested
}

void
Conductor::load( uint64_t nelem )
//******************************************************************************
// Reduction target indicating that all Partitioner chare groups have finished
// reading their part of the contiguously-numbered computational mesh graph and
// we are ready to compute the computational load
//! \param[in] nelem Total number of mesh elements (summed across all PEs)
//! \author J. Bakosi
//******************************************************************************
{
  m_print.diagend( "done" );    // "Reading mesh graph ..."

  // Compute load distribution given total work (nelem) and user-specified
  // virtualization
  uint64_t chunksize, remainder;
  m_nchare = static_cast<int>(
               tk::linearLoadDistributor(
                 g_inputdeck.get< tag::cmd, tag::virtualization >(),
                 nelem, CkNumPes(), chunksize, remainder ) );

  // Print out mesh graph stats
  m_print.section( "Input mesh graph statistics" );
  m_print.item( "Number of tetrahedra", nelem );
  tk::ExodusIIMeshReader er(g_inputdeck.get< tag::cmd, tag::io, tag::input >());
  auto npoin = er.readHeader();
  m_print.item( "Number of nodes", npoin );

  // Print out info on load distribution
  m_print.section( "Load distribution" );
  m_print.item( "Virtualization [0.0...1.0]",
                g_inputdeck.get< tag::cmd, tag::virtualization >() );
  m_print.item( "Load (number of tetrahedra)", nelem );
  m_print.item( "Number of processing elements", CkNumPes() );
  m_print.item( "Number of work units",
                std::to_string( m_nchare ) + " (" +
                std::to_string( m_nchare-1 ) + "*" +
                std::to_string( chunksize ) + "+" +
                std::to_string( chunksize+remainder ) + ")" );

  // Print out mesh partitioning configuration
  m_print.section( "Initial mesh partitioning" );
  m_print.Item< tk::ctr::PartitioningAlgorithm,
                tag::selected, tag::partitioner >();
  m_print.endsubsection();
}

void
Conductor::partition()
//******************************************************************************
// Reduction target indicating that all Partitioner chare groups have finished
// setting up the necessary data structures for partitioning the computational
// mesh and we are ready for partitioning
//! \author J. Bakosi
//******************************************************************************
{
  m_print.diagstart( "Partitioning and distributing mesh ..." );
  m_partitioner.partition( m_nchare );
}

void
Conductor::aveCost( tk::real c )
//******************************************************************************
// Reduction target estimating the average communication cost of merging the
// linear system
//! \param[in] c Communication cost summed across all PEs. The cost associated
//!   to a PE is a real number between 0 and 1, defined as the number of mesh
//!   points the PE does not own, i.e., needs to send to some other PE, divided
//!   by the total number of points the PE contributes to. The lower the better.
//! \details The average, computed here, gives an idea of the average
//!   communication cost across all PEs, while the standard deviation, computed
//!   by stdCost(), gives an idea on the expected load imbalance.
//! \author J. Bakosi
//******************************************************************************
{
  m_print.diagend( "done" );    // "Partitioning and distributing mesh ...";
  // Compute average and broadcast it back to all partitioners (PEs)
  m_avcost = c / CkNumPes();
  m_partitioner.stdCost( m_avcost );
}

void
Conductor::stdCost( tk::real c )
//******************************************************************************
// Reduction target estimating the standard deviation of the communication cost
// of merging the linear system
//! \param[in] c Sum of the squares of the communication cost minus the average,
//!   summed across all PEs. The cost associated to a PE is a real number
//!   between 0 and 1, defined as the number of mesh points the PE does not own,
//!   i.e., needs to send to some other PE, divided by the total number of
//!   points the PE contributes to. The lower the better.
//! \details The average, computed by avCost(), gives an idea of the average
//!   communication cost across all PEs, while the standard deviation, computed
//!   here, gives an idea on the expected load imbalance.
//! \author J. Bakosi
//******************************************************************************
{
  m_print.diag( "Linear system communication cost: avg = " +
                std::to_string( m_avcost ) + ", std = " +
                std::to_string( std::sqrt( c/CkNumPes() ) ) );
}

void
Conductor::rowcomplete()
//******************************************************************************
// Reduction target indicating that all linear system merger branches have done
// their part of storing and exporting global row ids
//! \details This function is a Charm++ reduction target that is called when
//!   all linear system merger branches have done their part of storing and
//!   exporting global row ids. Once this is done, we issue a broadcast to
//!   all Spawners and thus implicitly all Performer chares to continue with
//!   the initialization step.
//! \author J. Bakosi
//******************************************************************************
{
  m_linsysmerger.rowsreceived();
  m_performer.init( m_dt );
}

void
Conductor::initcomplete()
//******************************************************************************
//  Reduction target indicating that all Performer chares have finished their
//  initialization step and have already continued with start time stepping
//! \author J. Bakosi
//******************************************************************************
{
  m_print.diag( "Starting time stepping ..." );
  header();   // print out time integration header
}

void
Conductor::evaluateTime()
//******************************************************************************
//  Evaluate time step: decide if it is time to quit
//! \author J. Bakosi
//******************************************************************************
{
  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();

  // if at final stage of time step, finish time step just taken
  if (m_stage == 1) {
    // Increase number of iterations taken
    ++m_it;
    // Advance physical time to include time step just finished
    m_t += m_dt;
    // Truncate the size of last time step
    if (m_t > term) m_t = term;
    // Echo one-liner info on time step just taken
    report();
  }

  // if not final stage of time step or if neither max iterations nor max time
  // reached, will continue (by telling all linear system merger group
  // elements to prepare for a new rhs), otherwise, finish
  if (m_stage < 1 || (std::fabs(m_t-term) > eps && m_it < nstep))
    m_linsysmerger.enable_wait4rhs();
  else
    finish();
}

void
Conductor::advance()
//******************************************************************************
//  Reduction target indicating that all Performer chares have finished their
//  initialization step
//! \author J. Bakosi
//******************************************************************************
{
  // If not final stage, continue with next stage, increasing the stage. If
  // final stage, continue with next time step zeroing stage and using new time
  // step size.
  if (m_stage < 1)
    m_performer.advance( ++m_stage, m_dt, m_it, m_t );
  else
    m_performer.advance( m_stage=0, m_dt=computedt(), m_it, m_t );
}

tk::real
Conductor::computedt()
//******************************************************************************
// Compute size of next time step
//! \return Size of dt for the next time step
//! \author  J. Bakosi
//******************************************************************************
{
  // Simply return the constant user-defined initial dt for now
  return g_inputdeck.get< tag::discr, tag::dt >();
}

void
Conductor::finish()
//******************************************************************************
// Normal finish of time stepping
//! \author J. Bakosi
//******************************************************************************
{
  // Print out reason for stopping
  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  m_print.endsubsection();
  if (m_it >= nstep)
     m_print.note( "Normal finish, maximum number of iterations reached: " +
                   std::to_string( nstep ) );
   else
     m_print.note( "Normal finish, maximum time reached: " +
                   std::to_string( g_inputdeck.get<tag::discr,tag::term>() ) );

  // Quit
  mainProxy.finalize();
}

void
Conductor::header()
//******************************************************************************
// Print out time integration header
//! \author J. Bakosi
//******************************************************************************
{
  m_print.inthead( "Time integration", "Unstructured-mesh PDE solver testbed",
    "Legend: it - iteration count\n"
    "         t - time\n"
    "        dt - time step size\n"
    "       ETE - estimated time elapsed (h:m:s)\n"
    "       ETA - estimated time for accomplishment (h:m:s)\n"
    "       out - output-saved flags (F: field)\n",
    "\n      it             t            dt        ETE        ETA   out\n"
      " ---------------------------------------------------------------\n" );
  m_timer[ TimerTag::TIMESTEP ];
}

void
Conductor::report()
//******************************************************************************
// Print out one-liner report on time step
//! \author J. Bakosi
//******************************************************************************
{
  if (!(m_it % g_inputdeck.get< tag::interval, tag::tty >())) {

    // estimate time elapsed and time for accomplishment
    tk::Timer::Watch ete, eta;
    const auto& timer = tk::cref_find( m_timer, TimerTag::TIMESTEP );
    timer.eta( g_inputdeck.get< tag::discr, tag::term >() -
                 g_inputdeck.get< tag::discr, tag::t0 >(),
               m_t - g_inputdeck.get< tag::discr, tag::t0 >(),
               g_inputdeck.get< tag::discr, tag::nstep >(),
               m_it,
               ete,
               eta );

    // Output one-liner
    m_print << std::setfill(' ') << std::setw(8) << m_it << "  "
            << std::scientific << std::setprecision(6)
            << std::setw(12) << m_t << "  "
            << m_dt << "  "
            << std::setfill('0')
            << std::setw(3) << ete.hrs.count() << ":"
            << std::setw(2) << ete.min.count() << ":"
            << std::setw(2) << ete.sec.count() << "  "
            << std::setw(3) << eta.hrs.count() << ":"
            << std::setw(2) << eta.min.count() << ":"
            << std::setw(2) << eta.sec.count() << "  ";

    // Augment one-liner with output indicators
    if (!(m_it % g_inputdeck.get<tag::interval,tag::field>())) m_print << 'F';

    m_print << '\n';
  }
}

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "conductor.def.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif
