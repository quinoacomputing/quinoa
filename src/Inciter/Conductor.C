// *****************************************************************************
/*!
  \file      src/Inciter/Conductor.C
  \author    J. Bakosi
  \date      Mon 01 Aug 2016 08:27:40 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Conductor drives the time integration of a PDE
  \details   Conductor drives the time integration of a PDE
    The implementation uses the Charm++ runtime system and is fully asynchronous,
    overlapping computation, communication as well as I/O. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/Inciter/conductor.ci.
*/
// *****************************************************************************

#include <string>
#include <iostream>
#include <cstddef>
#include <unordered_set>

#include "Macro.h"
#include "Conductor.h"
#include "MeshNodes.h"
#include "PDEStack.h"
#include "ContainerUtil.h"
#include "LoadDistributor.h"
#include "ExodusIIMeshReader.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "DiagWriter.h"
#include "ParticleWriter.h"

#include "NoWarning/inciter.decl.h"

// Force the compiler to not instantiate the template below as it is
// instantiated in LinSys/LinSysMerger.C (only required on mac)
extern template class tk::LinSysMerger< inciter::CProxy_Conductor,
                                        inciter::CProxy_Performer >;
//extern template class tk::ParticleWriter< inciter::CProxy_Conductor >;

extern CProxy_Main mainProxy;

using inciter::Conductor;

Conductor::Conductor() :
  __dep(),
  m_print( g_inputdeck.get<tag::cmd,tag::verbose>() ? std::cout : std::clog ),
  m_nchare( 0 ),
  m_it( 0 ),
  m_t( g_inputdeck.get< tag::discr, tag::t0 >() ),
  m_dt( computedt() ),
  m_stage( 0 ),
  m_linsysmerger(),
  m_performer(),
  m_tracker(),
  m_particlewriter(),
  m_partitioner(),
  m_avcost( 0.0 ),
  m_npoin( 0 ),
  m_timer(),
  m_linsysbc(),
  m_diag( 5, 0.0 )   // <- WRONG FOR EVERYTHING BUT A SINGLE COMPNS
// *****************************************************************************
//  Constructor
//! \author J. Bakosi
// *****************************************************************************
{
  m_print.part( "Factory" );

  // Print out info data layout
  m_print.list( "Unknowns data layout (CMake: MESHNODE_DATA_LAYOUT)",
                std::list< std::string >{ tk::MeshNodes::layout() } );

  // Re-create partial differential equations stack for output
  PDEStack stack;

  // Print out information on factory
  m_print.eqlist( "Registered partial differential equations",
                  stack.factory(), stack.ntypes() );
  m_print.endpart();

  // Print out information on problem
  m_print.part( "Problem" );

  // Print out info on problem title
  if ( !g_inputdeck.get< tag::title >().empty() )
    m_print.title( g_inputdeck.get< tag::title >() );

  // Print out info on settings of selected partial differential equations
  m_print.pdes( "Partial differential equations integrated", stack.info() );

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

    // Enable SDAG waits
    wait4init();
    wait4report();

    // Print I/O filenames
    m_print.section( "Output filenames" );
    m_print.item( "Field", g_inputdeck.get< tag::cmd, tag::io, tag::output >()
                           + ".<chareid>" );
    m_print.item( "Diagnostics",
                  g_inputdeck.get< tag::cmd, tag::io, tag::diag >() );

    // Print output intervals
    m_print.section( "Output intervals" );
    m_print.item( "TTY", g_inputdeck.get< tag::interval, tag::tty>() );
    m_print.item( "Field", g_inputdeck.get< tag::interval, tag::field >() );
    m_print.item( "Diagnostics", g_inputdeck.get< tag::interval, tag::diag >() );
    m_print.endsubsection();

    // Output header for diagnostics output file
    tk::DiagWriter dw( g_inputdeck.get< tag::cmd, tag::io, tag::diag >(),
                       g_inputdeck.get< tag::flformat, tag::diag >(),
                       g_inputdeck.get< tag::prec, tag::diag >() );
    dw.header( { "r", "ru", "rv", "rw", "re" } );

    // Create (empty) worker array
    m_performer = PerformerProxy::ckNew();

    // Create (empty) particle tracker array
    m_tracker = TrackerProxy::ckNew();

    // Create ExodusII reader for reading side sets from file. When creating
    // LinSysMerger, er.readSideSets() reads all side sets from file, which is
    // a serial read, then send the same copy to all PEs. Performers then will
    // query the side sets from their local LinSysMerger branch.
    tk::ExodusIIMeshReader
      er( g_inputdeck.get< tag::cmd, tag::io, tag::input >() );

    // Create linear system merger chare group
    m_linsysmerger = LinSysMergerProxy::ckNew( thisProxy, m_performer,
                       er.readSidesets(),
                       g_inputdeck.get< tag::component >().nprop() );

    // Create particle writer Charm++ chare group
    m_particlewriter = ParticleWriterProxy::ckNew(
                         thisProxy,
                         g_inputdeck.get< tag::cmd, tag::io, tag::part >() );

    // Create mesh partitioner Charm++ chare group and start partitioning mesh
    m_print.diagstart( "Reading mesh graph ..." );
    m_partitioner = PartitionerProxy::ckNew( thisProxy, m_performer,
                                             m_linsysmerger, m_tracker,
                                             m_particlewriter );

  } else finish();      // stop if no time stepping requested
}

void
Conductor::load( uint64_t nelem )
// *****************************************************************************
// Reduction target indicating that all Partitioner chare groups have finished
// reading their part of the contiguously-numbered computational mesh graph and
// we are ready to compute the computational load
//! \param[in] nelem Total number of mesh elements (summed across all PEs)
//! \author J. Bakosi
// *****************************************************************************
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
  m_npoin = er.readHeader();
  m_print.item( "Number of nodes", m_npoin );

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
// *****************************************************************************
// Reduction target indicating that all Partitioner chare groups have finished
// setting up the necessary data structures for partitioning the computational
// mesh and we are ready for partitioning
//! \author J. Bakosi
// *****************************************************************************
{
  m_print.diagstart( "Partitioning and distributing mesh ..." );
  m_partitioner.partition( m_nchare );
}

void
Conductor::aveCost( tk::real c )
// *****************************************************************************
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
// *****************************************************************************
{
  m_print.diagend( "done" );    // "Partitioning and distributing mesh ...";
  // Compute average and broadcast it back to all partitioners (PEs)
  m_avcost = c / CkNumPes();
  m_partitioner.stdCost( m_avcost );
}

void
Conductor::stdCost( tk::real c )
// *****************************************************************************
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
// *****************************************************************************
{
  m_print.diag( "Linear system communication cost: avg = " +
                std::to_string( m_avcost ) + ", std = " +
                std::to_string( std::sqrt( c/CkNumPes() ) ) );
}

void
Conductor::rowcomplete()
// *****************************************************************************
// Reduction target indicating that all linear system merger branches have done
// their part of storing and exporting global row ids
//! \details This function is a Charm++ reduction target that is called when
//!   all linear system merger branches have done their part of storing and
//!   exporting global row ids. This is a necessary precondition to be done
//!   before we can issue a broadcast to all Performer chares to continue with
//!   the initialization step. The other, also necessary but by itself not
//!   sufficient, one is parcomplete(). Together rowcomplete() and
//!   parcomplete() are sufficient for continuing with the initialization. See
//!   also conductor.ci.
//! \author J. Bakosi
// *****************************************************************************
{
  m_linsysmerger.rowsreceived();
  trigger_row_complete();
}

void
Conductor::verifybc( CkReductionMsg* msg )
// *****************************************************************************
// Reduction target initiating verification of the boundary conditions set
//! \param[in] msg Serialized and concatenated vectors of BC nodelists
//! \details As a final step aggregating all node lists at which LinSysMerger
//!   will set boundary conditions, this function initiates verification that
//!   the aggregate boundary conditions to be set match the user's BCs. This
//!   function is a Charm++ reduction target that is called when all linear
//!   system merger branches have received from worker chares offering setting
//!   boundary conditions. This is the first step of the verification. First
//!   receives the aggregated node list at which LinSysMerger will set boundary
//!   conditions. Then we do a broadcast to all workers (Performers) to query
//!   the worker-owned node IDs on which a Dirichlet BC is set by the user (on
//!   any component of any PDEs integrated by the worker).
//! \author J. Bakosi
// *****************************************************************************
{
  // Deserialize final BC node list vector set by LinSysMerger
  PUP::fromMem creator( msg->getData() );
  creator | m_linsysbc;
  delete msg;

  // Issue broadcast querying the BCs set
  m_performer.requestBCs();
}

void
Conductor::doverifybc( CkReductionMsg* msg )
// *****************************************************************************
// Reduction target as a 2nd (final) of the verification of BCs
//! \param[in] msg Serialized and concatenated vectors of BC nodelists
//! \details This function finishes off the verification that the aggregate
//!   boundary conditions to be set by LinSysMerger objects match the user's
//!   BCs. This function is a Charm++ reduction target that is called by
//!   Performer::requestBCs() This is the last step of the verification, doing
//!   the actual verification of the BCs after receiving the aggregate BC node
//!   list from Performers querying user BCs from the input file.
//! \see Conductor::verifybc()
//! \author J. Bakosi
// *****************************************************************************
{
  // Deserialize final user BC node list vector
  PUP::fromMem creator( msg->getData() );
  std::vector< std::size_t > bc;
  creator | bc;
  delete msg;

  Assert( m_linsysbc == bc, "The boundary conditions node list to be set does "
          "not match the boundary condition node list based on the side sets "
          "given in the input file." );

  m_print.diag( "Boundary conditions verified (DEBUG only)" );

  m_linsysmerger.trigger_ver_complete();
}

void
Conductor::initcomplete()
// *****************************************************************************
//  Reduction target indicating that all Performer chares have finished their
//  initialization step and have already continued with start time stepping
//! \author J. Bakosi
// *****************************************************************************
{
  m_print.diag( "Starting time stepping ..." );
  header();   // print out time integration header
}

void
Conductor::diagnostics( tk::real* d, std::size_t n )
// *****************************************************************************
// Reduction target optionally collecting diagnostics, e.g., residuals, from all
// Performer chares
//! \param[in] d Diagnostics (sums) collected over all chares
//! \param[in] n Number of diagnostics in array d
//! \author J. Bakosi
// *****************************************************************************
{
  #ifdef NDEBUG
  IGNORE(n);
  #endif

  Assert( n == m_diag.size(),
          "Number of diagnostics contributed not equal to expected" );

  // Finish computing diagnostics, i.e., divide sums by the number of samples
  for (std::size_t i=0; i<m_diag.size(); ++i) m_diag[i] = d[i] / m_npoin;

  trigger_diag_complete();
}

void
Conductor::evaluateTime()
// *****************************************************************************
//  Evaluate time step: decide if it is time to quit
//! \author J. Bakosi
// *****************************************************************************
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
  }

  trigger_eval_complete();

  // if not final stage of time step or if neither max iterations nor max time
  // reached, will continue (by telling all linear system merger group
  // elements to prepare for a new rhs), otherwise finish
  if (m_stage < 1 || (std::fabs(m_t-term) > eps && m_it < nstep)) {
    m_linsysmerger.enable_wait4rhs();
    wait4report();      // re-enable SDAG wait for report
  } else finish();
}

void
Conductor::advance()
// *****************************************************************************
//  Reduction target indicating that all Performer chares have finished their
//  initialization step
//! \author J. Bakosi
// *****************************************************************************
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
// *****************************************************************************
// Compute size of next time step
//! \return Size of dt for the next time step
//! \author  J. Bakosi
// *****************************************************************************
{
  // Simply return the constant user-defined initial dt for now
  return g_inputdeck.get< tag::discr, tag::dt >();
}

void
Conductor::finish()
// *****************************************************************************
// Normal finish of time stepping
//! \author J. Bakosi
// *****************************************************************************
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
// *****************************************************************************
// Print out time integration header
//! \author J. Bakosi
// *****************************************************************************
{
  m_print.inthead( "Time integration", "Unstructured-mesh PDE solver testbed",
    "Legend: it - iteration count\n"
    "         t - time\n"
    "        dt - time step size\n"
    "       ETE - estimated time elapsed (h:m:s)\n"
    "       ETA - estimated time for accomplishment (h:m:s)\n"
    "       out - output-saved flags (F: field, D: diagnostics)\n",
    "\n      it             t            dt        ETE        ETA   out\n"
      " ---------------------------------------------------------------\n" );
  m_timer[ TimerTag::TIMESTEP ];
}

void
Conductor::report()
// *****************************************************************************
// Print out one-liner report on time step
//! \author J. Bakosi
// *****************************************************************************
{
  bool diag = false;

  // Append diagnostics file at selected times
  if (!(m_it % g_inputdeck.get< tag::interval, tag::diag >())) {
    tk::DiagWriter dw( g_inputdeck.get< tag::cmd, tag::io, tag::diag >(),
                       g_inputdeck.get< tag::flformat, tag::diag >(),
                       g_inputdeck.get< tag::prec, tag::diag >(),
                       std::ios_base::app );
    if (dw.diag( m_it, m_t, m_diag )) diag = true;
  }

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
    if (diag) m_print << 'D';

    m_print << '\n';
  }
}

#include "NoWarning/conductor.def.h"
