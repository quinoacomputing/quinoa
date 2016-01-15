//******************************************************************************
/*!
  \file      src/Inciter/Conductor.C
  \author    J. Bakosi
  \date      Fri 15 Jan 2016 08:26:02 AM MST
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
  m_stage( 0 ),
  m_gid( static_cast<std::size_t>(CkNumPes()) ),
  m_communication( m_gid.size() ),
  m_commbuilt( m_gid.size(), false ),
  m_nrecv( m_gid.size(), 0 ),
  m_start( m_gid.size() )
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

    // Activate SDAG waits
    wait4setup();

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
Conductor::flatten()
//******************************************************************************
// Reduction target indicating that all Partitioner chare groups have finished
// distributing mesh node IDs after partitioning and we are ready to start
// reordering mesh node IDs.
//! \author J. Bakosi
//******************************************************************************
{
  m_print.diagend( "done" );    // "Partitioning and distributing mesh ..."
  m_print.diagstart( "Prepare for reordering mesh nodes ..." );
  m_partitioner.flatten();
}

void
Conductor::addNodes( int pe, const std::vector< std::size_t >& gid )
//******************************************************************************
// Add global mesh node IDs from PE
//! \param[in] pe PE contribution coming from
//! \param[in] gid Global node indices resulting from the contributing PE
//!   reading its contiguously-numbered mesh elements from file
//! \author J. Bakosi
//******************************************************************************
{
  // Store node ids from pe
  for (auto p : gid) m_gid[ static_cast<std::size_t>(pe) ].insert( p );

  // Whenever a new list of node IDs has been received, attempt to build
  // communication maps for all PEs, even if not all PEs have contributed their
  // mesh node IDs. Since the node IDs arrive in arbitrary order, we cannot be
  // sure at this point that the node IDs for all PEs < PE have been received.
  // However, we don't wait for the nodes from all PEs to have been received to
  // attempt this. Building the map for a PE will only happen if the node IDs
  // for all PEs below the PE have already been received. See also buildComm().
  buildComm();

  // When the communication maps for all PEs have been computed, continue. Note
  // that PE 0 does not communicate, since there is no PEs with indices lower
  // than 0, hence we only test PEs > 0.
  if (std::all_of( std::next(m_commbuilt.cbegin()),
                   m_commbuilt.cend(),
                   [](const decltype(m_commbuilt)::value_type& m)
                   { return m; } ))
  {
    // Note that this assert should automatically be satisfied, since we can
    // only get here if the communication maps for all PEs have already been
    // built. This is ensured by the logic in buildComm(), which only sets an
    // entry in m_commbuilt to true if the corresponding entry in m_gid is not
    // empty. However, this was useful during development, so we leave it here.
    Assert( std::none_of( m_gid.cbegin(), m_gid.cend(),
                          [](const decltype(m_gid)::value_type& m)
                          { return m.empty(); } ),
            "Not all node IDs have been received" );

    // Compute node ID offset for all PEs. The offsets will be stored in vector
    // m_start. We use m_nrecv which, for each PE, contains the total number of
    // node IDs that a PE needs to receive from PEs with lower indices. Here we
    // compute the offset each PE will need to start assigning its new node IDs
    // from (for those nodes that are not assigned new IDs by any PEs with lower
    // indices). The offset for a PE is the offset for the previous PE plus the
    // number of node IDs the previous PE (uniquely) assigns new IDs for minus
    // the number of node IDs the previous PE receives for others.
    Assert( m_nrecv[0] == 0, "PE 0 must not receive any node IDs" );
    Assert( m_gid.size() == m_nrecv.size(), "Number of PE node IDs and number "
            "of to-be-received node IDs must be equal" );
    Assert( m_gid.size() == m_start.size(), "Number of PE node IDs and number "
            "of offsets must be equal" );
    m_start[0] = 0;
    for (std::size_t p=1; p<m_gid.size(); ++p)
      m_start[p] = m_start[p-1] + m_gid[p-1].size() - m_nrecv[p-1];

    m_print.diagend( "done" );  // "Prepare for reordering mesh nodes ..."
    // Continue with reordering, see also conductor.ci
    trigger_commmap_complete();
  }
}

void
Conductor::buildComm()
//******************************************************************************
// Attempt to build communication maps for for all PEs
//! \details This function attempts to build communication maps for all PEs. The
//!   container we store the maps is a vector of maps. The length of the vector
//!   is already set by constructor to be the number of PEs. Each vector entry
//!   is filled in by a map associating a vector global mesh node indices to a
//!   PE from which the PE (whose communication map we are building) will
//!   request new global node IDs from due to reordering. Note that only PEs
//!   with lower IDs than PE need to be considered here, since the global node
//!   reordering will start with PE 0 assigning new indices, followed by PE 1,
//!   assigning only to the nodes have not been encountered, and so on. Thus new
//!   IDs assigned by PEs lower than a given PE need to be communicated
//!   upstream. In other words, what we are doing here is collecting all mesh
//!   node IDs that a PE will contribute to that are also contributed to by PEs
//!   with lower IDs than PE.
//! \author J. Bakosi
//******************************************************************************
{
  // Only attempt to build the map for a PE if
  //  (1) we have not yet built the communication map for PE, and
  //  (2) we have received the node IDs for PE, and
  //  (4) we have the nodes from all PEs below PE.
  for (int pe=0; pe<static_cast<int>(m_gid.size()); ++pe) {
    const auto PE = static_cast< std::size_t >( pe );
    auto& cpe = m_communication[ PE ];  // attempt to build PE's comm map
    if (!m_commbuilt[PE] && !m_gid[PE].empty() && nodesComplete(pe)) {
      m_commbuilt[PE] = true;           // mark PE's communication map as built
      const auto& peid = m_gid[ PE ];   // global node IDs we are looing for
      for (auto i : peid) {             // try to find them all
        for (int p=0; p<pe; ++p) {      // on PEs < pe
          const auto& id = m_gid[ static_cast<std::size_t>(p) ];
          const auto it = id.find(i);
          if (it != end(id)) {          // node ID will be assigned by PE p
            cpe[ p ].insert( *it );
            ++m_nrecv[ PE ];            // count up to-be-received nodes for PE
            p = pe;      // get out, i.e., favor the lowest PE that has the node
          }
        }
      }
    }
  }
}

void
Conductor::reorder()
//******************************************************************************
// Start reordering mesh node IDs on all PEs
//! \author J. Bakosi
//******************************************************************************
{
  m_print.diagstart( "Reordering mesh nodes ..." );
  for (int p=0; p<CkNumPes(); ++p) {
    const auto P = static_cast< std::size_t >( p );
    m_partitioner[p].reorder( m_start[P], m_communication[P] );
  }
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
  m_print.diagend( "done" );    // "Reordering mesh nodes ...";
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
