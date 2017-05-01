// *****************************************************************************
/*!
  \file      src/Inciter/Transporter.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Transporter drives the time integration of transport equations
  \details   Transporter drives the time integration of transport equations.
    The implementation uses the Charm++ runtime system and is fully asynchronous,
    overlapping computation, communication as well as I/O. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/Inciter/transporter.ci.
*/
// *****************************************************************************

#include <string>
#include <iostream>
#include <cstddef>
#include <unordered_set>
#include <limits>
#include <cmath>

#include "Macro.h"
#include "Transporter.h"
#include "Fields.h"
#include "PDEStack.h"
#include "ContainerUtil.h"
#include "LoadDistributor.h"
#include "ExodusIIMeshReader.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "DiagWriter.h"

#include "NoWarning/inciter.decl.h"
#include "NoWarning/partitioner.decl.h"

// Force the compiler to not instantiate the template below as it is
// instantiated in LinSys/LinSysMerger.C (only required on mac)
extern template class tk::LinSysMerger< inciter::CProxy_Transporter,
                                        inciter::CProxy_Carrier,
                                        inciter::AuxSolverLumpMassDiff >;

// Force the compiler to not instantiate the template below as it is
// instantiated in Inciterer/Partitioner.C (only required with gcc 4.8.5)
extern template class
  inciter::Partitioner<
    inciter::CProxy_Transporter,
    inciter::CProxy_Carrier,
    tk::CProxy_LinSysMerger< inciter::CProxy_Transporter,
                             inciter::CProxy_Carrier,
                             inciter::AuxSolverLumpMassDiff >,
    tk::CProxy_ParticleWriter< inciter::CProxy_Transporter > >;

extern CProxy_Main mainProxy;

namespace inciter {

extern ctr::InputDeck g_inputdeck_defaults;
extern std::vector< PDE > g_pdes;

}

using inciter::Transporter;

Transporter::Transporter() :
  __dep(),
  m_print( g_inputdeck.get<tag::cmd,tag::verbose>() ? std::cout : std::clog ),
  m_nchare( 0 ),
  m_it( 0 ),
  m_t( g_inputdeck.get< tag::discr, tag::t0 >() ),
  m_dt( g_inputdeck.get< tag::discr, tag::dt >() ),
  m_stage( 0 ),
  m_linsysmerger(),
  m_carrier(),
  m_particlewriter(),
  m_partitioner(),
  m_avcost( 0.0 ),
  m_npoin( 0 ),
  m_timer(),
  m_linsysbc(),
  m_diag( g_inputdeck.get< tag::component >().nprop() * 2, 0.0 ),
  m_progPart( m_print, g_inputdeck.get< tag::cmd, tag::feedback >(),
              {{ "p", "d" }}, {{ CkNumPes(), CkNumPes() }} ),
  m_progGraph( m_print, g_inputdeck.get< tag::cmd, tag::feedback >(),
               {{ "g" }}, {{ CkNumPes() }} ),
  m_progReorder( m_print, g_inputdeck.get< tag::cmd, tag::feedback >(),
                 {{ "f", "m", "r", "b" }},
                 {{ CkNumPes(), CkNumPes(), CkNumPes(), CkNumPes() }} ),
  m_progSetup( m_print, g_inputdeck.get< tag::cmd, tag::feedback >(),
               {{ "r", "m", "b" }} ),
  m_progInit( m_print, g_inputdeck.get< tag::cmd, tag::feedback >(),
              {{ "i", "f", "l" }} ),
  m_progStep( m_print, g_inputdeck.get< tag::cmd, tag::feedback >(),
              {{ "r", "s", "l", "p" }} )
// *****************************************************************************
//  Constructor
//! \author J. Bakosi
// *****************************************************************************
{
  m_print.part( "Factory" );

  // Print out info data layout
  m_print.list( "Unknowns data layout (CMake: FIELD_DATA_LAYOUT)",
                std::list< std::string >{ tk::Fields::layout() } );

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
  const auto constdt = g_inputdeck.get< tag::discr, tag::dt >();
  const auto cfl = g_inputdeck.get< tag::discr, tag::cfl >();
  m_print.item( "Number of time steps", nstep );
  m_print.item( "Start time", t0 );
  m_print.item( "Terminate time", term );

  if (std::abs(constdt - g_inputdeck_defaults.get< tag::discr, tag::dt >()) >
        std::numeric_limits< tk::real >::epsilon())
    m_print.item( "Constant time step size", constdt );
  else if (std::abs(cfl - g_inputdeck_defaults.get< tag::discr, tag::cfl >()) >
             std::numeric_limits< tk::real >::epsilon())
    m_print.item( "CFL coefficient", cfl );

  m_print.item( "Mass diffusion coeff",
                g_inputdeck.get< tag::discr, tag::ctau >() );

  // If the desired max number of time steps is larger than zero, and the
  // termination time is larger than the initial time, and the constant time
  // step size (if that is used) is smaller than the duration of the time to be
  // simulated, we have work to do, otherwise, finish right away. If a constant
  // dt is not used, that part of the logic is always true as the default
  // constdt is zero, see inciter::ctr::InputDeck::InputDeck().
  if ( nstep != 0 && term > t0 && constdt < term-t0 ) {

    // Enable SDAG waits
    wait4eval();

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
    // This list is hard-coded here for CompFlow, so these are wrong for, e.g.,
    // scalar transport. This is punt for now and should be fixed so that these
    // column labels can be queried from the underlying PDEs, probably via a
    // a new polymorphic function of inciter::PDE.
    dw.header( { "r", "ru", "rv", "rw", "re",
                 "err(r)", "err(ru)", "err(rv)", "err(rw)", "err(re)" } );

    // Create (empty) worker array
    m_carrier = CarrierProxy::ckNew();

    // Create ExodusII reader for reading side sets from file. When creating
    // LinSysMerger, er.readSideSets() reads all side sets from file, which is
    // a serial read, then send the same copy to all PEs. Carriers then will
    // query the side sets from their local LinSysMerger branch.
    tk::ExodusIIMeshReader
      er( g_inputdeck.get< tag::cmd, tag::io, tag::input >() );

    // Read in side sets from file
    m_print.diagstart( "Reading side sets ..." );
    auto ss = er.readSidesets();
    m_print.diagend( "done" );

    // Verify that side sets to which boundary conditions are assigned by user
    // exist in mesh file
    std::unordered_set< int > conf;
    for (const auto& eq : g_pdes) eq.side( conf );
    for (auto i : conf)
      if (ss.find(i) == end(ss)) {
        m_print.diag( "WARNING: Boundary conditions specified on side set " +
          std::to_string(i) + " which does not exist in mesh file" );
        break;
      }

    // Create linear system merger chare group
    m_print.diag( "Creating linear system mergers" );
    m_linsysmerger = LinSysMergerProxy::ckNew( thisProxy, m_carrier, ss,
                       g_inputdeck.get< tag::component >().nprop(),
                       g_inputdeck.get< tag::cmd, tag::feedback >() );

    // Create particle writer Charm++ chare group. Note that by passing an empty
    // filename argument to the constructor, we tell the writer not to open a
    // file and not to perform I/O. To enable particle I/O, put in the filename
    // argument, commented out, instead of the empty string, and change the
    // number of particles (the constructor argument to m_particles) in the
    // initializer list of Carrier::Carrier(). This is basically a punt to
    // enable skipping H5Part I/O. Particles are a highly experimental feature
    // at this point.
    m_print.diag( "Creating particle writers" );
    m_particlewriter = ParticleWriterProxy::ckNew( thisProxy, "" );
                         //g_inputdeck.get< tag::cmd, tag::io, tag::part >() );

    // Create mesh partitioner Charm++ chare group and start partitioning mesh
    m_progGraph.start( "Creating partitioners and reading mesh graph ..." );
    m_partitioner = PartitionerProxy::ckNew( thisProxy, m_carrier,
                                             m_linsysmerger,
                                             m_particlewriter );

  } else finish();      // stop if no time stepping requested
}

void
Transporter::load( uint64_t nelem )
// *****************************************************************************
// Reduction target indicating that all Partitioner chare groups have finished
// reading their part of the contiguously-numbered computational mesh graph and
// we are ready to compute the computational load
//! \param[in] nelem Total number of mesh elements (summed across all PEs)
//! \author J. Bakosi
// *****************************************************************************
{
  m_progGraph.end();

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
  const auto ir = g_inputdeck.get< tag::selected, tag::initialamr >();
  if (ir == tk::ctr::InitialAMRType::UNIFORM)
    m_print.section( "Load distribution (before initial mesh refinement)" );
  else
    m_print.section( "Load distribution" );

  m_print.item( "Virtualization [0.0...1.0]",
                g_inputdeck.get< tag::cmd, tag::virtualization >() );
  m_print.item( "Load (number of tetrahedra)", nelem );
  m_print.item( "Number of processing elements", CkNumPes() );
  m_print.item( "Number of work units",
                std::to_string( m_nchare ) + " (" +
                std::to_string( m_nchare-1 ) + "*" +
                std::to_string( chunksize ) + "+" +
                std::to_string( chunksize+remainder ) + ')' );

  // Print out mesh partitioning configuration
  m_print.section( "Initial mesh partitioning" );
  m_print.Item< tk::ctr::PartitioningAlgorithm,
                tag::selected, tag::partitioner >();

  // Print out mesh refinement configuration and new mesh statistics
  if (ir == tk::ctr::InitialAMRType::UNIFORM) {
    m_print.section( "Mesh refinement" );
    m_print.Item< tk::ctr::InitialAMR, tag::selected, tag::initialamr >();
    m_print.item( "Final number of tetrahedra",
      std::to_string(nelem*8) + " (8*" + std::to_string(nelem) + ')' );
  }

  m_print.endsubsection();
}

void
Transporter::partition()
// *****************************************************************************
// Reduction target indicating that all Partitioner chare groups have finished
// setting up the necessary data structures for partitioning the computational
// mesh and we are ready for partitioning
//! \author J. Bakosi
// *****************************************************************************
{
  m_progPart.start( "Partitioning and distributing mesh ..." );
  m_partitioner.partition( m_nchare );
}

void
Transporter::distributed()
// *****************************************************************************
// Reduction target indicating that all Partitioner chare groups have finished
// distributing its global mesh node IDs and they are ready for preparing
// (flattening) their owned mesh node IDs for reordering
//! \author J. Bakosi
// *****************************************************************************
{
  m_progPart.end();
  m_progReorder.start( "Reordering mesh ..." );
  m_partitioner.flatten();
}

void
Transporter::aveCost( tk::real c )
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
  m_progReorder.end();
  m_print.diag( "Creating workers" );
  // Compute average and broadcast it back to all partitioners (PEs)
  m_avcost = c / CkNumPes();
  m_partitioner.stdCost( m_avcost );
}

void
Transporter::stdCost( tk::real c )
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
Transporter::setup()
// *****************************************************************************
// Reduction target indicating that all chare groups are ready for workers
//! \author J. Bakosi
// *****************************************************************************
{
  m_print.diag( "Reading mesh node coordinates, computing nodal volumes" );
  m_carrier.vol();
}

void
Transporter::volcomplete()
// *****************************************************************************
// Reduction target indicating that all Carriers have finished
// computing/receiving their part of the nodal volumes
//! \author J. Bakosi
// *****************************************************************************
{
  m_progSetup.start( "Computing row IDs, querying BCs, outputting mesh ...",
                     {{ CkNumPes(), m_nchare, CkNumPes() }} );
  m_carrier.setup();
}

void
Transporter::rowcomplete()
// *****************************************************************************
// Reduction target indicating that all linear system merger branches have done
// their part of storing and exporting global row ids
//! \details This function is a Charm++ reduction target that is called when
//!   all linear system merger branches have done their part of storing and
//!   exporting global row ids. This is a necessary precondition to be done
//!   before we can issue a broadcast to all Carrier chares to continue with
//!   the initialization step. The other, also necessary but by itself not
//!   sufficient, one is parcomplete(). Together rowcomplete() and
//!   parcomplete() are sufficient for continuing with the initialization. See
//!   also transporter.ci.
//! \author J. Bakosi
// *****************************************************************************
{
  m_progSetup.end();
  m_progInit.start( "Setting and outputting ICs, computing initial dt, "
                    "computing LHS ...",
                    {{ CkNumPes(), m_nchare, m_nchare }} );
  m_linsysmerger.rowsreceived();
  m_carrier.init();
}

void
Transporter::initcomplete()
// *****************************************************************************
//  Reduction target indicating that all Carrier chares have finished their
//  initialization step and have already continued with start time stepping
//! \author J. Bakosi
// *****************************************************************************
{
  m_progInit.end();
  m_print.diag( "Starting time stepping ..." );
  header();   // print out time integration header
  if ( g_inputdeck.get< tag::cmd, tag::feedback >() )
    m_progStep.start( "Time step stage ...",
      {{ m_nchare, CkNumPes(), m_nchare, m_nchare }} );
}

void
Transporter::diagnostics( tk::real* d, std::size_t n )
// *****************************************************************************
// Reduction target optionally collecting diagnostics, e.g., residuals, from all
// Carrier chares
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

  diag_complete();
}

void
Transporter::dt( tk::real* d, std::size_t n )
// *****************************************************************************
// Reduction target yielding a single minimum time step size across all workers
//! \param[in] d Minimum time step size collected over all chares
//! \param[in] n Size of data behind d
//! \author  J. Bakosi
// *****************************************************************************
{
  #ifdef NDEBUG
  IGNORE(n);
  #endif

  Assert( n == 1, "Size of min(dt) must be 1" );

  if (m_stage == 1 || m_it == 0) {
    // Use newly computed time step size
    m_dt = *d;
    // Truncate the size of last time step
    const auto term = g_inputdeck.get< tag::discr, tag::term >();
    if (m_t+m_dt > term) m_dt = term - m_t;;
  }

  // Advance to next time step stage
  m_carrier.advance( m_stage, m_dt, m_it, m_t );
}

void
Transporter::finish()
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
Transporter::header()
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
Transporter::evaluateTime()
// *****************************************************************************
// Evaluate time step and output one-liner report
//! \author J. Bakosi
// *****************************************************************************
{
  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();

  if ( g_inputdeck.get< tag::cmd, tag::feedback >() )
    m_progStep.end();

  if (m_stage < 1) {    // if at half-stage, simply go to next one

    ++m_stage;

  } else {      // if final stage of time step, finish time step just taken

    m_stage = 0;
    // Increase number of iterations taken
    ++m_it;
    // Advance physical time to include time step just finished
    m_t += m_dt;

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

      m_print << std::endl;
    }
  }

  wait4eval();

  // if neither max iterations nor max time reached, will continue (by telling
  // all linear system merger group elements to prepare for a new rhs),
  // otherwise finish
  if (std::fabs(m_t-term) > eps && m_it < nstep) {
    m_linsysmerger.enable_wait4rhs();
    if ( g_inputdeck.get< tag::cmd, tag::feedback >() )
      m_progStep.start( "Time step stage ..." );
  } else
    finish();
}

#include "NoWarning/transporter.def.h"
