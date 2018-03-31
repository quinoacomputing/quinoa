// *****************************************************************************
/*!
  \file      src/Inciter/Transporter.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
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
#include "UniPDF.h"
#include "PDFWriter.h"
#include "ContainerUtil.h"
#include "LoadDistributor.h"
#include "ExodusIIMeshReader.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "NodeDiagnostics.h"
#include "ElemDiagnostics.h"
#include "DiagWriter.h"

#include "NoWarning/inciter.decl.h"
#include "NoWarning/partitioner.decl.h"

extern CProxy_Main mainProxy;

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern ctr::InputDeck g_inputdeck_defaults;
extern std::vector< CGPDE > g_cgpde;
extern std::vector< DGPDE > g_dgpde;

}

using inciter::Transporter;

Transporter::Transporter() :
  m_print( g_inputdeck.get<tag::cmd,tag::verbose>() ? std::cout : std::clog ),
  m_nchare( 0 ),
  m_nelem( 0 ),
  m_chunksize( 0 ),
  m_remainder( 0 ),
  m_solver(),
  m_bc(),
  m_scheme( g_inputdeck.get< tag::selected, tag::scheme >() ),
  m_partitioner(),
  m_avcost( 0.0 ),
  m_V( 0.0 ),
  m_npoin( 0 ),
  m_minstat( {{ 0.0, 0.0 }} ),
  m_maxstat( {{ 0.0, 0.0 }} ),
  m_avgstat( {{ 0.0, 0.0 }} ),
  m_timer(),
  m_linsysbc(),
  m_progMesh( m_print, g_inputdeck.get< tag::cmd, tag::feedback >(),
               {{ "r", "f", "c" }}, {{ CkNumPes(), CkNumPes(), CkNumPes() }} ),
  m_progPart( m_print, g_inputdeck.get< tag::cmd, tag::feedback >(),
              {{ "p", "d" }}, {{ CkNumPes(), CkNumPes() }} ),
  m_progReorder( m_print, g_inputdeck.get< tag::cmd, tag::feedback >(),
                 {{ "f", "g", "q", "m", "r", "b" }},
                 {{ CkNumPes(), CkNumPes(), CkNumPes(), CkNumPes(),
                    CkNumPes(), CkNumPes() }} )
// *****************************************************************************
//  Constructor
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

  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  const auto t0 = g_inputdeck.get< tag::discr, tag::t0 >();
  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto constdt = g_inputdeck.get< tag::discr, tag::dt >();
  const auto cfl = g_inputdeck.get< tag::discr, tag::cfl >();
  const auto scheme = g_inputdeck.get< tag::selected, tag::scheme >();

  // Print discretization parameters
  m_print.section( "Discretization parameters" );
  m_print.Item< ctr::Scheme, tag::selected, tag::scheme >();
  if (scheme == ctr::SchemeType::MatCG || scheme == ctr::SchemeType::DiagCG) {
    auto fct = g_inputdeck.get< tag::discr, tag::fct >();
    m_print.item( "Flux-corrected transport (FCT)", fct );
    if (fct)
      m_print.item( "FCT mass diffusion coeff",
                    g_inputdeck.get< tag::discr, tag::ctau >() );
  }
  m_print.item( "Number of time steps", nstep );
  m_print.item( "Start time", t0 );
  m_print.item( "Terminate time", term );

  if (std::abs(constdt - g_inputdeck_defaults.get< tag::discr, tag::dt >()) >
        std::numeric_limits< tk::real >::epsilon())
    m_print.item( "Constant time step size", constdt );
  else if (std::abs(cfl - g_inputdeck_defaults.get< tag::discr, tag::cfl >()) >
             std::numeric_limits< tk::real >::epsilon())
    m_print.item( "CFL coefficient", cfl );

  // If the desired max number of time steps is larger than zero, and the
  // termination time is larger than the initial time, and the constant time
  // step size (if that is used) is smaller than the duration of the time to be
  // simulated, we have work to do, otherwise, finish right away. If a constant
  // dt is not used, that part of the logic is always true as the default
  // constdt is zero, see inciter::ctr::InputDeck::InputDeck().
  if ( nstep != 0 && term > t0 && constdt < term-t0 ) {

    // Enable SDAG waits
    wait4mesh();
    wait4stat();

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

    // Configure and write diagnostics file header
    diagHeader();

    // Create linear system solver group
    createSolver();

    // Create mesh partitioner AND boundary condition object group
    createPartitioner();

  } else finish();      // stop if no time stepping requested
}

void
Transporter::createSolver()
// *****************************************************************************
// Create linear solver
// *****************************************************************************
{
  // Create linear system solver callbacks
  std::vector< CkCallback > cbs {{
      CkCallback( CkReductionTarget(Transporter,comfinal), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,coord), thisProxy )
    , CkCallback( CkIndex_Transporter::diagnostics(nullptr), thisProxy )
  }};

  // Create linear system solver Charm++ chare group
  m_solver = tk::CProxy_Solver::
               ckNew( tk::CProxy_SolverShadow::ckNew(),
                      cbs,
                      g_inputdeck.get< tag::component >().nprop(),
                      g_inputdeck.get< tag::cmd, tag::feedback >() );
}

void
Transporter::createPartitioner()
// *****************************************************************************
// Create mesh partitioner AND boundary conditions group
// *****************************************************************************
{
  // Create ExodusII reader for reading side sets from file.
  tk::ExodusIIMeshReader er(g_inputdeck.get< tag::cmd, tag::io, tag::input >());

  // Read in side sets associated to mesh node IDs from file
  auto sidenodes = er.readSidesets();

  // Read side sets for boundary faces
  std::map< int, std::vector< std::size_t > > bface;

  std::vector< std::size_t > triinpoel;
  const auto scheme = g_inputdeck.get< tag::selected, tag::scheme >();

  // Read triangle boundary-face connectivity
  if (scheme == ctr::SchemeType::DG) {
    auto nbfac = er.readSidesetFaces( bface );
    er.readFaces( nbfac, triinpoel );
  }

  // Verify that side sets to which boundary conditions are assigned by user
  // exist in mesh file
  std::unordered_set< int > conf;
  for (const auto& eq : g_cgpde) eq.side( conf );
  for (auto i : conf)
  {
    if (sidenodes.find(i) == end(sidenodes)) {
      m_print.diag( "WARNING: Boundary conditions specified on side set " +
        std::to_string(i) + " which does not exist in mesh file" );
      break;
    }
  }

  // Create boundary conditions Charm++ chare group
  m_bc = inciter::CProxy_BoundaryConditions::ckNew( sidenodes );

  // Create partitioner callbacks
  std::vector< CkCallback > cbp {{
      CkCallback( CkReductionTarget(Transporter,refined), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,centroid), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,distributed), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,flattened), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,load), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,aveCost), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,stdCost), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,coord), thisProxy )
  }};

  // Start timer measuring preparation of the mesh for partitioning
  m_timer[ TimerTag::MESH_PREP ];

  // Create mesh partitioner Charm++ chare group and start preparing mesh
  m_progMesh.start( "Preparing mesh (read, optional refine, centroids) ..." );

  // Create mesh partitioner Charm++ chare group
  m_partitioner =
    CProxy_Partitioner::ckNew( cbp, thisProxy, m_solver, m_bc, m_scheme,
                               bface, triinpoel );
}

void
Transporter::diagHeader()
// *****************************************************************************
// Configure and write diagnostics file header
// *****************************************************************************
{
  // Output header for diagnostics output file
  tk::DiagWriter dw( g_inputdeck.get< tag::cmd, tag::io, tag::diag >(),
                     g_inputdeck.get< tag::flformat, tag::diag >(),
                     g_inputdeck.get< tag::prec, tag::diag >() );

  // Collect variables names for integral/diagnostics output
  std::vector< std::string > var;
  const auto scheme = g_inputdeck.get< tag::selected, tag::scheme >();
  if (scheme == ctr::SchemeType::MatCG || scheme == ctr::SchemeType::DiagCG)
    for (const auto& eq : g_cgpde) varnames( eq, var );
  else if (scheme == ctr::SchemeType::DG)
    for (const auto& eq : g_dgpde) varnames( eq, var );
  else Throw( "Diagnostics header not handled for discretization scheme" );

  const tk::ctr::Error opt;
  auto nv = var.size();
  std::vector< std::string > d;

  // Add 'L2(var)' for all variables as those are always computed
  const auto& l2name = opt.name( tk::ctr::ErrorType::L2 );
  for (std::size_t i=0; i<nv; ++i) d.push_back( l2name + '(' + var[i] + ')' );

  // Query user-requested diagnostics and augment diagnostics file header by
  // 'err(var)', where 'err' is the error type  configured, and var is the
  // variable computed, for all variables and all error types configured.
  const auto& err = g_inputdeck.get< tag::diag, tag::error >();
  for (const auto& e : err) {
    const auto& errname = opt.name( e );
    for (std::size_t i=0; i<nv; ++i)
      d.push_back( errname + '(' + var[i] + "-IC)" );
  }

  // Write diagnostics header
  dw.header( d );
}

void
Transporter::load( uint64_t nelem )
// *****************************************************************************
// Reduction target indicating that the mesh has been read from file
//! \details At this point all Partitioner chare groups have finished reading
//!   their part of the computational mesh and we are ready to compute the
//!   computational load
//! \param[in] nelem Total number of mesh elements (summed across all PEs)
// *****************************************************************************
{
  m_nelem = nelem;

  // Compute load distribution given total work (nelem) and user-specified
  // virtualization
  m_nchare = static_cast<int>(
               tk::linearLoadDistributor(
                 g_inputdeck.get< tag::cmd, tag::virtualization >(),
                 nelem, CkNumPes(), m_chunksize, m_remainder ) );

  // signal to runtime system that m_nchare is set
  load_complete();

  // Send total number of chares to all linear solver PEs, if they exist
  const auto scheme = g_inputdeck.get< tag::selected, tag::scheme >();
  if (scheme == ctr::SchemeType::MatCG || scheme == ctr::SchemeType::DiagCG)
    m_solver.nchare( m_nchare );
}

void
Transporter::centroid()
// *****************************************************************************
// Reduction target indicating that centroids have been computed all PEs
//! \details At this point all Partitioner chare groups have finished computing
//!   the cell centroids (f that was required for the mesh partitioner)
// *****************************************************************************
{
  centroid_complete();
}

void
Transporter::refined()
// *****************************************************************************
// Reduction target indicating that optional initial mesh refinement has been
// completed on all PEs
//! \details At this point all Partitioner chare groups have finished refining
//!   their mesh if that was requested by the user
// *****************************************************************************
{
  refine_complete();
}

void
Transporter::partition()
// *****************************************************************************
// Start partitioning the mesh
// *****************************************************************************
{
  m_progMesh.end();

  // Start timer measuring preparation of the mesh for partitioning
  const auto& timer = tk::cref_find( m_timer, TimerTag::MESH_PREP );
  m_print.diag( "Mesh preparation time: " + std::to_string( timer.dsec() ) +
                " sec" );

  // Print out mesh graph stats
  m_print.section( "Input mesh graph statistics" );
  m_print.item( "Number of tetrahedra", m_nelem );
  tk::ExodusIIMeshReader er(g_inputdeck.get< tag::cmd, tag::io, tag::input >());
  m_npoin = er.readHeader();
  m_print.item( "Number of nodes", m_npoin );

  // Print out info on load distribution
  const auto ir = g_inputdeck.get< tag::amr, tag::init >();
  if (!ir.empty())
    m_print.section( "Load distribution (before initial mesh refinement)" );
  else
    m_print.section( "Load distribution" );

  m_print.item( "Virtualization [0.0...1.0]",
                g_inputdeck.get< tag::cmd, tag::virtualization >() );
  m_print.item( "Load (number of tetrahedra)", m_nelem );
  m_print.item( "Number of processing elements", CkNumPes() );
  m_print.item( "Number of work units",
                std::to_string( m_nchare ) + " (" +
                std::to_string( m_nchare-1 ) + "*" +
                std::to_string( m_chunksize ) + "+" +
                std::to_string( m_chunksize+m_remainder ) + ')' );

  // Print out mesh partitioning configuration
  m_print.section( "Initial mesh partitioning" );
  m_print.Item< tk::ctr::PartitioningAlgorithm,
                tag::selected, tag::partitioner >();

  // Print out adaptive mesh refinement configuration
  const auto amr = g_inputdeck.get< tag::amr, tag::amr >();
  if (amr) {
    m_print.section( "Adaptive mesh refinement (AMR)" );
    m_print.ItemVec< ctr::AMRInitial >
                   ( g_inputdeck.get< tag::amr, tag::init >() );
    m_print.item( "Initial uniform levels",
                  g_inputdeck.get< tag::amr, tag::levels >() );
    m_print.Item< ctr::AMRError, tag::amr, tag::error >();
    // Print out initially refined  mesh statistics
    if (!ir.empty()) {
      m_print.section( "Initial mesh refinement" );
      m_print.item( "Final number of tetrahedra", "..." );
    }
  }

  m_print.endsubsection();

  m_progPart.start( "Partitioning and distributing mesh ..." );
  m_partitioner.partition( m_nchare );
}

void
Transporter::distributed()
// *****************************************************************************
// Reduction target indicating that all Partitioner chare groups have finished
// distributing its global mesh node IDs and they are ready for preparing
// (flattening) their owned mesh node IDs for reordering
// *****************************************************************************
{
  m_progPart.end();
  m_progReorder.start( "Reordering mesh (flatten, gather, query, mask, "
                       "reorder, bounds) ... " );
  m_partitioner.flatten();
}

void
Transporter::flattened()
// *****************************************************************************
// Reduction target indicating that all Partitioner chare groups have finished
// flattening its global mesh node IDs and they are ready for computing the
// communication maps required for node ID reordering
// *****************************************************************************
{
  m_partitioner.gather();
}

void
Transporter::aveCost( tk::real c )
// *****************************************************************************
// Reduction target estimating the average communication among all PEs
//! \param[in] c Communication cost summed across all PEs. The cost associated
//!   to a PE is a real number between 0 and 1, defined as the number of mesh
//!   points the PE does not own, i.e., needs to send to some other PE, divided
//!   by the total number of points the PE contributes to. The lower the better.
//! \details The average, computed here, gives an idea of the average
//!   communication cost across all PEs, while the standard deviation, computed
//!   by stdCost(), gives an idea on the expected load imbalance.
// *****************************************************************************
{
  m_progReorder.end();

  // Compute average and broadcast it back to all partitioners (PEs)
  m_avcost = c / CkNumPes();
  m_partitioner.stdCost( m_avcost );
}

void
Transporter::stdCost( tk::real c )
// *****************************************************************************
// Reduction target estimating the standard deviation of the communication cost
// acrosss all PEs
//! \param[in] c Sum of the squares of the communication cost minus the average,
//!   summed across all PEs. The cost associated to a PE is a real number
//!   between 0 and 1, defined as the number of mesh points the PE does not own,
//!   i.e., needs to send to some other PE, divided by the total number of
//!   points the PE contributes to. The lower the better.
//! \details The average, computed by avCost(), gives an idea of the average
//!   communication cost across all PEs, while the standard deviation, computed
//!   here, gives an idea on the expected load imbalance.
// *****************************************************************************
{
  m_print.diag( "Communication cost: avg = " + std::to_string( m_avcost ) +
                ", std = " + std::to_string( std::sqrt( c/CkNumPes() ) ) );
}

void
Transporter::coord()
// *****************************************************************************
// Reduction target indicating that all chare groups are ready for workers to
// start reading their mesh node coordinates
// *****************************************************************************
{
  // Tell the runtime system that every PE is done with dynamically inserting
  // Discretization chare array elements
  m_scheme.doneDiscInserting< tag::bcast >();

  // Tell the runtime system that every PE is done with dynamically inserting
  // Discretization chare array elements
  auto sch = g_inputdeck.get< tag::selected, tag::scheme >();
  if (sch == ctr::SchemeType::MatCG || sch == ctr::SchemeType::DiagCG)
    m_scheme.doneDistFCTInserting< tag::bcast >();

  m_scheme.coord< tag::bcast >();
}

void
Transporter::comfinal()
// *****************************************************************************
// Reduction target indicating that the communication has been established among
// PEs
// *****************************************************************************
{
  // Tell the runtime system that every PE is done with dynamically inserting
  // Discretization worker (MatCG, DiagCG, DG, ...) chare array elements
  m_scheme.doneInserting< tag::bcast >();
  com_complete();
}

void
Transporter::vol()
// *****************************************************************************
// Reduction target indicating that all workers have finished
// computing/receiving their part of the nodal volumes
// *****************************************************************************
{
  m_scheme.totalvol< tag::bcast >();
}

void
Transporter::totalvol( tk::real v )
// *****************************************************************************
// Reduction target summing total mesh volume across all workers
//! \param[in] v mesh volume
// *****************************************************************************
{
  m_V = v;
  m_partitioner.createWorkers();  // create "derived" workers (e.g., DG)
  m_scheme.stat< tag::bcast >();
}

void
Transporter::minstat( tk::real* d, std::size_t n )
// *****************************************************************************
// Reduction target yielding minimum mesh statistcs across all workers
//! \param[in] d Minimum mesh statistics collected over all chares
//! \param[in] n Size of data behind d
// *****************************************************************************
{
  #ifdef NDEBUG
  IGNORE(n);
  #endif

  Assert( n == m_minstat.size(),
          "Size of min(stat) must be " + std::to_string(m_minstat.size()) );

  m_minstat[0] = d[0];  // minimum edge length
  m_minstat[1] = d[1];  // minimum cell volume cubic root

  minstat_complete();
}

void
Transporter::maxstat( tk::real* d, std::size_t n )
// *****************************************************************************
// Reduction target yielding the maximum mesh statistics across all workers
//! \param[in] d Maximum mesh statistics collected over all chares
//! \param[in] n Size of data behind d
// *****************************************************************************
{
  #ifdef NDEBUG
  IGNORE(n);
  #endif

  Assert( n == m_maxstat.size(),
          "Size of max(stat) must be " + std::to_string(m_maxstat.size()) );

  m_maxstat[0] = d[0];  // maximum edge length
  m_maxstat[1] = d[1];  // maximum cell volume cubic root

  maxstat_complete();
}

void
Transporter::sumstat( tk::real* d, std::size_t n )
// *****************************************************************************
// Reduction target yielding the sum mesh statistics across all workers
//! \param[in] d Sum mesh statistics collected over all chares
//! \param[in] n Size of data behind d
// *****************************************************************************
{
  #ifdef NDEBUG
  IGNORE(n);
  #endif

  Assert( n == 2*m_avgstat.size(),
          "Size of sum(stat) must be " + std::to_string(2*m_avgstat.size()) );

  m_avgstat[0] = d[1] / d[0];      // average edge length
  m_avgstat[1] = d[3] / d[2];      // average cell volume cubic root

  sumstat_complete();
}

void
Transporter::pdfstat( CkReductionMsg* msg )
// *****************************************************************************
// Reduction target yielding PDF of mesh statistics across all workers
//! \param[in] msg Serialized PDF
// *****************************************************************************
{
  std::vector< tk::UniPDF > pdf;

  // Deserialize final PDF
  PUP::fromMem creator( msg->getData() );
  creator | pdf;
  delete msg;

  // Create new PDF file (overwrite if exists)
  tk::PDFWriter pdfe( "mesh_edge_pdf.txt" );
  // Output edgelength PDF
  pdfe.writeTxt( pdf[0], tk::ctr::PDFInfo{ {"PDF"}, {}, {"edgelength"} } );

  // Create new PDF file (overwrite if exists)
  tk::PDFWriter pdfv( "mesh_vol_pdf.txt" );
  // Output cell volume cubic root PDF
  pdfv.writeTxt( pdf[1], tk::ctr::PDFInfo{ {"PDF"}, {}, {"V^{1/3}"} } );

  pdfstat_complete();
}

void
Transporter::stat()
// *****************************************************************************
// Echo diagnostics mesh statistics
// *****************************************************************************
{
  m_print.diag( "Mesh statistics: min/max/avg(edgelength) = " +
                std::to_string( m_minstat[0] ) + " / " +
                std::to_string( m_maxstat[0] ) + " / " +
                std::to_string( m_avgstat[0] ) );
  m_print.diag( "Mesh statistics: min/max/avg(V^{1/3}) = " +
                std::to_string( m_minstat[1] ) + " / " +
                std::to_string( m_maxstat[1] ) + " / " +
                std::to_string( m_avgstat[1] ) );

  m_print.inthead( "Time integration", "Unstructured-mesh PDE solver testbed",
     "Legend: it - iteration count\n"
     "         t - time\n"
     "        dt - time step size\n"
     "       ETE - estimated time elapsed (h:m:s)\n"
     "       ETA - estimated time for accomplishment (h:m:s)\n"
     "       out - output-saved flags (F: field, D: diagnostics)\n",
     "\n      it             t            dt        ETE        ETA   out\n"
       " ---------------------------------------------------------------\n" );

  m_scheme.setup( m_V );
}

void
Transporter::start()
// *****************************************************************************
// Start time stepping
//! \note Only called if MatCG/DiagG is used
// *****************************************************************************
{
  m_scheme.dt< tag::bcast >();
}

void
Transporter::diagnostics( CkReductionMsg* msg )
// *****************************************************************************
// Reduction target optionally collecting diagnostics, e.g., residuals
//! \param[in] msg Serialized diagnostics vector aggregated across all PEs
// *****************************************************************************
{
  std::vector< std::vector< tk::real > > d;

  // Deserialize diagnostics vector
  PUP::fromMem creator( msg->getData() );
  creator | d;
  delete msg;

  auto ncomp = g_inputdeck.get< tag::component >().nprop();

  Assert( ncomp > 0, "Number of scalar components must be positive");
  Assert( d.size() == NUMDIAG, "Diagnostics vector size mismatch" );

  for (std::size_t i=0; i<d.size(); ++i)
     Assert( d[i].size() == ncomp,
             "Size mismatch at final stage of diagnostics aggregation" );

  // Allocate storage for those diagnostics that are always computed
  std::vector< tk::real > diag( ncomp, 0.0 );

  // Finish computing diagnostics
  for (std::size_t i=0; i<d[L2SOL].size(); ++i)
    diag[i] = sqrt( d[L2SOL][i] / m_V );
  
  // Query user-requested error types to output
  const auto& error = g_inputdeck.get< tag::diag, tag::error >();

  decltype(ncomp) n = 0;
  for (const auto& e : error) {
    n += ncomp;
    if (e == tk::ctr::ErrorType::L2) {
      // Finish computing the L2 norm of the numerical - analytical solution
     for (std::size_t i=0; i<d[L2ERR].size(); ++i)
       diag.push_back( sqrt( d[L2ERR][i] / m_V ) );
    } else if (e == tk::ctr::ErrorType::LINF) {
      // Finish computing the Linf norm of the numerical - analytical solution
      for (std::size_t i=0; i<d[LINFERR].size(); ++i)
        diag.push_back( d[LINFERR][i] );
    }
  }

  // Append diagnostics file at selected times
  tk::DiagWriter dw( g_inputdeck.get< tag::cmd, tag::io, tag::diag >(),
                     g_inputdeck.get< tag::flformat, tag::diag >(),
                     g_inputdeck.get< tag::prec, tag::diag >(),
                     std::ios_base::app );
  dw.diag( static_cast<uint64_t>(d[ITER][0]), d[TIME][0], d[DT][0], diag );

  // Evaluate whther to continue with next step
  m_scheme.eval< tag::bcast >();
}

void
Transporter::finish()
// *****************************************************************************
// Normal finish of time stepping
// *****************************************************************************
{
  mainProxy.finalize();
}

#include "NoWarning/transporter.def.h"
