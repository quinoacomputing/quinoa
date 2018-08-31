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
#include "MeshReader.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "NodeDiagnostics.h"
#include "ElemDiagnostics.h"
#include "DiagWriter.h"
#include "Callback.h"

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
  m_solver(),
  m_scheme( g_inputdeck.get< tag::discr, tag::scheme >() ),
  m_partitioner(),
  m_refiner(),
  m_sorter(),
  m_nelem( 0 ),
  m_npoin( 0 ),
  m_V( 0.0 ),
  m_minstat( {{ 0.0, 0.0, 0.0 }} ),
  m_maxstat( {{ 0.0, 0.0, 0.0 }} ),
  m_avgstat( {{ 0.0, 0.0, 0.0 }} ),
  m_timer(),
  m_linsysbc(),
  m_progMesh( m_print, g_inputdeck.get< tag::cmd, tag::feedback >(),
              {{ "p", "d", "r", "b", "c", "m", "r", "b" }},
              {{ "partition", "distribute", "refine", "bnd", "comm", "mask",
                  "reorder", "bounds"}} ),
  m_progWork( m_print, g_inputdeck.get< tag::cmd, tag::feedback >(),
              {{ "c", "b", "f", "g", "a" }},
              {{ "create", "bndface", "comfac", "ghost", "adj" }} )
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
  const auto scheme = g_inputdeck.get< tag::discr, tag::scheme >();

  // Print discretization parameters
  m_print.section( "Discretization parameters" );
  m_print.Item< ctr::Scheme, tag::discr, tag::scheme >();
  if (scheme == ctr::SchemeType::MatCG || scheme == ctr::SchemeType::DiagCG) {
    auto fct = g_inputdeck.get< tag::discr, tag::fct >();
    m_print.item( "Flux-corrected transport (FCT)", fct );
    if (fct)
      m_print.item( "FCT mass diffusion coeff",
                    g_inputdeck.get< tag::discr, tag::ctau >() );
  } else if (scheme == ctr::SchemeType::DG) {
    m_print.Item< ctr::Flux, tag::discr, tag::flux >();
  }
  m_print.item( "PE-locality mesh reordering",
                g_inputdeck.get< tag::discr, tag::reorder >() );
  m_print.item( "Number of time steps", nstep );
  m_print.item( "Start time", t0 );
  m_print.item( "Terminate time", term );

  if (std::abs(constdt - g_inputdeck_defaults.get< tag::discr, tag::dt >()) >
        std::numeric_limits< tk::real >::epsilon())
    m_print.item( "Constant time step size", constdt );
  else if (std::abs(cfl - g_inputdeck_defaults.get< tag::discr, tag::cfl >()) >
             std::numeric_limits< tk::real >::epsilon())
    m_print.item( "CFL coefficient", cfl );

  // Print out adaptive mesh refinement configuration
  const auto amr = g_inputdeck.get< tag::amr, tag::amr >();
  if (amr) {
    m_print.section( "Adaptive mesh refinement (AMR)" );
    m_print.refvar( g_inputdeck.get< tag::amr, tag::refvar >(),
                    g_inputdeck.get< tag::amr, tag::id >() );
    m_print.Item< ctr::AMRError, tag::amr, tag::error >();
    auto initamr = g_inputdeck.get< tag::amr, tag::initamr >();
    m_print.item( "Initial refinement", initamr );
    if (initamr) {
      const auto& initref = g_inputdeck.get< tag::amr, tag::init >();
      m_print.item( "Initial refinement steps", initref.size() );
      m_print.ItemVec< ctr::AMRInitial >( initref );
      m_print.edgeref( g_inputdeck.get< tag::amr, tag::edge >() );
    }
  }

  // If the desired max number of time steps is larger than zero, and the
  // termination time is larger than the initial time, and the constant time
  // step size (if that is used) is smaller than the duration of the time to be
  // simulated, we have work to do, otherwise, finish right away. If a constant
  // dt is not used, that part of the logic is always true as the default
  // constdt is zero, see inciter::ctr::InputDeck::InputDeck().
  if ( nstep != 0 && term > t0 && constdt < term-t0 ) {

    // Enable SDAG waits
    thisProxy.wait4stat();

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
  tk::SolverCallback cbs{
      CkCallback( CkReductionTarget(Transporter,partition), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,bounds), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,comfinal), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,disccreated), thisProxy )
  };

  // Create linear system solver Charm++ chare group
  m_solver = tk::CProxy_Solver::
               ckNew( tk::CProxy_SolverShadow::ckNew(),
                      cbs,
                      g_inputdeck.get< tag::component >().nprop() );
}

void
Transporter::createPartitioner()
// *****************************************************************************
// Create mesh partitioner AND boundary conditions group
// *****************************************************************************
{
  // Create mesh reader for reading side sets from file
  tk::MeshReader mr( g_inputdeck.get< tag::cmd, tag::io, tag::input >() );

  std::map< int, std::vector< std::size_t > > belem;
  std::map< int, std::vector< std::size_t > > faces;
  std::map< int, std::vector< std::size_t > > bnode;

  // Read boundary (side set) data from input file
  const auto scheme = g_inputdeck.get< tag::discr, tag::scheme >();
  if (scheme == ctr::SchemeType::DG) {
    // Read boundary-face connectivity on side sets
    mr.readSidesetFaces( belem, faces );
    // Verify boundarty condition (BC) side sets used exist in mesh file
    matchBCs( g_dgpde, belem );
  } else {
    // Read node lists on side sets
    bnode = mr.readSidesetNodes();
    // Verify boundarty condition (BC) side sets used exist in mesh file
    matchBCs( g_cgpde, bnode );
  }

  // Create partitioner callbacks (order matters)
  tk::PartitionerCallback cbp {
      CkCallback( CkReductionTarget(Transporter,load), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,distributed), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,refinserted), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,refined), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,flattened), thisProxy )
  };

  // Create refiner callbacks (order matters)
  tk::RefinerCallback cbr {
      CkCallback( CkReductionTarget(Transporter,matched), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,refined), thisProxy )
  };

  // Create sorter callbacks (order matters)
  tk::SorterCallback cbs {
      CkCallback( CkReductionTarget(Transporter,flattened), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,discinserted), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,workinserted), thisProxy )
  };

  // Start timer measuring preparation of the mesh for partitioning
  m_timer[ TimerTag::MESH_READ ];

  // Create mesh partitioner Charm++ chare group and start preparing mesh
  m_print.diag( "Reading mesh" );

  // Create empty mesh sorter Charm++ chare array
  m_sorter = CProxy_Sorter::ckNew();

  // Create empty mesh refiner Charm++ chare array
  m_refiner = CProxy_Refiner::ckNew();

  // Create mesh partitioner Charm++ chare group
  m_partitioner =
    CProxy_Partitioner::ckNew( cbp, cbr, cbs, thisProxy, m_solver, m_refiner,
                               m_sorter, m_scheme, belem, faces, bnode );
}

void
Transporter::load( uint64_t nelem, uint64_t npoin )
// *****************************************************************************
// Reduction target: the mesh has been read from file on all PEs
//! \param[in] nelem Total number of mesh elements (summed across all PEs)
//! \param[in] npoin Total number of mesh points (summed across all PEs)
// *****************************************************************************
{
  // Compute load distribution given total work (nelem) and user-specified
  // virtualization
  uint64_t chunksize, remainder;
  m_nchare = static_cast<int>(
               tk::linearLoadDistributor(
                 g_inputdeck.get< tag::cmd, tag::virtualization >(),
                 nelem, CkNumPes(), chunksize, remainder ) );

  // Send total number of chares to all linear solver PEs
  m_solver.nchare( m_nchare );

  // Start timer measuring preparation of the mesh for partitioning
  const auto& timer = tk::cref_find( m_timer, TimerTag::MESH_READ );
  m_print.diag( "Mesh read time: " + std::to_string( timer.dsec() ) + " sec" );

  // Print out mesh graph stats
  m_print.section( "Input mesh graph statistics" );
  m_print.item( "Number of tetrahedra", nelem );
  m_print.item( "Number of nodes", npoin );

  // Print out mesh partitioning configuration
  m_print.section( "Mesh partitioning" );
  m_print.Item< tk::ctr::PartitioningAlgorithm,
                tag::selected, tag::partitioner >();

  // Print out info on load distribution
  m_print.section( "Initial load distribution" );
  m_print.item( "Virtualization [0.0...1.0]",
                g_inputdeck.get< tag::cmd, tag::virtualization >() );
  m_print.item( "Number of tetrahedra", nelem );
  m_print.item( "Number of processing elements",
                std::to_string( CkNumPes() ) + " (" +
                std::to_string( CkNumNodes() ) + 'x' +
                std::to_string( CkNumPes()/CkNumNodes() ) + ')' );
  m_print.item( "Number of work units", m_nchare );

  m_print.endsubsection();

  // Query number of initial mesh refinement steps
  int nref = 0;
  if (g_inputdeck.get< tag::amr, tag::initamr >())
    nref = static_cast<int>( g_inputdeck.get< tag::amr, tag::init >().size() );

  m_progMesh.start( "Preparing mesh", {{ CkNumPes(), CkNumPes(), nref,
    m_nchare, m_nchare, m_nchare, m_nchare, m_nchare }} );
}

void
Transporter::partition()
// *****************************************************************************
// Reduction target: Reduction target: all Solver (PEs) have computed the number
// of chares they will recieve contributions from during linear solution
// *****************************************************************************
{
  m_partitioner.partition( m_nchare );
}

void
Transporter::distributed()
// *****************************************************************************
// Reduction target: all PEs have distrbuted their mesh after partitioning
// *****************************************************************************
{
  m_partitioner.refine();
}

void
Transporter::refinserted( int error )
// *****************************************************************************
// Reduction target: all PEs have created the mesh refiners
//! \param[in] Error aggregated across all PEs with operator max
// *****************************************************************************
{
  if (error) {
    m_print << "\n>>> ERROR: A worker chare was not assigned any mesh "
               "elements. This can happen in SMP-mode with a large +ppn "
               "parameter (number of worker threads per logical node) and is "
               "most likely the fault of the mesh partitioning algorithm not "
               "tolerating the case when it is asked to divide the "
               "computational domain into a number of partitions different "
               "than the number of ranks it is called on, i.e., in case of "
               "overdecomposition and/or calling the partitioner in SMP mode "
               "with +ppn larger than 1. Solution 1: Try a different "
               "partitioning algorithm (e.g., rcb instead of mj). Solution 2: "
               "Decrease +ppn.";
    finish();
  } else {
     m_refiner.doneInserting();
  }
}

void
Transporter::matched( std::size_t extra )
// *****************************************************************************
// Reduction target: all mesh refiner chares have distributed their newly added
// node IDs that are shared among chares
//! \param[in] extra Max number of edges/chare collected across all chares that
//!   correction due to refinement along chare boundaries
// *****************************************************************************
{
//std::cout << "max extra: " << extra << '\n';
  // If at least a single edge on a chare still needs correction, do correction,
  // otherwise, this initial mesh refinement step is complete
  if (extra > 0)
    m_refiner.correctref();
  else {
    m_progMesh.inc< REFINE >();
    m_refiner.nextref();
  }
}

void
Transporter::refined( std::size_t nelem, std::size_t npoin )
// *****************************************************************************
// Reduction target: all PEs have refined their mesh
//! \param[in] nelem Total number of elements in mesh across the whole problem
//! \param[in] npoin Total number of points in mesh across the whole problem
// *****************************************************************************
{
  m_sorter.doneInserting();

  m_nelem = nelem;
  m_npoin = npoin;
}

void
Transporter::bounds()
// *****************************************************************************
// Reduction target: all Solver (PEs) have computed their row bounds
// *****************************************************************************
{
  m_sorter.createDiscWorkers();
}

void
Transporter::discinserted()
// *****************************************************************************
// Reduction target: all Discretization chares have been inserted
// *****************************************************************************
{
  m_scheme.doneDiscInserting< tag::bcast >();
}

void
Transporter::disccreated()
// *****************************************************************************
// Reduction target: all Discretization constructors have been called
// *****************************************************************************
{
  m_progMesh.end();

  if (g_inputdeck.get< tag::amr, tag::initamr >()) {
    m_print.section( "Refined mesh graph statistics" );
    m_print.item( "Number of tetrahedra", m_nelem );
    m_print.item( "Number of nodes", m_npoin );
    m_print.endsubsection();
  }

  m_progWork.start( "Preparing workers",
                    {{ m_nchare, m_nchare, m_nchare, m_nchare, m_nchare }} );

  m_sorter.createWorkers();

  auto sch = g_inputdeck.get< tag::discr, tag::scheme >();
  if (sch == ctr::SchemeType::MatCG || sch == ctr::SchemeType::DiagCG)
    m_scheme.doneDistFCTInserting< tag::bcast >();

  m_scheme.vol< tag::bcast >();
}

void
Transporter::workinserted()
// *****************************************************************************
// Reduction target: all worker (derived discretization) chares have been
// inserted
// *****************************************************************************
{
  m_scheme.doneInserting< tag::bcast >();
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
  const auto scheme = g_inputdeck.get< tag::discr, tag::scheme >();
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
Transporter::flattened()
// *****************************************************************************
// Reduction target indicating that all Partitioner chare groups have finished
// flattening its global mesh node IDs and they are ready for computing the
// communication maps required for node ID reordering
// *****************************************************************************
{
  //m_sorter.gather();
}

void
Transporter::comfinal()
// *****************************************************************************
// Reduction target indicating that the communication has been established among
// PEs
// *****************************************************************************
{
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
  m_scheme.stat< tag::bcast >();
}

void
Transporter::minstat( tk::real d0, tk::real d1, tk::real d2 )
// *****************************************************************************
// Reduction target yielding minimum mesh statistcs across all workers
//! \param[in] d0 Minimum mesh statistics collected over all chares
//! \param[in] d1 Minimum mesh statistics collected over all chares
//! \param[in] d2 Minimum mesh statistics collected over all chares
// *****************************************************************************
{
  m_minstat[0] = d0;  // minimum edge length
  m_minstat[1] = d1;  // minimum cell volume cubic root
  m_minstat[2] = d2;  // minimum number of cells on chare

  minstat_complete();
}

void
Transporter::maxstat( tk::real d0, tk::real d1, tk::real d2 )
// *****************************************************************************
// Reduction target yielding the maximum mesh statistics across all workers
//! \param[in] d0 Maximum mesh statistics collected over all chares
//! \param[in] d1 Maximum mesh statistics collected over all chares
//! \param[in] d2 Maximum mesh statistics collected over all chares
// *****************************************************************************
{
  m_maxstat[0] = d0;  // maximum edge length
  m_maxstat[1] = d1;  // maximum cell volume cubic root
  m_maxstat[2] = d2;  // maximum number of cells on chare

  maxstat_complete();
}

void
Transporter::sumstat( tk::real d0, tk::real d1, tk::real d2, tk::real d3,
                      tk::real d4, tk::real d5 )
// *****************************************************************************
// Reduction target yielding the sum mesh statistics across all workers
//! \param[in] d0 Sum mesh statistics collected over all chares
//! \param[in] d1 Sum mesh statistics collected over all chares
//! \param[in] d2 Sum mesh statistics collected over all chares
//! \param[in] d3 Sum mesh statistics collected over all chares
//! \param[in] d4 Sum mesh statistics collected over all chares
//! \param[in] d5 Sum mesh statistics collected over all chares
// *****************************************************************************
{
  m_avgstat[0] = d1 / d0;      // average edge length
  m_avgstat[1] = d3 / d2;      // average cell volume cubic root
  m_avgstat[2] = d5 / d4;      // average number of cells per chare

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

  // Create new PDF file (overwrite if exists)
  tk::PDFWriter pdfn( "mesh_ntet_pdf.txt" );
  // Output number of cells PDF
  pdfn.writeTxt( pdf[2], tk::ctr::PDFInfo{ {"PDF"}, {}, {"ntets"} } );

  pdfstat_complete();
}

void
Transporter::stat()
// *****************************************************************************
// Echo diagnostics mesh statistics
// *****************************************************************************
{
  m_progWork.end();

  m_print.diag( "Mesh statistics: min/max/avg(edgelength) = " +
                std::to_string( m_minstat[0] ) + " / " +
                std::to_string( m_maxstat[0] ) + " / " +
                std::to_string( m_avgstat[0] ) );
  m_print.diag( "Mesh statistics: min/max/avg(V^{1/3}) = " +
                std::to_string( m_minstat[1] ) + " / " +
                std::to_string( m_maxstat[1] ) + " / " +
                std::to_string( m_avgstat[1] ) );
  m_print.diag( "Mesh statistics: min/max/avg(ntets) = " +
              std::to_string( static_cast<std::size_t>(m_minstat[2]) ) + " / " +
              std::to_string( static_cast<std::size_t>(m_maxstat[2]) ) + " / " +
              std::to_string( static_cast<std::size_t>(m_avgstat[2]) ) );

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
