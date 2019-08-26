// *****************************************************************************
/*!
  \file      src/Inciter/Transporter.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
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

#include "Macro.hpp"
#include "Transporter.hpp"
#include "Fields.hpp"
#include "PDEStack.hpp"
#include "UniPDF.hpp"
#include "PDFWriter.hpp"
#include "ContainerUtil.hpp"
#include "LoadDistributor.hpp"
#include "MeshReader.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "NodeDiagnostics.hpp"
#include "ElemDiagnostics.hpp"
#include "DiagWriter.hpp"
#include "Callback.hpp"

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
  m_ncit( 0 ),
  m_nt0refit( 0 ),
  m_ndtrefit( 0 ),
  m_scheme( g_inputdeck.get< tag::discr, tag::scheme >() ),
  m_partitioner(),
  m_refiner(),
  m_meshwriter(),
  m_sorter(),
  m_nelem( 0 ),
  m_npoin_larger( 0 ),
  m_meshvol( 0.0 ),
  m_minstat( {{ 0.0, 0.0, 0.0 }} ),
  m_maxstat( {{ 0.0, 0.0, 0.0 }} ),
  m_avgstat( {{ 0.0, 0.0, 0.0 }} ),
  m_timer(),
  m_progMesh( m_print, g_inputdeck.get< tag::cmd, tag::feedback >(),
              ProgMeshPrefix, ProgMeshLegend ),
  m_progWork( m_print, g_inputdeck.get< tag::cmd, tag::feedback >(),
              ProgWorkPrefix, ProgWorkLegend )
// *****************************************************************************
//  Constructor
// *****************************************************************************
{
  // Echo configuration to screen
  info();

  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  const auto t0 = g_inputdeck.get< tag::discr, tag::t0 >();
  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto constdt = g_inputdeck.get< tag::discr, tag::dt >();

  // If the desired max number of time steps is larger than zero, and the
  // termination time is larger than the initial time, and the constant time
  // step size (if that is used) is smaller than the duration of the time to be
  // simulated, we have work to do, otherwise, finish right away. If a constant
  // dt is not used, that part of the logic is always true as the default
  // constdt is zero, see inciter::ctr::InputDeck::InputDeck().
  if ( nstep != 0 && term > t0 && constdt < term-t0 ) {

    // Enable SDAG waits for collecting mesh statistics
    thisProxy.wait4stat();

    // Configure and write diagnostics file header
    diagHeader();

    // Create mesh partitioner AND boundary condition object group
    createPartitioner();

  } else finish( 0, t0 );      // stop if no time stepping requested
}

Transporter::Transporter( CkMigrateMessage* m ) :
  CBase_Transporter( m ),
  m_print( g_inputdeck.get<tag::cmd,tag::verbose>() ? std::cout : std::clog ),
  m_progMesh( m_print, g_inputdeck.get< tag::cmd, tag::feedback >(),
              ProgMeshPrefix, ProgMeshLegend ),
    m_progWork( m_print, g_inputdeck.get< tag::cmd, tag::feedback >(),
              ProgWorkPrefix, ProgWorkLegend )
// *****************************************************************************
//  Migrate constructor: returning from a checkpoint
//! \param[in] m Charm++ migrate message
// *****************************************************************************
{
   m_print.diag( "Restarted from checkpoint" );
   info();
   inthead();
}

void
Transporter::info()
// *****************************************************************************
// Echo configuration to screen
// *****************************************************************************
{
  m_print.part( "Factory" );

  // Print out info data layout
  m_print.list( "Unknowns data layout (CMake: FIELD_DATA_LAYOUT)",
                std::list< std::string >{ tk::Fields::layout() } );

  // Re-create partial differential equations stack for output
  PDEStack stack;

  // Print out information on PDE factories
  m_print.eqlegend();
  m_print.eqlist( "Registered PDEs using continuous Galerkin (CG) methods",
                  stack.cgfactory(), stack.cgntypes() );
  m_print.eqlist( "Registered PDEs using discontinuous Galerkin (DG) methods",
                  stack.dgfactory(), stack.dgntypes() );
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

  if (scheme == ctr::SchemeType::DiagCG) {
    auto fct = g_inputdeck.get< tag::discr, tag::fct >();
    m_print.item( "Flux-corrected transport (FCT)", fct );
    if (fct)
      m_print.item( "FCT mass diffusion coeff",
                    g_inputdeck.get< tag::discr, tag::ctau >() );
  } else if (scheme == ctr::SchemeType::DG ||
             scheme == ctr::SchemeType::P0P1 || scheme == ctr::SchemeType::DGP1 ||
             scheme == ctr::SchemeType::DGP2 || scheme == ctr::SchemeType::PDG)
  {
    m_print.Item< ctr::Flux, tag::discr, tag::flux >();
    m_print.Item< ctr::Limiter, tag::discr, tag::limiter >();
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

  // Print out adaptive polynomial refinement configuration
  if (scheme == ctr::SchemeType::PDG) {
    m_print.section( "Polynomial refinement (p-ref)" );
    m_print.item( "p-refinement",
                  g_inputdeck.get< tag::pref, tag::pref >() );
    m_print.Item< ctr::PrefIndicator, tag::pref, tag::indicator >();
    m_print.item( "Max degrees of freedom",
                  g_inputdeck.get< tag::pref, tag::ndofmax >() );
    m_print.item( "Tolerance",
                  g_inputdeck.get< tag::pref, tag::tolref >() );
  }

  // Print out adaptive mesh refinement configuration
  const auto amr = g_inputdeck.get< tag::amr, tag::amr >();
  if (amr) {
    m_print.section( "Mesh refinement (h-ref)" );
    m_print.refvar( g_inputdeck.get< tag::amr, tag::refvar >(),
                    g_inputdeck.get< tag::amr, tag::id >() );
    m_print.Item< ctr::AMRError, tag::amr, tag::error >();
    auto t0ref = g_inputdeck.get< tag::amr, tag::t0ref >();
    m_print.item( "Refinement at t<0 (t0ref)", t0ref );
    if (t0ref) {
      const auto& initref = g_inputdeck.get< tag::amr, tag::init >();
      m_print.item( "Initial refinement steps", initref.size() );
      m_print.ItemVec< ctr::AMRInitial >( initref );
      m_print.edgeref( g_inputdeck.get< tag::amr, tag::edge >() );

      auto rmax =
        std::numeric_limits< kw::amr_xminus::info::expect::type >::max();
      auto eps =
        std::numeric_limits< kw::amr_xminus::info::expect::type >::epsilon();
     
      auto xminus = g_inputdeck.get< tag::amr, tag::xminus >();
      if (std::abs( xminus - rmax ) > eps)
        m_print.item( "Initial refinement x-", xminus );
      auto xplus = g_inputdeck.get< tag::amr, tag::xplus >();
      if (std::abs( xplus - rmax ) > eps)
        m_print.item( "Initial refinement x+", xplus );

      auto yminus = g_inputdeck.get< tag::amr, tag::yminus >();
      if (std::abs( yminus - rmax ) > eps)
        m_print.item( "Initial refinement y-", yminus );
      auto yplus = g_inputdeck.get< tag::amr, tag::yplus >();
      if (std::abs( yplus - rmax ) > eps)
        m_print.item( "Initial refinement y+", yplus );

      auto zminus = g_inputdeck.get< tag::amr, tag::zminus >();
      if (std::abs( zminus - rmax ) > eps)
        m_print.item( "Initial refinement z-", zminus );
      auto zplus = g_inputdeck.get< tag::amr, tag::zplus >();
      if (std::abs( zplus - rmax ) > eps)
        m_print.item( "Initial refinement z+", zplus );
    }
    auto dtref = g_inputdeck.get< tag::amr, tag::dtref >();
    m_print.item( "Refinement at t>0 (dtref)", dtref );
    if (dtref) {
      auto dtfreq = g_inputdeck.get< tag::amr, tag::dtfreq >();
      m_print.item( "Mesh refinement frequency, t>0", dtfreq );
      m_print.item( "Uniform-only mesh refinement, t>0",
                    g_inputdeck.get< tag::amr, tag::dtref_uniform >() );
    }
    m_print.item( "Refinement tolerance",
                  g_inputdeck.get< tag::amr, tag::tolref >() );
    m_print.item( "De-refinement tolerance",
                  g_inputdeck.get< tag::amr, tag::tolderef >() );
  }

  // Print I/O filenames
  m_print.section( "Output filenames and directories" );
  m_print.item( "Field output file(s)",
    g_inputdeck.get< tag::cmd, tag::io, tag::output >() +
    ".e-s.<meshid>.<numchares>.<chareid>" );
  m_print.item( "Diagnostics file",
                g_inputdeck.get< tag::cmd, tag::io, tag::diag >() );
  m_print.item( "Checkpoint/restart directory",
                g_inputdeck.get< tag::cmd, tag::io, tag::restart >() + '/' );

  // Print output intervals
  m_print.section( "Output intervals" );
  m_print.item( "TTY", g_inputdeck.get< tag::interval, tag::tty>() );
  m_print.item( "Field", g_inputdeck.get< tag::interval, tag::field >() );
  m_print.item( "Diagnostics",
                g_inputdeck.get< tag::interval, tag::diag >() );
  m_print.item( "Checkpoint/restart",
                g_inputdeck.get< tag::cmd, tag::rsfreq >() );
  m_print.endsubsection();
}

void
Transporter::createPartitioner()
// *****************************************************************************
// Create mesh partitioner AND boundary conditions group
// *****************************************************************************
{
  // Create mesh reader for reading side sets from file
  tk::MeshReader mr( g_inputdeck.get< tag::cmd, tag::io, tag::input >() );

  std::map< int, std::vector< std::size_t > > bface;
  std::map< int, std::vector< std::size_t > > faces;
  std::map< int, std::vector< std::size_t > > bnode;

  // Read boundary (side set) data from input file
  const auto scheme = g_inputdeck.get< tag::discr, tag::scheme >();
  const auto centering = ctr::Scheme().centering( scheme );
  if (centering == tk::Centering::ELEM) {
    // Read boundary-face connectivity on side sets
    mr.readSidesetFaces( bface, faces );
    // Verify boundarty condition (BC) side sets used exist in mesh file
    matchBCs( g_dgpde, bface );
  } else if (centering == tk::Centering::NODE) {
    // Read node lists on side sets
    bnode = mr.readSidesetNodes();
    // Verify boundarty condition (BC) side sets used exist in mesh file
    matchBCs( g_cgpde, bnode );
  }

  // Create partitioner callbacks (order matters)
  tk::PartitionerCallback cbp {{
      CkCallback( CkReductionTarget(Transporter,load), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,distributed), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,refinserted), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,refined), thisProxy )
  }};

  // Create refiner callbacks (order matters)
  tk::RefinerCallback cbr {{
      CkCallback( CkReductionTarget(Transporter,edges), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,compatibility), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,bndint), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,matched), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,refined), thisProxy )
  }};

  // Create sorter callbacks (order matters)
  tk::SorterCallback cbs {{
      CkCallback( CkReductionTarget(Transporter,queried), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,responded), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,discinserted), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,workinserted), thisProxy )
  }};

  // Start timer measuring preparation of the mesh for partitioning
  m_timer[ TimerTag::MESH_READ ];

  // Create mesh partitioner Charm++ chare group and start preparing mesh
  m_print.diag( "Reading mesh" );

  // Create empty mesh sorter Charm++ chare array (bound to workers)
  m_sorter = CProxy_Sorter::ckNew( m_scheme.arrayoptions() );

  // Create empty mesh refiner chare array (bound to workers)
  m_refiner = CProxy_Refiner::ckNew( m_scheme.arrayoptions() );

  // Create MeshWriter chare group
  m_meshwriter = tk::CProxy_MeshWriter::ckNew(
                    g_inputdeck.get< tag::selected, tag::filetype >(),
                    centering,
                    g_inputdeck.get< tag::cmd, tag::benchmark >() );

  // Create mesh partitioner Charm++ chare nodegroup
  m_partitioner =
    CProxy_Partitioner::ckNew( cbp, cbr, cbs, thisProxy, m_refiner, m_sorter,
                               m_meshwriter, m_scheme, bface, faces, bnode );
}

void
Transporter::load( std::size_t nelem, std::size_t npoin )
// *****************************************************************************
// Reduction target: the mesh has been read from file on all PEs
//! \param[in] nelem Total number of mesh elements (summed across all PEs)
//! \param[in] npoin Total number of mesh nodes (summed across all PEs). Note
//!    that in parallel this is larger than the total number of points in the
//!    mesh, because the boundary nodes are double-counted.
// *****************************************************************************
{
  m_npoin_larger = npoin;

  // Compute load distribution given total work (nelem) and user-specified
  // virtualization
  uint64_t chunksize, remainder;
  m_nchare = static_cast<int>(
               tk::linearLoadDistributor(
                 g_inputdeck.get< tag::cmd, tag::virtualization >(),
                 nelem, CkNumPes(), chunksize, remainder ) );

  // Start timer measuring preparation of the mesh for partitioning
  const auto& timer = tk::cref_find( m_timer, TimerTag::MESH_READ );
  m_print.diag( "Mesh read time: " + std::to_string( timer.dsec() ) + " sec" );

  // Print out mesh graph stats
  m_print.section( "Input mesh graph statistics" );
  m_print.item( "Number of tetrahedra", nelem );

  // Print out mesh partitioning configuration
  m_print.section( "Mesh partitioning" );
  m_print.Item< tk::ctr::PartitioningAlgorithm,
                tag::selected, tag::partitioner >();

  // Print out info on load distribution
  m_print.section( "Initial load distribution" );
  m_print.item( "Virtualization [0.0...1.0]",
                g_inputdeck.get< tag::cmd, tag::virtualization >() );
  m_print.item( "Number of tetrahedra", nelem );
  m_print.item( "Number of work units", m_nchare );

  m_print.endsubsection();

  // Tell meshwriter the total number of chares
  m_meshwriter.nchare( m_nchare );

  // Query number of initial mesh refinement steps
  int nref = 0;
  if (g_inputdeck.get< tag::amr, tag::t0ref >())
    nref = static_cast<int>( g_inputdeck.get< tag::amr, tag::init >().size() );

  m_progMesh.start( "Preparing mesh", {{ CkNumPes(), CkNumPes(), nref,
    m_nchare, m_nchare, m_nchare, m_nchare }} );

  // Partition the mesh
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
//! \param[in] error aggregated across all PEs with operator max
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

    finish( 0, g_inputdeck.get< tag::discr, tag::t0 >() );

  } else {

     m_refiner.doneInserting();

  }
}

void
Transporter::edges()
// *****************************************************************************
// Reduction target: all mesh refiner chares have setup their boundary edges
// *****************************************************************************
{
  m_refiner.refine();
}

void
Transporter::compatibility( int modified )
// *****************************************************************************
// Reduction target: all mesh refiner chares have received a round of edges,
// and ran their compatibility algorithm
//! \param[in] modified Sum acorss all workers, if nonzero, mesh is modified
//! \details This is called iteratively, until convergence by Refiner. At this
//!   point all Refiner chares have received a round of edge data (tags whether
//!   an edge needs to be refined, etc.), and applied the compatibility
//!   algorithm independent of other Refiner chares. We keep going until the
//!   mesh is no longer modified by the compatibility algorithm (based on a new
//!   round of edge data communication started in Refiner::comExtra().
// *****************************************************************************
{
  if (modified)
    m_refiner.comExtra();
  else
    m_refiner.correctref();
}

void
Transporter::matched( std::size_t nextra,
                      std::size_t nref,
                      std::size_t nderef,
                      std::size_t initial )
// *****************************************************************************
// Reduction target: all mesh refiner chares have matched/corrected the tagging
// of chare-boundary edges, all chares are ready to perform refinement.
//! \param[in] nextra Sum (across all chares) of the number of edges on each
//!   chare that need correction along chare boundaries
//! \param[in] nref Sum of number of refined tetrahedra across all chares.
//! \param[in] nderef Sum of number of derefined tetrahedra across all chares.
//! \param[in] initial Sum of contributions from all chares. If larger than
//!    zero, we are during time stepping and if zero we are during setup.
// *****************************************************************************
{
  // If at least a single edge on a chare still needs correction, do correction,
  // otherwise, this mesh refinement step is complete
  if (nextra > 0) {

    ++m_ncit;
    m_refiner.comExtra();

  } else {

    if (initial > 0) {

      if (!g_inputdeck.get< tag::cmd, tag::feedback >()) {
        m_print.diag( { "t0ref", "nref", "nderef", "ncorr" },
                      { ++m_nt0refit, nref, nderef, m_ncit } );
      }
      m_progMesh.inc< REFINE >();

    } else {

      m_print.diag( { "dtref", "nref", "nderef", "ncorr" },
                    { ++m_ndtrefit, nref, nderef, m_ncit }, false );

    }

    m_ncit = 0;
    m_refiner.perform();

  }
}

void
Transporter::bndint( tk::real sx, tk::real sy, tk::real sz, tk::real cb )
// *****************************************************************************
// Compute surface integral across the whole problem and perform leak-test
//! \param[in] sx X component of vector summed
//! \param[in] sy Y component of vector summed
//! \param[in] sz Z component of vector summed
//! \param[in] cb Invoke callback if positive
//! \details This function aggregates partial surface integrals across the
//!   boundary faces of the whole problem. After this global sum a
//!   non-zero vector result indicates a leak, e.g., a hole in the boundary,
//!   which indicates an error in the boundary face data structures used to
//!   compute the partial surface integrals.
// *****************************************************************************
{
  std::stringstream err;
  if (cb < 0.0) {  // called from Refiner
    err << "Mesh boundary leaky after mesh refinement step; this is due to a "
     "problem with updating the side sets used to specify boundary conditions "
     "on faces, required for DG methods: ";
  } else if (cb > 0.0) {  // called from DG
    err << "Mesh boundary leaky during initialization of the DG algorithm; this "
    "is due to incorrect or incompletely specified boundary conditions for a "
    "given input mesh: ";
  }

  auto eps = std::numeric_limits< tk::real >::epsilon() * 1.0e+3; // ~ 2.0e-13

  if (std::abs(sx) > eps || std::abs(sy) > eps || std::abs(sz) > eps) {
    err << "Integral result must be a zero vector: " << std::setprecision(12) <<
           std::abs(sx) << ", " << std::abs(sy) << ", " << std::abs(sz) <<
           ", eps = " << eps;
    Throw( err.str() );
  }

  if (cb > 0.0) m_scheme.bcast< Scheme::resizeComm >();
}

void
Transporter::refined( std::size_t nelem, std::size_t npoin )
// *****************************************************************************
// Reduction target: all PEs have refined their mesh
//! \param[in] nelem Total number of elements in mesh across the whole problem
//! \param[in] npoin Total number of mesh nodes (summed across all PEs). Note
//!    that in parallel this is larger than the total number of points in the
//!    mesh, because the boundary nodes are double-counted.
// *****************************************************************************
{
  m_sorter.doneInserting();

  m_nelem = nelem;
  m_npoin_larger = npoin;

  m_sorter.setup( m_npoin_larger );
}

void
Transporter::queried()
// *****************************************************************************
// Reduction target: all Sorter chares have queried their boundary nodes
// *****************************************************************************
{
  m_sorter.response();
}

void
Transporter::responded()
// *****************************************************************************
// Reduction target: all Sorter chares have responded with their boundary nodes
// *****************************************************************************
{
  m_sorter.start();
}

void
Transporter::resized()
// *****************************************************************************
// Reduction target: all worker chares have resized their own data after
// mesh refinement
//! \note Only used for nodal schemes
// *****************************************************************************
{
  m_scheme.disc().vol();
  m_scheme.bcast< Scheme::lhs >();
}

void
Transporter::discinserted()
// *****************************************************************************
// Reduction target: all Discretization chares have been inserted
// *****************************************************************************
{
  m_scheme.disc().doneInserting();
}

void
Transporter::disccreated()
// *****************************************************************************
// Reduction target: all Discretization constructors have been called
// *****************************************************************************
{
  m_progMesh.end();

  if (g_inputdeck.get< tag::amr, tag::t0ref >()) {
    m_print.section( "Initially (t<0) refined mesh graph statistics" );
    m_print.item( "Number of tetrahedra", m_nelem );
    m_print.endsubsection();
  }

  m_refiner.sendProxy();

  auto sch = g_inputdeck.get< tag::discr, tag::scheme >();
  if (sch == ctr::SchemeType::DiagCG) m_scheme.fct().doneInserting();

  m_scheme.disc().vol();
}

void
Transporter::workinserted()
// *****************************************************************************
// Reduction target: all worker (derived discretization) chares have been
// inserted
// *****************************************************************************
{
  m_scheme.bcast< Scheme::doneInserting >();
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
  if (scheme == ctr::SchemeType::DiagCG || scheme == ctr::SchemeType::ALECG)
    for (const auto& eq : g_cgpde) varnames( eq, var );
  else if (scheme == ctr::SchemeType::DG ||
           scheme == ctr::SchemeType::P0P1 || scheme == ctr::SchemeType::DGP1 ||
           scheme == ctr::SchemeType::DGP2 || scheme == ctr::SchemeType::PDG)
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
Transporter::comfinal( int initial )
// *****************************************************************************
// Reduction target indicating that communication maps have been setup
//! \param[in] initial Sum of contributions from all chares. If larger than
//!    zero, we are during time stepping and if zero we are during setup.
// *****************************************************************************
// [Discretization-specific communication maps]
{
  if (initial > 0) {
    m_progWork.end();
    m_scheme.bcast< Scheme::setup >();
    // Turn on automatic load balancing
    tk::CProxy_LBSwitch::ckNew( g_inputdeck.get<tag::cmd,tag::verbose>() );
  } else {
    m_scheme.bcast< Scheme::lhs >();
  }
}
// [Discretization-specific communication maps]

void
Transporter::totalvol( tk::real v, tk::real initial )
// *****************************************************************************
// Reduction target summing total mesh volume across all workers
//! \param[in] v Mesh volume summed across the whole problem
//! \param[in] initial Sum of contributions from all chares. If larger than
//!    zero, we are during time stepping and if zero we are during setup.
// *****************************************************************************
{
  m_meshvol = v;

  if (initial > 0.0)
    m_scheme.disc().stat( m_meshvol );
  else
    m_scheme.bcast< Scheme::resized >();
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
  // cppcheck-suppress containerOutOfBounds
  pdfe.writeTxt( pdf[0],
                 tk::ctr::PDFInfo{ {"PDF"}, {}, {"edgelength"}, 0, 0.0 } );

  // Create new PDF file (overwrite if exists)
  tk::PDFWriter pdfv( "mesh_vol_pdf.txt" );
  // Output cell volume cubic root PDF
  pdfv.writeTxt( pdf[1],
                 tk::ctr::PDFInfo{ {"PDF"}, {}, {"V^{1/3}"}, 0, 0.0 } );

  // Create new PDF file (overwrite if exists)
  tk::PDFWriter pdfn( "mesh_ntet_pdf.txt" );
  // Output number of cells PDF
  pdfn.writeTxt( pdf[2],
                 tk::ctr::PDFInfo{ {"PDF"}, {}, {"ntets"}, 0, 0.0 } );

  pdfstat_complete();
}

void
Transporter::stat()
// *****************************************************************************
// Echo diagnostics on mesh statistics
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
  m_print.diag( "Mesh statistics: min/max/avg(ntets) = " +
              std::to_string( static_cast<std::size_t>(m_minstat[2]) ) + " / " +
              std::to_string( static_cast<std::size_t>(m_maxstat[2]) ) + " / " +
              std::to_string( static_cast<std::size_t>(m_avgstat[2]) ) );

  // Print out time integration header to screen
  inthead();

  m_progWork.start( "Preparing workers",
                    {{ m_nchare, m_nchare, m_nchare, m_nchare, m_nchare }} );
  // Create "derived-class" workers
  m_sorter.createWorkers();
}

void
Transporter::inthead()
// *****************************************************************************
// Print out time integration header to screen
// *****************************************************************************
{
  m_print.inthead( "Time integration", "Navier-Stokes solver",
  "Legend: it - iteration count\n"
  "         t - time\n"
  "        dt - time step size\n"
  "       ETE - estimated time elapsed (h:m:s)\n"
  "       ETA - estimated time for accomplishment (h:m:s)\n"
  "       EGT - estimated grind time (ms/status)\n"
  "       out - output-saved flags\n"
  "             f - field\n"
  "             d - diagnostics\n"
  "             h - h-refinement\n"
  "             l - load balancing\n"
  "             r - checkpoint/restart\n",
  "\n      it             t            dt        ETE        ETA        EGT  out\n"
    " -------------------------------------------------------------------------\n" );
}

void
Transporter::diagnostics( CkReductionMsg* msg )
// *****************************************************************************
// Reduction target optionally collecting diagnostics, e.g., residuals
//! \param[in] msg Serialized diagnostics vector aggregated across all PEs
//! \note Only used for nodal schemes
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
    diag[i] = sqrt( d[L2SOL][i] / m_meshvol );
  
  // Query user-requested error types to output
  const auto& error = g_inputdeck.get< tag::diag, tag::error >();

  decltype(ncomp) n = 0;
  for (const auto& e : error) {
    n += ncomp;
    if (e == tk::ctr::ErrorType::L2) {
      // Finish computing the L2 norm of the numerical - analytical solution
     for (std::size_t i=0; i<d[L2ERR].size(); ++i)
       diag.push_back( sqrt( d[L2ERR][i] / m_meshvol ) );
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

  // Evaluate whether to continue with next step
  m_scheme.bcast< Scheme::refine >();
}

void
Transporter::resume()
// *****************************************************************************
// Resume execution from checkpoint/restart files
//! \details This is invoked by Charm++ after the checkpoint is done, as well as
//!   when the restart (returning from a checkpoint) is complete
// *****************************************************************************
{
  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto eps = std::numeric_limits< tk::real >::epsilon();

  // If neither max iterations nor max time reached, continue, otherwise finish
  if (std::fabs(m_t-term) > eps && m_it < nstep)
    m_scheme.bcast< Scheme::evalLB >();
  else
    mainProxy.finalize();
}

void
Transporter::checkpoint( tk::real it, tk::real t )
// *****************************************************************************
// Save checkpoint/restart files
//! \param[in] it Iteration count
//! \param[in] t Physical time
// *****************************************************************************
{
  m_it = static_cast< uint64_t >( it );
  m_t = t;

  const auto& restart = g_inputdeck.get< tag::cmd, tag::io, tag::restart >();
  CkCallback res( CkIndex_Transporter::resume(), thisProxy );
  CkStartCheckpoint( restart.c_str(), res );
}

void
Transporter::finish( tk::real it, tk::real t )
// *****************************************************************************
// Normal finish of time stepping
//! \param[in] it Iteration count
//! \param[in] t Physical time
// *****************************************************************************
{
  checkpoint( it, t );
}

#include "NoWarning/transporter.def.h"
