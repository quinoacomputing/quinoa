// *****************************************************************************
/*!
  \file      src/Inciter/Transporter.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
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

#include <brigand/algorithms/for_each.hpp>

#include "Macro.hpp"
#include "Transporter.hpp"
#include "Fields.hpp"
#include "PDEStack.hpp"
#include "UniPDF.hpp"
#include "PDFWriter.hpp"
#include "ContainerUtil.hpp"
#include "LoadDistributor.hpp"
#include "MeshReader.hpp"
#include "Inciter/Types.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "NodeDiagnostics.hpp"
#include "ElemDiagnostics.hpp"
#include "DiagWriter.hpp"
#include "Callback.hpp"
#include "CartesianProduct.hpp"

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
  m_npoin( 0 ),
  m_finished( 0 ),
  m_meshvol( 0.0 ),
  m_minstat( {{ 0.0, 0.0, 0.0 }} ),
  m_maxstat( {{ 0.0, 0.0, 0.0 }} ),
  m_avgstat( {{ 0.0, 0.0, 0.0 }} ),
  m_timer(),
  m_progMesh( g_inputdeck.get< tag::cmd, tag::feedback >(),
              ProgMeshPrefix, ProgMeshLegend ),
  m_progWork( g_inputdeck.get< tag::cmd, tag::feedback >(),
              ProgWorkPrefix, ProgWorkLegend )
// *****************************************************************************
//  Constructor
// *****************************************************************************
{
  // Echo configuration to screen
  info( printer() );

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

  } else finish();      // stop if no time stepping requested
}

Transporter::Transporter( CkMigrateMessage* m ) :
  CBase_Transporter( m ),
  m_progMesh( g_inputdeck.get< tag::cmd, tag::feedback >(),
              ProgMeshPrefix, ProgMeshLegend ),
  m_progWork( g_inputdeck.get< tag::cmd, tag::feedback >(),
              ProgWorkPrefix, ProgWorkLegend )
// *****************************************************************************
//  Migrate constructor: returning from a checkpoint
//! \param[in] m Charm++ migrate message
// *****************************************************************************
{
   auto print = printer();
   print.diag( "Restarted from checkpoint" );
   info( print );
   inthead( print );
}

void
Transporter::info( const InciterPrint& print )
// *****************************************************************************
// Echo configuration to screen
//! \param[in] print Pretty printer object to use for printing
// *****************************************************************************
{
  print.part( "Factory" );

  // Print out info data layout
  print.list( "Unknowns data layout (CMake: FIELD_DATA_LAYOUT)",
              std::list< std::string >{ tk::Fields::layout() } );

  // Re-create partial differential equations stack for output
  PDEStack stack;

  // Print out information on PDE factories
  print.eqlegend();
  print.eqlist( "Registered PDEs using continuous Galerkin (CG) methods",
                stack.cgfactory(), stack.cgntypes() );
  print.eqlist( "Registered PDEs using discontinuous Galerkin (DG) methods",
                stack.dgfactory(), stack.dgntypes() );
  print.endpart();

  // Print out information on problem
  print.part( "Problem" );

  // Print out info on problem title
  if ( !g_inputdeck.get< tag::title >().empty() )
    print.title( g_inputdeck.get< tag::title >() );

  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  const auto t0 = g_inputdeck.get< tag::discr, tag::t0 >();
  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto constdt = g_inputdeck.get< tag::discr, tag::dt >();
  const auto cfl = g_inputdeck.get< tag::discr, tag::cfl >();
  const auto scheme = g_inputdeck.get< tag::discr, tag::scheme >();

  // Print discretization parameters
  print.section( "Discretization parameters" );
  print.Item< ctr::Scheme, tag::discr, tag::scheme >();

  if (scheme == ctr::SchemeType::DiagCG) {
    auto fct = g_inputdeck.get< tag::discr, tag::fct >();
    print.item( "Flux-corrected transport (FCT)", fct );
    if (fct) {
      print.item( "FCT mass diffusion coeff",
                  g_inputdeck.get< tag::discr, tag::ctau >() );
      print.item( "FCT small number",
                  g_inputdeck.get< tag::discr, tag::fcteps >() );
      print.item( "Clipping FCT",
                  g_inputdeck.get< tag::discr, tag::fctclip >() );
    }
  } else if (scheme == ctr::SchemeType::DG ||
             scheme == ctr::SchemeType::P0P1 || scheme == ctr::SchemeType::DGP1 ||
             scheme == ctr::SchemeType::DGP2 || scheme == ctr::SchemeType::PDG)
  {
    print.Item< ctr::Limiter, tag::discr, tag::limiter >();
  }
  print.item( "PE-locality mesh reordering",
              g_inputdeck.get< tag::discr, tag::pelocal_reorder >() );
  print.item( "Operator-access mesh reordering",
              g_inputdeck.get< tag::discr, tag::operator_reorder >() );
  auto steady = g_inputdeck.get< tag::discr, tag::steady_state >();
  print.item( "Local time stepping", steady );
  if (steady) {
    print.item( "L2-norm residual convergence criterion",
                g_inputdeck.get< tag::discr, tag::residual >() );
    print.item( "Convergence criterion component index",
                g_inputdeck.get< tag::discr, tag::rescomp >() );
  }
  print.item( "Number of time steps", nstep );
  print.item( "Start time", t0 );
  print.item( "Terminate time", term );

  if (std::abs(constdt - g_inputdeck_defaults.get< tag::discr, tag::dt >()) >
        std::numeric_limits< tk::real >::epsilon())
    print.item( "Constant time step size", constdt );
  else if (std::abs(cfl - g_inputdeck_defaults.get< tag::discr, tag::cfl >()) >
             std::numeric_limits< tk::real >::epsilon())
    print.item( "CFL coefficient", cfl );

  // Print out info on settings of selected partial differential equations
  print.pdes( "Partial differential equations integrated", stack.info() );

  // Print out adaptive polynomial refinement configuration
  if (scheme == ctr::SchemeType::PDG) {
    print.section( "Polynomial refinement (p-ref)" );
    print.item( "p-refinement",
                g_inputdeck.get< tag::pref, tag::pref >() );
    print.Item< ctr::PrefIndicator, tag::pref, tag::indicator >();
    print.item( "Max degrees of freedom",
                g_inputdeck.get< tag::pref, tag::ndofmax >() );
    print.item( "Tolerance",
                g_inputdeck.get< tag::pref, tag::tolref >() );
  }

  // Print out adaptive mesh refinement configuration
  const auto amr = g_inputdeck.get< tag::amr, tag::amr >();
  if (amr) {
    print.section( "Mesh refinement (h-ref)" );
    print.refvar( g_inputdeck.get< tag::amr, tag::refvar >(),
                  g_inputdeck.get< tag::amr, tag::id >() );
    print.Item< ctr::AMRError, tag::amr, tag::error >();
    auto t0ref = g_inputdeck.get< tag::amr, tag::t0ref >();
    print.item( "Refinement at t<0 (t0ref)", t0ref );
    if (t0ref) {
      const auto& initref = g_inputdeck.get< tag::amr, tag::init >();
      print.item( "Initial refinement steps", initref.size() );
      print.ItemVec< ctr::AMRInitial >( initref );
      print.ItemVecLegend< ctr::AMRInitial >();
      print.edgeref( g_inputdeck.get< tag::amr, tag::edge >() );

      auto eps =
        std::numeric_limits< kw::amr_xminus::info::expect::type >::epsilon();
     
      auto xminus = g_inputdeck.get< tag::amr, tag::xminus >();
      auto xminus_default = g_inputdeck_defaults.get< tag::amr, tag::xminus >();
      if (std::abs( xminus - xminus_default ) > eps)
        print.item( "Initial refinement x-", xminus );
      auto xplus = g_inputdeck.get< tag::amr, tag::xplus >();
      auto xplus_default = g_inputdeck_defaults.get< tag::amr, tag::xplus >();
      if (std::abs( xplus - xplus_default ) > eps)
        print.item( "Initial refinement x+", xplus );

      auto yminus = g_inputdeck.get< tag::amr, tag::yminus >();
      auto yminus_default = g_inputdeck_defaults.get< tag::amr, tag::yminus >();
      if (std::abs( yminus - yminus_default ) > eps)
        print.item( "Initial refinement y-", yminus );
      auto yplus = g_inputdeck.get< tag::amr, tag::yplus >();
      auto yplus_default = g_inputdeck_defaults.get< tag::amr, tag::yplus >();
      if (std::abs( yplus - yplus_default ) > eps)
        print.item( "Initial refinement y+", yplus );

      auto zminus = g_inputdeck.get< tag::amr, tag::zminus >();
      auto zminus_default = g_inputdeck_defaults.get< tag::amr, tag::zminus >();
      if (std::abs( zminus - zminus_default ) > eps)
        print.item( "Initial refinement z-", zminus );
      auto zplus = g_inputdeck.get< tag::amr, tag::zplus >();
      auto zplus_default = g_inputdeck_defaults.get< tag::amr, tag::zplus >();
      if (std::abs( zplus - zplus_default ) > eps)
        print.item( "Initial refinement z+", zplus );
    }
    auto dtref = g_inputdeck.get< tag::amr, tag::dtref >();
    print.item( "Refinement at t>0 (dtref)", dtref );
    if (dtref) {
      auto dtfreq = g_inputdeck.get< tag::amr, tag::dtfreq >();
      print.item( "Mesh refinement frequency, t>0", dtfreq );
      print.item( "Uniform-only mesh refinement, t>0",
                  g_inputdeck.get< tag::amr, tag::dtref_uniform >() );
    }
    print.item( "Refinement tolerance",
                g_inputdeck.get< tag::amr, tag::tolref >() );
    print.item( "De-refinement tolerance",
                g_inputdeck.get< tag::amr, tag::tolderef >() );
  }

  // Print I/O filenames
  print.section( "Output filenames and directories" );
  const auto& of = g_inputdeck.get< tag::cmd, tag::io, tag::output >();
  print.item( "Volume field output file(s)",
              of + ".e-s.<meshid>.<numchares>.<chareid>" );
  print.item( "Surface field output file(s)",
              of + "-surf.<surfid>.e-s.<meshid>.<numchares>.<chareid>" );
  print.item( "History output file(s)", of + ".hist.{pointid}" );
  print.item( "Diagnostics file",
              g_inputdeck.get< tag::cmd, tag::io, tag::diag >() );
  print.item( "Checkpoint/restart directory",
              g_inputdeck.get< tag::cmd, tag::io, tag::restart >() + '/' );

  // Print output intervals
  print.section( "Output intervals" );
  print.item( "TTY", g_inputdeck.get< tag::interval, tag::tty>() );
  print.item( "Field", g_inputdeck.get< tag::interval, tag::field >() );
  print.item( "Diagnostics",
              g_inputdeck.get< tag::interval, tag::diag >() );
  print.item( "Checkpoint/restart",
              g_inputdeck.get< tag::cmd, tag::rsfreq >() );

  const auto outsets = g_inputdeck.outsets();
  if (!outsets.empty()) {
    print.section( "Output fields" );
    print.item( "Surface side set(s)", tk::parameters( outsets ) );
  }

  const auto& pt = g_inputdeck.get< tag::history, tag::point >();
  const auto& id = g_inputdeck.get< tag::history, tag::id >();
  if (!pt.empty()) {
    print.section( "Output time history" );
    for (std::size_t p=0; p<pt.size(); ++p) {
      std::stringstream ss;
      auto prec = g_inputdeck.get< tag::prec, tag::history >();
      ss << std::setprecision( static_cast<int>(prec) );
      ss << of << ".hist." << id[p];
      print.longitem( "At point " + id[p] + ' ' + tk::parameters(pt[p]),
                      ss.str() );
    }
  }

  print.endsubsection();
}

bool
Transporter::matchBCs( std::map< int, std::vector< std::size_t > >& bnd )
// *****************************************************************************
 // Verify that side sets specified in the control file exist in mesh file
 //! \details This function does two things: (1) it verifies that the side
 //!   sets used in the input file (either to which boundary conditions (BC)
 //!   are assigned or listed as field output by the user in the
 //!   input file) all exist among the side sets read from the input mesh
 //!   file and errors out if at least one does not, and (2) it matches the
 //!   side set ids at which the user has configured BCs (or listed as an output
 //!   surface) to side set ids read from the mesh file and removes those face
 //!   and node lists associated to side sets that the user did not set BCs or
 //!   listed as field output on (as they will not need processing further since
 //!   they will not be used).
 //! \param[in,out] bnd Node or face lists mapped to side set ids
 //! \return True if sidesets have been used and found in mesh
// *****************************************************************************
 {
   // Query side set ids at which BCs assigned for all BC types for all PDEs
   using PDEsBCs =
     tk::cartesian_product< ctr::parameters::Keys, ctr::bc::Keys >;
   std::unordered_set< int > usedsets;
   brigand::for_each< PDEsBCs >( UserBC( g_inputdeck, usedsets ) );

   // Add sidesets requested for field output
   const auto& ss = g_inputdeck.get< tag::cmd, tag::io, tag::surface >();
   for (const auto& s : ss) {
     std::stringstream conv( s );
     int num;
     conv >> num;
     usedsets.insert( num );
   }

   // Find user-configured side set ids among side sets read from mesh file
   std::unordered_set< int > sidesets_used;
   for (auto i : usedsets) {       // for all side sets used in control file
     if (bnd.find(i) != end(bnd))  // used set found among side sets in file
       sidesets_used.insert( i );  // store side set id configured as BC
     else {
       Throw( "Boundary conditions specified on side set " +
         std::to_string(i) + " which does not exist in mesh file" );
     }
   }

   // Remove sidesets not used (will not process those further)
   tk::erase_if( bnd, [&]( auto& item ) {
     return sidesets_used.find( item.first ) == end(sidesets_used);
   });

   return !bnd.empty();
 }

void
Transporter::createPartitioner()
// *****************************************************************************
// Create mesh partitioner AND boundary conditions group
// *****************************************************************************
{
  // Create mesh reader for reading side sets from file
  tk::MeshReader mr( g_inputdeck.get< tag::cmd, tag::io, tag::input >() );

  // Read out total number of mesh points from mesh file
  m_npoin = mr.npoin();

  std::map< int, std::vector< std::size_t > > bface;
  std::map< int, std::vector< std::size_t > > faces;
  std::map< int, std::vector< std::size_t > > bnode;

  // Read boundary (side set) data from input file
  const auto scheme = g_inputdeck.get< tag::discr, tag::scheme >();
  const auto centering = ctr::Scheme().centering( scheme );

  // Read boundary-face connectivity on side sets
  mr.readSidesetFaces( bface, faces );

  bool bcs_set = false;
  if (centering == tk::Centering::ELEM) {

    // Verify boundarty condition (BC) side sets used exist in mesh file
    bcs_set = matchBCs( bface );

  } else if (centering == tk::Centering::NODE) {

    // Read node lists on side sets
    bnode = mr.readSidesetNodes();
    // Verify boundarty condition (BC) side sets used exist in mesh file
    bcs_set = matchBCs( bnode );
    bcs_set = bcs_set || matchBCs( bface );
  }

  auto print = printer();

  // Warn on no BCs
  if (!bcs_set) print << "\n>>> WARNING: No boundary conditions set\n\n";

  // Create partitioner callbacks (order matters)
  tk::PartitionerCallback cbp {{
      CkCallback( CkReductionTarget(Transporter,load), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,distributed), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,refinserted), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,refined), thisProxy )
  }};

  // Create refiner callbacks (order matters)
  tk::RefinerCallback cbr {{
      CkCallback( CkReductionTarget(Transporter,queriedRef), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,respondedRef), thisProxy )
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
  print.diag( "Reading mesh" );

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
Transporter::load( std::size_t nelem )
// *****************************************************************************
// Reduction target: the mesh has been read from file on all PEs
//! \param[in] nelem Total number of mesh elements (summed across all PEs)
// *****************************************************************************
{
  // Compute load distribution given total work (nelem) and user-specified
  // virtualization
  uint64_t chunksize, remainder;
  m_nchare = static_cast<int>(
               tk::linearLoadDistributor(
                 g_inputdeck.get< tag::cmd, tag::virtualization >(),
                 nelem, CkNumPes(), chunksize, remainder ) );

  auto print = printer();

  // Start timer measuring preparation of the mesh for partitioning
  const auto& timer = tk::cref_find( m_timer, TimerTag::MESH_READ );
  print.diag( "Mesh read time: " + std::to_string( timer.dsec() ) + " sec" );

  // Print out mesh partitioning configuration
  print.section( "Mesh partitioning" );
  print.Item< tk::ctr::PartitioningAlgorithm,
              tag::selected, tag::partitioner >();

  // Print out info on load distribution
  print.section( "Initial load distribution" );
  print.item( "Virtualization [0.0...1.0]",
              g_inputdeck.get< tag::cmd, tag::virtualization >() );
  print.item( "Number of tetrahedra", nelem );
  print.item( "Number of points", m_npoin );
  print.item( "Number of work units", m_nchare );

  print.endsubsection();

  // Tell meshwriter the total number of chares
  m_meshwriter.nchare( m_nchare );

  // Query number of initial mesh refinement steps
  int nref = 0;
  if (g_inputdeck.get< tag::amr, tag::t0ref >())
    nref = static_cast<int>( g_inputdeck.get< tag::amr, tag::init >().size() );

  m_progMesh.start( print, "Preparing mesh", {{ CkNumPes(), CkNumPes(), nref,
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

    printer() << "\n>>> ERROR: A worker chare was not assigned any mesh "
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
Transporter::queriedRef()
// *****************************************************************************
// Reduction target: all Sorter chares have queried their boundary nodes
// *****************************************************************************
{
  m_refiner.response();
}

void
Transporter::respondedRef()
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
//!    zero, we are before time stepping, if zero we are during time stepping.
// *****************************************************************************
{
  // If at least a single edge on a chare still needs correction, do correction,
  // otherwise, this mesh refinement step is complete
  if (nextra > 0) {

    ++m_ncit;
    m_refiner.comExtra();

  } else {

    auto print = printer();

    if (initial > 0) {

      if (!g_inputdeck.get< tag::cmd, tag::feedback >()) {
        const auto& initref = g_inputdeck.get< tag::amr, tag::init >();
        ctr::AMRInitial opt;
        print.diag( { "t0ref", "type", "nref", "nderef", "ncorr" },
                    { std::to_string(m_nt0refit),
                      opt.code(initref[m_nt0refit]),
                      std::to_string(nref),
                      std::to_string(nderef),
                      std::to_string(m_ncit) } );
        ++m_nt0refit;
      }
      m_progMesh.inc< REFINE >( print );

    } else {

      auto dtref_uni = g_inputdeck.get< tag::amr, tag::dtref_uniform >();
      print.diag( { "dtref", "type", "nref", "nderef", "ncorr" },
                  { std::to_string(++m_ndtrefit),
                    (dtref_uni?"uniform":"error"),
                    std::to_string(nref),
                    std::to_string(nderef),
                    std::to_string(m_ncit) },
                  false );

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
  if (cb < 0.0) {
    err << "Mesh boundary leaky after mesh refinement step; this is due to a "
     "problem with updating the side sets used to specify boundary conditions "
     "on faces: ";
  } else if (cb > 0.0) {
    err << "Mesh boundary leaky during initialization; this is due to "
    "incorrect or incompletely specified boundary conditions for a given input "
    "mesh: ";
  }

  auto eps = 1.0e-10;
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
//! \param[in] npoin Total number of mesh points (summed across all PEs). Note
//!    that in parallel this is larger than the total number of points in the
//!    mesh, because the boundary nodes are multi-counted.
// *****************************************************************************
{
  m_sorter.doneInserting();

  m_nelem = nelem;

  m_sorter.setup( npoin );
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
Transporter::startEsup()
// *****************************************************************************
// Reduction target: all worker chares have generated their own esup
//! \note Only used for cell-centered schemes
// *****************************************************************************
{
  m_scheme.bcast< Scheme::nodeNeighSetup >();
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
Transporter::disccreated( std::size_t npoin )
// *****************************************************************************
// Reduction target: all Discretization constructors have been called
//! \param[in] npoin Total number of mesh points (summed across all PEs). Note
//!  that as opposed to npoin in refined(), this npoin is not multi-counted, and
//!  thus should be correct in parallel.
// *****************************************************************************
{
  m_npoin = npoin;

  auto print = printer();

  m_progMesh.end( print );

  if (g_inputdeck.get< tag::amr, tag::t0ref >()) {
    print.section( "Initially (t<0) refined mesh graph statistics" );
    print.item( "Number of tetrahedra", m_nelem );
    print.item( "Number of points", m_npoin );
    print.endsubsection();
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

  // Augment diagnostics variables by L2-norm of the residual
  if (scheme == ctr::SchemeType::DiagCG || scheme == ctr::SchemeType::ALECG) {
    for (std::size_t i=0; i<nv; ++i) d.push_back( "L2(d" + var[i] + ')' );
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
    auto print = printer();
    m_progWork.end( print );
    m_scheme.bcast< Scheme::setup >();
    // Turn on automatic load balancing
    tk::CProxy_LBSwitch::ckNew();
    print.diag( "Load balancing on (if enabled in Charm++)" );
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
//!    zero, we are during setup, if zero, during time stepping.
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
  auto print = printer();

  print.diag( "Mesh statistics: min/max/avg(edgelength) = " +
              std::to_string( m_minstat[0] ) + " / " +
              std::to_string( m_maxstat[0] ) + " / " +
              std::to_string( m_avgstat[0] ) );
  print.diag( "Mesh statistics: min/max/avg(V^{1/3}) = " +
              std::to_string( m_minstat[1] ) + " / " +
              std::to_string( m_maxstat[1] ) + " / " +
              std::to_string( m_avgstat[1] ) );
  print.diag( "Mesh statistics: min/max/avg(ntets) = " +
              std::to_string( static_cast<std::size_t>(m_minstat[2]) ) + " / " +
              std::to_string( static_cast<std::size_t>(m_maxstat[2]) ) + " / " +
              std::to_string( static_cast<std::size_t>(m_avgstat[2]) ) );

  // Print out time integration header to screen
  inthead( print );

  m_progWork.start( print, "Preparing workers",
                    {{ m_nchare, m_nchare, m_nchare, m_nchare, m_nchare }} );

  // Create "derived-class" workers
  m_sorter.createWorkers();
}

void
Transporter::boxvol( tk::real v )
// *****************************************************************************
// Reduction target computing total volume of IC box
//! \param[in] v Total volume within user-specified box IC
// *****************************************************************************
{
  if (v > 0.0) printer().diag( "Box IC volume: " + std::to_string(v) );
  m_scheme.bcast< Scheme::box >( v );
}

void
Transporter::inthead( const InciterPrint& print )
// *****************************************************************************
// Print out time integration header to screen
//! \param[in] print Pretty printer object to use for printing
// *****************************************************************************
{
  print.inthead( "Time integration", "Navier-Stokes solver",
  "Legend: it - iteration count\n"
  "         t - physics time\n"
  "        dt - physics time step size\n"
  "       ETE - estimated wall-clock time elapsed (h:m:s)\n"
  "       ETA - estimated wall-clock time for accomplishment (h:m:s)\n"
  "       EGT - estimated grind wall-clock time (ms/timestep)\n"
  "       flg - status flags, legend:\n"
  "             f - field (volume and surface)\n"
  "             d - diagnostics\n"
  "             t - physics time history\n"
  "             h - h-refinement\n"
  "             l - load balancing\n"
  "             r - checkpoint\n",
  "\n      it             t            dt        ETE        ETA        EGT  flg\n"
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

  // Finish computing the L2 norm of the residual and append
  const auto scheme = g_inputdeck.get< tag::discr, tag::scheme >();
  std::vector< tk::real > l2res( d[L2RES].size(), 0.0 );
  if (scheme == ctr::SchemeType::DiagCG || scheme == ctr::SchemeType::ALECG)
    for (std::size_t i=0; i<d[L2RES].size(); ++i) {
      l2res[i] = std::sqrt( d[L2RES][i] / m_meshvol );
      diag.push_back( l2res[i] );
    }

  // Append diagnostics file at selected times
  tk::DiagWriter dw( g_inputdeck.get< tag::cmd, tag::io, tag::diag >(),
                     g_inputdeck.get< tag::flformat, tag::diag >(),
                     g_inputdeck.get< tag::prec, tag::diag >(),
                     std::ios_base::app );
  dw.diag( static_cast<uint64_t>(d[ITER][0]), d[TIME][0], d[DT][0], diag );

  // Evaluate whether to continue with next step
  m_scheme.bcast< Scheme::refine >( l2res );
}

void
Transporter::resume()
// *****************************************************************************
// Resume execution from checkpoint/restart files
//! \details This is invoked by Charm++ after the checkpoint is done, as well as
//!   when the restart (returning from a checkpoint) is complete
// *****************************************************************************
{
  if (not m_finished) {
    // If just restarted from a checkpoint, Main( CkMigrateMessage* msg ) has
    // increased nrestart in g_inputdeck, but only on PE 0, so broadcast.
    auto nrestart = g_inputdeck.get< tag::cmd, tag::io, tag::nrestart >();
    m_scheme.bcast< Scheme::evalLB >( nrestart );
  } else
    mainProxy.finalize();
}

void
Transporter::checkpoint( int finished )
// *****************************************************************************
// Save checkpoint/restart files
//! \param[in] finished True if finished with time stepping
// *****************************************************************************
{
  m_finished = finished;

  const auto benchmark = g_inputdeck.get< tag::cmd, tag::benchmark >();

  if (!benchmark) {
    const auto& restart = g_inputdeck.get< tag::cmd, tag::io, tag::restart >();
    CkCallback res( CkIndex_Transporter::resume(), thisProxy );
    CkStartCheckpoint( restart.c_str(), res );
  } else {
    resume();
  }
}

void
Transporter::finish()
// *****************************************************************************
// Normal finish of time stepping
// *****************************************************************************
{
  checkpoint( /* finished = */ 1 );
}

#include "NoWarning/transporter.def.h"
