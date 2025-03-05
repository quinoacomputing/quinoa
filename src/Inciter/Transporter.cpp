// *****************************************************************************
/*!
  \file      src/Inciter/Transporter.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
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

extern ctr::InputDeck g_inputdeck_defaults;
extern ctr::InputDeck g_inputdeck;
extern std::vector< CGPDE > g_cgpde;
extern std::vector< DGPDE > g_dgpde;
extern std::vector< FVPDE > g_fvpde;

}

using inciter::Transporter;

Transporter::Transporter() :
  m_input( input() ),
  m_nchare( m_input.size() ),
  m_meshid(),
  m_ncit( m_nchare.size(), 0 ),
  m_nload( 0 ),
  m_ntrans( 0 ),
  m_ndtmsh( 0 ),
  m_dtmsh(),
  m_npart( 0 ),
  m_nstat( 0 ),
  m_ndisc( 0 ),
  m_nchk( 0 ),
  m_ncom( 0 ),
  m_nt0refit( m_nchare.size(), 0 ),
  m_ndtrefit( m_nchare.size(), 0 ),
  m_noutrefit( m_nchare.size(), 0 ),
  m_noutderefit( m_nchare.size(), 0 ),
  m_scheme(),
  m_partitioner(),
  m_refiner(),
  m_meshwriter(),
  m_sorter(),
  m_nelem( m_nchare.size() ),
  m_npoin(),
  m_finished( m_nchare.size(), 0 ),
  m_meshvol( m_nchare.size() ),
  m_minstat( m_nchare.size() ),
  m_maxstat( m_nchare.size() ),
  m_avgstat( m_nchare.size() ),
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

  const auto nstep = g_inputdeck.get< tag::nstep >();
  const auto t0 = g_inputdeck.get< tag::t0 >();
  const auto term = g_inputdeck.get< tag::term >();
  const auto constdt = g_inputdeck.get< tag::dt >();

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

std::vector< std::string >
Transporter::input()
// *****************************************************************************
// Generate list of input mesh filenames configured by the user
//! \return List of input mesh filenames configured by the user
//! \details If the input file is given on the command line, a single solver
//!   will be instantiated on the single mesh, solving potentially multiple
//!   systems of (potentially coupled) equations. If the input file is not given
//!   on the command line, the mesh files are expected to be configured in the
//!   control/input file, associating a potentially different mesh to each
//!   solver. Both configurations allow the solution of coupled systems, but the
//!   first one solves all equations on the same mesh, while the latter can
//!   couple solutions computed on multiple different meshes.
// *****************************************************************************
{
  // Query input mesh filename specified on the command line
  const auto& cmdinput = g_inputdeck.get< tag::cmd, tag::io, tag::input >();

  // Extract mesh filenames specified in the control file (assigned to solvers)
  std::vector< std::string > ctrinput;
  for (const auto& im : g_inputdeck.get< tag::mesh >()) {
    ctrinput.push_back(im.get< tag::filename >());
  }

  ErrChk( not cmdinput.empty() or not ctrinput.empty(),
    "Either a single input mesh must be given on the command line or multiple "
    "meshes must be configured in the control file." );

   // Prepend control file path to mesh filenames in given in control file
  if (not ctrinput.empty()) {
     const auto& ctr = g_inputdeck.get< tag::cmd, tag::io, tag::control >();
     auto path = ctr.substr( 0, ctr.find_last_of("/")+1 );
     for (auto& f : ctrinput) f = path + f;
  }

  if (cmdinput.empty()) return ctrinput; else return { cmdinput };
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
  print.eqlist( "Registered PDEs using continuous Galerkin (CG) methods",
                stack.cgfactory(), stack.cgntypes() );
  print.eqlist( "Registered PDEs using discontinuous Galerkin (DG) methods",
                stack.dgfactory(), stack.dgntypes() );
  print.eqlist( "Registered PDEs using finite volume (DG) methods",
                stack.fvfactory(), stack.fvntypes() );
  print.endpart();

  // Print out information on problem
  print.part( "Problem" );

  // Print out info on problem title
  if ( !g_inputdeck.get< tag::title >().empty() )
    print.title( g_inputdeck.get< tag::title >() );

  const auto nstep = g_inputdeck.get< tag::nstep >();
  const auto t0 = g_inputdeck.get< tag::t0 >();
  const auto term = g_inputdeck.get< tag::term >();
  const auto constdt = g_inputdeck.get< tag::dt >();
  const auto cfl = g_inputdeck.get< tag::cfl >();
  const auto scheme = g_inputdeck.get< tag::scheme >();

  // Print discretization parameters
  print.section( "Discretization parameters" );
  print.Item< ctr::Scheme, tag::scheme >();
  print.item( "Implicit-Explicit Runge-Kutta",
              g_inputdeck.get< tag::imex_runge_kutta >() );

  if (g_inputdeck.centering() == tk::Centering::ELEM)
  {
    print.Item< ctr::Limiter, tag::limiter >();

    print.item("Shock detection based limiting",
      g_inputdeck.get< tag::shock_detector_coeff >());

    if (g_inputdeck.get< tag::accuracy_test >())
    {
      print.item("WARNING: order-of-accuracy testing enabled",
        "robustness corrections inactive");
    }

    print.item("Limited solution projection",
      g_inputdeck.get< tag::limsol_projection >());
  }
  print.item( "PE-locality mesh reordering",
              g_inputdeck.get< tag::pelocal_reorder >() );
  print.item( "Operator-access mesh reordering",
              g_inputdeck.get< tag::operator_reorder >() );
  auto steady = g_inputdeck.get< tag::steady_state >();
  print.item( "Local time stepping", steady );
  if (steady) {
    print.item( "L2-norm residual convergence criterion",
                g_inputdeck.get< tag::residual >() );
    print.item( "Convergence criterion component index",
                g_inputdeck.get< tag::rescomp >() );
  }
  print.item( "Number of time steps", nstep );
  print.item( "Start time", t0 );
  print.item( "Terminate time", term );

  if (constdt > std::numeric_limits< tk::real >::epsilon())
    print.item( "Constant time step size", constdt );
  else if (cfl > std::numeric_limits< tk::real >::epsilon())
  {
    print.item( "CFL coefficient", cfl );
  }

  const auto& meshes = g_inputdeck.get< tag::mesh >();
  const auto& depvars = g_inputdeck.get< tag::depvar >();
  for (std::size_t i=0; i<meshes.size(); ++i) {
    print.item( "Dependent var name and assoc. mesh", std::string{depvars[i]}
      + " - " + meshes[i].get< tag::filename >() );
  }

  const auto& rbmotion = g_inputdeck.get< tag::rigid_body_motion >();
  if (rbmotion.get< tag::rigid_body_movt >()) {
    const auto& rbdof = rbmotion.get< tag::rigid_body_dof >();
    print.item( "Rigid body motion DOF", rbdof );
    if (rbdof == 3)
      print.item( "Rigid body 3-DOF symmetry plane",
        rbmotion.get< tag::symmetry_plane >() );
  }

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
    auto maxlevels = g_inputdeck.get< tag::amr, tag::maxlevels >();
    print.item( "Maximum mesh refinement levels", maxlevels );
    print.Item< ctr::AMRError, tag::amr, tag::error >();
    auto t0ref = g_inputdeck.get< tag::amr, tag::t0ref >();
    print.item( "Refinement at t<0 (t0ref)", t0ref );
    if (t0ref) {
      const auto& initref = g_inputdeck.get< tag::amr, tag::initial >();
      print.item( "Initial refinement steps", initref.size() );

      auto eps = std::numeric_limits< tk::real >::epsilon();

      const auto& amr_coord = g_inputdeck.get< tag::amr, tag::coords >();
      const auto& amr_defcoord = g_inputdeck_defaults.get< tag::amr, tag::coords >();

      auto xminus = amr_coord.get< tag::xminus >();
      auto xminus_default = amr_defcoord.get< tag::xminus >();
      if (std::abs( xminus - xminus_default ) > eps)
        print.item( "Initial refinement x-", xminus );
      auto xplus = amr_coord.get< tag::xplus >();
      auto xplus_default = amr_defcoord.get< tag::xplus >();
      if (std::abs( xplus - xplus_default ) > eps)
        print.item( "Initial refinement x+", xplus );

      auto yminus = amr_coord.get< tag::yminus >();
      auto yminus_default = amr_defcoord.get< tag::yminus >();
      if (std::abs( yminus - yminus_default ) > eps)
        print.item( "Initial refinement y-", yminus );
      auto yplus = amr_coord.get< tag::yplus >();
      auto yplus_default = amr_defcoord.get< tag::yplus >();
      if (std::abs( yplus - yplus_default ) > eps)
        print.item( "Initial refinement y+", yplus );

      auto zminus = amr_coord.get< tag::zminus >();
      auto zminus_default = amr_defcoord.get< tag::zminus >();
      if (std::abs( zminus - zminus_default ) > eps)
        print.item( "Initial refinement z-", zminus );
      auto zplus = amr_coord.get< tag::zplus >();
      auto zplus_default = amr_defcoord.get< tag::zplus >();
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
                g_inputdeck.get< tag::amr, tag::tol_refine >() );
    print.item( "De-refinement tolerance",
                g_inputdeck.get< tag::amr, tag::tol_derefine >() );
  }

  // Print out ALE configuration
  const auto ale = g_inputdeck.get< tag::ale, tag::ale >();
  if (ale) {
    print.section( "Arbitrary Lagrangian-Eulerian (ALE) mesh motion" );
    auto dvcfl = g_inputdeck.get< tag::ale, tag::dvcfl >();
    print.item( "Volume-change CFL coefficient", dvcfl );
    print.Item< ctr::MeshVelocity, tag::ale, tag::mesh_velocity >();
    print.Item< ctr::MeshVelocitySmoother, tag::ale, tag::smoother >();
    print.item( "Mesh motion dimensions", tk::parameters(
                g_inputdeck.get< tag::ale, tag::mesh_motion >() ) );
    const auto& meshforce = g_inputdeck.get< tag::ale, tag::meshforce >();
    print.item( "Mesh velocity force coefficients", tk::parameters(meshforce) );
    print.item( "Vorticity multiplier",
                g_inputdeck.get< tag::ale, tag::vortmult >() );
    print.item( "Mesh velocity linear solver tolerance",
                g_inputdeck.get< tag::ale, tag::tolerance >() );
    print.item( "Mesh velocity linear solver maxit",
                g_inputdeck.get< tag::ale, tag::maxit >() );
    const auto& dir = g_inputdeck.get< tag::ale, tag::dirichlet >();
    if (not dir.empty())
      print.item( "Mesh velocity Dirichlet BC sideset(s)",
                  tk::parameters( dir ) );
    const auto& sym = g_inputdeck.get< tag::ale, tag::symmetry >();
    if (not sym.empty())
      print.item( "Mesh velocity symmetry BC sideset(s)",
                  tk::parameters( sym ) );
    std::size_t i = 1;
    for (const auto& m : g_inputdeck.get< tag::ale, tag::move >()) {
       tk::ctr::UserTable opt;
       print.item( opt.group() + ' ' + std::to_string(i) + " interpreted as",
                   opt.name( m.get< tag::fntype >() ) );
       const auto& s = m.get< tag::sideset >();
       if (not s.empty())
         print.item( "Moving sideset(s) with table " + std::to_string(i),
                     tk::parameters(s));
       ++i;
    }
  }

  // Print I/O filenames
  print.section( "Input/Output filenames and directories" );
  print.item( "Input mesh(es)", tk::parameters( m_input ) );
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
  print.section( "Output intervals (in units of iteration count)" );
  print.item( "TTY", g_inputdeck.get< tag::ttyi>() );
  print.item( "Field and surface",
              g_inputdeck.get< tag::field_output, tag::interval >() );
  print.item( "History",
              g_inputdeck.get< tag::history_output, tag::interval >() );
  print.item( "Diagnostics",
              g_inputdeck.get< tag::diagnostics, tag::interval >() );
  print.item( "Checkpoint/restart",
              g_inputdeck.get< tag::cmd, tag::rsfreq >() );
  auto tf = g_inputdeck.get< tag::field_output, tag::time_interval >();
  auto th = g_inputdeck.get< tag::history_output, tag::time_interval >();
  if (tf>0.0 || th>0.0) {
    print.section( "Output intervals (in units of physics time)" );
    if (tf > 0.0) print.item( "Field and surface", tf );
    if (th > 0.0) print.item( "History", th );
  }
  const auto& rf = g_inputdeck.get< tag::field_output, tag::time_range >();
  const auto& rh = g_inputdeck.get< tag::history_output, tag::time_range >();
  if (not rf.empty() or not rh.empty()) {
    print.section( "Output time ranges (in units of physics time)" );
    print.item("Field output { mintime, maxtime, dt }", tk::parameters(rf));
    print.item("History output { mintime, maxtime, dt }", tk::parameters(rh));
  }

  //// Print output variables: fields and surfaces
  //const auto nodeoutvars = g_inputdeck.outvars( tk::Centering::NODE );
  //const auto elemoutvars = g_inputdeck.outvars( tk::Centering::ELEM );
  //const auto outsets = g_inputdeck.outsets();
  //if (!nodeoutvars.empty() || !elemoutvars.empty() || !outsets.empty())
  //   print.section( "Output fields" );
  //if (!nodeoutvars.empty())
  //  print.item( "Node field(s)", tk::parameters(nodeoutvars) );
  //if (!elemoutvars.empty())
  //  print.item( "Elem field(s)", tk::parameters(elemoutvars) );
  //if (!aliases.empty())
  //  print.item( "Alias(es)", tk::parameters(aliases) );
  //if (!outsets.empty())
  //  print.item( "Surface side set(s)", tk::parameters(outsets) );

  // Print output variables: history
  const auto& pt = g_inputdeck.get< tag::history_output, tag::point >();
  if (!pt.empty()) {
    print.section( "Output time history" );
    for (std::size_t p=0; p<pt.size(); ++p) {
      const auto& id = pt[p].get< tag::id >();
      std::stringstream ss;
      auto prec = g_inputdeck.get< tag::history_output, tag::precision >();
      ss << std::setprecision( static_cast<int>(prec) );
      ss << of << ".hist." << id;
      print.longitem( "At point " + id + ' ' +
        tk::parameters(pt[p].get<tag::coord>()), ss.str() );
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
  std::unordered_set< int > usedsets;

  using bclist = ctr::bclist::Keys;
  const auto& bcs = g_inputdeck.get< tag::bc >();
  brigand::for_each< bclist >( UserBC( g_inputdeck, usedsets ) );

  // Query side sets of time dependent and inlet BCs (since these are not a part
  // of tag::bc)
  for (const auto& bci : bcs) {
    for (const auto& b : bci.get< tag::inlet >()) {
      for (auto i : b.get< tag::sideset >())
        usedsets.insert(static_cast<int>(i));
    }
    for (const auto& b : bci.get< tag::timedep >()) {
      for (auto i : b.get< tag::sideset >())
        usedsets.insert(static_cast<int>(i));
    }
  }

  // Query side sets of boundaries prescribed as moving with ALE
  for (const auto& move : g_inputdeck.get< tag::ale, tag::move >())
    for (auto i : move.get< tag::sideset >())
      usedsets.insert(static_cast<int>(i));

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
    //else {
    //  Throw( "Boundary conditions specified on side set " +
    //    std::to_string(i) + " which does not exist in mesh file" );
    //}
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
  auto scheme = g_inputdeck.get< tag::scheme >();
  auto centering = ctr::Scheme().centering( scheme );
  auto print = printer();

  // Create partitioner callbacks (order important)
  tk::PartitionerCallback cbp {{
      CkCallback( CkReductionTarget(Transporter,load), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,partitioned), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,distributed), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,refinserted), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,refined), thisProxy )
  }};

  // Create refiner callbacks (order important)
  tk::RefinerCallback cbr {{
      CkCallback( CkReductionTarget(Transporter,queriedRef), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,respondedRef), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,compatibility), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,bndint), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,matched), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,refined), thisProxy )
  }};

  // Create sorter callbacks (order important)
  tk::SorterCallback cbs {{
      CkCallback( CkReductionTarget(Transporter,queried), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,responded), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,discinserted), thisProxy )
    , CkCallback( CkReductionTarget(Transporter,workinserted), thisProxy )
  }};

  // Start timer measuring preparation of mesh(es) for partitioning
  m_timer[ TimerTag::MESH_READ ];

  // Start preparing mesh(es)
  print.diag( "Reading mesh(es)" );

  // Create (discretization) Scheme chare worker arrays for all meshes
  for ([[maybe_unused]] const auto& filename : m_input)
    m_scheme.emplace_back( g_inputdeck.get< tag::scheme >(),
                           g_inputdeck.get< tag::ale, tag::ale >(),
                           need_linearsolver(),
                           centering );

  ErrChk( !m_input.empty(), "No input mesh" );

  // Read boundary (side set) data from a list of input mesh files
  std::size_t meshid = 0;
  for (const auto& filename : m_input) {
    // Create mesh reader for reading side sets from file
    tk::MeshReader mr( filename );

    // Read out total number of mesh points from mesh file
    m_npoin.push_back( mr.npoin() );

    std::map< int, std::vector< std::size_t > > bface;
    std::map< int, std::vector< std::size_t > > faces;
    std::map< int, std::vector< std::size_t > > bnode;

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
      bool bcnode_set = matchBCs( bnode );
      bool bcface_set = matchBCs( bface );
      bcs_set = bcface_set or bcnode_set;
    }

    // Warn on no BCs
    if (!bcs_set) print << "\n>>> WARNING: No boundary conditions set\n\n";

    auto opt = m_scheme[meshid].arrayoptions();
    // Create empty mesh refiner chare array (bound to workers)
    m_refiner.push_back( CProxy_Refiner::ckNew(opt) );
    // Create empty mesh sorter Charm++ chare array (bound to workers)
    m_sorter.push_back( CProxy_Sorter::ckNew(opt) );

    // Create MeshWriter chare group for mesh
    m_meshwriter.push_back(
      tk::CProxy_MeshWriter::ckNew(
        g_inputdeck.get< tag::field_output, tag::filetype >(),
        centering,
        g_inputdeck.get< tag::cmd, tag::benchmark >(),
        m_input.size() ) );

    // Create mesh partitioner Charm++ chare nodegroup for all meshes
    m_partitioner.push_back(
      CProxy_Partitioner::ckNew( meshid, filename, cbp, cbr, cbs,
        thisProxy, m_refiner.back(), m_sorter.back(), m_meshwriter.back(),
        m_scheme, bface, faces, bnode ) );

    ++meshid;
  }
}

void
Transporter::load( std::size_t meshid, std::size_t nelem )
// *****************************************************************************
// Reduction target: the mesh has been read from file on all PEs
//! \param[in] meshid Mesh id (summed accross all compute nodes)
//! \param[in] nelem Number of mesh elements per mesh (summed across all
//!    compute nodes)
// *****************************************************************************
{
  meshid /= static_cast< std::size_t >( CkNumNodes() );
  Assert( meshid < m_nelem.size(), "MeshId indexing out" );
  m_nelem[meshid] = nelem;

  // Compute load distribution given total work (nelem) and user-specified
  // virtualization
  uint64_t chunksize, remainder;
  m_nchare[meshid] = static_cast<int>(
    tk::linearLoadDistributor(
       g_inputdeck.get< tag::cmd, tag::virtualization >(),
       m_nelem[meshid], CkNumPes(), chunksize, remainder ) );

  // Store sum of meshids (across all chares, key) for each meshid (value).
  // This is used to look up the mesh id after collectives that sum their data.
  m_meshid[ static_cast<std::size_t>(m_nchare[meshid])*meshid ] = meshid;
  Assert( meshid < m_nelem.size(), "MeshId indexing out" );

  // Tell the meshwriter for this mesh the total number of its chares
  m_meshwriter[meshid].nchare( meshid, m_nchare[meshid] );

  if (++m_nload == m_nelem.size()) {     // all meshes have been loaded
    m_nload = 0;
    auto print = printer();

    // Start timer measuring preparation of the mesh for partitioning
    const auto& timer = tk::cref_find( m_timer, TimerTag::MESH_READ );
    print.diag( "Mesh read time: " + std::to_string( timer.dsec() ) + " sec" );

    // Print out mesh partitioning configuration
    print.section( "Mesh partitioning" );
    print.Item< tk::ctr::PartitioningAlgorithm, tag::partitioning >();
    print.item( "Virtualization [0.0...1.0]",
                g_inputdeck.get< tag::cmd, tag::virtualization >() );
    // Print out initial mesh statistics
    meshstat( "Initial load distribution" );

    // Query number of initial mesh refinement steps
    int nref = 0;
    if (g_inputdeck.get< tag::amr, tag::t0ref >())
      nref = static_cast<int>(g_inputdeck.get< tag::amr,
        tag::initial >().size());

    m_progMesh.start( print, "Preparing mesh", {{ CkNumPes(), CkNumPes(), nref,
      m_nchare[0], m_nchare[0], m_nchare[0], m_nchare[0] }} );

    // Partition first mesh
    m_partitioner[0].partition( m_nchare[0] );
  }
}

void
Transporter::partitioned( std::size_t meshid )
// *****************************************************************************
// Reduction target: a mesh has been partitioned
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  if (++m_npart == m_nelem.size()) {     // all meshes have been partitioned
    m_npart = 0;
  } else { // partition next mesh
    m_partitioner[meshid+1].partition( m_nchare[meshid+1] );
  }
}

void
Transporter::distributed( std::size_t meshid )
// *****************************************************************************
// Reduction target: all compute nodes have distributed their mesh after
// partitioning
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  m_partitioner[meshid].refine();
}

void
Transporter::refinserted( std::size_t meshid, std::size_t error )
// *****************************************************************************
// Reduction target: all compute nodes have created the mesh refiners
//! \param[in] meshid Mesh id (aggregated across all compute nodes with operator
//!   max)
//! \param[in] error Error code (aggregated across all compute nodes with
//!   operator max)
// *****************************************************************************
{
  if (error) {

    printer() << "\n>>> ERROR: A worker chare was not assigned any mesh "
              "elements after distributing mesh " + std::to_string(meshid) +
              ". This can happen in SMP-mode with a large +ppn "
              "parameter (number of worker threads per logical node) and is "
              "most likely the fault of the mesh partitioning algorithm not "
              "tolerating the case when it is asked to divide the "
              "computational domain into a number of partitions different "
              "than the number of ranks it is called on, i.e., in case of "
              "overdecomposition and/or calling the partitioner in SMP mode "
              "with +ppn larger than 1. Solution 1: Try a different "
              "partitioning algorithm (e.g., rcb instead of mj). Solution 2: "
              "Decrease +ppn.";
    finish( meshid );

  } else {

     m_refiner[meshid].doneInserting();

  }
}

void
Transporter::queriedRef( std::size_t meshid )
// *****************************************************************************
// Reduction target: all Refiner chares have queried their boundary edges
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  m_refiner[meshid].response();
}

void
Transporter::respondedRef( std::size_t meshid )
// *****************************************************************************
// Reduction target: all Refiner chares have setup their boundary edges
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  m_refiner[meshid].refine();
}

void
Transporter::compatibility( std::size_t meshid )
// *****************************************************************************
// Reduction target: all Refiner chares have received a round of edges,
// and have run their compatibility algorithm
//! \param[in] meshid Mesh id (aggregated across all chares using operator max)
//! \details This is called iteratively, until convergence by Refiner. At this
//!   point all Refiner chares have received a round of edge data (tags whether
//!   an edge needs to be refined, etc.), and applied the compatibility
//!   algorithm independent of other Refiner chares. We keep going until the
//!   mesh is no longer modified by the compatibility algorithm, based on a new
//!   round of edge data communication started in Refiner::comExtra().
// *****************************************************************************
{
  m_refiner[meshid].correctref();
}

void
Transporter::matched( std::size_t summeshid,
                      std::size_t nextra,
                      std::size_t nref,
                      std::size_t nderef,
                      std::size_t sumrefmode )
// *****************************************************************************
// Reduction target: all Refiner chares have matched/corrected the tagging
// of chare-boundary edges, all chares are ready to perform refinement.
//! \param[in] summeshid Mesh id (summed across all chares)
//! \param[in] nextra Sum (across all chares) of the number of edges on each
//!   chare that need correction along chare boundaries
//! \param[in] nref Sum of number of refined tetrahedra across all chares.
//! \param[in] nderef Sum of number of derefined tetrahedra across all chares.
//! \param[in] sumrefmode Sum of contributions from all chares, encoding
//!   refinement mode of operation.
// *****************************************************************************
{
  auto meshid = tk::cref_find( m_meshid, summeshid );

  // If at least a single edge on a chare still needs correction, do correction,
  // otherwise, this mesh refinement step is complete
  if (nextra > 0) {

    ++m_ncit[meshid];
    m_refiner[meshid].comExtra();

  } else {

    auto print = printer();

    // decode refmode
    auto refmode = static_cast< Refiner::RefMode >(
                     sumrefmode / static_cast<std::size_t>(m_nchare[meshid]) );

    if (refmode == Refiner::RefMode::T0REF) {

      if (!g_inputdeck.get< tag::cmd, tag::feedback >()) {
        ctr::AMRInitial opt;
        print.diag( { "meshid", "t0ref", "type", "nref", "nderef", "ncorr" },
                    { std::to_string(meshid),
                      std::to_string(m_nt0refit[meshid]),
                      "initial",
                      std::to_string(nref),
                      std::to_string(nderef),
                      std::to_string(m_ncit[meshid]) } );
        ++m_nt0refit[meshid];
      }
      m_progMesh.inc< REFINE >( print );

    } else if (refmode == Refiner::RefMode::DTREF) {

      auto dtref_uni = g_inputdeck.get< tag::amr, tag::dtref_uniform >();
      print.diag( { "meshid", "dtref", "type", "nref", "nderef", "ncorr" },
                  { std::to_string(meshid),
                    std::to_string(++m_ndtrefit[meshid]),
                    (dtref_uni?"uniform":"error"),
                    std::to_string(nref),
                    std::to_string(nderef),
                    std::to_string(m_ncit[meshid]) },
                  false );

    } else if (refmode == Refiner::RefMode::OUTREF) {

      print.diag( { "meshid", "outref", "nref", "nderef", "ncorr" },
                  { std::to_string(meshid),
                    std::to_string(++m_noutrefit[meshid]),
                    std::to_string(nref),
                    std::to_string(nderef),
                    std::to_string(m_ncit[meshid]) }, false );

    } else if (refmode == Refiner::RefMode::OUTDEREF) {

      print.diag( { "meshid", "outderef", "nref", "nderef", "ncorr" },
                  { std::to_string(meshid),
                    std::to_string(++m_noutderefit[meshid]),
                    std::to_string(nref),
                    std::to_string(nderef),
                    std::to_string(m_ncit[meshid]) },
                  false );

    } else Throw( "RefMode not implemented" );

    m_ncit[meshid] = 0;
    m_refiner[meshid].perform();

  }
}

void
Transporter::bndint( tk::real sx, tk::real sy, tk::real sz, tk::real cb,
                     tk::real summeshid )
// *****************************************************************************
// Compute surface integral across the whole problem and perform leak-test
//! \param[in] sx X component of vector summed
//! \param[in] sy Y component of vector summed
//! \param[in] sz Z component of vector summed
//! \param[in] cb Invoke callback if positive
//! \param[in] summeshid Mesh id (summed accross all chares)
//! \details This function aggregates partial surface integrals across the
//!   boundary faces of the whole problem. After this global sum a
//!   non-zero vector result indicates a leak, e.g., a hole in the boundary,
//!   which indicates an error in the boundary face data structures used to
//!   compute the partial surface integrals.
// *****************************************************************************
{
  auto meshid = tk::cref_find( m_meshid, static_cast<std::size_t>(summeshid) );

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

  if (cb > 0.0) m_scheme[meshid].ghosts().resizeComm();
}

void
Transporter::refined( std::size_t summeshid,
                      std::size_t nelem,
                      std::size_t npoin )
// *****************************************************************************
// Reduction target: all chares have refined their mesh
//! \param[in] summeshid Mesh id (summed accross all Refiner chares)
//! \param[in] nelem Total number of elements in mesh summed across the
//!   distributed mesh
//! \param[in] npoin Total number of mesh points summed across the distributed
//!   mesh. Note that in parallel this is larger than the number of points in
//!   the mesh, because the boundary nodes are multi-counted. But we only need
//!   an equal or larger than npoin for Sorter::setup, so this is okay.
// *****************************************************************************
{
  auto meshid = tk::cref_find( m_meshid, summeshid );

  // Store new number of elements for initially refined mesh
  m_nelem[meshid] = nelem;

  m_sorter[meshid].doneInserting();
  m_sorter[meshid].setup( npoin );
}

void
Transporter::queried( std::size_t meshid )
// *****************************************************************************
// Reduction target: all Sorter chares have queried their boundary edges
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  m_sorter[meshid].response();
}

void
Transporter::responded( std::size_t meshid )
// *****************************************************************************
// Reduction target: all Sorter chares have responded with their boundary edges
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  m_sorter[meshid].start();
}

void
Transporter::resized( std::size_t meshid )
// *****************************************************************************
// Reduction target: all worker chares have resized their own mesh data after
// AMR or ALE
//! \param[in] meshid Mesh id
//! \note Only used for nodal schemes
// *****************************************************************************
{
  m_scheme[meshid].disc().vol();
  m_scheme[meshid].bcast< Scheme::lhs >();
}

void
Transporter::startEsup( std::size_t meshid )
// *****************************************************************************
// Reduction target: all worker chares have generated their own esup
//! \param[in] meshid Mesh id
//! \note Only used for cell-centered schemes
// *****************************************************************************
{
  m_scheme[meshid].ghosts().nodeNeighSetup();
}

void
Transporter::discinserted( std::size_t meshid )
// *****************************************************************************
// Reduction target: all Discretization chares have been inserted
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  m_scheme[meshid].disc().doneInserting();
}

void
Transporter::meshstat( const std::string& header ) const
// *****************************************************************************
// Print out mesh statistics
//! \param[in] header Section header
// *****************************************************************************
{
  auto print = printer();

  print.section( header );

  if (m_nelem.size() > 1) {
    print.item( "Number of tetrahedra (per mesh)",tk::parameters(m_nelem) );
    print.item( "Number of points (per mesh)", tk::parameters(m_npoin) );
    print.item( "Number of work units (per mesh)", tk::parameters(m_nchare) );
  }

  print.item( "Total number of tetrahedra",
              std::accumulate( begin(m_nelem), end(m_nelem), 0UL ) );
  print.item( "Total number of points",
              std::accumulate( begin(m_npoin), end(m_npoin), 0UL ) );
  print.item( "Total number of work units",
              std::accumulate( begin(m_nchare), end(m_nchare), 0 ) );

  print.endsubsection();
}

bool
Transporter::need_linearsolver() const
// *****************************************************************************
//  Decide if we need a linear solver for ALE
//! \return True if ALE will neeed a linear solver
// *****************************************************************************
{
  auto smoother = g_inputdeck.get< tag::ale, tag::smoother >();

  if ( g_inputdeck.get< tag::ale, tag::ale >() and
       (smoother == ctr::MeshVelocitySmootherType::LAPLACE ||
        smoother == ctr::MeshVelocitySmootherType::HELMHOLTZ) )
  {
     return true;
  } else {
     return false;
  }
}

void
Transporter::disccreated( std::size_t summeshid, std::size_t npoin )
// *****************************************************************************
// Reduction target: all Discretization constructors have been called
//! \param[in] summeshid Mesh id (summed accross all chares)
//! \param[in] npoin Total number of mesh points (summed across all chares)
//!  Note that as opposed to npoin in refined(), this npoin is not
//!  multi-counted, and thus should be correct in parallel.
// *****************************************************************************
{
  auto meshid = tk::cref_find( m_meshid, summeshid );
  //std::cout << "Trans: " << meshid << " Transporter::disccreated()\n";

  // Update number of mesh points for mesh, since it may have been refined
  if (g_inputdeck.get< tag::amr, tag::t0ref >()) m_npoin[meshid] = npoin;

  if (++m_ndisc == m_nelem.size()) { // all Disc arrays have been created
    m_ndisc = 0;
    auto print = printer();
    m_progMesh.end( print );
    if (g_inputdeck.get< tag::amr, tag::t0ref >())
      meshstat( "Initially (t<0) refined mesh graph statistics" );
  }

  m_refiner[meshid].sendProxy();

  if (g_inputdeck.get< tag::ale, tag::ale >())
    m_scheme[meshid].ale().doneInserting();

  if (need_linearsolver())
    m_scheme[meshid].conjugategradients().doneInserting();

  m_scheme[meshid].disc().vol();
}

void
Transporter::workinserted( std::size_t meshid )
// *****************************************************************************
// Reduction target: all worker (derived discretization) chares have been
// inserted
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  m_scheme[meshid].bcast< Scheme::doneInserting >();
}

void
Transporter::diagHeader()
// *****************************************************************************
// Configure and write diagnostics file header
// *****************************************************************************
{
  for (std::size_t m=0; m<m_input.size(); ++m) {

   // Output header for diagnostics output file
    auto mid = m_input.size() > 1 ? std::string( '.' + std::to_string(m) ) : "";
    tk::DiagWriter dw( g_inputdeck.get< tag::cmd, tag::io, tag::diag >() + mid,
      g_inputdeck.get< tag::diagnostics, tag::format >(),
      g_inputdeck.get< tag::diagnostics, tag::precision >() );

    // Collect variables names for integral/diagnostics output
    std::vector< std::string > var;
    const auto scheme = g_inputdeck.get< tag::scheme >();
    if ( scheme == ctr::SchemeType::ALECG ||
         scheme == ctr::SchemeType::OversetFE )
      for (const auto& eq : g_cgpde) varnames( eq, var );
    else if ( scheme == ctr::SchemeType::DG ||
              scheme == ctr::SchemeType::P0P1 ||
              scheme == ctr::SchemeType::DGP1 ||
              scheme == ctr::SchemeType::DGP2 ||
              scheme == ctr::SchemeType::PDG )
      for (const auto& eq : g_dgpde) varnames( eq, var );
    else if (scheme == ctr::SchemeType::FV)
      for (const auto& eq : g_fvpde) varnames( eq, var );
    else Throw( "Diagnostics header not handled for discretization scheme" );

    const tk::ctr::Error opt;
    auto nv = var.size() / m_input.size();
    std::vector< std::string > d;

    // Add 'L2(var)' for all variables as those are always computed
    const auto& l2name = opt.name( tk::ctr::ErrorType::L2 );
    for (std::size_t i=0; i<nv; ++i) d.push_back( l2name + '(' + var[i] + ')' );

    // Query user-requested diagnostics and augment diagnostics file header by
    // 'err(var)', where 'err' is the error type  configured, and var is the
    // variable computed, for all variables and all error types configured.
    const auto& err = g_inputdeck.get< tag::diagnostics, tag::error >();
    const auto& errname = opt.name( err );
    for (std::size_t i=0; i<nv; ++i)
      d.push_back( errname + '(' + var[i] + "-IC)" );

    // Augment diagnostics variables by L2-norm of the residual and total energy
    if ( scheme == ctr::SchemeType::ALECG ||
         scheme == ctr::SchemeType::OversetFE ||
         scheme == ctr::SchemeType::FV )
    {
      for (std::size_t i=0; i<nv; ++i) d.push_back( "L2(d" + var[i] + ')' );
    }
    d.push_back( "mE" );

    // Write diagnostics header
    dw.header( d );

  }
}

void
Transporter::doneInsertingGhosts(std::size_t meshid)
// *****************************************************************************
// Reduction target indicating all "ghosts" insertions are done
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  m_scheme[meshid].ghosts().doneInserting();
  m_scheme[meshid].ghosts().startCommSetup();
}

void
Transporter::comfinal( std::size_t initial, std::size_t summeshid )
// *****************************************************************************
// Reduction target indicating that communication maps have been setup
//! \param[in] initial Sum of contributions from all chares. If larger than
//!    zero, we are during time stepping and if zero we are during setup.
//! \param[in] summeshid Mesh id (summed accross the distributed mesh)
// *****************************************************************************
// [Discretization-specific communication maps]
{
  auto meshid = tk::cref_find( m_meshid, static_cast<std::size_t>(summeshid) );

  if (initial > 0) {
    m_scheme[meshid].bcast< Scheme::setup >();
    // Turn on automatic load balancing
    if (++m_ncom == m_nelem.size()) { // all worker arrays have finished
      m_ncom = 0;
      auto print = printer();
      m_progWork.end( print );
      tk::CProxy_LBSwitch::ckNew();
      print.diag( "Load balancing on (if enabled in Charm++)" );
    }
  } else {
    m_scheme[meshid].bcast< Scheme::lhs >();
  }
}
// [Discretization-specific communication maps]

void
Transporter::totalvol( tk::real v, tk::real initial, tk::real summeshid )
// *****************************************************************************
// Reduction target summing total mesh volume across all workers
//! \param[in] v Mesh volume summed across the distributed mesh
//! \param[in] initial Sum of contributions from all chares. If larger than
//!    zero, we are during setup, if zero, during time stepping.
//! \param[in] summeshid Mesh id (summed accross the distributed mesh)
// *****************************************************************************
{
  auto meshid = tk::cref_find( m_meshid, static_cast<std::size_t>(summeshid) );

  m_meshvol[meshid] = v;

  if (initial > 0.0)   // during initialization
    m_scheme[meshid].disc().stat( v );
  else                  // during ALE or AMR
    m_scheme[meshid].bcast< Scheme::resized >();
}

void
Transporter::minstat( tk::real d0, tk::real d1, tk::real d2, tk::real rmeshid )
// *****************************************************************************
// Reduction target yielding minimum mesh statistcs across all workers
//! \param[in] d0 Minimum mesh statistics collected over all chares
//! \param[in] d1 Minimum mesh statistics collected over all chares
//! \param[in] d2 Minimum mesh statistics collected over all chares
//! \param[in] rmeshid Mesh id as a real
// *****************************************************************************
{
  auto meshid = static_cast<std::size_t>(rmeshid);

  m_minstat[meshid][0] = d0;  // minimum edge length
  m_minstat[meshid][1] = d1;  // minimum cell volume cubic root
  m_minstat[meshid][2] = d2;  // minimum number of cells on chare

  minstat_complete(meshid);
}

void
Transporter::maxstat( tk::real d0, tk::real d1, tk::real d2, tk::real rmeshid )
// *****************************************************************************
// Reduction target yielding the maximum mesh statistics across all workers
//! \param[in] d0 Maximum mesh statistics collected over all chares
//! \param[in] d1 Maximum mesh statistics collected over all chares
//! \param[in] d2 Maximum mesh statistics collected over all chares
//! \param[in] rmeshid Mesh id as a real
// *****************************************************************************
{
  auto meshid = static_cast<std::size_t>(rmeshid);

  m_maxstat[meshid][0] = d0;  // maximum edge length
  m_maxstat[meshid][1] = d1;  // maximum cell volume cubic root
  m_maxstat[meshid][2] = d2;  // maximum number of cells on chare

  maxstat_complete(meshid);
}

void
Transporter::sumstat( tk::real d0, tk::real d1, tk::real d2, tk::real d3,
                      tk::real d4, tk::real d5, tk::real summeshid )
// *****************************************************************************
// Reduction target yielding the sum mesh statistics across all workers
//! \param[in] d0 Sum mesh statistics collected over all chares
//! \param[in] d1 Sum mesh statistics collected over all chares
//! \param[in] d2 Sum mesh statistics collected over all chares
//! \param[in] d3 Sum mesh statistics collected over all chares
//! \param[in] d4 Sum mesh statistics collected over all chares
//! \param[in] d5 Sum mesh statistics collected over all chares
//! \param[in] summeshid Mesh id (summed accross the distributed mesh)
// *****************************************************************************
{
  auto meshid = tk::cref_find( m_meshid, static_cast<std::size_t>(summeshid) );

  m_avgstat[meshid][0] = d1 / d0;      // average edge length
  m_avgstat[meshid][1] = d3 / d2;      // average cell volume cubic root
  m_avgstat[meshid][2] = d5 / d4;      // average number of cells per chare

  sumstat_complete(meshid);
}

void
Transporter::pdfstat( CkReductionMsg* msg )
// *****************************************************************************
// Reduction target yielding PDF of mesh statistics across all workers
//! \param[in] msg Serialized PDF
// *****************************************************************************
{
  std::size_t meshid;
  std::vector< tk::UniPDF > pdf;

  // Deserialize final PDF
  PUP::fromMem creator( msg->getData() );
  creator | meshid;
  creator | pdf;
  delete msg;

  auto id = std::to_string(meshid);

  // Create new PDF file (overwrite if exists)
  tk::PDFWriter pdfe( "mesh_edge_pdf." + id + ".txt" );
  // Output edgelength PDF
  // cppcheck-suppress containerOutOfBounds
  pdfe.writeTxt( pdf[0],
                 tk::ctr::PDFInfo{ {"PDF"}, {}, {"edgelength"}, 0, 0.0 } );

  // Create new PDF file (overwrite if exists)
  tk::PDFWriter pdfv( "mesh_vol_pdf." + id + ".txt" );
  // Output cell volume cubic root PDF
  pdfv.writeTxt( pdf[1],
                 tk::ctr::PDFInfo{ {"PDF"}, {}, {"V^{1/3}"}, 0, 0.0 } );

  // Create new PDF file (overwrite if exists)
  tk::PDFWriter pdfn( "mesh_ntet_pdf." + id + ".txt" );
  // Output number of cells PDF
  pdfn.writeTxt( pdf[2],
                 tk::ctr::PDFInfo{ {"PDF"}, {}, {"ntets"}, 0, 0.0 } );

  pdfstat_complete(meshid);
}

void
Transporter::stat()
// *****************************************************************************
// Echo diagnostics on mesh statistics
// *****************************************************************************
{
  auto print = printer();

  if (++m_nstat == m_nelem.size()) {     // stats from all meshes have arrived
    m_nstat = 0;
    for (std::size_t i=0; i<m_nelem.size(); ++i) {
      print.diag(
        "Mesh " + std::to_string(i) +
        " distribution statistics: min/max/avg(edgelength) = " +
        std::to_string( m_minstat[i][0] ) + " / " +
        std::to_string( m_maxstat[i][0] ) + " / " +
        std::to_string( m_avgstat[i][0] ) + ", " +
        "min/max/avg(V^{1/3}) = " +
        std::to_string( m_minstat[i][1] ) + " / " +
        std::to_string( m_maxstat[i][1] ) + " / " +
        std::to_string( m_avgstat[i][1] ) + ", " +
        "min/max/avg(ntets) = " +
        std::to_string( static_cast<std::size_t>(m_minstat[i][2]) ) + " / " +
        std::to_string( static_cast<std::size_t>(m_maxstat[i][2]) ) + " / " +
        std::to_string( static_cast<std::size_t>(m_avgstat[i][2]) ) );
    }

    // Print out time integration header to screen
    inthead( print );

    m_progWork.start( print, "Preparing workers",
      {{ m_nchare[0], m_nchare[0], m_nchare[0], m_nchare[0], m_nchare[0] }} );

    // Create "derived-class" workers
    for (std::size_t i=0; i<m_nelem.size(); ++i) m_sorter[i].createWorkers();
  }
}

void
Transporter::boxvol( tk::real* meshdata, int n )
// *****************************************************************************
// Reduction target computing total volume of IC mesh blocks and box
//! \param[in] meshdata Vector containing volumes of all IC mesh blocks,
//!   volume of IC box, and mesh id as a real summed across the distributed mesh
//! \param[in] n Size of vector, automatically computed by Charm
// *****************************************************************************
{
  Assert(n>=2, "mesh data size incorrect");

  // extract summed mesh id from vector
  tk::real summeshid = meshdata[n-1];
  auto meshid = tk::cref_find( m_meshid, static_cast<std::size_t>(summeshid) );

  // extract summed box volume from vector
  tk::real v = meshdata[n-2];
  if (v > 0.0) printer().diag( "Box IC volume: " + std::to_string(v) );

  // extract summed mesh block volumes from the vector
  std::vector< tk::real > blockvols;
  for (std::size_t blid=0; blid<(static_cast<std::size_t>(n)-2); ++blid) {
    blockvols.push_back(meshdata[blid]);
    if (blockvols[blid] > 0.0)
      printer().diag( "Mesh block " + std::to_string(blid) +
        " discrete volume: " + std::to_string(blockvols[blid]) );
  }

  m_scheme[meshid].bcast< Scheme::box >( v, blockvols );
}

void
Transporter::solutionTransferred()
// *****************************************************************************
// Reduction target broadcasting to Schemes after mesh transfer
// *****************************************************************************
{
  if (++m_ntrans == m_nelem.size()) {    // all meshes have been loaded
    m_ntrans = 0;
    for (auto& m : m_scheme) m.bcast< Scheme::transferSol >();
  }
}

void
Transporter::minDtAcrossMeshes( CkReductionMsg* advMsg )
// *****************************************************************************
// Reduction target that computes minimum timestep across all meshes
//! \param[in] advMsg Reduction msg containing minimum timestep and total
//!   surface force information
// *****************************************************************************
{
  // obtain results of reduction from reduction-msg
  CkReduction::tupleElement* results = nullptr;
  int num_reductions = 0;
  advMsg->toTuple(&results, &num_reductions);

// ignore the old-style-cast warning from clang for this code
#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wold-style-cast"
#endif

  tk::real mindt = *(tk::real*)results[0].data;
  std::array< tk::real, 3 > F;
  F[0] = *(tk::real*)results[1].data;
  F[1] = *(tk::real*)results[2].data;
  F[2] = *(tk::real*)results[3].data;

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

  m_dtmsh.push_back(mindt);

  if (++m_ndtmsh == m_nelem.size()) {    // all meshes have been loaded
    Assert(m_dtmsh.size() == m_nelem.size(), "Incorrect size of dtmsh");

    // compute minimum dt across meshes
    tk::real dt = std::numeric_limits< tk::real >::max();
    for (auto idt : m_dtmsh) dt = std::min(dt, idt);

    // clear dt-vector and counter
    m_dtmsh.clear();
    m_ndtmsh = 0;

    // broadcast to advance time step
    for (auto& m : m_scheme) {
      m.bcast< Scheme::advance >( dt, F );
    }
  }
}

void
Transporter::inthead( const InciterPrint& print )
// *****************************************************************************
// Print out time integration header to screen
//! \param[in] print Pretty printer object to use for printing
// *****************************************************************************
{
  auto refined = g_inputdeck.get< tag::field_output, tag::refined >();
  const auto scheme = g_inputdeck.get< tag::scheme >();
  if (refined && scheme == ctr::SchemeType::DG) {
    printer() << "\n>>> WARNING: Ignoring refined field output for DG(P0)\n\n";
    refined = false;
  }

  print.inthead( "Time integration", "Euler/Navier-Stokes solver",
  "Legend: it - iteration count\n"
  "         t - physics time\n"
  "        dt - physics time step size\n"
  "       ETE - estimated wall-clock time elapsed (h:m:s)\n"
  "       ETA - estimated wall-clock time for accomplishment (h:m:s)\n"
  "       EGT - estimated grind wall-clock time (ms/timestep)\n"
  "       flg - status flags, legend:\n"
  "             f - " + std::string(refined ? "refined " : "")
                      + "field (volume and surface)\n"
  "             d - diagnostics\n"
  "             t - physics time history\n"
  "             h - h-refinement\n"
  "             l - load balancing\n"
  "             r - checkpoint\n"
  "             a - ALE mesh velocity linear solver did not converge\n",
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
  std::size_t meshid, ncomp;
  std::vector< std::vector< tk::real > > d;

  // Deserialize diagnostics vector
  PUP::fromMem creator( msg->getData() );
  creator | meshid;
  creator | ncomp;
  creator | d;
  delete msg;

  auto id = std::to_string(meshid);

  Assert( ncomp > 0, "Number of scalar components must be positive");
  Assert( d.size() == NUMDIAG, "Diagnostics vector size mismatch" );

  for (std::size_t i=0; i<d.size(); ++i)
     Assert( d[i].size() == ncomp, "Size mismatch at final stage of "
             "diagnostics aggregation for mesh " + id );

  // Allocate storage for those diagnostics that are always computed
  std::vector< tk::real > diag( ncomp, 0.0 );

  // Finish computing diagnostics
  for (std::size_t i=0; i<d[L2SOL].size(); ++i)
    diag[i] = sqrt( d[L2SOL][i] / m_meshvol[meshid] );
  
  // Query user-requested error types to output
  const auto& error = g_inputdeck.get< tag::diagnostics, tag::error >();

  decltype(ncomp) n = 0;
  n += ncomp;
  if (error == tk::ctr::ErrorType::L2) {
   // Finish computing the L2 norm of the numerical - analytical solution
   for (std::size_t i=0; i<d[L2ERR].size(); ++i)
     diag.push_back( sqrt( d[L2ERR][i] / m_meshvol[meshid] ) );
  } else if (error == tk::ctr::ErrorType::LINF) {
    // Finish computing the Linf norm of the numerical - analytical solution
    for (std::size_t i=0; i<d[LINFERR].size(); ++i)
      diag.push_back( d[LINFERR][i] );
  }

  // Finish computing the L2 norm of the residual and append
  const auto scheme = g_inputdeck.get< tag::scheme >();
  std::vector< tk::real > l2res( d[L2RES].size(), 0.0 );
  if (scheme == ctr::SchemeType::ALECG || scheme == ctr::SchemeType::OversetFE) {
    for (std::size_t i=0; i<d[L2RES].size(); ++i) {
      l2res[i] = std::sqrt( d[L2RES][i] / m_meshvol[meshid] );
      diag.push_back( l2res[i] );
    }
  }
  else if (scheme == ctr::SchemeType::FV) {
    for (std::size_t i=0; i<d[L2RES].size(); ++i) {
      l2res[i] = std::sqrt( d[L2RES][i] );
      diag.push_back( l2res[i] );
    }
  }

  // Append total energy
  diag.push_back( d[TOTALSOL][0] );

  // Append diagnostics file at selected times
  auto filename = g_inputdeck.get< tag::cmd, tag::io, tag::diag >();
  if (m_nelem.size() > 1) filename += '.' + id;
  tk::DiagWriter dw( filename,
    g_inputdeck.get< tag::diagnostics, tag::format >(),
    g_inputdeck.get< tag::diagnostics, tag::precision >(),
    std::ios_base::app );
  dw.diag( static_cast<uint64_t>(d[ITER][0]), d[TIME][0], d[DT][0], diag );

  // Continue time step
  m_scheme[meshid].bcast< Scheme::refine >( l2res );
}

void
Transporter::resume()
// *****************************************************************************
// Resume execution from checkpoint/restart files
//! \details This is invoked by Charm++ after the checkpoint is done, as well as
//!   when the restart (returning from a checkpoint) is complete
// *****************************************************************************
{
  if (std::any_of(begin(m_finished), end(m_finished), [](auto f){return !f;})) {
    // If just restarted from a checkpoint, Main( CkMigrateMessage* msg ) has
    // increased nrestart in g_inputdeck, but only on PE 0, so broadcast.
    auto nrestart = g_inputdeck.get< tag::cmd, tag::io, tag::nrestart >();
    for (std::size_t i=0; i<m_nelem.size(); ++i)
      m_scheme[i].bcast< Scheme::evalLB >( nrestart );
  } else
    mainProxy.finalize();
}

void
Transporter::checkpoint( std::size_t finished, std::size_t meshid )
// *****************************************************************************
// Save checkpoint/restart files
//! \param[in] finished Nonzero if finished with time stepping
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  m_finished[meshid] = finished;

  if (++m_nchk == m_nelem.size()) { // all worker arrays have checkpointed
    m_nchk = 0;
    if (not g_inputdeck.get< tag::cmd, tag::benchmark >()) {
      const auto& restart = g_inputdeck.get< tag::cmd, tag::io, tag::restart >();
      CkCallback res( CkIndex_Transporter::resume(), thisProxy );
      CkStartCheckpoint( restart.c_str(), res );
    } else {
      resume();
    }
  }
}

void
Transporter::finish( std::size_t meshid )
// *****************************************************************************
// Normal finish of time stepping
//! \param[in] meshid Mesh id
// *****************************************************************************
{
  checkpoint( /* finished = */ 1, meshid );
}

#include "NoWarning/transporter.def.h"
