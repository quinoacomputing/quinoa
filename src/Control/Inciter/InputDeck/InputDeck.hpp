// *****************************************************************************
/*!
  \file      src/Control/Inciter/InputDeck/InputDeck.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter's new input deck definition
  \details   This file defines the heterogeneous struct that is used for storing
     the data from user input during the control file parsing of the
     computational shock hydrodynamics tool, Inciter.
*/
// *****************************************************************************
#ifndef InputDeck_h
#define InputDeck_h

#include <getopt.h>
#include "Types.hpp"
#include "TaggedTuple.hpp"
#include "Inciter/CmdLine/CmdLine.hpp"
#include "Transfer.hpp"
#include "Inciter/OutVar.hpp"
#include "Inciter/Options/PDE.hpp"
#include "Inciter/Options/Problem.hpp"
#include "Inciter/Options/Scheme.hpp"
#include "Inciter/Options/Limiter.hpp"
#include "Inciter/Options/Flux.hpp"
#include "Inciter/Options/Initiate.hpp"
#include "Inciter/Options/AMRInitial.hpp"
#include "Inciter/Options/AMRError.hpp"
#include "Inciter/Options/PrefIndicator.hpp"
#include "Inciter/Options/MeshVelocity.hpp"
#include "Inciter/Options/MeshVelocitySmoother.hpp"
#include "Inciter/Options/Material.hpp"
#include "Options/PartitioningAlgorithm.hpp"
#include "Options/TxtFloatFormat.hpp"
#include "Options/FieldFile.hpp"
#include "Options/Error.hpp"
#include "Options/UserTable.hpp"

namespace inciter {

namespace ctr {

using ncomp_t = std::size_t;

using bclist = tk::TaggedTuple< brigand::list<
  tag::dirichlet,   std::vector< std::size_t >,
  tag::symmetry,    std::vector< std::size_t >,
  tag::outlet,      std::vector< std::size_t >,
  tag::farfield,    std::vector< std::size_t >,
  tag::extrapolate, std::vector< std::size_t >,
  tag::noslipwall,  std::vector< std::size_t >
> >;

// Transport
using transportList = tk::TaggedTuple< brigand::list<
    tag::physics,        PhysicsType,
    tag::ncomp,          std::size_t,
    tag::intsharp,       int,
    tag::intsharp_param, tk::real,
    tag::problem,        ProblemType,
    tag::diffusivity,    std::vector< tk::real >,
    tag::lambda,         std::vector< tk::real >,
    tag::u0,             std::vector< tk::real >
> >;

// CompFlow
using compflowList = tk::TaggedTuple< brigand::list<
    tag::physics, PhysicsType,
    tag::problem, ProblemType,
    tag::alpha,   tk::real,
    tag::beta,    tk::real,
    tag::betax,   tk::real,
    tag::betay,   tk::real,
    tag::betaz,   tk::real,
    tag::r0,      tk::real,
    tag::p0,      tk::real,
    tag::ce,      tk::real,
    tag::kappa,   tk::real
> >;

// MultiMat
using multimatList = tk::TaggedTuple< brigand::list<
  tag::physics,          PhysicsType,
  tag::nmat,             std::size_t,
  tag::prelax,           uint64_t,
  tag::prelax_timescale, tk::real,
  tag::intsharp,         int,
  tag::intsharp_param,   tk::real,
  tag::rho0constraint,   uint64_t,
  tag::dt_sos_massavg,   int,
  tag::problem,          ProblemType,
  tag::viscous,          bool
> >;

// MultiSpecies
using multispeciesList = tk::TaggedTuple< brigand::list<
  tag::physics,          PhysicsType,
  tag::nspec,            std::size_t,
  tag::problem,          ProblemType,
  tag::viscous,          bool
> >;

// Material/EOS object
using materialList = tk::TaggedTuple< brigand::list<
  tag::eos,          MaterialType,
  tag::id,           std::vector< uint64_t >,
  tag::gamma,        std::vector< tk::real >,
  tag::pstiff,       std::vector< tk::real >,
  tag::w_gru,        std::vector< tk::real >,
  tag::A_jwl,        std::vector< tk::real >,
  tag::B_jwl,        std::vector< tk::real >,
  tag::C_jwl,        std::vector< tk::real >,
  tag::R1_jwl,       std::vector< tk::real >,
  tag::R2_jwl,       std::vector< tk::real >,
  tag::rho0_jwl,     std::vector< tk::real >,
  tag::de_jwl,       std::vector< tk::real >,
  tag::rhor_jwl,     std::vector< tk::real >,
  tag::Tr_jwl,       std::vector< tk::real >,
  tag::Pr_jwl,       std::vector< tk::real >,
  tag::mu,           std::vector< tk::real >,
  tag::yield_stress, std::vector< tk::real >,
  tag::cv,           std::vector< tk::real >,
  tag::k,            std::vector< tk::real >
> >;

// Boundary conditions block
using bcList = tk::TaggedTuple< brigand::list<
  tag::mesh,        std::vector< std::size_t >,
  tag::dirichlet,   std::vector< std::size_t >,
  tag::symmetry,    std::vector< std::size_t >,
  tag::outlet,      std::vector< std::size_t >,
  tag::farfield,    std::vector< std::size_t >,
  tag::extrapolate, std::vector< std::size_t >,
  tag::noslipwall,  std::vector< std::size_t >,
  tag::stag_point,  std::vector< tk::real >,
  tag::radius,      tk::real,
  tag::velocity,    std::vector< tk::real >,
  tag::pressure,    tk::real,
  tag::density,     tk::real,
  tag::temperature, tk::real,
  tag::mass_fractions, std::vector< tk::real >,
  tag::materialid,  std::size_t,
  tag::inlet,       std::vector<
    tk::TaggedTuple< brigand::list<
      tag::sideset,      std::vector< uint64_t >,
      tag::velocity,     std::vector< tk::real >,
      tag::materialid,   std::size_t
    > >
  >,
  tag::timedep,     std::vector<
    tk::TaggedTuple< brigand::list<
      tag::sideset,    std::vector< uint64_t >,
      tag::fn,         std::vector< tk::real >
    > >
  >
> >;

// IC box
using boxList = tk::TaggedTuple< brigand::list<
  tag::materialid,     std::size_t,
  tag::volume,         tk::real,
  tag::mass,           tk::real,
  tag::density,        tk::real,
  tag::velocity,       std::vector< tk::real >,
  tag::pressure,       tk::real,
  tag::energy,         tk::real,
  tag::energy_content, tk::real,
  tag::temperature,    tk::real,
  tag::mass_fractions, std::vector< tk::real >,
  tag::xmin,           tk::real,
  tag::xmax,           tk::real,
  tag::ymin,           tk::real,
  tag::ymax,           tk::real,
  tag::zmin,           tk::real,
  tag::zmax,           tk::real,
  tag::orientation,    std::vector< tk::real >,
  tag::initiate,       inciter::ctr::InitiateType,
  tag::point,          std::vector< tk::real >,
  tag::init_time,      tk::real,
  tag::front_width,    tk::real,
  tag::front_speed,    tk::real
> >;

// IC meshblock
using meshblockList = tk::TaggedTuple< brigand::list<
  tag::blockid,        std::uint64_t,
  tag::materialid,     std::size_t,
  tag::volume,         tk::real,
  tag::mass,           tk::real,
  tag::density,        tk::real,
  tag::velocity,       std::vector< tk::real >,
  tag::pressure,       tk::real,
  tag::energy,         tk::real,
  tag::energy_content, tk::real,
  tag::temperature,    tk::real,
  tag::mass_fractions, std::vector< tk::real >,
  tag::initiate,       inciter::ctr::InitiateType,
  tag::point,          std::vector< tk::real >,
  tag::init_time,      tk::real,
  tag::front_width,    tk::real,
  tag::front_speed,    tk::real
> >;

// Initial conditions (ic) block
using icList = tk::TaggedTuple< brigand::list<
  tag::materialid,     std::size_t,
  tag::pressure,       tk::real,
  tag::temperature,    tk::real,
  tag::mass_fractions, std::vector< tk::real >,
  tag::density,        tk::real,
  tag::energy,         tk::real,
  tag::velocity,       std::vector< tk::real >,
  tag::box,            std::vector< boxList >,
  tag::meshblock,      std::vector< meshblockList >
> >;

// Overset mesh block
using meshList = tk::TaggedTuple< brigand::list<
  tag::filename,    std::string,
  tag::location,    std::vector< tk::real >,
  tag::orientation, std::vector< tk::real >,
  tag::velocity,    std::vector< tk::real >
> >;

// Field output block
using fieldOutputList = tk::TaggedTuple< brigand::list<
  tag::interval, uint32_t,
  tag::time_interval, tk::real,
  tag::time_range,    std::vector< tk::real >,
  tag::refined,       bool,
  tag::filetype,      tk::ctr::FieldFileType,
  tag::sideset,       std::vector< uint64_t >,
  tag::outvar,        std::vector< OutVar >,
  tag::elemalias,     std::vector< std::string >,  // only for error checking
  tag::elemvar,       std::vector< std::string >,  // only for error checking
  tag::nodealias,     std::vector< std::string >,  // only for error checking
  tag::nodevar,       std::vector< std::string >   // only for error checking
> >;

// Diagnostics block
using diagnosticsList = tk::TaggedTuple< brigand::list<
  tag::interval, uint32_t,
  tag::error,         tk::ctr::ErrorType,
  tag::format,        tk::ctr::TxtFloatFormatType,
  tag::precision,     std::streamsize
> >;

// History output block
using historyOutputList = tk::TaggedTuple< brigand::list<
  tag::interval, uint32_t,
  tag::time_interval, tk::real,
  tag::time_range,    std::vector< tk::real >,
  tag::format,        tk::ctr::TxtFloatFormatType,
  tag::precision,     std::streamsize,
  tag::point,         std::vector<
    tk::TaggedTuple< brigand::list<
      tag::id,    std::string,
      tag::coord, std::vector< tk::real >
    > >
  >
> >;

using ConfigMembers = brigand::list<

  tag::title, std::string,

  // Command line parameters
  tag::cmd, CmdLine,

  // time stepping options
  tag::nstep,            uint64_t,
  tag::term,             tk::real,
  tag::t0,               tk::real,
  tag::dt,               tk::real,
  tag::cfl,              tk::real,
  tag::ttyi,             uint32_t,
  tag::imex_runge_kutta, uint32_t,
  tag::imex_maxiter,     uint32_t,
  tag::imex_reltol,      tk::real,
  tag::imex_abstol,      tk::real,

  // steady-state solver options
  tag::steady_state, bool,
  tag::residual,     tk::real,
  tag::rescomp,      uint32_t,

  // mesh partitioning and reordering/sorting choices
  tag::partitioning,     tk::ctr::PartitioningAlgorithmType,
  tag::pelocal_reorder,  bool,
  tag::operator_reorder, bool,

  // discretization scheme choices
  tag::scheme,      SchemeType,
  tag::ndof,        std::size_t,
  tag::rdof,        std::size_t,
  tag::flux,        FluxType,
  tag::lowspeed_kp, tk::real,

  // limiter options
  tag::limiter,              LimiterType,
  tag::cweight,              tk::real,
  tag::shock_detector_coeff, tk::real,
  tag::accuracy_test,        bool,
  tag::limsol_projection,    bool,

  // PDE options
  tag::ncomp,        std::size_t,
  tag::pde,          PDEType,
  tag::transport,    transportList,
  tag::compflow,     compflowList,
  tag::multimat,     multimatList,
  tag::multispecies, multispeciesList,

  // Dependent variable name
  tag::depvar, std::vector< char >,

  tag::sys, std::map< std::size_t, std::size_t >,

  tag::material, std::vector< materialList >,

  tag::matidxmap, tk::TaggedTuple< brigand::list<
    tag::eosidx, std::vector< std::size_t >,
    tag::matidx, std::vector< std::size_t >,
    tag::solidx, std::vector< std::size_t >
  > >,

  // Conditions
  tag::bc, std::vector< bcList >,
  tag::ic, icList,
  tag::mesh, std::vector< meshList >,
  tag::transfer, std::vector< Transfer >,

  // ALE block
  // ---------------------------------------------------------------------------
  tag::ale, tk::TaggedTuple< brigand::list<
    tag::ale,           bool,
    tag::smoother,      MeshVelocitySmootherType,
    tag::mesh_velocity, MeshVelocityType,
    tag::mesh_motion,   std::vector< std::size_t >,
    tag::meshforce,     std::vector< tk::real >,
    tag::dvcfl,         tk::real,
    tag::vortmult,      tk::real,
    tag::maxit,         std::size_t,
    tag::tolerance,     tk::real,
    tag::dirichlet,     std::vector< std::size_t >,
    tag::symmetry,      std::vector< std::size_t >,
    tag::move,          std::vector<
      tk::TaggedTuple< brigand::list<
        tag::sideset, std::vector< uint64_t >,
        tag::fntype,  tk::ctr::UserTableType,
        tag::fn,      std::vector< tk::real >
      > >
    >
  > >,

  // p-refinement block
  // ---------------------------------------------------------------------------
  tag::pref, tk::TaggedTuple< brigand::list<
    tag::pref,      bool,
    tag::indicator, PrefIndicatorType,
    tag::ndofmax,   std::size_t,
    tag::tolref,    tk::real
  > >,

  // AMR block
  // ---------------------------------------------------------------------------
  tag::amr, tk::TaggedTuple< brigand::list<
    tag::amr,           bool,
    tag::t0ref,         bool,
    tag::dtref,         bool,
    tag::dtref_uniform, bool,
    tag::dtfreq,        std::size_t,
    tag::maxlevels,     std::size_t,
    tag::initial,       std::vector< AMRInitialType >,
    tag::edgelist,      std::vector< std::size_t >,
    tag::coords,        tk::TaggedTuple< brigand::list<
      tag::xminus, tk::real,
      tag::xplus,  tk::real,
      tag::yminus, tk::real,
      tag::yplus,  tk::real,
      tag::zminus, tk::real,
      tag::zplus,  tk::real
    > >,
    tag::error,         AMRErrorType,
    tag::refvar,        std::vector< char >,
    tag::tol_refine,    tk::real,
    tag::tol_derefine,  tk::real
  > >,

  // Output options
  tag::field_output,   fieldOutputList,
  tag::diagnostics,    diagnosticsList,
  tag::history_output, historyOutputList
>;

// Class storing the Config params
class InputDeck : public tk::TaggedTuple< ConfigMembers > {

  public:
    //! Set of tags to ignore when printing this InputDeck
    using ignore = CmdLine::ignore;

    //! \brief Constructor: set defaults
    //! \param[in] cl Previously parsed and store command line
    //! \details Anything not set here is initialized by the compiler using the
    //!   default constructor for the corresponding type.
    explicit InputDeck( const CmdLine& cl = {} ) {
      // Set previously parsed command line
      get< tag::cmd >() = cl;
      // Default time stepping params
      get< tag::dt >() = 0.0;
      get< tag::cfl >() = 0.0;
      // Default AMR settings
      auto rmax =
        std::numeric_limits< tk::real >::max() / 100;
      get< tag::amr, tag::coords, tag::xminus >() = rmax;
      get< tag::amr, tag::coords, tag::xplus >() = -rmax;
      get< tag::amr, tag::coords, tag::yminus >() = rmax;
      get< tag::amr, tag::coords, tag::yplus >() = -rmax;
      get< tag::amr, tag::coords, tag::zminus >() = rmax;
      get< tag::amr, tag::coords, tag::zplus >() = -rmax;

      // -----------------------------------------------------------------------
      /*
                                  Keyword vector
         The following code generates a vector of keywords, for the sole purpose
         of documentation and user-help (accessible via the --helpctr and
         --helpkw cmdline arguments. If entries for a keyword are not added to
         this vector, help will not be output for it. The entries follow the
         function signature of the tk::entry_t constructor (defined in
         Base/Types.hpp).
      */
      // -----------------------------------------------------------------------
      std::set< tk::entry_t > keywords;

      keywords.insert({"inciter",
        "Start configuration block for inciter",
        R"(This keyword is used to select inciter. Inciter, is a continuum-realm
        shock hydrodynamics tool, solving a system of PDEs. The entire control
        file must be enclosed within the inciter block)", "block-title"});

      keywords.insert({"title", "Title", R"(The title may be specified in
        the input file. It is optional.)", "string"});

      // -----------------------------------------------------------------------
      // time stepping options
      // -----------------------------------------------------------------------

      keywords.insert({"nstep", "Set number of time steps to take",
        R"(This keyword is used to specify the number of time steps to take in a
        simulation. The number of time steps are used in conjunction with the
        maximmum time specified by keyword 'term': the simulation stops whichever
        is reached first. Both 'nstep' and 'term' can be left unspecified, in
        which case their default values are used. See also 'term'.)", "uint"});

      keywords.insert({"term", "Set maximum physical time to simulate",
        R"(This keyword is used to specify the termination time in a simulation.
        The termination time and number of time steps, specified by 'nstep', are
        used in conjunction to determine when to stop a simulation: whichever is
        reached first. Both 'nstep' and 'term' can be left unspecified, in which
        case their default values are used. See also 'nstep'.)", "real"});

      keywords.insert({"t0", "Set starting non-dimensional time",
        R"(This keyword is used to specify the starting time in a simulation.)",
        "real"});

      keywords.insert({"dt", "Select constant time step size",
        R"(This keyword is used to specify the time step size that used as a
        constant during simulation. Setting 'cfl' and 'dt' are mutually
        exclusive. If both 'cfl' and 'dt' are set, 'dt' wins.)", "real"});

      keywords.insert({"cfl",
      "Set the Courant-Friedrichs-Lewy (CFL) coefficient",
      R"(This keyword is used to specify the CFL coefficient for
      variable-time-step-size simulations. Setting 'cfl' and 'dt' are mutually
      exclusive. If both 'cfl' and 'dt' are set, 'dt' wins.)", "real"});

      keywords.insert({"ttyi", "Set screen output interval",
        R"(This keyword is used to specify the interval in time steps for screen
        output during a simulation.)", "uint"});

      keywords.insert({"imex_runge_kutta",
        "Toggle use of IMplicit-EXplicit Runge-Kutta scheme",
        R"(This keywords is used to turn IMEX integrator on/off for solid materials
        in a multimat run. Plastic terms are integrated implicitly in time. This
        flag will activate an Implicit-Explicit Runge-Kutta scheme to replace the
        explicit one that is usually used. Scheme taken from Cavaglieri, D., &
        Bewley, T. (2015). Low-storage implicit/explicit Runge–Kutta schemes for
        the simulation of stiff high-dimensional ODE systems. Journal of
        Computational Physics, 286, 172-193.)", "uint 0/1"});

      keywords.insert({"imex_maxiter",
        "Set maximum number of iterations for non-linear solver with IMEX-RK scheme",
        R"(This keywords is used to specify the maximum number of iterations that
        the non-linear solver uses to obtain the implicit unknowns within the
        Implicit-Explicit Runge-Kutta scheme.)", "uint"});

      keywords.insert({"imex_reltol",
        "Set relative tolerance for non-linear solver with IMEX-RK scheme",
        R"(This keywords is used to specify the relative tolerance that
        the non-linear solver uses to obtain the implicit unknowns within the
        Implicit-Explicit Runge-Kutta scheme.)", "real"});

      keywords.insert({"imex_abstol",
        "Set absolute tolerance for non-linear solver with IMEX-RK scheme",
        R"(This keywords is used to specify the absolute tolerance that
        the non-linear solver uses to obtain the implicit unknowns within the
        Implicit-Explicit Runge-Kutta scheme.)", "real"});

      // -----------------------------------------------------------------------
      // steady-state solver options
      // -----------------------------------------------------------------------

      keywords.insert({"steady_state", "March to steady state",
        R"(This keyword is used indicate that local time stepping should be used
        to march towards a stationary solution.)", "bool"});

      keywords.insert({"residual",
        "Set the convergence criterion for the residual to reach",
        R"(This keyword is used to specify a convergence criterion for the local
        time stepping marching to steady state, below which the simulation is
        considered converged.)", "real"});

      keywords.insert({"rescomp",
        "Equation system component index for convergence",
        R"(This keyword is used to specify a single integer that is used to denote
        the equation component index in the complete system of equations
        configured, to use for the convergence criterion for local
        time stepping marching towards steady state.)", "uint"});

      // -----------------------------------------------------------------------
      // mesh partitioning and reordering/sorting choices
      // -----------------------------------------------------------------------

      keywords.insert({"partitioning",
        "Select mesh partitioning algorithm",
        R"(This keyword is used to select a mesh partitioning algorithm. See
        Control/Options/PartitioningAlgorithm.hpp for valid options.)",
        "string"});

      keywords.insert({"rcb",
        "Select recursive coordinate bisection mesh partitioner",
        R"(This keyword is used to select the recursive coordinate bisection (RCB)
        mesh partitioner. RCB is a geometry-based partitioner used to distribute
        an input mesh among processing elements. See
        Control/Options/PartitioningAlgorithm.hpp for other valid options.)"});

      keywords.insert({"rib",
        "Select recursive inertial bisection mesh partitioner",
        R"(This keyword is used to select the recursive inertial bisection (RIB)
        mesh partitioner. RIB is a geometry-based partitioner used to distribute
        an input mesh among processing elements. See
        Control/Options/PartitioningAlgorithm.hpp for other valid options.)"});

      keywords.insert({"hsfc",
        "Select Hilbert Space Filling Curve (HSFC) mesh partitioner",
        R"(This keyword is used to select the Hilbert Space Filling Curve (HSFC)
        mesh partitioner. HSFC is a geometry-based partitioner used to distribute
        an input mesh among processing elements. See
        Control/Options/PartitioningAlgorithm.hpp for other valid options.)"});

      keywords.insert({"phg",
        "Select parallel hypergraph mesh partitioner",
        R"(This keyword is used to select the parallel hypergraph (PHG)
        mesh partitioner. PHG is a graph-based partitioner used to distribute an
        input mesh among processing elements. See
        Control/Options/PartitioningAlgorithm.hpp for other valid options.)"});

      keywords.insert({"mj",
        "Select multi-jagged (MJ) mesh partitioner",
        R"(This keyword is used to select the multi-jagged (MJ) mesh partitioner.
        MJ is a geometry-based partitioner used to distribute an input mesh among
        processing elements. See
        Control/Options/PartitioningAlgorithm.hpp for other valid options.)"});

      keywords.insert({"pelocal_reorder",
        "PE-local reorder",
        R"(This keyword is used in inciter as a keyword in the inciter...end block
        as "pelocal_reorder true" (or false) to do (or not do) a global
        distributed mesh reordering across all PEs that yields an approximately
        continuous mesh node ID order as mesh partitions are assigned to PEs after
        mesh partitioning. This reordering is optional.)", "bool"});

      keywords.insert({"operator_reorder",
        "Operator-access reorder",
        R"(This keyword is used in inciter as a keyword in the inciter...end block
        as "operator_reorder on" (or off) to do (or not do) a local mesh node
        reordering based on the PDE operator access pattern. This reordering is
        optional.)", "bool"});

      // -----------------------------------------------------------------------
      // discretization scheme choices
      // -----------------------------------------------------------------------

      keywords.insert({"scheme", "Select discretization scheme",
        R"(This keyword is used to select a spatial discretization scheme,
        necessarily connected to the temporal discretization scheme. See
        Control/Inciter/Options/Scheme.hpp for valid options.)", "string"});

      keywords.insert({"alecg",
        "Select continuous Galerkin with ALE + Runge-Kutta",
        R"(This keyword is used to select the continuous Galerkin finite element
        scheme in the arbitrary Lagrangian-Eulerian (ALE) reference frame combined
        with Runge-Kutta (RK) time stepping.
        See Control/Inciter/Options/Scheme.hpp for other valid options.)"});

      keywords.insert({"oversetfe",
        "Select continuous Galerkin finite element with overset meshes + "
        "Runge-Kutta",
        R"(This keyword is used to select the continuous Galerkin finite element
        scheme with Runge-Kutta (RK) time stepping, combined with overset grids.
        See Control/Inciter/Options/Scheme.hpp for other valid options.)"});

      keywords.insert({"dg",
        "Select 1st-order discontinuous Galerkin discretization + Runge-Kutta",
        R"(This keyword is used to select the first-order accurate discontinuous
        Galerkin, DG(P0), spatial discretiztaion used in Inciter. As this is first
        order accurate, it is intended for testing and debugging purposes only.
        Selecting this spatial discretization also selects the Runge-Kutta scheme
        for time discretization. See Control/Inciter/Options/Scheme.hpp for other
        valid options.)"});

      keywords.insert({"p0p1",
        "Select 2nd-order finite volume discretization + Runge-Kutta",
        R"(This keyword is used to select the second-order accurate finite volume,
        P0P1, spatial discretiztaion used in Inciter. This method uses a
        least-squares procedure to reconstruct the second-order solution from the
        first-order one. Selecting this spatial discretization also selects the
        Runge-Kutta scheme for time discretization.
        See Control/Inciter/Options/Scheme.hpp for other valid options.)"});

      keywords.insert({"dgp1",
        "Select 2nd-order discontinuous Galerkin discretization + Runge-Kutta",
        R"(This keyword is used to select the second-order accurate discontinuous
        Galerkin, DG(P1), spatial discretiztaion used in Inciter. Selecting this
        spatial discretization also selects the Runge-Kutta scheme for time
        discretization. See Control/Inciter/Options/Scheme.hpp for other
        valid options.)"});

      keywords.insert({"dgp2",
        "Select 3nd-order discontinuous Galerkin discretization + Runge-Kutta",
        R"(This keyword is used to select the third-order accurate discontinuous
        Galerkin, DG(P2), spatial discretiztaion used in Inciter. Selecting this
        spatial discretization also selects the Runge-Kutta scheme for time
        discretization. See Control/Inciter/Options/Scheme.hpp for other
        valid options.)"});

      keywords.insert({"pdg",
        "Select p-adaptive discontinuous Galerkin discretization + Runge-Kutta",
        R"(This keyword is used to select the polynomial adaptive discontinuous
        Galerkin spatial discretizaion used in Inciter. Selecting this spatial
        discretization also selects the Runge-Kutta scheme for time
        discretization. See Control/Inciter/Options/Scheme.hpp for other valid
        options.)"});

      keywords.insert({"fv",
        "Select 2nd-order finite volume discretization + Runge-Kutta",
        R"(This keyword is used to select the second-order accurate finite volume,
        P0P1, spatial discretiztaion used in Inciter. This method uses a
        least-squares procedure to reconstruct the second-order solution from the
        first-order one. See Control/Inciter/Options/Scheme.hpp for other valid
        options.)"});

      keywords.insert({"ndof", "Number of evolved solution DOFs",
        R"(The number of solution DOFs that are evolved.)", "uint"});

      keywords.insert({"rdof", "Total number of solution DOFs",
        R"(The total number of solution DOFs, including the reconstructed and the
        evolved ones.)", "uint"});

      // -----------------------------------------------------------------------
      // limiter options
      // -----------------------------------------------------------------------

      keywords.insert({"limiter", "Select limiter function",
        R"(This keyword is used to select a limiter function, used for
        discontinuous Galerkin (DG) spatial discretization used in inciter. See
        Control/Inciter/Options/Limiter.hpp for valid options.)", "string"});

      keywords.insert({"nolimiter", "No limiter used",
        R"(This keyword is used for discontinuous Galerkin (DG) spatial
        discretization without any limiter in inciter. See
        Control/Inciter/Options/Limiter.hpp for other valid options.)"});

      keywords.insert({"wenop1",
        "Select the Weighted Essentially Non-Oscillatory (WENO) limiter for DGP1",
        R"(This keyword is used to select the Weighted Essentially Non-Oscillatory
        limiter used for discontinuous Galerkin (DG) P1 spatial discretization
        used in inciter. See Control/Inciter/Options/Limiter.hpp for other valid
        options.)"});

      keywords.insert({"cweight",
        "Set value for central linear weight used by WENO, cweight",
        R"(This keyword is used to set the central linear weight used for the
        central stencil in the Weighted Essentially Non-Oscillatory (WENO) limiter
        for discontinuous Galerkin (DG) methods.)", "real"});

      keywords.insert({"superbeep1",
        "Select the Superbee limiter for DGP1",
        R"(This keyword is used to select the Superbee limiter used for
        discontinuous Galerkin (DG) P1 spatial discretization used in inciter.
        See Control/Inciter/Options/Limiter.hpp for other valid options.)"});

      keywords.insert({"shock_detector_coeff",
        "Configure the coefficient used in shock indicator",
        R"(This keyword can be used to configure the coefficient used in the
        threshold calculation for the shock indicator.)", "real"});

      keywords.insert({"vertexbasedp1",
        "Select the vertex-based limiter for DGP1",
        R"(This keyword is used to select the vertex-based limiter used for
        discontinuous Galerkin (DG) P1 spatial discretization used in inciter.
        Ref. Kuzmin, D. (2010). A vertex-based hierarchical slope limiter for
        p-adaptive discontinuous Galerkin methods. Journal of computational and
        applied mathematics, 233(12), 3077-3085.
        See Control/Inciter/Options/Limiter.hpp for other valid options.)"});

      keywords.insert({"accuracy_test", "Toggle accuracy test setup",
        R"(This keyword is used to specify if the current setup is for an
        order-of-accuracy testing, used for discontinuous Galerkin (DG) spatial
        discretization in inciter. This deactivates certain robustness corrections
        which might impact order-of-accuracy. Only intended for simple test
        problems and not for real problems.)", "bool"});

      keywords.insert({"limsol_projection",
        "Toggle limited solution projection",
        R"(This keyword is used to specify limited solution projection.
        This is used for discontinuous Galerkin (DG) spatial discretization in
        inciter, for multi-material hydrodynamics. This uses a projection to
        obtain bulk momentum and material energies from the limited primitive
        quantities. This step is essential to obtain closure-law obeying limited
        quantities. See Pandare et al. (2023). On the Design of Stable,
        Consistent, and Conservative High-Order Methods for Multi-Material
        Hydrodynamics. J Comp Phys (490).)", "bool"});

      // -----------------------------------------------------------------------
      // flux options
      // -----------------------------------------------------------------------

      keywords.insert({"flux", "Select flux function",
        R"(This keyword is used to select a flux function, used for
        discontinuous Galerkin (DG) spatial discretization used in inciter. See
        Control/Inciter/Options/Flux.hpp for valid options.)", "string"});

      keywords.insert({"laxfriedrichs",
        "Select Lax-Friedrichs flux function",
        R"(This keyword is used to select the Lax-Friedrichs flux function used
        for discontinuous Galerkin (DG) spatial discretization used in inciter.
        See Control/Inciter/Options/Flux.hpp for other valid options.)"});

      keywords.insert({"hllc",
        "Select the Harten-Lax-van Leer-Contact (HLLC) flux function",
        R"(This keyword is used to select the Harten-Lax-van Leer-Contact flux
        function used for discontinuous Galerkin (DG) spatial discretization
        used in inciter. See Control/Inciter/Options/Flux.hpp for other valid
        options.)"});

      keywords.insert({"upwind", "Select the upwind flux function",
        R"(This keyword is used to select the upwind flux
        function used for discontinuous Galerkin (DG) spatial discretization
        used in inciter. It is only usable for scalar transport.
        See Control/Inciter/Options/Flux.hpp for other valid options.)"});

      keywords.insert({"ausm",
        "Select the Advection Upstream Splitting Method (AUSM) flux function",
        R"(This keyword is used to select the AUSM flux
        function used for discontinuous Galerkin (DG) spatial discretization
        used in inciter. It is only set up for for multi-material hydro, and
        not selectable for anything else.)"});

      keywords.insert({"ldfss",
        "Select the Low Diffusion Flux Splitting Scheme (LDFSS)",
        R"(This keyword is used to select the LDFSS flux
        function used for discontinuous Galerkin (DG) spatial discretization
        used in inciter. It is only set up for for multi-material hydro, and
        not selectable for anything else.)"});

      keywords.insert({"lowspeed_kp",
        "Select the low-speed coefficient K_p in the AUSM+up flux function",
        R"(This keyword is used to select the low-speed coefficient K_p in the
        AUSM+up flux function used for the DG or FV spatial discretization for
        multi-material hydro, and not used for anything else. The default
        value is 0, and recommended value for low speed flows (Mach < 0.1) is
        1.)"});

      keywords.insert({"hll",
        "Select the Harten-Lax-vanLeer (HLL) flux function",
        R"(This keyword is used to select the HLL flux
        function used for discontinuous Galerkin (DG) spatial discretization
        used in inciter. It is only set up for for multi-material hydro, and
        not selectable for anything else.)"});

      // -----------------------------------------------------------------------
      // PDE keywords
      // -----------------------------------------------------------------------

      keywords.insert({"problem",
        "Specify problem configuration for partial differential equation solver",
        R"(This keyword is used to specify the problem configuration for the
        partial differential equation solver.)", "string"});

      keywords.insert({"transport",
        "Start configuration block for an transport equation",
        R"(This keyword is used to introduce a transport block, used to
        specify the configuration for a transport equation type.)",
        "block-title"});

      keywords.insert({"ncomp",
        "Set number of scalar components for a system of transport equations",
        R"(This keyword is used to specify the number of scalar
        components of transport (linear advection) equations.)", "uint"});

      keywords.insert({"compflow",
        "Start configuration block for the compressible flow equations",
        R"(This keyword is used to introduce the compflow block, used to
        specify the configuration for a system of partial differential equations,
        governing single material compressible fluid flow.)", "block-title"});

      keywords.insert({"multimat",
        "Start configuration block for the compressible multi-material equations",
        R"(This keyword is used to introduce the multimat block,
        used to specify the configuration for a system of partial differential
        equations, governing compressible multi-material hydrodynamics assuming
        velocity equilibrium (single velocity).)", "block-title"});

      keywords.insert({"multispecies",
        "Start configuration block for the compressible multi-species equations",
        R"(This keyword is used to introduce the multispecies block,
        used to specify the configuration for a system of partial differential
        equations, governing compressible multi-species fluid dynamics.)",
        "block-title"});

      keywords.insert({"nmat",
        "Set number of materials for the multi-material system",
        R"(This keyword is used to specify the number of materials for
        multi-material flow, see also the keyword 'multimat'.)", "uint"});

      keywords.insert({"nspec",
        "Set number of species for the multi-species system",
        R"(This keyword is used to specify the number of species for
        multi-species flow, see also the keyword 'multispecies'.)", "uint"});

      keywords.insert({"prelax",
        "Toggle multi-material finite pressure relaxation",
        R"(This keyword is used to turn finite pressure relaxation between
        multiple materials on/off. It is used only for the multi-material solver,
        and has no effect when used for the other PDE types.)", "uint 0/1"});

      keywords.insert({"prelax_timescale",
        "Time-scale for multi-material finite pressure relaxation",
        R"(This keyword is used to specify the time-scale at which finite pressure
        relaxation between multiple materials occurs. The default value of 0.25
        corresponds to a relaxation time that is 4 times the time required for a
        sound wave to pass through a computational element. It is used only for
        multimat, and has no effect for the other PDE types.)", "real"});

      keywords.insert({"intsharp",
        "Toggle multi-material interface sharpening",
        R"(This keyword is used to turn interface sharpening on/off. It uses the
        multi-material THINC interface reconstruction.
        Ref. Pandare A. K., Waltz J., & Bakosi J. (2021) Multi-Material
        Hydrodynamics with Algebraic Sharp Interface Capturing. Computers &
        Fluids, doi: https://doi.org/10.1016/j.compfluid.2020.104804. It is used
        for the multi-material and the transport solver, and has no effect when
        used for the other PDE types.)", "uint 0/1"});

      keywords.insert({"intsharp_param",
        "Parameter for multi-material interface sharpening",
        R"(This keyword is used to specify the parameter for the interface
        sharpening. This parameter affects how many cells the material interfaces
        span, after the use of sharpening. It is used for multimat and transport,
        and has no effect for the other PDE types.)", "real" });

      keywords.insert({"rho0constraint",
        "Toggle the density constraint correction",
        R"(This keyword is used to toggle the density constraint in solid
        dynamics on/off. It is used only for the multi-material solver in the
        presence of solids. The default is 1 (on).)", "uint 0/1"});

      keywords.insert({"dt_sos_massavg",
        "Toggle method for calculating speed of sound used for time step in a cell",
        R"(This keyword is used to specify the method to calculate the speed of
        sound in a cell used for the time step. If set to 1, the speed of sound
        will be calculated using the mass average, rather than the maximum value
        across materials. It is used for multimat, and has no effect for the
        other PDE types.)", "uint 0/1" });

      // Dependent variable name
      keywords.insert({"depvar",
        "Select dependent variable name for PDE.",
        R"(Select dependent variable name for PDE.)", "string"});

      // -----------------------------------------------------------------------
      // physics choices
      // -----------------------------------------------------------------------

      keywords.insert({"physics",
        "Specify the physics configuration for a system of PDEs",
        R"(This keyword is used to select the physics configuration for a
        particular PDE system. Valid options depend on the system of PDEs in
        which the keyword is used.)", "string"});

      keywords.insert({"advection",
        "Specify the advection physics",
        R"(This keyword is used to select the advection physics for the transport
        PDE system. Only usable for 'transport'.)"});

      keywords.insert({"advdiff",
        "Specify the advection + diffusion physics",
        R"(This keyword is used to select the advection + diffusion physics
        for transport PDEs. Only usable for 'transport'.)"});

      keywords.insert({"euler",
        "Specify the Euler (inviscid) compressible flow physics",
        R"(This keyword is used to select the Euler (inviscid) compressible
        flow physics configuration. Usable for 'compflow' and 'multimat')"});

      keywords.insert({"energy_pill",
        "Specify the energy pill physics",
        R"(This keyword is used to select an energy pill initialization as physics
        configuration for multiple material compressible flow. Parameters for the
        linearly traveling front are required to be specified when energy_pill is
        selected. See 'linear' for more details. Currently setup only for
        'multimat')"});

      // -----------------------------------------------------------------------
      // material/eos object
      // -----------------------------------------------------------------------

      keywords.insert({"material",
        "Start configuration block for material (eos) properties",
        R"(This keyword is used to introduce a material block, used to
        specify material properties.)", "vector block-title"});

      keywords.insert({"id", "ID",
        R"(This keyword is used to specify an ID, a positive integer. Usage is
        context specific, i.e. what block it is specified in. E.g. Inside the
        material block, it is used to specify a block consisting of IDs of
        materials of that EOS type)", "vector of uints"});

      keywords.insert({"eos", "Select equation of state (type)",
        R"(This keyword is used to select an equation of state for a material.)",
        "string"});

      keywords.insert({"gamma", "ratio of specific heats",
        R"(This keyword is used to specify the material property, ratio of
        specific heats.)", "vector of reals"});

      keywords.insert({"pstiff", "EoS stiffness parameter",
        R"(This keyword is used to specify the material property, stiffness
        parameter in the stiffened gas equation of state.)", "vector of reals"});

      keywords.insert({"w_gru", "Grueneisen coefficient",
        R"(This keyword is used to specify the material property, Gruneisen
        coefficient for the Jones-Wilkins-Lee equation of state.)",
        "vector of reals"});

      keywords.insert({"A_jwl", "JWL EoS A parameter",
        R"(This keyword is used to specify the material property A (units: Pa)
        for the Jones-Wilkins-Lee equation of state.)", "vector of reals"});

      keywords.insert({"B_jwl", "JWL EoS B parameter",
        R"(This keyword is used to specify the material property B (units: Pa)
        for the Jones-Wilkins-Lee equation of state.)", "vector of reals"});

      keywords.insert({"C_jwl", "JWL EoS C parameter",
        R"(This keyword is used to specify the material property C (units: Pa)
        for the Jones-Wilkins-Lee equation of state.)", "vector of reals"});

      keywords.insert({"R1_jwl", "JWL EoS R1 parameter",
        R"(This keyword is used to specify the material property R1 for the
        Jones-Wilkins-Lee equation of state.)", "vector of reals"});

      keywords.insert({"R2_jwl", "JWL EoS R2 parameter",
        R"(This keyword is used to specify the material property R2 for the
        Jones-Wilkins-Lee equation of state.)", "vector of reals"});

      keywords.insert({"rho0_jwl", "JWL EoS rho0 parameter",
        R"(This keyword is used to specify the material property rho0, which is
        the density of initial state (units: kg/m3) for the Jones-Wilkins-Lee
        equation of state.)", "vector of reals"});

      keywords.insert({"de_jwl", "JWL EoS de parameter",
        R"(This keyword is used to specify the material property de, which is the
        heat of detonation for products; and for reactants, it is chosen such that
        the ambient internal energy (e0) is 0 (units: J/kg). Used for the
        Jones-Wilkins-Lee equation of state.)", "vector of reals"});

      keywords.insert({"rhor_jwl", "JWL EoS rhor parameter",
        R"(This keyword is used to specify the material property rhor, which is
        the density of reference state (units: kg/m3) for the Jones-Wilkins-Lee
        equation of state.)", "vector of reals"});

      keywords.insert({"Tr_jwl", "JWL EoS Tr parameter",
        R"(This keyword is used to specify the material property Tr, which is the
        temperature of reference state (units: K) for the Jones-Wilkins-Lee
        equation of state.)", "vector of reals"});

      keywords.insert({"Pr_jwl", "JWL EoS er parameter",
        R"(This keyword is used to specify the material property Pr, which is the
        pressure at the reference state (units: Pa) for the Jones-Wilkins-Lee
        equation of state. It is used to calculate the reference temperature for
        the EoS.)", "vector of reals"});

      keywords.insert({"mu", "shear modulus/dynamic viscosity",
        R"(This keyword is used to specify the material property, shear modulus
        for solids, or dynamic viscosity for fluids.)", "vector of reals"});

      keywords.insert({"yield_stress", "Yield stress of solid material",
        R"(This keyword is used to specify the material property yield stress,
        which indicates the stress (units: Pa) after which the material begins
        plastic flow.)", "vector of reals"});

      keywords.insert({"cv", "specific heat at constant volume",
        R"(This keyword is used to specify the material property, specific heat at
        constant volume.)", "vector of reals"});

      keywords.insert({"k", "heat conductivity",
        R"(This keyword is used to specify the material property, heat
        conductivity.)", "vector of reals"});

      keywords.insert({"stiffenedgas",
        "Select the stiffened gas equation of state",
        R"(This keyword is used to select the stiffened gas equation of state.)"});

      keywords.insert({"jwl", "Select the JWL equation of state",
        R"(This keyword is used to select the Jones, Wilkins, Lee equation of
        state.)"});

      keywords.insert({"smallshearsolid",
        "Select the SMALLSHEARSOLID equation of state",
        R"(This keyword is used to select the small shear strain equation of state
        for solids. This EOS uses a small-shear approximation for the elastic
        contribution, and a stiffened gas EOS for the hydrodynamic contribution of
        the internal energy See Plohr, J. N., & Plohr, B. J. (2005). Linearized
        analysis of Richtmyer–Meshkov flow for elastic materials. Journal of Fluid
        Mechanics, 537, 55-89 for further details.)"});

      keywords.insert({"godunovromenski_aluminum",
        "Select the GODUNOVROMENSKIALUMINUM equation of state",
        R"(This keyword is used to select the Godunov-Romenski equation of
        state for solids and a hydro EoS for aluminum. These function were
        taken from Barton, Philip T. "An interface-capturing Godunov method
        for the simulation of compressible solid-fluid problems." Journal
        of Computational Physics 390 (2019): 25-50.)"});

      keywords.insert({"matidxmap",
      "AUTO-GENERATED Material index map for EOS",
      R"(The following AUTO-GENERATED data structure is used to index into the
      correct material vector entry. This is done using the following three maps:
      1. eosidx: This vector provides the eos-index (value) in the
      vector<tag::material> for the given user-spec material id (index).
      2. matidx: This vector provides the material-index (value) inside the
      vector<tag::material>[eosidx] block for the given user-specified
      material id (index).
      3. solidx: This vector provides the solid-index (value) assigned to
      the given user-specified material id (index). It is 0 for fluids.)"});

      // -----------------------------------------------------------------------
      // output object
      // -----------------------------------------------------------------------

      keywords.insert({"field_output",
        "Start of field_output input block",
        R"(This keyword is used to start a block in the input file containing the
        list and settings of requested field output.)", "block-title"});

      keywords.insert({"interval",
        "Set interval (in units of iteration count)",
        R"(This keyword is used to specify an interval in units of iteration count
        (i.e., number of time steps). This must be used within a relevant
        block.)", "uint"});

      keywords.insert({"time_interval",
        "Set interval (in units of physics time)",
        R"(This keyword is used to specify an interval in units of physics time.
        This must be used within a relevant block.)", "real"});

      keywords.insert({"time_range",
        "Configure physics time range for output (in units of physics time)",
        R"(This keyword is used to configure field-, or history-output, specifying
        a start time, a stop time, and an output frequency in physics time units.
        Example: 'time_range = {0.2, 0.3, 0.001}', which specifies that from t=0.2 to
        t=0.3 output should happen at physics time units of dt=0.001. This must be
        used within a relevant block.)", "vector of 3 reals"});

      keywords.insert({"refined", "Toggle refined field output on/off",
        R"(This keyword can be used to turn on/off refined field output, which
        refines the mesh and evaluates the solution on the refined mesh for saving
        the solution.)", "bool"});

      keywords.insert({"filetype", "Select output file type",
        R"(This keyword is used to specify the output file type of
        mesh-based field output in a field_output block.)", "string" });

      keywords.insert({"elemvar",
        "Specify list of elem-centered variables for output",
        R"(This keyword is used to specify elem-centered variables for output to
        file. It is used in field_output blocks.)", "vector of string"});

      keywords.insert({"nodevar",
        "Specify list of node-centered variables for output",
        R"(This keyword is used to specify node-centered variables for output to
        file. It is used in field_output blocks.)", "vector of string"});

      keywords.insert({"sideset",
        "Specify list of side sets",
        R"(This keyword is used to specify side sets. Usage is context specific,
        i.e. depends on what block it is specified in. Eg. in the field_output
        block it specifies the sidesets on which field output is desired.)",
        "vector of uints"});

      keywords.insert({"diagnostics",
        "Specify the diagnostics block",
        R"(This keyword is used to introduce the dagnostics block, used to
        configure diagnostics output.)", "block-title"});

      keywords.insert({"interval",
        "Set interval (in units of iteration count)",
        R"(This keyword is used to specify an interval in units of iteration count
        (i.e., number of time steps). Usage is context specific, and this must
        be used within a relevant block.)", "uint"});

      keywords.insert({"error", "Select an error type",
        R"(This keyword is used to select the error type. Used either in the
        diagnostics block to specify error norm, or in the AMR block to specify
        the error for solution-adaptive mesh refinement.)", "string"});

      keywords.insert({"format",
        "Specify the ASCII floating-point output format",
        R"(This keyword is used to select the
        floating-point output format for ASCII floating-point number output.
        Valid options are 'default', 'fixed', and 'scientific'. For more info on
        these various formats, see
        http://en.cppreference.com/w/cpp/io/manip/fixed.)", "string"});

      keywords.insert({"precision",
        "Precision in digits for ASCII floating-point output",
        R"(This keyword is used to select
        the precision in digits for ASCII floating-point real number output.
        Example: "precision=10", which selects ten digits for floating-point
        output, e.g., 3.141592654. The number of digits must be larger than zero
        and lower than the maximum representable digits for the given
        floating-point type. For more info on setting the precision in C++, see
        http://en.cppreference.com/w/cpp/io/manip/setprecision, and
        http://en.cppreference.com/w/cpp/types/numeric_limits/digits10)", "uint"});

      keywords.insert({"history_output",
        "Start of history_output input block",
        R"(This keyword is used to start a block in the input file containing the
        descriptions and settings of requested history output.)", "string"});

      keywords.insert({"point",
        "Start configuration block for history point, or a single point in IC blocks",
        R"(This keyword is used to either introduce a vector block used to
        specify probes for history output, or in the IC/BC block to specify
        a single point. When used in history output, it takes sub-entries of
        'id' and 'coord'. When used in IC/BC, directly takes three reals as
        coordinates.)", "vector block-title"});

      keywords.insert({"coord", "Specify point coordinates",
        R"(This keyword is used to specify coordinates of the history-point.)",
        "3 reals"});

      keywords.insert({"exodusii", "Select ExodusII output",
        R"(This keyword is used to select the
        ExodusII output file type readable by, e.g., ParaView of either a requested
        probability density function (PDF) within a pdfs ... end block or for
        mesh-based field output in a field_output ... end block. Example:
        "filetype exodusii", which selects ExodusII file output. For more info on
        ExodusII, see http://sourceforge.net/projects/exodusii.)"});

      keywords.insert({"density", "Request/specify density",
        R"(This keyword is used to request/specify the density. Usage is
        context specific. When specifed as 'elemvar' or 'nodevar' inside
        field_output, requests density as an output quantity. Otherwise,
        specifies density at IC/BC.)", "real"});

      keywords.insert({"x-momentum",
        "Request x-momentum",
        R"(This keyword is used to request the fluid x-momentum as an output
        variable.)"});

      keywords.insert({"y-momentum",
        "Request y-momentum",
        R"(This keyword is used to request the fluid y-momentum as an output
        variable.)"});

      keywords.insert({"z-momentum",
        "Request z-momentum",
        R"(This keyword is used to request the fluid z-momentum as an output
        variable.)"});

      keywords.insert({
        "specific_total_energy", "Request specific total energy",
        R"(This keyword is used to request the specific total energy as an output
        variable.)"});

      keywords.insert({
        "volumetric_total_energy", "Request total volumetric energy",
        R"(This keyword is used to request the volumetric total energy as an output
        variable.)"});

      keywords.insert({"x-velocity",
        "Request x-velocity",
        R"(This keyword is used to request the fluid x-velocity as an output
        variable.)"});

      keywords.insert({"y-velocity",
        "Request y-velocity",
        R"(This keyword is used to request the fluid y-velocity as an output
        variable.)"});

      keywords.insert({"z-velocity",
        "Request z-velocity",
        R"(This keyword is used to request the fluid z-velocity as an output
        variable.)"});

      keywords.insert({"pressure", "Request/specify pressure",
        R"(This keyword is used to request/specify the pressure. Usage is
        context specific. When specifed as 'elemvar' or 'nodevar' inside
        field_output, requests pressure as an output quantity. Otherwise,
        specifies pressure at IC/BC.)", "real"});

      keywords.insert({"material_indicator",
        "Request material_indicator",
        R"(This keyword is used to request the material indicator function as an
        output variable.)"});

      keywords.insert({"analytic",
        "Request analytic solution",
        R"(This keyword is used to request the analytic solution (if exist) as an
        output variable.)"});

      keywords.insert({"l2", "Select the L2 norm",
        R"(This keyword is used to enable computing the L2 norm.)"});

      keywords.insert({"linf", "Select the L_{infinity} norm",
        R"(This keyword is used to enable computing the L-infinity norm.)"});

      keywords.insert({"default",
        "Select the default ASCII floating-point output",
        R"(This keyword is used to select the
        'default' floating-point output format for ASCII floating-point real
        number output. For more info on these various formats, see
        http://en.cppreference.com/w/cpp/io/manip/fixed.)"});

      keywords.insert({"fixed",
        "Select the fixed ASCII floating-point output",
        R"(This keyword is used to select the
        'fixed' floating-point output format for ASCII floating-point real
        number output. For more info on these various formats, see
        http://en.cppreference.com/w/cpp/io/manip/fixed.)"});

      keywords.insert({"scientific",
        "Select the scientific ASCII floating-point output",
        R"(This keyword is used to select the
        'scientific' floating-point output format for ASCII floating-point real
        number output. For more info on these various formats, see
        http://en.cppreference.com/w/cpp/io/manip/fixed.)"});

      // -----------------------------------------------------------------------
      // ALE options
      // -----------------------------------------------------------------------

      keywords.insert({"ale", "Start block configuring ALE",
        R"(This keyword is used to introduce the ale block, used to
        configure arbitrary Lagrangian-Eulerian (ALE) mesh movement.)",
        "block-title"});

      keywords.insert({"smoother", "Select mesh velocity smoother",
        R"(This keyword is used to select a mesh velocity smoother option, used
        for Arbitrary-Lagrangian-Eulerian (ALE) mesh motion. Valid options are
        'laplace', 'helmholtz', and 'none')", "string"});

      keywords.insert({"mesh_velocity", "Select mesh velocity",
        R"(This keyword is used to select a mesh velocity option, used for
        Arbitrary-Lagrangian-Eulerian (ALE) mesh motion. Valid options are
        'sine', 'fluid', and 'user_defined".)", "string"});

      keywords.insert({"mesh_motion",
        "List of dimension indices that are allowed to move in ALE calculations",
        R"(This keyword is used to specify a list of integers (0, 1, or 2) whose
        coordinate directions corresponding to x, y, or z are allowed to move with
        the mesh velocity in ALE calculations. Useful for 1D/2D problems.)",
        "vector of uints"});

      keywords.insert({"meshforce", "Set ALE mesh force model parameter(s)",
        R"(This keyword is used to specify a vector of real numbers used to
        parameterize a mesh force model for ALE. The length of the vector must
        exactly be 4.)", "vector of 4 reals"});

      keywords.insert({"dvcfl",
        "Set the volume-change Courant-Friedrichs-Lewy (CFL) coefficient",
        R"(This keyword is used to specify the volume-change (dV/dt) CFL coefficient
        for variable-time-step-size simulations due to volume change in time in
        arbitrary-Lagrangian-Eulerian (ALE) calculations. Setting 'dvcfl' only has
        effect in ALE calculations and used together with 'cfl'. See also J. Waltz,
        N.R. Morgan, T.R. Canfield, M.R.J. Charest, L.D. Risinger, J.G. Wohlbier, A
        three-dimensional finite element arbitrary Lagrangian–Eulerian method for
        shock hydrodynamics on unstructured grids, Computers & Fluids, 92: 172-187,
        2014.)", "real"});

      keywords.insert({"vortmult",
        "Configure vorticity multiplier for ALE mesh velocity",
        R"(This keyword is used to configure the multiplier for the vorticity term
        in the mesh velocity smoother (mesh_velocity=fluid) or for the potential
        gradient for the Helmholtz mesh velocity (mesh_velocity=helmholtz) for ALE
        mesh motion. For 'fluid' this is coefficient c2 in Eq.(36) of Waltz,
        Morgan, Canfield, Charest, Risinger, Wohlbier, A three-dimensional finite
        element arbitrary Lagrangian–Eulerian method for shock hydrodynamics on
        unstructured grids, Computers & Fluids, 2014, and for 'helmholtz', this
        is coefficient a1 in Eq.(23) of Bakosi, Waltz, Morgan, Improved ALE mesh
        velocities for complex flows, International Journal for Numerical Methods
        in Fluids, 2017. )", "real"});

      keywords.insert({"maxit",
        "Set the max number of iterations for the ALE mesh velocity linear solve",
        R"(This keyword is used to specify the maximum number of linear solver
        iterations taken to converge the mesh velocity linear solve in
        arbitrary-Lagrangian-Eulerian (ALE) calculations. See also J. Waltz,
        N.R. Morgan, T.R. Canfield, M.R.J. Charest, L.D. Risinger, J.G. Wohlbier, A
        three-dimensional finite element arbitrary Lagrangian–Eulerian method for
        shock hydrodynamics on unstructured grids, Computers & Fluids, 92: 172-187,
        2014.)", "uint"});

      keywords.insert({"tolerance",
        "Set the tolerance for the ALE mesh velocity linear solve",
        R"(This keyword is used to specify the tolerance to converge the mesh
        velocity linear solve for in
        arbitrary-Lagrangian-Eulerian (ALE) calculations. See also J. Waltz,
        N.R. Morgan, T.R. Canfield, M.R.J. Charest, L.D. Risinger, J.G. Wohlbier, A
        three-dimensional finite element arbitrary Lagrangian–Eulerian method for
        shock hydrodynamics on unstructured grids, Computers & Fluids, 92: 172-187,
        2014.)", "real"});

      keywords.insert({"move",
        "Start configuration block configuring surface movement in ALE",
        R"(This keyword is used to introduce a move block, used to
        configure surface movement for ALE simulations.)", "vector block-title"});

      keywords.insert({"fntype",
        "Select how a user-defined function is interpreted",
        R"(This keyword is used to select how a user-defined function should be
        interpreted.)", "string"});

      keywords.insert({"fn", "Specify a discrete user-defined function",
        R"(This keyword is used to specify a user-defined function block with
        discrete points, listed inside the fn block. Used in ale mesh motion and
        time-dependent BC specification)", "reals" });

      keywords.insert({"laplace",
        "Select the Laplace mesh velocity smoother for ALE",
        R"(This keyword is used to select the 'Laplace' mesh velocity smoother
        for Arbitrary-Lagrangian-Eulerian (ALE) mesh motion.)"});

      keywords.insert({"helmholtz",
        "Select the Helmholtz velocity for ALE",
        R"(This keyword is used to select the a velocity, computed from the
        Helmholtz-decomposition as the mesh velocity for
        Arbitrary-Lagrangian-Eulerian (ALE) mesh motion. See J. Bakosi, J. Waltz,
        N. Morgan, Improved ALE mesh velocities for complex flows, Int. J. Numer.
        Meth. Fl., 1-10, 2017, https://doi.org/10.1002/fld.4403.)"});

      keywords.insert({"none", "Select none option",
        R"(This keyword is used to select the 'none' option from a list of
        configuration options.)"});

      keywords.insert({"sine",
        "Prescribe sinusoidal mesh velocity for ALE",
        R"(This keyword is used to prescribe a sinusoidal mesh velocity
        for Arbitrary-Lagrangian-Eulerian (ALE) mesh motion.)"});

      keywords.insert({"fluid", "Select the fluid velocity for ALE",
        R"(This keyword is used to select the 'fluid' velocity as the mesh velocity
        for Arbitrary-Lagrangian-Eulerian (ALE) mesh motion.)"});

      // -----------------------------------------------------------------------
      // h/p adaptation objects
      // -----------------------------------------------------------------------

      keywords.insert({"amr",
        "Start configuration block configuring adaptive mesh refinement",
        R"(This keyword is used to introduce the amr block, used to
        configure adaptive mesh refinement.)", "block-title"});

      keywords.insert({"t0ref", "Enable mesh refinement at t<0",
        R"(This keyword is used to enable initial mesh refinement, which can be
        configured to perform multiple levels of mesh refinement based on various
        refinement criteria and configuration settings.)", "bool"});

      keywords.insert({"dtref", "Enable mesh refinement at t>0",
        R"(This keyword is used to enable solution-adaptive mesh refinement during
        time stepping.)", "bool"});

      keywords.insert({"dtref_uniform",
        "Enable mesh refinement at t>0 but only perform uniform refinement",
        R"(This keyword is used to force uniform-only solution-adaptive mesh
        refinement during time stepping.)", "bool"});

      keywords.insert({"dtfreq",
        "Set mesh refinement frequency during time stepping",
        R"(This keyword is used to configure the frequency of mesh refinement
        during time stepping.)", "uint"});

      keywords.insert({"maxlevels",
        "Set maximum allowed mesh refinement levels",
        R"(This keyword is used to configure the maximum allowed mesh refinement
        levels.)", "uint"});

      keywords.insert({"initial",
        "Configure initial mesh refinement (before time stepping)",
        R"(This keyword is used to add to a list of initial mesh refinement types
        that happens before t = 0. Allowed options are 'uniform',
        'uniform_derefine', 'initial_conditions', 'coords', 'edgelist')",
        "vector of strings"});

      keywords.insert({"coords",
        "Configure initial refinement using coordinate planes",
        R"(This keyword can be used to configure entire volumes on a given side of
        a plane in 3D space. The keyword introduces an coords block within
        an amr block. All edges of the input mesh will be tagged for refinement
        whose end-points lie within the given ranges.
        Example: 'xminus 0.5' refines all edges whose end-point coordinates are
        less than 0.5. Multiple specifications are understood by combining with
        a logical AND. That is: 'xminus 0.5 yplus 0.3' refines all edges whose
        end-point x coordinates are less than 0.5 AND y coordinates are larger than
        0.3.)", "block-title"});

      keywords.insert({"xminus",
        "Configure initial refinement for coordinates lower than an x-normal plane",
        R"(This keyword can be used to configure a mesh refinement volume for edges
        whose end-points are less than the x coordinate of a plane perpendicular
        to coordinate x in 3D space. The keyword must be used in a coords-block
        within an amr-block with syntax 'xminus <real>'. All edges of the
        input mesh will be tagged for refinement whose end-points lie less than (-)
        the real number given. Example: 'xminus 0.5' refines all edges whose end-point
        x-coordinates are less than 0.5.)", "real"});

      keywords.insert({"xplus",
        "Configure initial refinement for coordinates larger than an x-normal plane",
        R"(This keyword can be used to configure a mesh refinement volume for edges
        whose end-points are larger than the x coordinate of a plane perpendicular
        to coordinate x in 3D space. The keyword must be used in a coords-block
        within an amr-block with syntax 'xplus <real>'. All edges of the
        input mesh will be tagged for refinement whose end-points lie larger than
        (+) the real number given. Example: 'xplus 0.5' refines all edges whose
        end-point coordinates are larger than 0.5.)", "real"});

      keywords.insert({"yminus",
        "Configure initial refinement for coordinates lower than an y-normal plane",
        R"(This keyword can be used to configure a mesh refinement volume for edges
        whose end-points are less than the y coordinate of a plane perpendicular
        to coordinate y in 3D space. The keyword must be used in a coords-block
        within an amr-block with syntax 'yminus <real>'. All edges of the
        input mesh will be tagged for refinement whose end-points lie less than (-)
        the real number given. Example: 'yminus 0.5' refines all edges whose end-point
        coordinates are less than 0.5.)", "real"});

      keywords.insert({"yplus",
        "Configure initial refinement for coordinates larger than an y-normal plane",
        R"(This keyword can be used to configure a mesh refinement volume for edges
        whose end-points are larger than the y coordinate of a plane perpendicular
        to coordinate y in 3D space. The keyword must be used in a coords-block
        within an amr-block with syntax 'yplus <real>'. All edges of the
        input mesh will be tagged for refinement whose end-points lie larger than
        (+) the real number given. Example: 'yplus 0.5' refines all edges whose
        end-point coordinates are larger than 0.5.)", "real"});

      keywords.insert({"zminus",
        "Configure initial refinement for coordinates lower than an z-normal plane",
        R"(This keyword can be used to configure a mesh refinement volume for edges
        whose end-points are less than the z coordinate of a plane perpendicular
        to coordinate z in 3D space. The keyword must be used in a coords-block
        within an amr-block with syntax 'zminus <real>'. All edges of the
        input mesh will be tagged for refinement whose end-points lie less than (-)
        the real number given. Example: 'zminus 0.5' refines all edges whose end-point
        coordinates are less than 0.5.)", "real"});

      keywords.insert({"zplus",
        "Configure initial refinement for coordinates larger than an z-normal plane",
        R"(This keyword can be used to configure a mesh refinement volume for edges
        whose end-points are larger than the z coordinate of a plane perpendicular
        to coordinate z in 3D space. The keyword must be used in a coords-block
        within an amr-block with syntax 'zplus <real>'. All edges of the
        input mesh will be tagged for refinement whose end-points lie larger than
        (+) the real number given. Example: 'zplus 0.5' refines all edges whose
        end-point coordinates are larger than 0.5.)", "real"});

      keywords.insert({"edgelist",
        "Configure edge-node pairs for initial refinement",
        R"(This keyword can be used to configure a list of edges that are explicitly
        tagged for initial refinement during setup in inciter. The keyword
        introduces an edgelist block within an amr block and must
        contain a list of integer pairs, i.e., the number of ids must be even,
        denoting the end-points of the nodes (=edge) which should be tagged for
        refinement.)", "vector of uints"});

      keywords.insert({"error",
        "Configure the error type for solution-adaptive mesh refinement",
        R"(This keyword is used to select the algorithm used to estimate the error
        for solution-adaptive mesh refinement. Available options are 'jump' and
        'hessian')", "string"});

      keywords.insert({"refvar",
        "Configure dependent variables used for adaptive mesh refinement",
        R"(This keyword is used to configured a list of dependent variables that
        trigger adaptive mesh refinement based on estimating their numerical error.
        These refinement variables are used for both initial (i.e., before time
        stepping) mesh refinement as well as during time stepping. Only previously
        (i.e., earlier in the input file) selected dependent variables can be
        configured as refinement variables. Dependent variables are required to be
        defined in all equation system configuration blocks, e.g., transport ...
        end, by using the 'depvar' keyword. Example: transport depvar c end amr
        refvar c end end. Selecting a particular scalar component in a system is
        done by appending the equation number to the refvar: Example: transport
        depvar q ncomp 3 end amr refvar q1 q2 end end, which configures two
        refinement variables: the first and third scalar component of the previously
        configured transport equation system.)", "vector of char"});

      keywords.insert({"tol_refine", "Configure refine tolerance",
        R"(This keyword is used to set the tolerance used to tag an edge for
        refinement if the relative error exceeds this value.)", "real"});

      keywords.insert({"tol_derefine",
        "Configure derefine tolerance",
        R"(This keyword is used to set the tolerance used to tag an edge for
        derefinement if the relative error is below this value.)", "real"});

      keywords.insert({"uniform",
        "Select uniform initial mesh refinement",
        R"(This keyword is used to select uniform initial mesh refinement.)",
        "string"});

      keywords.insert({"uniform_derefine",
        "Select uniform initial mesh de-refinement",
        R"(This keyword is used to select uniform initial mesh de-refinement.)",
        "string"});

      keywords.insert({"initial_conditions",
        "Select initial-conditions-based initial mesh refinement",
        R"(This keyword is used to select initial-conditions-based initial mesh
        refinement.)", "string"});

      keywords.insert({"jump",
        "Error estimation based on the solution jump normalized by solution value",
        R"(This keyword is used to select the jump-based error indicator for
        solution-adaptive mesh refinement. The error is estimated by computing the
        magnitude of the jump in the solution value normalized by the solution
        value.)", "string"});

      keywords.insert({"hessian",
        "Error estimation based on the Hessian normalized by solution value",
        R"(This keyword is used to select the Hessian-based error indicator for
        solution-adaptive mesh refinement. The error is estimated by computing the
        Hessian (2nd derivative matrix) of the solution normalized by sum of the
        absolute values of the gradients at edges-end points.)", "string"});

      keywords.insert({"pref",
        "Start configuration block configuring p-adaptive refinement",
        R"(This keyword is used to introduce the pref block, to
        configure p-adaptive refinement)", "block-title"});

      keywords.insert({"indicator",
        "Configure the specific adaptive indicator for p-adaptive DG scheme",
        R"(This keyword can be used to configure a specific type of adaptive
        indicator for p-adaptive refinement  of the DG scheme. The keyword must
        be used in a pref block. Available options are 'pref_spectral_decay' and
        'pref_non_conformity'.)", "string"});

      keywords.insert({"ndofmax",
        "Configure the maximum number of degree of freedom for p-adaptive DG",
        R"(This keyword can be used to configure a maximum number of degree of
        freedom for p-adaptive refinement  of the DG scheme. The keyword must
        be used in a pref block.)", "uint either 4 or 10"});

      keywords.insert({"tolref",
        "Configure the tolerance for p-refinement for p-adaptive DG",
        R"(This keyword can be used to configure a tolerance for p-adaptive
        refinement  for the DG scheme. The keyword must be used in a pref
        block. All elements with a refinement indicator larger than this
        tolerance will be p-refined.)", "real between 0 and 1"});

      keywords.insert({"spectral_decay",
        "Select the spectral-decay indicator for p-adaptive DG scheme",
        R"(This keyword is used to select the spectral-decay indicator used for
        p-adaptive discontinuous Galerkin (DG) discretization used in inciter.
        See Control/Inciter/Options/PrefIndicator.hpp for other valid options.)",
        "string"});

      keywords.insert({"non_conformity",
        "Select the non-conformity indicator for p-adaptive DG scheme",
        R"(This keyword is used to select the non-conformity indicator used for
        p-adaptive discontinuous Galerkin (DG) discretization used in inciter.
        See Control/Inciter/Options/PrefIndicator.hpp for other valid options.)",
        "string"});

      // -----------------------------------------------------------------------
      // boundary condition options
      // -----------------------------------------------------------------------

      keywords.insert({"bc",
        "Start configuration block for boundary conditions",
        R"(This keyword is used to introduce the bc block, used for
        boundary conditions. This is a vector block, where each vector entry
        specifies BCs for a particular mesh)", "vector block-title"});

      keywords.insert({"mesh",
        "List meshes on which the following BCs apply",
        R"(This keyword is used to list multiple meshes on which the boundary
        conditions listed in this particular bc-block apply.)",
        "vector of uints"});

      keywords.insert({"dirichlet",
        "List sidesets with Dirichlet boundary conditions",
        R"(This keyword is used to list Dirichlet sidesets.
        This keyword is used to list multiple sidesets on
        which a prescribed Dirichlet BC is then applied. Such prescribed BCs
        at each point in space and time are evaluated using a built-in function,
        e.g., using the method of manufactured solutions.)", "vector of uint(s)"});

      keywords.insert({"symmetry",
        "List sidesets with symmetry boundary conditions",
        R"(This keyword is used to list (multiple) symmetry BC sidesets.)",
        "vector of uint(s)"});

      keywords.insert({"inlet",
        "List sidesets with inlet boundary conditions",
        R"(This keyword is used to list (multiple) inlet BC sidesets.)",
        "vector of uint(s)"});

      keywords.insert({"outlet",
        "List sidesets with outlet boundary conditions",
        R"(This keyword is used to list (multiple) outlet BC sidesets.)",
        "vector of uint(s)"});

      keywords.insert({"farfield",
        "List sidesets with farfield boundary conditions",
        R"(This keyword is used to list (multiple) farfield BC sidesets.
        Keywords allowed in a bc_farfield block are 'density', 'velocity',
        'pressure')", "vector of uint(s)"});

      keywords.insert({"extrapolate",
        "List sidesets with Extrapolation boundary conditions",
        R"(This keyword is used to list (multiple) extrapolate BC sidesets.)",
        "vector of uint(s)"});

      keywords.insert({"noslipwall",
        "List sidesets with no-slip wall boundary conditions",
        R"(This keyword is used to list (multiple) no-slip wall BC sidesets.)",
        "vector of uint(s)"});

      keywords.insert({"stag",
        "List sidesets with stagnation boundary conditions",
        R"(This keyword is used to list (multiple) stagnation BC sidesets.)",
        "vector of uint(s)"});

      keywords.insert({"timedep",
        "Start configuration block describing time dependent boundary conditions",
        R"(This keyword is used to introduce a bc_timedep block, used to
        specify the configuration of time dependent boundary conditions for a
        partial differential equation. A discrete function in time t in the form
        of a table with 6 columns (t, pressure(t), density(t), vx(t), vy(t), vz(t))
        is expected inside a fn ... end block, specified within the bc_timedep
        block. Multiple such bc_timedep blocks can be specified for different
        time dependent BCs on different groups of side sets.)", "block-title"});

      keywords.insert({"radius", "Specify a radius",
        R"(This keyword is used to specify a radius, used, e.g., in specifying a
        point in 3D space for setting a stagnation (velocity vector = 0).)",
        "real"});

      keywords.insert({"velocity", "Specify velocity",
        R"(This keyword is used to configure a velocity vector used in a
        context-specific way, e.g., for boundary or initial conditions, or
        specifying overset mesh velocity.)", "vector of 3 reals"});

      // -----------------------------------------------------------------------
      // IC object
      // -----------------------------------------------------------------------

      keywords.insert({"ic",
        "Introduce an ic block used to configure initial conditions",
        R"(This keyword is used to introduce an ic block used to set initial
        conditions.)", "block-title"});

      keywords.insert({"materialid", "Specify material id",
        R"(This keyword is used to configure the material id within an IC box,
        IC mesh-block, farfield BC, or in the background as a part of the
        initialization.)", "uint"});

      keywords.insert({"temperature", "Specify temperature",
        R"(This keyword is used to configure temperature, used for, e.g.,
        boundary or initial conditions.)" , "real"});

      keywords.insert({"mass_fractions", "Specify species mass fractions",
        R"(This keyword is used to configure species mass fractions, used for,
        e.g., boundary or initial conditions.)" , "vector of reals"});

      keywords.insert({"box",
        "Introduce a box block used to assign initial conditions",
        R"(This keyword is used to introduce a IC box block used to assign
        initial conditions within a box given by spatial coordinates.)",
        "vector block-title"});

      keywords.insert({"meshblock",
        "Introduce a meshblock block used to assign initial conditions",
        R"(This keyword is used to introduce a IC meshblock block used to
        assign initial conditions within a mesh block specified in the mesh file.)",
        "vector block-title"});

      keywords.insert({"blockid", "Specify mesh block id",
        R"(This keyword is used to configure the mesh block id within the
        meshblock-block as a part of the initialization. It is strongly
        recommended to use contiguous block ids in mesh file starting from 1.)",
        "uint"});

      keywords.insert({"volume", "Specify volume",
        R"(This keyword is used to configure the volume of a meshblock.)",
        "real"});

      keywords.insert({"mass", "Specify mass",
        R"(This keyword is used to configure the mass within a box/meshblock.)",
        "real"});

      keywords.insert({"energy", "Specify energy per unit mass",
        R"(This keyword is used to configure energy per unit mass, used for, e.g.,
        boundary or initial conditions.)", "real"});

      keywords.insert({"energy_content", "Specify energy per unit volume",
        R"(This keyword is used to configure energy per unit volume, used for
        initial conditions.)", "real"});

      keywords.insert({"xmin", "Minimum x coordinate",
        R"(This keyword used to configure a minimum x coordinate to specify
        a box.)", "real"});

      keywords.insert({"xmax", "Maximum x coordinate",
        R"(This keyword used to configure a maximum x coordinate to specify
        a box.)", "real"});

      keywords.insert({"ymin", "Minimum y coordinate",
        R"(This keyword used to configure a minimum y coordinate to specify
        a box.)", "real"});

      keywords.insert({"ymax", "Maximum y coordinate",
        R"(This keyword used to configure a maximum y coordinate to specify
        a box.)", "real"});

      keywords.insert({"zmin", "Minimum z coordinate",
        R"(This keyword used to configure a minimum z coordinate to specify
        a box.)", "real"});

      keywords.insert({"zmax", "Maximum z coordinate",
        R"(This keyword used to configure a maximum z coordinate to specify
        a box.)", "real"});

      keywords.insert({"orientation", "Configure orientation",
        R"(Configure orientation of an IC box for rotation about centroid of box.)",
        "vector of 3 reals"});

      keywords.insert({"initiate", "Initiation type",
        R"(This keyword is used to select an initiation type to configure how
        values are assigned for a box/meshblock initialization. This can be used
        to specify, how the values are assigned to mesh nodes within a box. Uses:
        (1) impulse: assign the full values at t=0 for all points in a box,
        (2) linear: use a linear function in time and space, configured with an
        initiation point in space, a constant velocity of the growing spherical
        front in time (and space) linearly, and width of the front and assigns
        values to mesh points falling within the growing spherical shell inside
        a configured box.)", "string"});

      keywords.insert({"init_time","Specify the initialization time",
        R"(This keyword is used to specify the time at which the propagating front
        is initialized for a mesh block or box IC, with 'initiate linear' type.
        Delays in initializing separate mesh blocks or boxes can be achieved using
        different initialization times.)", "real"});

      keywords.insert({"front_width", "Specify a front width",
        R"(This keyword is used to specify the width of the propagating front for
        a mesh block or box IC, with 'initiate linear' type. The suggested value
        of the front width is about 4-5 times the mesh size inside the mesh block
        or box.)", "real"});

      keywords.insert({"front_speed", "Specify a front speed",
        R"(This keyword is used to specify the speed at which a front propagates
        for a mesh block or box IC, with 'initiate linear' type.)", "real"});

      keywords.insert({"impulse",
        "Select the impulse initiation type, for a box/meshblock IC",
        R"(This keyword can be used to select the 'impulse' initiation/assignment
        type for box initial conditions. It simply assigns the prescribed values
        to all mesh points within a configured box at t=0.)"});

      keywords.insert({"linear",
        "Select the linear initiation type, for a box/meshblock IC",
        R"(This keyword is be used to specify the 'linear' initiation parameters
        for a particular box or meshblock, as a part of the 'energy_pill'
        initialization. Linear initiation uses a linear function in time and space,
        configured with an initiation point in space, a constant velocity of the
        growing spherical front in time (and space) linearly, and width of the front
        and assigns values to mesh points falling within the growing spherical shell
        inside a configured box or meshblock. The following keywords are required
        in the box/meshblock block if 'linear' is used: 'init_time',
        'front_width', 'front_speed')"});

      // -----------------------------------------------------------------------
      // Overset mesh object
      // -----------------------------------------------------------------------

      keywords.insert({"mesh",
        "Start configuration block assigning a mesh to a solver",
        R"(This keyword is used to introduce a mesh block, used to
        assign and configure a mesh to a solver.)", "vector block-title"});

      keywords.insert({"filename", "Set filename",
        R"(Set filename, e.g., mesh filename for solver coupling.)", "string"});

      keywords.insert({"location", "Configure location",
        R"(Configure location of a mesh relative to another.)",
        "vector of 3 reals"});

      keywords.insert({"orientation", "Configure orientation",
        R"(Configure orientation of a mesh relative to another.)",
        "vector of 3 reals"});

      // -----------------------------------------------------------------------
      // pre-configured problems
      // -----------------------------------------------------------------------

      keywords.insert({"user_defined",
        "Select user-defined specification for a problem",
        R"(This keyword is used to select the user-defined specification for an
        option. This could be a 'problem' to be solved by a partial differential
        equation, but can also be a 'user-defined' mesh velocity specification for
        ALE mesh motion.)"});

      keywords.insert({"shear_diff",
        "Select the shear + diffusion test problem ",
        R"(This keyword is used to select the shear diffusion test problem. The
        initial and boundary conditions are specified to set up the test problem
        suitable to exercise and test the advection and diffusion terms of the
        scalar transport equation.)" });

      keywords.insert({"slot_cyl",
        "Select Zalesak's slotted cylinder test problem",
        R"(This keyword is used to select Zalesak's slotted cylinder test
        problem. The initial and boundary conditions are specified to set up the
        test problem suitable to exercise and test the advection and diffusion
        terms of the scalar transport equation.)"});

      keywords.insert({"gauss_hump",
        "Select advection of 2D Gaussian hump test problem",
        R"(This keyword is used to select the advection of 2D Gaussian hump test
        problem. The initial and boundary conditions are specified to set up the
        test problem suitable to exercise and test the advection
        terms of the scalar transport equation.)"});

      keywords.insert({"cyl_advect",
        "Select advection of cylinder test problem",
        R"(This keyword is used to select the advection of cylinder test
        problem. The initial and boundary conditions are specified to set up the
        test problem suitable to exercise and test the advection
        terms of the scalar transport equation.)"});

      keywords.insert({"cyl_vortex",
        "Select deformation of cylinder in a vortex test problem",
        R"(This keyword is used to select the test problem which deforms a cylinder
        in a vortical velocity field. The initial and boundary conditions are
        specified to set up the test problem suitable to exercise and test the
        advection terms of the scalar transport equation.)"});

      keywords.insert({"vortical_flow",
        "Select the vortical flow test problem ",
        R"(This keyword is used to select the vortical flow test problem. The
        purpose of this test problem is to test velocity errors generated by spatial
        operators in the presence of 3D vorticity and in particluar the
        superposition of planar and vortical flows, analogous to voritcity
        stretching. For more details, see Waltz,
        et. al, "Manufactured solutions for the three-dimensional Euler equations
        with relevance to Inertial Confinement Fusion", Journal of Computational
        Physics 267 (2014) 196-209.)"});

      keywords.insert({"nl_energy_growth",
        "Select the nonlinear energy growth test problem",
        R"(This keyword is used to select the nonlinear energy growth test problem.
        The purpose of this test problem is to test nonlinear, time dependent energy
        growth and the subsequent development of pressure gradients due to coupling
        between the internal energy and the equation of state. For more details,
        see Waltz, et. al, "Manufactured
        solutions for the three-dimensional Euler equations with relevance to
        Inertial Confinement Fusion", Journal of Computational Physics 267 (2014)
        196-209.)"});

      keywords.insert({"rayleigh_taylor",
        "Select the Rayleigh-Taylor test problem ",
        R"(This keyword is used to select the Rayleigh-Taylor unstable configuration
        test problem. The purpose of this test problem is to assess time dependent
        fluid motion in the presence of Rayleigh-Taylor unstable conditions, i.e.
        opposing density and pressure gradients.
        For more details, see Waltz, et. al, "Manufactured solutions for the
        three-dimensional Euler equations with relevance to Inertial Confinement
        Fusion", Journal of Computational Physics 267 (2014) 196-209.)"});

      keywords.insert({"taylor_green",
        "Select the Taylor-Green test problem ",
        R"(This keyword is used to select the Taylor-Green vortex test problem. The
        purpose of this problem is to test time accuracy and the correctness of the
        discretization of the viscous term in the Navier-Stokes equation. For more
        details on the flow, see G.I. Taylor, A.E.
        Green, "Mechanism of the Production of Small Eddies from Large Ones", Proc.
        R. Soc. Lond. A 1937 158 499-521; DOI: 10.1098/rspa.1937.0036. Published 3
        February 1937.)"});

      keywords.insert({"sod_shocktube",
        "Select the Sod shock-tube test problem ",
        R"(This keyword is used to select the Sod shock-tube test problem. The
        purpose of this test problem is to test the correctness of the
        approximate Riemann solver and its shock and interface capturing
        capabilities. For more details, see
        G. A. Sod, "A Survey of Several Finite Difference Methods for Systems of
        Nonlinear Hyperbolic Conservation Laws", J. Comput. Phys., 27 (1978)
        1–31.)"});

      keywords.insert({"rotated_sod_shocktube",
        "Select the rotated Sod shock-tube test problem ",
        R"(This keyword is used to select the rotated Sod shock-tube test problem.
        This the same as Sod shocktube but the geometry is rotated about X, Y, Z
        each by 45 degrees (in that order) so that none of the domain boundary align
        with any of the coordinate directions. The purpose of this test problem is
        to test the correctness of the approximate Riemann solver and its shock and
        interface capturing capabilities in an arbitrarily oriented geometry.
        For more details on the Sod
        problem, see G. A. Sod, "A Survey of Several Finite Difference Methods for
        Systems of Nonlinear Hyperbolic Conservation Laws", J. Comput. Phys., 27
        (1978) 1–31.)"});

      keywords.insert({"shedding_flow",
        "Select the Shedding flow test problem ",
        R"(This keyword is used to select the Shedding flow test problem. It
        describe a quasi-2D inviscid flow over a triangular wedge in tetrahedron
        grid. The purpose of this test problem is to test the capability of DG
        scheme for retaining the shape of vortices and also different error
        indicator behavior for this external flow problem when p-adaptive DG scheme
        is applied.)"});

      keywords.insert({"sedov_blastwave",
        "Select the Sedov blast-wave test problem ",
        R"(This keyword is used to select the Sedov blast-wave test problem. The
        purpose of this test problem is to test the correctness of the
        approximate Riemann solver and its strong shock and interface capturing
        capabilities.)"});

      keywords.insert({"interface_advection",
        "Select the interface advection test problem ",
        R"(This keyword is used to select the interface advection test problem.
        The purpose of this test problem is to test the well-balancedness of the
        multi-material discretization and its interface capturing
        capabilities.)"});

      keywords.insert({"gauss_hump_compflow",
        "Select advection of 2D Gaussian hump test problem",
        R"(This keyword is used to select the advection of 2D Gaussian hump test
        problem. The initial and boundary conditions are specified to set up the
        test problem suitable to exercise and test the advection terms of the
        Euler equations. The baseline of the density distribution in this testcase
        is 1 instead of 0 in gauss_hump_transport which enables it to be the
        regression testcase for p-adaptive DG scheme.)"});

      keywords.insert({"waterair_shocktube",
        "Select the water-air shock-tube test problem ",
        R"(This keyword is used to select the Water-air shock-tube test problem.
        The purpose of this test problem is to test the correctness of the
        multi-material pressure relaxation procedure and its interface capturing
        capabilities. For more details, see
        Chiapolino, A., Saurel, R., & Nkonga, B. (2017). Sharpening diffuse
        interfaces with compressible fluids on unstructured meshes. Journal of
        Computational Physics, 340, 389-417.)"});

      keywords.insert({"shock_hebubble",
        "Select the shock He-bubble test problem ",
        R"(This keyword is used to select the shock He-bubble test problem. The
        purpose of this test problem is to test the correctness of the
        multi-material algorithm and its shock-interface interaction
        capabilities. For more details, see
        Quirk, J. J., & Karni, S. (1996). On the dynamics of a shock–bubble
        interaction. Journal of Fluid Mechanics, 318, 129-163.)"});

      keywords.insert({"underwater_ex",
        "Select the underwater explosion test problem ",
        R"(This keyword is used to select the underwater explosion test problem.
        The purpose of this test problem is to test the correctness of the
        multi-material algorithm and its interface capturing capabilities in the
        presence of strong shocks and large deformations.
        For more details, see
        Chiapolino, A., Saurel, R., & Nkonga, B. (2017). Sharpening diffuse
        interfaces with compressible fluids on unstructured meshes. Journal of
        Computational Physics, 340, 389-417.)"});

      keywords.insert({"shockdensity_wave",
        "Select the shock-density wave test problem ",
        R"(This keyword is used to select the shock-density wave test problem.
        THe purpose of this test problem is to assess the accuracy of high order
        method in predicting the interaction of a density wave with a shock front.
        For more details, see Yu, L., Matthias
        I. (2014). Discontinuous Galerkin method for multicomponent chemically
        reacting flows and combustion. Journal of Computational Physics, 270,
        105-137.)"});

      keywords.insert({"equilinterface_advect",
        "Select the advection of equilibrium interface problem ",
        R"(This keyword is used to select the advection of equilibrium interface
        problem. This is a manufactured problem with source terms with nonlinear
        solutions near the material interface. Source terms are used to ensure
        that the conservation laws are satisfied by the manufactured solution.)"});

      keywords.insert({"sinewave_packet",
        "Select the advection of sinewave packet problem ",
        R"(This keyword is used to select the advection of sinewave packet
        problem.)"});

      keywords.insert({"richtmyer_meshkov",
        "Select the Richtmyer-Meshkov instability problem ",
        R"(This keyword is used to select the Richtmyer-Meshkov instability
        problem. In this problem, a shock hits a perturbed material interface.)"});

      // -----------------------------------------------------------------------

      // Initialize help: fill own keywords
      tk::ctr::Info ctrinfoFill(get< tag::cmd, tag::ctrinfo >());
      for (const auto& i : keywords) {
        ctrinfoFill.fill(i);
      }
    }

    //! Query scheme centering
    //! \return Scheme centering
    tk::Centering centering() const
    { return ctr::Scheme().centering( get< tag::scheme >() ); }

    /** @name Pack/Unpack: Serialize InputDeck object for Charm++ */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er& p ) { tk::TaggedTuple< ConfigMembers >::pup(p); }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] c InputDeck object reference
    friend void operator|( PUP::er& p, InputDeck& c ) { c.pup(p); }
    //@}

};

} // ctr::
} // inciter::

#endif // InputDeck_h
