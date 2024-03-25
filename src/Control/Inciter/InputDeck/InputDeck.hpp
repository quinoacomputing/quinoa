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
#include "SimpTaggedTuple.hpp"
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

using bclist = tk::SimpTaggedTuple< brigand::list<
  tag::dirichlet,   std::vector< std::size_t >,
  tag::symmetry,    std::vector< std::size_t >,
  tag::inlet,       std::vector< std::size_t >,
  tag::outlet,      std::vector< std::size_t >,
  tag::farfield,    std::vector< std::size_t >,
  tag::extrapolate, std::vector< std::size_t >
> >;

using newbox = tk::SimpTaggedTuple< brigand::list<
  tag::materialid,     std::size_t,
  tag::volume,         tk::real,
  tag::mass,           tk::real,
  tag::density,        tk::real,
  tag::velocity,       std::vector< tk::real >,
  tag::pressure,       tk::real,
  tag::energy,         tk::real,
  tag::energy_content, tk::real,
  tag::temperature,    tk::real,
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

using newmeshblock = tk::SimpTaggedTuple< brigand::list<
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
  tag::initiate,       inciter::ctr::InitiateType,
  tag::point,          std::vector< tk::real >,
  tag::init_time,      tk::real,
  tag::front_width,    tk::real,
  tag::front_speed,    tk::real
> >;

using ConfigMembers = brigand::list<

  tag::title, std::string,

  // Command line parameters
  tag::cmd, CmdLine,

  // Control file help object
  tag::ctrinfo, tk::ctr::HelpFactory,

  // time stepping options
  // ---------------------------------------------------------------------------
  tag::nstep, uint64_t,
  tag::term,  tk::real,
  tag::t0,    tk::real,
  tag::dt,    tk::real,
  tag::cfl,   tk::real,
  tag::ttyi,  uint32_t,

  // steady-state solver options
  // ---------------------------------------------------------------------------
  tag::steady_state, bool,
  tag::residual,     tk::real,
  tag::rescomp,      uint32_t,

  // mesh partitioning and reordering/sorting choices
  // ---------------------------------------------------------------------------
  tag::partitioning,     tk::ctr::PartitioningAlgorithmType,
  tag::pelocal_reorder,  bool,
  tag::operator_reorder, bool,

  // discretization scheme choices
  // ---------------------------------------------------------------------------
  tag::scheme, SchemeType,
  tag::ndof,   std::size_t,
  tag::rdof,   std::size_t,
  tag::flux,   FluxType,

  // limiter options
  // ---------------------------------------------------------------------------
  tag::limiter,              LimiterType,
  tag::cweight,              tk::real,
  tag::shock_detector_coeff, tk::real,
  tag::accuracy_test,        bool,
  tag::limsol_projection,    bool,
  tag::fct,                  bool,
  tag::fctclip,              bool,
  tag::fcteps,               tk::real,
  tag::ctau,                 tk::real,
  tag::sysfct,               bool,
  tag::sysfctvar,            std::vector< std::size_t >,

  // PDE options
  // ---------------------------------------------------------------------------
  tag::ncomp, std::size_t,
  tag::pde,   PDEType,

  // Transport
  // ---------------------------------------------------------------------------
  tag::transport, tk::SimpTaggedTuple< brigand::list<
    tag::ncomp,          std::size_t,
    tag::intsharp,       int,
    tag::intsharp_param, tk::real,
    tag::problem,        ProblemType,
    tag::diffusivity,    std::vector< tk::real >,
    tag::lambda,         std::vector< tk::real >,
    tag::u0,             std::vector< tk::real >
  > >,

  // CompFlow
  // ---------------------------------------------------------------------------
  tag::compflow, tk::SimpTaggedTuple< brigand::list<
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
  > >,

  // MultiMat
  // ---------------------------------------------------------------------------
  tag::multimat, tk::SimpTaggedTuple< brigand::list<
    tag::nmat,             std::size_t,
    tag::prelax,           uint64_t,
    tag::prelax_timescale, tk::real,
    tag::intsharp,         int,
    tag::intsharp_param,   tk::real,
    tag::problem,          ProblemType
  > >,

  // Dependent variable name
  tag::depvar, std::vector< char >,

  tag::sys, std::map< std::size_t, std::size_t >,

  // physics choices
  tag::physics, PhysicsType,

  // Material/EOS object
  // ---------------------------------------------------------------------------
  tag::material, std::vector<
    tk::SimpTaggedTuple< brigand::list<
      tag::eos,      MaterialType,
      tag::id,       std::vector< uint64_t >,
      tag::gamma,    std::vector< tk::real >,
      tag::pstiff,   std::vector< tk::real >,
      tag::w_gru,    std::vector< tk::real >,
      tag::A_jwl,    std::vector< tk::real >,
      tag::B_jwl,    std::vector< tk::real >,
      tag::C_jwl,    std::vector< tk::real >,
      tag::R1_jwl,   std::vector< tk::real >,
      tag::R2_jwl,   std::vector< tk::real >,
      tag::rho0_jwl, std::vector< tk::real >,
      tag::de_jwl,   std::vector< tk::real >,
      tag::rhor_jwl, std::vector< tk::real >,
      tag::Tr_jwl,   std::vector< tk::real >,
      tag::Pr_jwl,   std::vector< tk::real >,
      tag::mu,       std::vector< tk::real >,
      tag::cv,       std::vector< tk::real >,
      tag::k,        std::vector< tk::real >
    > >
  >,

  tag::matidxmap, tk::SimpTaggedTuple< brigand::list<
    tag::eosidx, std::vector< std::size_t >,
    tag::matidx, std::vector< std::size_t >,
    tag::solidx, std::vector< std::size_t >
  > >,

  // Field output block
  // ---------------------------------------------------------------------------
  tag::field_output, tk::SimpTaggedTuple< brigand::list<
    tag::iter_interval, uint32_t,
    tag::time_interval, tk::real,
    tag::time_range,    std::vector< tk::real >,
    tag::refined,       bool,
    tag::filetype,      tk::ctr::FieldFileType,
    tag::sideset,       std::vector< uint64_t >,
    tag::outvar,        std::vector< OutVar >
  > >,

  // Diagnostics block
  // ---------------------------------------------------------------------------
  tag::diagnostics, tk::SimpTaggedTuple< brigand::list<
    tag::iter_interval, uint32_t,
    tag::error,         tk::ctr::ErrorType,
    tag::format,        tk::ctr::TxtFloatFormatType,
    tag::precision,     uint32_t
  > >,

  // History output block
  // ---------------------------------------------------------------------------
  tag::history_output, tk::SimpTaggedTuple< brigand::list<
    tag::iter_interval, uint32_t,
    tag::time_interval, tk::real,
    tag::time_range,    std::vector< tk::real >,
    tag::format,        tk::ctr::TxtFloatFormatType,
    tag::precision,     uint32_t,
    tag::point,         std::vector<
      tk::SimpTaggedTuple< brigand::list<
        tag::id,    std::string,
        tag::coord, std::vector< tk::real >
      > >
    >
  > >,

  // ALE block
  // ---------------------------------------------------------------------------
  tag::ale, tk::SimpTaggedTuple< brigand::list<
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
      tk::SimpTaggedTuple< brigand::list<
        tag::sideset, std::vector< uint64_t >,
        tag::fntype,  tk::ctr::UserTableType,
        tag::fn,      std::vector< tk::real >
      > >
    >
  > >,

  // p-refinement block
  // ---------------------------------------------------------------------------
  tag::pref, tk::SimpTaggedTuple< brigand::list<
    tag::pref,      bool,
    tag::indicator, PrefIndicatorType,
    tag::ndofmax,   std::size_t,
    tag::tolref,    tk::real
  > >,

  // AMR block
  // ---------------------------------------------------------------------------
  tag::amr, tk::SimpTaggedTuple< brigand::list<
    tag::amr,           bool,
    tag::t0ref,         bool,
    tag::dtref,         bool,
    tag::dtref_uniform, bool,
    tag::dtfreq,        std::size_t,
    tag::maxlevels,     std::size_t,
    tag::initial,       std::vector< AMRInitialType >,
    tag::edgelist,      std::vector< std::size_t >,
    tag::coords,        tk::SimpTaggedTuple< brigand::list<
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

  // Boundary conditions block
  // ---------------------------------------------------------------------------
  tag::bc, std::vector<
    tk::SimpTaggedTuple< brigand::list<
      tag::mesh,        std::vector< std::size_t >,
      tag::dirichlet,   std::vector< std::size_t >,
      tag::symmetry,    std::vector< std::size_t >,
      tag::inlet,       std::vector< std::size_t >,
      tag::outlet,      std::vector< std::size_t >,
      tag::farfield,    std::vector< std::size_t >,
      tag::extrapolate, std::vector< std::size_t >,
      tag::stag_point,  std::vector< tk::real >,
      tag::radius,      tk::real,
      tag::velocity,    std::vector< tk::real >,
      tag::pressure,    tk::real,
      tag::density,     tk::real,
      tag::sponge,      tk::SimpTaggedTuple< brigand::list<
        tag::sideset,     std::vector< uint64_t >,
        tag::vparam,      std::vector< tk::real >,
        tag::pparam,      tk::real
      > >,
      tag::timedep,     std::vector<
        tk::SimpTaggedTuple< brigand::list<
          tag::sideset,   std::vector< uint64_t >,
          tag::fn,        std::vector< tk::real >
        > >
      >
    > >
  >,

  // Initial conditions (ic) block
  // ---------------------------------------------------------------------------
  tag::ic, tk::SimpTaggedTuple< brigand::list<
    tag::materialid,  std::size_t,
    tag::pressure,    tk::real,
    tag::temperature, tk::real,
    tag::density,     tk::real,
    tag::energy,      tk::real,
    tag::velocity,    std::vector< tk::real >,
    tag::box,         std::vector< newbox >,
    tag::meshblock,   std::vector< newmeshblock >
  > >,

  // Overset mesh block
  // ---------------------------------------------------------------------------
  tag::mesh, std::vector<
    tk::SimpTaggedTuple< brigand::list<
      tag::filename,    std::string,
      tag::location,    std::vector< tk::real >,
      tag::orientation, std::vector< tk::real >,
      tag::velocity,    std::vector< tk::real >
    > >
  >,

  tag::transfer, std::vector< Transfer >
>;

// Class storing the Config params
class InputDeck : public tk::SimpTaggedTuple< ConfigMembers > {

  public:

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
        std::numeric_limits< kw::amr_xminus::info::expect::type >::max() / 100;
      get< tag::amr, tag::coords, tag::xminus >() = rmax;
      get< tag::amr, tag::coords, tag::xplus >() = -rmax;
      get< tag::amr, tag::coords, tag::yminus >() = rmax;
      get< tag::amr, tag::coords, tag::yplus >() = -rmax;
      get< tag::amr, tag::coords, tag::zminus >() = rmax;
      get< tag::amr, tag::coords, tag::zplus >() = -rmax;
      //// Initialize help: fill own keywords
      //const auto& ctrinfoFill = tk::ctr::Info( get< tag::ctrinfo >() );
      //brigand::for_each< keywords >( ctrinfoFill );
    }

    //! Extract list of mesh filenames (each assigned to a solver)
    std::vector< std::string > mesh() const {
      std::vector< std::string > meshes;
      const auto& mdeck = this->get< tag::mesh >();
      for (const auto& im : mdeck) {
        meshes.push_back(im.get< tag::filename >());
      }
      return meshes;
    }

    //! Query scheme centering
    //! \return Scheme centering
    tk::Centering centering() const
    { return ctr::Scheme().centering( get< tag::scheme >() ); }

    /** @name Pack/Unpack: Serialize InputDeck object for Charm++ */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er& p ) { tk::SimpTaggedTuple< ConfigMembers >::pup(p); }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] c InputDeck object reference
    friend void operator|( PUP::er& p, InputDeck& c ) { c.pup(p); }
    //@}

};

} // ctr::
} // inciter::

#endif // InputDeck_h
