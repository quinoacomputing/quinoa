// *****************************************************************************
/*!
  \file      src/Control/Inciter/InputDeck/New2InputDeck.hpp
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
#ifndef New2InputDeck_h
#define New2InputDeck_h

#include <getopt.h>
#include "SimpTaggedTuple.hpp"
#include "Inciter/CmdLine/CmdLine.hpp"
#include "Inciter/Options/PDE.hpp"
#include "Inciter/Options/Scheme.hpp"
#include "Transfer.hpp"
#include "Options/PartitioningAlgorithm.hpp"
#include "Inciter/NewOutVar.hpp"

namespace newtag {
  DEFTAG(cmd);
  DEFTAG(title);
  DEFTAG(nstep);
  DEFTAG(term);
  DEFTAG(t0);
  DEFTAG(dt);
  DEFTAG(cfl);
  DEFTAG(ttyi);
  DEFTAG(steady_state);
  DEFTAG(residual);
  DEFTAG(rescomp);
  DEFTAG(partitioning);
  DEFTAG(pelocal_reorder);
  DEFTAG(operator_reorder);
  DEFTAG(scheme);
  DEFTAG(ndof);
  DEFTAG(rdof);
  DEFTAG(flux);
  DEFTAG(limiter);
  DEFTAG(cweight);
  DEFTAG(shock_detector_coeff);
  DEFTAG(accuracy_test);
  DEFTAG(limsol_projection);
  DEFTAG(fct);
  DEFTAG(fctclip);
  DEFTAG(fcteps);
  DEFTAG(ctau);
  DEFTAG(sysfct);
  DEFTAG(sysfctvar);

  DEFTAG(ncomp);
  DEFTAG(pde);
  DEFTAG(problem);
  DEFTAG(transport);
  DEFTAG(compflow);
  DEFTAG(multimat);
  DEFTAG(diffusivity);
  DEFTAG(lambda);
  DEFTAG(u0);
  DEFTAG(alpha);
  DEFTAG(beta);
  DEFTAG(betax);
  DEFTAG(betay);
  DEFTAG(betaz);
  DEFTAG(r0);
  DEFTAG(p0);
  DEFTAG(ce);
  DEFTAG(kappa);
  DEFTAG(nmat);
  DEFTAG(prelax);
  DEFTAG(prelax_timescale);
  DEFTAG(intsharp);
  DEFTAG(intsharp_param);
  DEFTAG(depvar);
  DEFTAG(sys);
  DEFTAG(physics);

  DEFTAG(material);
  DEFTAG(eos);
  DEFTAG(id);
  DEFTAG(gamma);
  DEFTAG(pstiff);
  DEFTAG(w_gru);
  DEFTAG(A_jwl);
  DEFTAG(B_jwl);
  DEFTAG(C_jwl);
  DEFTAG(R1_jwl);
  DEFTAG(R2_jwl);
  DEFTAG(rho0_jwl);
  DEFTAG(de_jwl);
  DEFTAG(rhor_jwl);
  DEFTAG(Tr_jwl);
  DEFTAG(Pr_jwl);
  DEFTAG(mu);
  DEFTAG(cv);
  DEFTAG(k);
  DEFTAG(matidxmap);
  DEFTAG(eosidx);
  DEFTAG(matidx);
  DEFTAG(solidx);

  DEFTAG(field_output);
  DEFTAG(iter_interval);
  DEFTAG(time_interval);
  DEFTAG(time_range);
  DEFTAG(refined);
  DEFTAG(filetype);
  DEFTAG(sideset);
  DEFTAG(outvar);
  DEFTAG(diagnostics);
  DEFTAG(error);
  DEFTAG(format);
  DEFTAG(precision);
  DEFTAG(history_output);
  DEFTAG(point);
  DEFTAG(coord);

  DEFTAG(ale);
  DEFTAG(smoother);
  DEFTAG(mesh_velocity);
  DEFTAG(mesh_motion);
  DEFTAG(meshforce);
  DEFTAG(dvcfl);
  DEFTAG(vortmult);
  DEFTAG(maxit);
  DEFTAG(tolerance);
  DEFTAG(move);
  DEFTAG(fntype);
  DEFTAG(fn);

  DEFTAG(amr);
  DEFTAG(t0ref);
  DEFTAG(dtref);
  DEFTAG(dtref_uniform);
  DEFTAG(dtfreq);
  DEFTAG(maxlevels);
  DEFTAG(initial);
  DEFTAG(edgelist);
  DEFTAG(coords);
  DEFTAG(xminus);
  DEFTAG(xplus);
  DEFTAG(yminus);
  DEFTAG(yplus);
  DEFTAG(zminus);
  DEFTAG(zplus);
  DEFTAG(refvar);
  DEFTAG(tol_refine);
  DEFTAG(tol_derefine);
  DEFTAG(pref);
  DEFTAG(indicator);
  DEFTAG(ndofmax);
  DEFTAG(tolref);

  DEFTAG(bc);
  DEFTAG(mesh);
  DEFTAG(dirichlet);
  DEFTAG(symmetry);
  DEFTAG(inlet);
  DEFTAG(outlet);
  DEFTAG(farfield);
  DEFTAG(extrapolate);
  DEFTAG(stag_point);
  DEFTAG(sponge);
  DEFTAG(vparam);
  DEFTAG(pparam);
  DEFTAG(radius);
  DEFTAG(timedep);

  DEFTAG(ic);
  DEFTAG(materialid);
  DEFTAG(pressure);
  DEFTAG(temperature);
  DEFTAG(velocity);
  DEFTAG(box);
  DEFTAG(meshblock);
  DEFTAG(blockid);
  DEFTAG(volume);
  DEFTAG(mass);
  DEFTAG(density);
  DEFTAG(energy);
  DEFTAG(energy_content);
  DEFTAG(xmin);
  DEFTAG(xmax);
  DEFTAG(ymin);
  DEFTAG(ymax);
  DEFTAG(zmin);
  DEFTAG(zmax);
  DEFTAG(orientation);
  DEFTAG(initiate);
  DEFTAG(init_time);
  DEFTAG(front_width);
  DEFTAG(front_speed);

  DEFTAG(filename);
  DEFTAG(location);
  DEFTAG(transfer);

  using bclist = tk::TaggedTuple< brigand::list<
    dirichlet,   std::vector< std::size_t >,
    symmetry,    std::vector< std::size_t >,
    inlet,       std::vector< std::size_t >,
    outlet,      std::vector< std::size_t >,
    farfield,    std::vector< std::size_t >,
    extrapolate, std::vector< std::size_t >
  > >;

  using newbox = tk::SimpTaggedTuple< brigand::list<
    newtag::materialid,     std::size_t,
    newtag::volume,         tk::real,
    newtag::mass,           tk::real,
    newtag::density,        tk::real,
    newtag::velocity,       std::vector< tk::real >,
    newtag::pressure,       tk::real,
    newtag::energy,         tk::real,
    newtag::energy_content, tk::real,
    newtag::temperature,    tk::real,
    newtag::xmin,           tk::real,
    newtag::xmax,           tk::real,
    newtag::ymin,           tk::real,
    newtag::ymax,           tk::real,
    newtag::zmin,           tk::real,
    newtag::zmax,           tk::real,
    newtag::orientation,    std::vector< tk::real >,
    newtag::initiate,       inciter::ctr::InitiateType,
    newtag::point,          std::vector< tk::real >,
    newtag::init_time,      tk::real,
    newtag::front_width,    tk::real,
    newtag::front_speed,    tk::real
  > >;

  using newmeshblock = tk::SimpTaggedTuple< brigand::list<
    newtag::blockid,        std::uint64_t,
    newtag::materialid,     std::size_t,
    newtag::volume,         tk::real,
    newtag::mass,           tk::real,
    newtag::density,        tk::real,
    newtag::velocity,       std::vector< tk::real >,
    newtag::pressure,       tk::real,
    newtag::energy,         tk::real,
    newtag::energy_content, tk::real,
    newtag::temperature,    tk::real,
    newtag::initiate,       inciter::ctr::InitiateType,
    newtag::point,          std::vector< tk::real >,
    newtag::init_time,      tk::real,
    newtag::front_width,    tk::real,
    newtag::front_speed,    tk::real
  > >;
} // newtag::

namespace inciter {

namespace ctr {

using ncomp_t = std::size_t;

using ConfigMembers = brigand::list<

  newtag::title, std::string,

  // Command line parameters
  newtag::cmd, CmdLine,

  // time stepping options
  // ---------------------------------------------------------------------------
  newtag::nstep, uint64_t,
  newtag::term,  tk::real,
  newtag::t0,    tk::real,
  newtag::dt,    tk::real,
  newtag::cfl,   tk::real,
  newtag::ttyi,  uint32_t,

  // steady-state solver options
  // ---------------------------------------------------------------------------
  newtag::steady_state, bool,
  newtag::residual,     tk::real,
  newtag::rescomp,      uint32_t,

  // mesh partitioning and reordering/sorting choices
  // ---------------------------------------------------------------------------
  newtag::partitioning,     tk::ctr::PartitioningAlgorithmType,
  newtag::pelocal_reorder,  bool,
  newtag::operator_reorder, bool,

  // discretization scheme choices
  // ---------------------------------------------------------------------------
  newtag::scheme, SchemeType,
  newtag::ndof,   std::size_t,
  newtag::rdof,   std::size_t,
  newtag::flux,   FluxType,

  // limiter options
  // ---------------------------------------------------------------------------
  newtag::limiter,              LimiterType,
  newtag::cweight,              tk::real,
  newtag::shock_detector_coeff, tk::real,
  newtag::accuracy_test,        bool,
  newtag::limsol_projection,    bool,
  newtag::fct,                  bool,
  newtag::fctclip,              bool,
  newtag::fcteps,               tk::real,
  newtag::ctau,                 tk::real,
  newtag::sysfct,               bool,
  newtag::sysfctvar,            std::vector< std::size_t >,

  // PDE options
  // ---------------------------------------------------------------------------
  newtag::ncomp, std::size_t,
  newtag::pde,   PDEType,

  // Transport
  // ---------------------------------------------------------------------------
  newtag::transport, tk::SimpTaggedTuple< brigand::list<
    newtag::ncomp,          std::size_t,
    newtag::intsharp,       int,
    newtag::intsharp_param, tk::real,
    newtag::problem,        ProblemType,
    newtag::diffusivity,    std::vector< tk::real >,
    newtag::lambda,         std::vector< tk::real >,
    newtag::u0,             std::vector< tk::real >
  > >,

  // CompFlow
  // ---------------------------------------------------------------------------
  newtag::compflow, tk::SimpTaggedTuple< brigand::list<
    newtag::problem, ProblemType,
    newtag::alpha,   tk::real,
    newtag::beta,    tk::real,
    newtag::betax,   tk::real,
    newtag::betay,   tk::real,
    newtag::betaz,   tk::real,
    newtag::r0,      tk::real,
    newtag::p0,      tk::real,
    newtag::ce,      tk::real,
    newtag::kappa,   tk::real
  > >,

  // MultiMat
  // ---------------------------------------------------------------------------
  newtag::multimat, tk::SimpTaggedTuple< brigand::list<
    newtag::nmat,             std::size_t,
    newtag::prelax,           uint64_t,
    newtag::prelax_timescale, tk::real,
    newtag::intsharp,         int,
    newtag::intsharp_param,   tk::real,
    newtag::problem,          ProblemType
  > >,

  // Dependent variable name
  newtag::depvar, std::vector< char >,

  newtag::sys, std::map< std::size_t, std::size_t >,

  // physics choices
  newtag::physics, PhysicsType,

  // Material/EOS object
  // ---------------------------------------------------------------------------
  newtag::material, std::vector<
    tk::SimpTaggedTuple< brigand::list<
      newtag::eos,      MaterialType,
      newtag::id,       std::vector< uint64_t >,
      newtag::gamma,    std::vector< tk::real >,
      newtag::pstiff,   std::vector< tk::real >,
      newtag::w_gru,    std::vector< tk::real >,
      newtag::A_jwl,    std::vector< tk::real >,
      newtag::B_jwl,    std::vector< tk::real >,
      newtag::C_jwl,    std::vector< tk::real >,
      newtag::R1_jwl,   std::vector< tk::real >,
      newtag::R2_jwl,   std::vector< tk::real >,
      newtag::rho0_jwl, std::vector< tk::real >,
      newtag::de_jwl,   std::vector< tk::real >,
      newtag::rhor_jwl, std::vector< tk::real >,
      newtag::Tr_jwl,   std::vector< tk::real >,
      newtag::Pr_jwl,   std::vector< tk::real >,
      newtag::mu,       std::vector< tk::real >,
      newtag::cv,       std::vector< tk::real >,
      newtag::k,        std::vector< tk::real >
    > >
  >,

  newtag::matidxmap, tk::SimpTaggedTuple< brigand::list<
    newtag::eosidx, std::vector< std::size_t >,
    newtag::matidx, std::vector< std::size_t >,
    newtag::solidx, std::vector< std::size_t >
  > >,

  // Field output block
  // ---------------------------------------------------------------------------
  newtag::field_output, tk::SimpTaggedTuple< brigand::list<
    newtag::iter_interval, uint32_t,
    newtag::time_interval, tk::real,
    newtag::time_range,    std::vector< tk::real >,
    newtag::refined,       bool,
    newtag::filetype,      tk::ctr::FieldFileType,
    newtag::sideset,       std::vector< uint64_t >,
    newtag::outvar,        std::vector< NewOutVar >
  > >,

  // Diagnostics block
  // ---------------------------------------------------------------------------
  newtag::diagnostics, tk::SimpTaggedTuple< brigand::list<
    newtag::iter_interval, uint32_t,
    newtag::error,         tk::ctr::ErrorType,
    newtag::format,        tk::ctr::TxtFloatFormatType,
    newtag::precision,     uint32_t
  > >,

  // History output block
  // ---------------------------------------------------------------------------
  newtag::history_output, tk::SimpTaggedTuple< brigand::list<
    newtag::iter_interval, uint32_t,
    newtag::time_interval, tk::real,
    newtag::time_range,    std::vector< tk::real >,
    newtag::format,        tk::ctr::TxtFloatFormatType,
    newtag::precision,     uint32_t,
    newtag::point,         std::vector<
      tk::SimpTaggedTuple< brigand::list<
        newtag::id,    std::string,
        newtag::coord, std::vector< tk::real >
      > >
    >
  > >,

  // ALE block
  // ---------------------------------------------------------------------------
  newtag::ale, tk::SimpTaggedTuple< brigand::list<
    newtag::ale,           bool,
    newtag::smoother,      MeshVelocitySmootherType,
    newtag::mesh_velocity, MeshVelocityType,
    newtag::mesh_motion,   std::vector< std::size_t >,
    newtag::meshforce,     std::vector< tk::real >,
    newtag::dvcfl,         tk::real,
    newtag::vortmult,      tk::real,
    newtag::maxit,         std::size_t,
    newtag::tolerance,     tk::real,
    newtag::dirichlet,     std::vector< std::size_t >,
    newtag::symmetry,      std::vector< std::size_t >,
    newtag::move,          std::vector<
      tk::SimpTaggedTuple< brigand::list<
        newtag::sideset, std::vector< uint64_t >,
        newtag::fntype,  tk::ctr::UserTableType,
        newtag::fn,      std::vector< tk::real >
      > >
    >
  > >,

  // p-refinement block
  // ---------------------------------------------------------------------------
  newtag::pref, tk::SimpTaggedTuple< brigand::list<
    newtag::pref,      bool,
    newtag::indicator, PrefIndicatorType,
    newtag::ndofmax,   std::size_t,
    newtag::tolref,    tk::real
  > >,

  // AMR block
  // ---------------------------------------------------------------------------
  newtag::amr, tk::SimpTaggedTuple< brigand::list<
    newtag::amr,           bool,
    newtag::t0ref,         bool,
    newtag::dtref,         bool,
    newtag::dtref_uniform, bool,
    newtag::dtfreq,        std::size_t,
    newtag::maxlevels,     std::size_t,
    newtag::initial,       std::vector< AMRInitialType >,
    newtag::edgelist,      std::vector< std::size_t >,
    newtag::coords,        tk::SimpTaggedTuple< brigand::list<
      newtag::xminus, tk::real,
      newtag::xplus,  tk::real,
      newtag::yminus, tk::real,
      newtag::yplus,  tk::real,
      newtag::zminus, tk::real,
      newtag::zplus,  tk::real
    > >,
    newtag::error,         AMRErrorType,
    newtag::refvar,        std::vector< char >,
    newtag::tol_refine,    tk::real,
    newtag::tol_derefine,  tk::real
  > >,

  // Boundary conditions block
  // ---------------------------------------------------------------------------
  newtag::bc, std::vector<
    tk::SimpTaggedTuple< brigand::list<
      newtag::mesh,        std::vector< std::size_t >,
      newtag::dirichlet,   std::vector< std::size_t >,
      newtag::symmetry,    std::vector< std::size_t >,
      newtag::inlet,       std::vector< std::size_t >,
      newtag::outlet,      std::vector< std::size_t >,
      newtag::farfield,    std::vector< std::size_t >,
      newtag::extrapolate, std::vector< std::size_t >,
      newtag::stag_point,  std::vector< tk::real >,
      newtag::radius,      tk::real,
      newtag::velocity,    std::vector< tk::real >,
      newtag::pressure,    tk::real,
      newtag::density,     tk::real,
      newtag::sponge,      tk::SimpTaggedTuple< brigand::list<
        newtag::sideset,     std::vector< uint64_t >,
        newtag::vparam,      std::vector< tk::real >,
        newtag::pparam,      tk::real
      > >,
      newtag::timedep,     std::vector<
        tk::SimpTaggedTuple< brigand::list<
          newtag::sideset,   std::vector< uint64_t >,
          newtag::fn,        std::vector< tk::real >
        > >
      >
    > >
  >,

  // Initial conditions (ic) block
  // ---------------------------------------------------------------------------
  newtag::ic, tk::SimpTaggedTuple< brigand::list<
    newtag::materialid,  std::size_t,
    newtag::pressure,    tk::real,
    newtag::temperature, tk::real,
    newtag::density,     tk::real,
    newtag::energy,      tk::real,
    newtag::velocity,    std::vector< tk::real >,
    newtag::box,         std::vector< newtag::newbox >,
    newtag::meshblock,   std::vector< newtag::newmeshblock >
  > >,

  // Overset mesh block
  // ---------------------------------------------------------------------------
  newtag::mesh, std::vector<
    tk::SimpTaggedTuple< brigand::list<
      newtag::filename,    std::string,
      newtag::location,    std::vector< tk::real >,
      newtag::orientation, std::vector< tk::real >,
      newtag::velocity,    std::vector< tk::real >
    > >
  >,

  newtag::transfer, std::vector< Transfer >
>;

// Class storing the Config params
class New2InputDeck : public tk::SimpTaggedTuple< ConfigMembers > {

  public:

    //! \brief Constructor: set defaults
    //! \param[in] cl Previously parsed and store command line
    //! \details Anything not set here is initialized by the compiler using the
    //!   default constructor for the corresponding type.
    explicit New2InputDeck( const CmdLine& cl = {} ) {
      // Set previously parsed command line
      get< newtag::cmd >() = cl;
      // Default time stepping params
      get< newtag::dt >() = 0.0;
      get< newtag::cfl >() = 0.0;
      // Default AMR settings
      auto rmax =
        std::numeric_limits< kw::amr_xminus::info::expect::type >::max() / 100;
      get< newtag::amr, newtag::coords, newtag::xminus >() = rmax;
      get< newtag::amr, newtag::coords, newtag::xplus >() = -rmax;
      get< newtag::amr, newtag::coords, newtag::yminus >() = rmax;
      get< newtag::amr, newtag::coords, newtag::yplus >() = -rmax;
      get< newtag::amr, newtag::coords, newtag::zminus >() = rmax;
      get< newtag::amr, newtag::coords, newtag::zplus >() = -rmax;
    }

    //! Extract list of mesh filenames (each assigned to a solver)
    std::vector< std::string > mesh() const {
      std::vector< std::string > meshes;
      const auto& mdeck = this->get< newtag::mesh >();
      for (const auto& im : mdeck) {
        meshes.push_back(im.get< newtag::filename >());
      }
      return meshes;
    }

    //! Query scheme centering
    //! \return Scheme centering
    tk::Centering centering() const
    { return ctr::Scheme().centering( get< newtag::scheme >() ); }

    /** @name Pack/Unpack: Serialize New2InputDeck object for Charm++ */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er& p ) { tk::SimpTaggedTuple< ConfigMembers >::pup(p); }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] c New2InputDeck object reference
    friend void operator|( PUP::er& p, New2InputDeck& c ) { c.pup(p); }
    //@}

};

} // ctr::
} // inciter::

#endif // New2InputDeck_h
