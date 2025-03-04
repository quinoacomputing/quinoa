// *****************************************************************************
/*!
  \file      src/Control/Inciter/InputDeck/LuaParser.cpp
  \copyright 2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter's lua input deck file parser
  \details   This file defines the input deck, i.e., control file, parser for
    the computational shock hydrodynamics tool, Inciter.
*/
// *****************************************************************************

#include <ostream>
#include <type_traits>

#include "QuinoaConfig.hpp"

#include "NoWarning/pegtl.hpp"

#include "Print.hpp"
#include "Exception.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "Inciter/InputDeck/LuaParser.hpp"
#include "PDE/MultiMat/MultiMatIndexing.hpp"
#include "PDE/MultiSpecies/MultiSpeciesIndexing.hpp"

namespace tk {
namespace grm {

  //! \brief Case-insensitive character comparison functor
  struct CaseInsensitiveCharLess {
    //! Function call operator
    //! \param[in] lhs Left character of the comparitor operand
    //! \param[in] rhs Right character of the comparitor operand
    //! \return Boolean indicating the result of the comparison
    bool operator() ( char lhs, char rhs ) const {
      return std::tolower( lhs ) < std::tolower( rhs );
    }
  };

  //! \brief Parser-lifetime storage for dependent variables selected.
  //! \details Used to track the dependent variable of differential equations
  //!   (i.e., models) assigned during parsing. It needs to be case insensitive
  //!   since we only care about whether the variable is selected or not and not
  //!   whether it denotes a full variable (upper case) or a fluctuation (lower
  //!   case). This is true for both inserting variables into the set as well as
  //!   at matching terms of products in parsing requested statistics.
  static std::set< char, CaseInsensitiveCharLess > depvars;

} // grm::
} // tk::

using inciter::LuaParser;

LuaParser::LuaParser( const tk::Print& /*print*/,
                      const ctr::CmdLine& cmdline,
                      ctr::InputDeck& inputdeck ) :
  m_filename( cmdline.get< tag::io, tag::control >() )
// *****************************************************************************
//  Constructor
// //! \param[in] print Pretty printer
//! \param[in] cmdline Command line stack
//! \param[in,out] inputdeck Input deck stack where data is stored during
//!    parsing
// *****************************************************************************
{
  // Make sure there is a filename
  Assert( !m_filename.empty(), "No filename specified" );

  // Local file stream handle
  std::ifstream q;

  // Check if file exists, throw exception if it does not
  q.open( m_filename, std::ifstream::in );
  ErrChk( q.good(), "Failed to open file: " + m_filename );

  // Attempt to read a character, throw if it fails
  // It is curious that on some systems opening a directory instead of a file
  // with the above ifstream::open() call does not set the failbit. Thus we get
  // here fine, so we try to read a character from it. If it is a directory or
  // an empty file the read will fail, so we throw. Read more at: http://
  // stackoverflow.com/questions/9591036/
  // ifstream-open-doesnt-set-error-bits-when-argument-is-a-directory.
  q.get();
  ErrChk( q.good(), "Failed to read from file: " + m_filename );

  // Close it
  q.close();
  ErrChk( !q.fail(), "Failed to close file: " + m_filename );

  // Create InputDeck (a tagged tuple) to store parsed input
  ctr::InputDeck ideck( cmdline );

  // Read lua file into sol object
  sol::state lua_deck;
  lua_deck.open_libraries(sol::lib::base);
  lua_deck.script_file(m_filename);

  // Store inputdeck parameters from sol object into tagged tuple
  storeInputDeck(lua_deck["inciter"], ideck);

  inputdeck = std::move( ideck );
}

void
LuaParser::storeInputDeck(
  const sol::table& lua_ideck,
  ctr::InputDeck& gideck )
// *****************************************************************************
//  Store lua inputdeck in custom struct
//! \param[in] lua_ideck Lua inputdeck parsed by sol2
//! \param[in,out] gideck Inciter's inputdeck storage
// *****************************************************************************
{
  // TODO: explore replacing storeIfSpecd() and storeOptIfSpecd with sol::get_or()

  storeIfSpecd< std::string >(
    lua_ideck, "title", gideck.get< tag::title >(), "No title");

  // time stepping options
  // ---------------------------------------------------------------------------
  storeIfSpecd< uint64_t >(
    lua_ideck, "nstep", gideck.get< tag::nstep >(),
    std::numeric_limits< uint64_t >::max());
  storeIfSpecd< tk::real >(
    lua_ideck, "term", gideck.get< tag::term >(),
    std::numeric_limits< tk::real >::max());
  storeIfSpecd< tk::real >(
    lua_ideck, "t0", gideck.get< tag::t0 >(), 0.0);
  storeIfSpecd< tk::real >(
    lua_ideck, "dt", gideck.get< tag::dt >(), 0.0);
  storeIfSpecd< tk::real >(
    lua_ideck, "cfl", gideck.get< tag::cfl >(), 0.0);
  storeIfSpecd< uint32_t >(
    lua_ideck, "ttyi", gideck.get< tag::ttyi >(), 1);
  storeIfSpecd< bool >(
    lua_ideck, "steady_state", gideck.get< tag::steady_state >(), false);
  storeIfSpecd< tk::real >(
    lua_ideck, "residual", gideck.get< tag::residual >(), 1.0e-8);
  storeIfSpecd< uint32_t >(
    lua_ideck, "rescomp", gideck.get< tag::rescomp >(), 1);
  storeIfSpecd< uint32_t >(
    lua_ideck, "imex_runge_kutta", gideck.get< tag::imex_runge_kutta >(), 0);
  storeIfSpecd< uint32_t >(
    lua_ideck, "imex_maxiter", gideck.get< tag::imex_maxiter >(), 50);
  storeIfSpecd< tk::real >(
    lua_ideck, "imex_reltol", gideck.get< tag::imex_reltol >(), 1.0e-02);
  storeIfSpecd< tk::real >(
    lua_ideck, "imex_abstol", gideck.get< tag::imex_abstol >(), 1.0e-04);

  if (gideck.get< tag::dt >() < 1e-12 && gideck.get< tag::cfl >() < 1e-12)
    Throw("No time step calculation policy has been selected in the "
      "preceeding block. Use keyword 'dt' to set a constant or 'cfl' to set an "
      "adaptive time step size calculation policy.");

  // partitioning/reordering options
  // ---------------------------------------------------------------------------
  storeOptIfSpecd< tk::ctr::PartitioningAlgorithmType,
    tk::ctr::PartitioningAlgorithm >(
    lua_ideck, "partitioning", gideck.get< tag::partitioning >(),
    tk::ctr::PartitioningAlgorithmType::RCB);
  storeIfSpecd< bool >(
    lua_ideck, "pelocal_reorder", gideck.get< tag::pelocal_reorder >(),
    false);
  storeIfSpecd< bool >(
    lua_ideck, "operator_reorder", gideck.get< tag::operator_reorder >(),
    false);

  // discretization scheme options
  // ---------------------------------------------------------------------------
  using inciter::ctr::SchemeType;
  storeOptIfSpecd< SchemeType, inciter::ctr::Scheme >(
    lua_ideck, "scheme", gideck.get< tag::scheme >(), SchemeType::ALECG);
  storeOptIfSpecd< inciter::ctr::LimiterType, inciter::ctr::Limiter >(
    lua_ideck, "limiter", gideck.get< tag::limiter >(),
    inciter::ctr::LimiterType::NOLIMITER);
  storeIfSpecd< tk::real >(
    lua_ideck, "cweight", gideck.get< tag::cweight >(), 1.0);
  storeIfSpecd< tk::real >(
    lua_ideck, "shock_detector_coeff",
    gideck.get< tag::shock_detector_coeff >(), 1.0);
  storeIfSpecd< bool >(
    lua_ideck, "accuracy_test", gideck.get< tag::accuracy_test >(), false);
  storeIfSpecd< bool >(
    lua_ideck, "limsol_projection", gideck.get< tag::limsol_projection >(),
    true);
  storeIfSpecd< tk::real >(
    lua_ideck, "lowspeed_kp", gideck.get< tag::lowspeed_kp >(), 0.0);

  // configure solutions DOFs
  auto scheme = gideck.get< tag::scheme >();
  auto& ndof = gideck.get< tag::ndof >();
  auto& rdof = gideck.get< tag::rdof >();
  ndof = rdof = 1;
  if (scheme == SchemeType::P0P1 || scheme == SchemeType::FV) {
    ndof = 1; rdof = 4;
  } else if (scheme == SchemeType::DGP1) {
    ndof = rdof = 4;
  } else if (scheme == SchemeType::DGP2) {
    ndof = rdof = 10;
  } else if (scheme == SchemeType::PDG) {
    ndof = rdof = 10;
    gideck.get< tag::pref, tag::pref >() = true;
  } else if (scheme != SchemeType::DG &&
      scheme != SchemeType::ALECG &&
      scheme != SchemeType::OversetFE) {
    Throw("Scheme type not configured in configure_scheme");
  }

  // PDE options
  // ---------------------------------------------------------------------------

  char depvar_cnt = 'a';
  gideck.get< tag::depvar >().resize(1);

  // check transport
  if (lua_ideck["transport"].valid()) {

    checkBlock< inciter::ctr::transportList::Keys >(lua_ideck["transport"],
      "transport");

    gideck.get< tag::pde >() = inciter::ctr::PDEType::TRANSPORT;
    storeIfSpecd< std::size_t >(
      lua_ideck["transport"], "ncomp",
      gideck.get< tag::transport, tag::ncomp >(), 1);
    storeIfSpecd< int >(
      lua_ideck["transport"], "intsharp",
      gideck.get< tag::transport, tag::intsharp >(), 0);
    storeIfSpecd< tk::real >(
      lua_ideck["transport"], "intsharp_param",
      gideck.get< tag::transport, tag::intsharp_param >(), 1.8);
    storeOptIfSpecd< inciter::ctr::ProblemType, inciter::ctr::Problem >(
      lua_ideck["transport"], "problem",
      gideck.get< tag::transport, tag::problem >(),
      inciter::ctr::ProblemType::USER_DEFINED);
    storeVecIfSpecd< tk::real >(
      lua_ideck["transport"], "diffusivity",
      gideck.get< tag::transport, tag::diffusivity >(), {0.0, 0.0, 0.0});
    storeVecIfSpecd< tk::real >(
      lua_ideck["transport"], "u0",
      gideck.get< tag::transport, tag::u0 >(), {0.0, 0.0, 0.0});
    storeVecIfSpecd< tk::real >(
      lua_ideck["transport"], "lambda",
      gideck.get< tag::transport, tag::lambda >(), {0.0, 0.0, 0.0});
    gideck.get< tag::depvar >()[0] = 'c';
    storeOptIfSpecd< inciter::ctr::FluxType, inciter::ctr::Flux >(
      lua_ideck, "flux", gideck.get< tag::flux >(),
      inciter::ctr::FluxType::UPWIND);
    storeOptIfSpecd< inciter::ctr::PhysicsType, inciter::ctr::Physics >(
      lua_ideck["transport"], "physics",
      gideck.get< tag::transport, tag::physics >(),
      inciter::ctr::PhysicsType::ADVECTION);

    // store number of equations in PDE system
    gideck.get< tag::ncomp >() =
      gideck.get< tag::transport, tag::ncomp >();
  }

  // check compflow
  if (lua_ideck["compflow"].valid()) {

    checkBlock< inciter::ctr::compflowList::Keys >(lua_ideck["compflow"],
      "compflow");

    gideck.get< tag::pde >() = inciter::ctr::PDEType::COMPFLOW;
    storeOptIfSpecd< inciter::ctr::ProblemType, inciter::ctr::Problem >(
      lua_ideck["compflow"], "problem",
      gideck.get< tag::compflow, tag::problem >(),
      inciter::ctr::ProblemType::USER_DEFINED);
    storeOptIfSpecd< inciter::ctr::PhysicsType, inciter::ctr::Physics >(
      lua_ideck["compflow"], "physics",
      gideck.get< tag::compflow, tag::physics >(),
      inciter::ctr::PhysicsType::EULER);

    // problem parameters for MMS
    storeIfSpecd< tk::real >(lua_ideck["compflow"], "alpha",
      gideck.get< tag::compflow, tag::alpha >(), 0.0);
    storeIfSpecd< tk::real >(lua_ideck["compflow"], "beta",
      gideck.get< tag::compflow, tag::beta >(), 0.0);
    storeIfSpecd< tk::real >(lua_ideck["compflow"], "betax",
      gideck.get< tag::compflow, tag::betax >(), 0.0);
    storeIfSpecd< tk::real >(lua_ideck["compflow"], "betay",
      gideck.get< tag::compflow, tag::betay >(), 0.0);
    storeIfSpecd< tk::real >(lua_ideck["compflow"], "betaz",
      gideck.get< tag::compflow, tag::betaz >(), 0.0);
    storeIfSpecd< tk::real >(lua_ideck["compflow"], "r0",
      gideck.get< tag::compflow, tag::r0 >(), 0.0);
    storeIfSpecd< tk::real >(lua_ideck["compflow"], "p0",
      gideck.get< tag::compflow, tag::p0 >(), 0.0);
    storeIfSpecd< tk::real >(lua_ideck["compflow"], "ce",
      gideck.get< tag::compflow, tag::ce >(), 0.0);
    storeIfSpecd< tk::real >(lua_ideck["compflow"], "kappa",
      gideck.get< tag::compflow, tag::kappa >(), 0.0);

    gideck.get< tag::depvar >()[0] = 'a';
    storeOptIfSpecd< inciter::ctr::FluxType, inciter::ctr::Flux >(
      lua_ideck, "flux", gideck.get< tag::flux >(),
      inciter::ctr::FluxType::HLLC);

    // store number of equations in PDE system
    gideck.get< tag::ncomp >() = 5;
  }

  // check multimat
  if (lua_ideck["multimat"].valid()) {

    checkBlock< inciter::ctr::multimatList::Keys >(lua_ideck["multimat"],
      "multimat");

    gideck.get< tag::pde >() = inciter::ctr::PDEType::MULTIMAT;
    storeIfSpecd< std::size_t >(
      lua_ideck["multimat"], "nmat",
      gideck.get< tag::multimat, tag::nmat >(), 2);
    storeIfSpecd< uint64_t >(
      lua_ideck["multimat"], "prelax",
      gideck.get< tag::multimat, tag::prelax >(), 1);
    storeIfSpecd< tk::real >(
      lua_ideck["multimat"], "prelax_timescale",
      gideck.get< tag::multimat, tag::prelax_timescale >(), 0.25);
    storeIfSpecd< int >(
      lua_ideck["multimat"], "intsharp",
      gideck.get< tag::multimat, tag::intsharp >(), 0);
    storeIfSpecd< tk::real >(
      lua_ideck["multimat"], "intsharp_param",
      gideck.get< tag::multimat, tag::intsharp_param >(), 1.8);
    storeIfSpecd< uint64_t >(
      lua_ideck["multimat"], "rho0constraint",
      gideck.get< tag::multimat, tag::rho0constraint >(), 1);
    storeIfSpecd< int >(
      lua_ideck["multimat"], "dt_sos_massavg",
      gideck.get< tag::multimat, tag::dt_sos_massavg >(), 0);
    storeOptIfSpecd< inciter::ctr::ProblemType, inciter::ctr::Problem >(
      lua_ideck["multimat"], "problem",
      gideck.get< tag::multimat, tag::problem >(),
      inciter::ctr::ProblemType::USER_DEFINED);
    storeOptIfSpecd< inciter::ctr::PhysicsType, inciter::ctr::Physics >(
      lua_ideck["multimat"], "physics",
      gideck.get< tag::multimat, tag::physics >(),
      inciter::ctr::PhysicsType::EULER);
    gideck.get< tag::depvar >()[0] = 'a';
    storeOptIfSpecd< inciter::ctr::FluxType, inciter::ctr::Flux >(
      lua_ideck, "flux", gideck.get< tag::flux >(),
      inciter::ctr::FluxType::AUSM);

    // number of equations in PDE system are determined based on materials
  }

  // check multispecies
  if (lua_ideck["multispecies"].valid()) {

    checkBlock< inciter::ctr::multispeciesList::Keys >(lua_ideck["multispecies"],
      "multispecies");

    gideck.get< tag::pde >() = inciter::ctr::PDEType::MULTISPECIES;
    storeIfSpecd< std::size_t >(
      lua_ideck["multispecies"], "nspec",
      gideck.get< tag::multispecies, tag::nspec >(), 1);
    storeOptIfSpecd< inciter::ctr::ProblemType, inciter::ctr::Problem >(
      lua_ideck["multispecies"], "problem",
      gideck.get< tag::multispecies, tag::problem >(),
      inciter::ctr::ProblemType::USER_DEFINED);
    storeOptIfSpecd< inciter::ctr::PhysicsType, inciter::ctr::Physics >(
      lua_ideck["multispecies"], "physics",
      gideck.get< tag::multispecies, tag::physics >(),
      inciter::ctr::PhysicsType::EULER);
    gideck.get< tag::depvar >()[0] = 'a';
    storeOptIfSpecd< inciter::ctr::FluxType, inciter::ctr::Flux >(
      lua_ideck, "flux", gideck.get< tag::flux >(),
      inciter::ctr::FluxType::AUSM);

    // store number of equations in PDE system
    // nspec: species mass conservation equations,
    // 3: momentum equations,
    // 1: total energy equation.
    gideck.get< tag::ncomp >() =
      gideck.get< tag::multispecies, tag::nspec >() + 3 + 1;
  }

  // number of species, for future use
  std::size_t nspec(1);
  if (gideck.get< tag::pde >() == inciter::ctr::PDEType::MULTISPECIES)
    nspec = gideck.get< tag::multispecies, tag::nspec >();

  // add depvar to deck::depvars so it can be selected as outvar later
  tk::grm::depvars.insert( gideck.get< tag::depvar >()[0] );

  // Assemble material blocks
  // ---------------------------------------------------------------------------

  // solid counters
  std::size_t tmat(0), imatcntr(0), mtypei(0), isolcntr(0);
  bool is_solid(false);
  std::set< std::size_t > matidset;

  // material vector
  if (lua_ideck["material"].valid()) {

    // size material map vectors
    std::size_t nmat(1);
    if (gideck.get< tag::pde >() == inciter::ctr::PDEType::MULTIMAT)
      nmat = gideck.get< tag::multimat, tag::nmat >();
    gideck.get< tag::matidxmap, tag::eosidx >().resize(nmat);
    gideck.get< tag::matidxmap, tag::matidx >().resize(nmat);
    gideck.get< tag::matidxmap, tag::solidx >().resize(nmat);

    // size material vector appropriately
    // size of the material vector is the number of distinct types of materials
    sol::table sol_mat = lua_ideck["material"];
    gideck.get< tag::material >().resize(sol_mat.size());
    // species vector size is one, since all species are only of one type for now
    gideck.get< tag::species >().resize(1);

    // store material properties
    for (std::size_t i=0; i<gideck.get< tag::material >().size(); ++i) {

      checkBlock< inciter::ctr::materialList::Keys >(sol_mat[i+1], "material");

      auto& mati_deck = gideck.get< tag::material >()[i];
      // eos
      storeOptIfSpecd< inciter::ctr::MaterialType, inciter::ctr::Material >(
        sol_mat[i+1], "eos", mati_deck.get< tag::eos >(),
        inciter::ctr::MaterialType::STIFFENEDGAS);

      // material ids in this eos (default is for compflow i.e. single mat)
      storeVecIfSpecd< uint64_t >(
        sol_mat[i+1], "id", mati_deck.get< tag::id >(),
        std::vector< uint64_t >(1,1));

      // Track total number of materials in multiple material blocks (eos's)
      tmat += mati_deck.get< tag::id >().size();

      // Check for repeating user specified material ids
      for (auto midx : mati_deck.get< tag::id >()) {
        if (!matidset.count(midx))
          matidset.insert(midx);
        else
          Throw("Repeating material id specified");
      }

      std::size_t ntype = mati_deck.get< tag::id >().size();
      // cv
      if (!sol_mat[i+1]["cv"].valid())
        sol_mat[i+1]["cv"] = std::vector< tk::real >(ntype, 717.5);
      checkStoreMatProp(sol_mat[i+1], "cv", ntype, mati_deck.get< tag::cv >());

      // reset solid-marker
      is_solid = false;

      // Stiffened-gas materials
      if (mati_deck.get< tag::eos >() ==
        inciter::ctr::MaterialType::STIFFENEDGAS) {
        // gamma
        checkStoreMatProp(sol_mat[i+1], "gamma", ntype,
          mati_deck.get< tag::gamma >());

        // pstiff
        if (!sol_mat[i+1]["pstiff"].valid())
          sol_mat[i+1]["pstiff"] = std::vector< tk::real >(ntype, 0.0);
        checkStoreMatProp(sol_mat[i+1], "pstiff", ntype,
          mati_deck.get< tag::pstiff >());

        // mu (dynamic viscosity) and 'viscous' keyword
        if (!sol_mat[i+1]["mu"].valid())
          sol_mat[i+1]["mu"] = std::vector< tk::real >(ntype, 0.0);
        else gideck.get< tag::multimat, tag::viscous >() = true;
        checkStoreMatProp(sol_mat[i+1], "mu", ntype, mati_deck.get< tag::mu >());
      }
      // Small-shear solid materials
      else if (mati_deck.get< tag::eos >() ==
        inciter::ctr::MaterialType::SMALLSHEARSOLID) {
        // gamma
        checkStoreMatProp(sol_mat[i+1], "gamma", ntype,
          mati_deck.get< tag::gamma >());

        // pstiff
        if (!sol_mat[i+1]["pstiff"].valid())
          sol_mat[i+1]["pstiff"] = std::vector< tk::real >(ntype, 0.0);
        checkStoreMatProp(sol_mat[i+1], "pstiff", ntype,
          mati_deck.get< tag::pstiff >());

        // mu
        checkStoreMatProp(sol_mat[i+1], "mu", ntype,
          mati_deck.get< tag::mu >());

        // yield_stress
        if (!sol_mat[i+1]["yield_stress"].valid())
          sol_mat[i+1]["yield_stress"] =
            std::vector< tk::real >(ntype, 300.0e+06);
        checkStoreMatProp(sol_mat[i+1], "yield_stress", ntype,
          mati_deck.get< tag::yield_stress >());

        // assign solid
        is_solid = true;
      }
      // Godunov-Romenski aluminum materials
      else if (mati_deck.get< tag::eos >() ==
        inciter::ctr::MaterialType::GODUNOVROMENSKIALUMINUM) {
        // gamma
        checkStoreMatProp(sol_mat[i+1], "gamma", ntype,
          mati_deck.get< tag::gamma >());

        // mu
        checkStoreMatProp(sol_mat[i+1], "mu", ntype,
          mati_deck.get< tag::mu >());

        // yield_stress
        if (!sol_mat[i+1]["yield_stress"].valid())
          sol_mat[i+1]["yield_stress"] =
            std::vector< tk::real >(ntype, 300.0e+06);
        checkStoreMatProp(sol_mat[i+1], "yield_stress", ntype,
          mati_deck.get< tag::yield_stress >());

        // assign solid
        is_solid = true;
      }
      // JWL materials
      else if (mati_deck.get< tag::eos >() == inciter::ctr::MaterialType::JWL) {
        // w_gru
        checkStoreMatProp(sol_mat[i+1], "w_gru", ntype,
          mati_deck.get< tag::w_gru >());

        // a_jwl
        checkStoreMatProp(sol_mat[i+1], "A_jwl", ntype,
          mati_deck.get< tag::A_jwl >());

        // b_jwl
        checkStoreMatProp(sol_mat[i+1], "B_jwl", ntype,
          mati_deck.get< tag::B_jwl >());

        // R1_jwl
        checkStoreMatProp(sol_mat[i+1], "R1_jwl", ntype,
          mati_deck.get< tag::R1_jwl >());

        // R2_jwl
        checkStoreMatProp(sol_mat[i+1], "R2_jwl", ntype,
          mati_deck.get< tag::R2_jwl >());

        // rho0_jwl
        checkStoreMatProp(sol_mat[i+1], "rho0_jwl", ntype,
          mati_deck.get< tag::rho0_jwl >());

        // de_jwl
        checkStoreMatProp(sol_mat[i+1], "de_jwl", ntype,
          mati_deck.get< tag::de_jwl >());

        // Pr_jwl
        checkStoreMatProp(sol_mat[i+1], "Pr_jwl", ntype,
          mati_deck.get< tag::Pr_jwl >());

        // rhor_jwl
        if (sol_mat[i+1]["rhor_jwl"].valid()) {
          checkStoreMatProp(sol_mat[i+1], "rhor_jwl", ntype,
            mati_deck.get< tag::rhor_jwl >());
        }
        // Tr_jwl
        else if (sol_mat[i+1]["Tr_jwl"].valid()) {
          checkStoreMatProp(sol_mat[i+1], "Tr_jwl", ntype,
            mati_deck.get< tag::Tr_jwl >());
        }
        else
          Throw("Either reference density or reference temperature must be "
            "specified for JWL equation of state (EOS).");
      }
      // Thermally-perfect gas materials
      else if (mati_deck.get< tag::eos >() ==
        inciter::ctr::MaterialType::THERMALLYPERFECTGAS) {

        if (!lua_ideck["species"].valid())
          Throw("Species block must be specified for thermally perfect gas");
        sol::table sol_spc = lua_ideck["species"];

        // We have assumed that nmat == 1 always for multi species, and that
        // all species are of a single type, so that the outer species vector
        // is of size one
        auto& spci_deck = gideck.get< tag::species >()[0];

        // species ids (default is for single species)
        storeVecIfSpecd< uint64_t >(
          sol_spc[i+1], "id", spci_deck.get< tag::id >(),
          std::vector< uint64_t >(1,1));

        Assert(nspec == spci_deck.get< tag::id >().size(),
          "Number of ids in species-block not equal to number of species");

        // gamma
        checkStoreMatProp(sol_spc[i+1], "gamma", nspec,
          spci_deck.get< tag::gamma >());
        // R
        checkStoreMatProp(sol_spc[i+1], "R", nspec,
          spci_deck.get< tag::R >());
        // cp_coeff
        checkStoreMatPropVecVec(sol_spc[i+1], "cp_coeff", nspec, 3, 8,
          spci_deck.get< tag::cp_coeff >());
        // t_range
        checkStoreMatPropVec(sol_spc[i+1], "t_range", nspec, 4,
          spci_deck.get< tag::t_range >());
        // dH_ref
        checkStoreMatProp(sol_spc[i+1], "dH_ref", nspec,
          spci_deck.get< tag::dH_ref >());
      }

      // Generate mapping between material index and eos parameter index
      auto& eosmap = gideck.get< tag::matidxmap, tag::eosidx >();
      auto& idxmap = gideck.get< tag::matidxmap, tag::matidx >();
      auto& solidxmap = gideck.get< tag::matidxmap, tag::solidx >();
      for (auto midx : mati_deck.get< tag::id >()) {
        midx -= 1;
        eosmap[midx] = mtypei;
        idxmap[midx] = imatcntr;
        if (is_solid) {
          // add to solid-counter
          ++isolcntr;
          solidxmap[midx] = isolcntr;
        }
        ++imatcntr;
      }
      // end of materials for this eos, thus reset index counter
      imatcntr = 0;
      // increment material-type/eos-type index counter
      ++mtypei;
    }

    // Error checking on material ids
    // -------------------------------------------------------------------------

    // Total number of materials
    if (tmat != nmat) 
      Throw("The total number of materials in all the material blocks (" +
        std::to_string(tmat) +
        ") is not equal to the number of materials specified 'nmat'.");

    // Contiguous and 1-based material ids
    if (*matidset.begin() != 1)
      Throw("Material ids specified in material blocks not one-based. "
        "Material ids must begin with one.");
    std::size_t icount(1);
    for (auto midx : matidset) {
      if (midx != icount)
        Throw("Material ids specified in material blocks have a gap. "
          "Material ids must be contiguous.");
      ++icount;
    }

    // Set up number of PDEs for multimat
    if (gideck.get< tag::pde >() == inciter::ctr::PDEType::MULTIMAT) {
      auto ntot = nmat + nmat + 3 + nmat;
      // if solid EOS, add components
      const auto& solidx = gideck.get< tag::matidxmap, tag::solidx >();
      for (std::size_t i=0; i<solidx.size(); ++i) {
        if (solidx[i] > 0)
          ntot += 9;
      }
      gideck.get< tag::ncomp >() = ntot;
    }
  }

  // Mesh specification block (for overset)
  // ---------------------------------------------------------------------------
  if (lua_ideck["mesh"].valid()) {
    const sol::table& lua_mesh = lua_ideck["mesh"];
    auto& mesh_deck = gideck.get< tag::mesh >();
    mesh_deck.resize(lua_mesh.size());

    for (std::size_t i=0; i<mesh_deck.size(); ++i) {

      checkBlock< inciter::ctr::meshList::Keys >(lua_mesh[i+1], "mesh");

      // filename
      storeIfSpecd< std::string >(lua_mesh[i+1], "filename",
        mesh_deck[i].get< tag::filename >(), "");

      // location
      storeVecIfSpecd< tk::real >(lua_mesh[i+1], "location",
        mesh_deck[i].get< tag::location >(), {0.0, 0.0, 0.0});
      if (mesh_deck[i].get< tag::location >().size() != 3)
        Throw("Mesh location requires 3 coordinates.");

      // orientation
      storeVecIfSpecd< tk::real >(lua_mesh[i+1], "orientation",
        mesh_deck[i].get< tag::orientation >(), {0.0, 0.0, 0.0});
      if (mesh_deck[i].get< tag::orientation >().size() != 3)
        Throw("Mesh orientation requires 3 rotation angles.");

      // mass
      storeIfSpecd< tk::real >(lua_mesh[i+1], "mass",
        mesh_deck[i].get< tag::mass >(), 0.0);

      // Transfer object
      if (i > 0) {
        gideck.get< tag::transfer >().emplace_back( 0, i );

        // assign depvar
        ++depvar_cnt;
        gideck.get< tag::depvar >().push_back(depvar_cnt);

        // add depvar to deck::depvars so it can be selected as outvar later
        tk::grm::depvars.insert(depvar_cnt);
      }
    }
  }
  else {
    // TODO: remove double-specification of defaults
    auto& mesh_deck = gideck.get< tag::mesh >();
    mesh_deck.resize(1);
    mesh_deck[0].get< tag::filename >() =
      gideck.get< tag::cmd, tag::io, tag::input >();
    mesh_deck[0].get< tag::location >() = {0.0, 0.0, 0.0};
    mesh_deck[0].get< tag::orientation >() = {0.0, 0.0, 0.0};
    mesh_deck[0].get< tag::mass >() = 0.0;
  }

  Assert(gideck.get< tag::mesh >().size() == gideck.get< tag::depvar >().size(),
    "Number of depvar not equal to the number of meshes.");

  // Rigid body motion block for overset meshes
  // ---------------------------------------------------------------------------
  if (lua_ideck["rigid_body_motion"].valid()) {

    Assert(gideck.get< tag::mesh >().size() > 1,
      "Multiple meshes (overset) needed for rigid body motion.");

    // check that rigid body mass is provided
    const auto mesh_deck = gideck.get< tag::mesh >();
    for (std::size_t i=1; i<mesh_deck.size(); ++i) {
      Assert(mesh_deck[i].get< tag::mass >() > 1e-10,
        "Mass of body required for overset meshes with rigid body motion.");
    }

    auto& rbm_deck = gideck.get< tag::rigid_body_motion >();

    rbm_deck.get< tag::rigid_body_movt >() = true;

    // degrees of freedom
    storeIfSpecd< std::size_t >(
      lua_ideck["rigid_body_motion"], "rigid_body_dof",
      rbm_deck.get< tag::rigid_body_dof >(), 3);
    if (rbm_deck.get< tag::rigid_body_dof >() != 3 &&
      rbm_deck.get< tag::rigid_body_dof >() != 6)
      Throw("Only 3 or 6 rigid body DOFs supported.");

    // symmetry plane
    storeIfSpecd< std::size_t >(
      lua_ideck["rigid_body_motion"], "symmetry_plane",
      rbm_deck.get< tag::symmetry_plane >(), 0);
    if (rbm_deck.get< tag::symmetry_plane >() > 3)
      Throw("Rigid body motion symmetry plane must be 1(x), 2(y), or 3(z).");
    if (rbm_deck.get< tag::symmetry_plane >() == 0 &&
      rbm_deck.get< tag::rigid_body_dof >() == 3)
      Throw(
        "Rigid body motion symmetry plane must be specified for 3 DOF motion.");
    // reset to 0-based indexing
    rbm_deck.get< tag::symmetry_plane >() -= 1;
  }
  else {
    // TODO: remove double-specification of defaults
    auto& rbm_deck = gideck.get< tag::rigid_body_motion >();
    rbm_deck.get< tag::rigid_body_movt >() = false;
    rbm_deck.get< tag::rigid_body_dof >() = 0;
    rbm_deck.get< tag::symmetry_plane >() = 0;
  }

  // Field output block
  // ---------------------------------------------------------------------------
  if (lua_ideck["field_output"].valid()) {

    checkBlock< inciter::ctr::fieldOutputList::Keys >(lua_ideck["field_output"],
      "field_output");

    auto& fo_deck = gideck.get< tag::field_output >();

    // interval iteration
    storeIfSpecd< uint32_t >(
      lua_ideck["field_output"], "interval",
      fo_deck.get< tag::interval >(),
      std::numeric_limits< uint32_t >::max());

    // interval physical time
    storeIfSpecd< tk::real >(
      lua_ideck["field_output"], "time_interval",
      fo_deck.get< tag::time_interval >(),
      std::numeric_limits< tk::real >::max());

    // interval time range
    storeVecIfSpecd< tk::real >(
      lua_ideck["field_output"], "time_range",
      fo_deck.get< tag::time_range >(), {});

    // refined mesh field output
    storeIfSpecd< bool >(
      lua_ideck["field_output"], "refined", fo_deck.get< tag::refined >(),
      false);
    gideck.get< tag::cmd, tag::io, tag::refined >() = fo_deck.get< tag::refined >();

    // filetype
    storeOptIfSpecd< tk::ctr::FieldFileType, tk::ctr::FieldFile >(
      lua_ideck["field_output"], "filetype", fo_deck.get< tag::filetype >(),
      tk::ctr::FieldFileType::EXODUSII);

    // sidesets for field output
    storeVecIfSpecd< uint64_t >(
      lua_ideck["field_output"], "sideset", fo_deck.get< tag::sideset >(), {});

    // Assign outvar
    auto& foutvar = fo_deck.get< tag::outvar >();
    std::size_t nevar(0), nnvar(0), tensorcompvar(0);
    std::size_t nmat(1);
    if (gideck.get< tag::pde >() == inciter::ctr::PDEType::MULTIMAT)
      nmat = gideck.get< tag::multimat, tag::nmat >();

    // elem aliases
    std::vector< std::string > elemalias;
    if (lua_ideck["field_output"]["elemvar"].valid())
      nevar = sol::table(lua_ideck["field_output"]["elemvar"]).size();
    storeVecIfSpecd< std::string >(
      lua_ideck["field_output"], "elemalias", elemalias,
      std::vector< std::string >(nevar, ""));

    // node aliases
    std::vector< std::string > nodealias;
    if (lua_ideck["field_output"]["nodevar"].valid())
      nnvar = sol::table(lua_ideck["field_output"]["nodevar"]).size();
    storeVecIfSpecd< std::string >(
      lua_ideck["field_output"], "nodealias", nodealias,
      std::vector< std::string >(nnvar, ""));

    // error check on aliases
    if (elemalias.size() != nevar)
      Throw("elemalias should have the same size as elemvar.");
    if (nodealias.size() != nnvar)
      Throw("nodealias should have the same size as nodevar.");

    // element variables
    if (lua_ideck["field_output"]["elemvar"].valid()) {
      for (std::size_t i=0; i<nevar; ++i) {
        std::string varname(lua_ideck["field_output"]["elemvar"][i+1]);
        std::string alias(elemalias[i]);
        // add extra outvars for tensor components
        if (varname.find("_tensor") != std::string::npos) tensorcompvar += 8;
        addOutVar(varname, alias, gideck.get< tag::depvar >(), nmat,
          nspec, gideck.get< tag::pde >(), tk::Centering::ELEM, foutvar);
      }
    }

    // node variables
    if (lua_ideck["field_output"]["nodevar"].valid()) {
      for (std::size_t i=0; i<nnvar; ++i) {
        std::string varname(lua_ideck["field_output"]["nodevar"][i+1]);
        std::string alias(nodealias[i]);
        // add extra outvars for tensor components
        if (varname.find("_tensor") != std::string::npos) tensorcompvar += 8;
        addOutVar(varname, alias, gideck.get< tag::depvar >(), nmat,
          nspec, gideck.get< tag::pde >(), tk::Centering::NODE, foutvar);
      }
    }

    Assert(foutvar.size() == (nevar + nnvar + tensorcompvar),
      "Incorrectly sized outvar vector.");
  }
  else {
    // TODO: remove double-specification of defaults
    auto& fo_deck = gideck.get< tag::field_output >();
    fo_deck.get< tag::interval >() =
      std::numeric_limits< uint32_t >::max();
    fo_deck.get< tag::time_interval >() =
      std::numeric_limits< tk::real >::max();
    fo_deck.get< tag::time_range >() = {};
    fo_deck.get< tag::refined >() = false;
    fo_deck.get< tag::filetype >() = tk::ctr::FieldFileType::EXODUSII;
    fo_deck.get< tag::sideset >() = {};
  }

  // Diagnostics output block
  // ---------------------------------------------------------------------------
  if (lua_ideck["diagnostics"].valid()) {

    checkBlock< inciter::ctr::diagnosticsList::Keys >(lua_ideck["diagnostics"],
      "diagnostics");

    auto& diag_deck = gideck.get< tag::diagnostics >();

    // interval iteration
    storeIfSpecd< uint32_t >(
      lua_ideck["diagnostics"], "interval",
      diag_deck.get< tag::interval >(), 1);

    // error norm
    storeOptIfSpecd< tk::ctr::ErrorType, tk::ctr::Error >(
      lua_ideck["diagnostics"], "error", diag_deck.get< tag::error >(),
      tk::ctr::ErrorType::L2);

    // float format
    storeOptIfSpecd< tk::ctr::TxtFloatFormatType, tk::ctr::TxtFloatFormat >(
      lua_ideck["diagnostics"], "format", diag_deck.get< tag::format >(),
      tk::ctr::TxtFloatFormatType::DEFAULT);

    // precision
    storeIfSpecd< std::streamsize >(
      lua_ideck["diagnostics"], "precision",
      diag_deck.get< tag::precision >(), std::cout.precision());
  }
  else {
    // TODO: remove double-specification of defaults
    auto& diag_deck = gideck.get< tag::diagnostics >();
    diag_deck.get< tag::interval >() = 1;
    diag_deck.get< tag::error >() = tk::ctr::ErrorType::L2;
    diag_deck.get< tag::format >() = tk::ctr::TxtFloatFormatType::DEFAULT;
    diag_deck.get< tag::precision >() = std::cout.precision();
  }

  // History output block
  // ---------------------------------------------------------------------------
  if (lua_ideck["history_output"].valid()) {

    checkBlock< inciter::ctr::historyOutputList::Keys >(
      lua_ideck["history_output"], "history_output");

    auto& hist_deck = gideck.get< tag::history_output >();

    // interval iteration
    storeIfSpecd< uint32_t >(
      lua_ideck["history_output"], "interval",
      hist_deck.get< tag::interval >(),
      std::numeric_limits< uint32_t >::max());

    // interval time
    storeIfSpecd< tk::real >(
      lua_ideck["history_output"], "time_interval",
      hist_deck.get< tag::time_interval >(),
      std::numeric_limits< tk::real >::max());

    // interval time range
    storeVecIfSpecd< tk::real >(
      lua_ideck["history_output"], "time_range",
      hist_deck.get< tag::time_range >(), {});

    // point probes
    if (lua_ideck["history_output"]["point"].valid()) {
      const sol::table& sol_pt = lua_ideck["history_output"]["point"];
      hist_deck.get< tag::point >().resize(sol_pt.size());

      for (std::size_t i=0; i<hist_deck.get< tag::point >().size(); ++i) {
        auto& pti = hist_deck.get< tag::point >()[i];
        storeIfSpecd< std::string >(
          sol_pt[i+1], "id", pti.get< tag::id >(), "p");
        storeVecIfSpecd< tk::real >(
          sol_pt[i+1], "coord", pti.get< tag::coord >(), {});
      }
    }

    // float format
    storeOptIfSpecd< tk::ctr::TxtFloatFormatType, tk::ctr::TxtFloatFormat >(
      lua_ideck["history_output"], "format", hist_deck.get< tag::format >(),
      tk::ctr::TxtFloatFormatType::DEFAULT);

    // precision
    storeIfSpecd< std::streamsize >(
      lua_ideck["history_output"], "precision",
      hist_deck.get< tag::precision >(), std::cout.precision());

    // error check point
    for (std::size_t i=0; i<hist_deck.get< tag::point >().size(); ++i) {
      if (hist_deck.get< tag::point >()[i].get< tag::coord >().size() != 3)
      Throw("Three reals required for point coordinates in history_output.");
    }
  }
  else {
    // TODO: remove double-specification of defaults
    auto& hist_deck = gideck.get< tag::history_output >();
    hist_deck.get< tag::interval >() =
      std::numeric_limits< uint32_t >::max();
    hist_deck.get< tag::time_interval >() =
      std::numeric_limits< tk::real >::max();
    hist_deck.get< tag::time_range >() = {};
    hist_deck.get< tag::precision >() = std::cout.precision();
    hist_deck.get< tag::point >().resize(0);
  }

  // ALE block
  // ---------------------------------------------------------------------------
  gideck.get< tag::ale, tag::ale >() = false;
  if (lua_ideck["ale"].valid()) {
    auto& ale_deck = gideck.get< tag::ale >();
    ale_deck.get< tag::ale >() = true;

    // Mesh velocity smoother
    storeOptIfSpecd< inciter::ctr::MeshVelocitySmootherType,
      inciter::ctr::MeshVelocitySmoother >(lua_ideck["ale"], "smoother",
      ale_deck.get< tag::smoother >(), inciter::ctr::MeshVelocitySmootherType::NONE);

    // Mesh velocity
    storeOptIfSpecd< inciter::ctr::MeshVelocityType,
      inciter::ctr::MeshVelocity >(lua_ideck["ale"], "mesh_velocity",
      ale_deck.get< tag::mesh_velocity >(), inciter::ctr::MeshVelocityType::SINE);

    // Mesh motion direction
    storeVecIfSpecd< std::size_t >(lua_ideck["ale"], "mesh_motion",
      ale_deck.get< tag::mesh_motion >(), { 0, 1, 2 });

    // Mesh force
    storeVecIfSpecd< tk::real >(lua_ideck["ale"], "meshforce",
      ale_deck.get< tag::meshforce >(), { 0, 0, 0, 0 });

    // Dirichlet
    storeVecIfSpecd< std::size_t >(lua_ideck["ale"], "dirichlet",
      ale_deck.get< tag::dirichlet >(), {});

    // Symmetry
    storeVecIfSpecd< std::size_t >(lua_ideck["ale"], "symmetry",
      ale_deck.get< tag::symmetry >(), {});

    // Move sidesets with user defined function
    if (lua_ideck["ale"]["move"].valid()) {
      const sol::table& sol_mv = lua_ideck["ale"]["move"];
      ale_deck.get< tag::move >().resize(sol_mv.size());

      for (std::size_t i=0; i<ale_deck.get< tag::move >().size(); ++i) {
        auto& mvi = ale_deck.get< tag::move >()[i];
        storeOptIfSpecd< tk::ctr::UserTableType, tk::ctr::UserTable >(
          sol_mv[i+1], "fntype", mvi.get< tag::fntype >(),
          tk::ctr::UserTableType::POSITION);
        storeVecIfSpecd< uint64_t >(
          sol_mv[i+1], "sideset", mvi.get< tag::sideset >(), {});
        storeVecIfSpecd< tk::real >(
          sol_mv[i+1], "fn", mvi.get< tag::fn >(), {});

        // error checking on user-def function
        if (mvi.get< tag::fn >().size() % 4 != 0)
          Throw("Incomplete user-defined function for ALE sideset movement. An "
          "R->R^3 function is expected, the number of descrete entries must be "
          "divisible by 4: one 'column' for the abscissa, and 3 for the "
          "ordinate.");
      }
    }

    // dv-CFL
    storeIfSpecd< tk::real >(lua_ideck["ale"], "dvcfl", ale_deck.get< tag::dvcfl >(),
      0.01);

    // Vorticity multiplier
    storeIfSpecd< tk::real >(lua_ideck["ale"], "vortmult",
      ale_deck.get< tag::vortmult >(), 0.0);

    // Mesh velocity max iterations
    storeIfSpecd< std::size_t >(lua_ideck["ale"], "maxit",
      ale_deck.get< tag::maxit >(), 5);

    // Mesh velocity max iterations
    storeIfSpecd< tk::real >(lua_ideck["ale"], "tolerance",
      ale_deck.get< tag::tolerance >(), 1e-2);
  }

  // AMR block
  // ---------------------------------------------------------------------------
  gideck.get< tag::amr, tag::amr >() = false;
  gideck.get< tag::amr, tag::maxlevels >() = 2; // this is needed for outref
  if (lua_ideck["amr"].valid()) {
    auto& amr_deck = gideck.get< tag::amr >();
    amr_deck.get< tag::amr >() = true;

    // Initial refinement toggle
    storeIfSpecd< bool >(lua_ideck["amr"], "t0ref", amr_deck.get< tag::t0ref >(),
      false);

    // Mesh refinement during time-stepping toggle
    storeIfSpecd< bool >(lua_ideck["amr"], "dtref", amr_deck.get< tag::dtref >(),
      false);

    // Uniform mesh refinement during time-stepping toggle
    storeIfSpecd< bool >(lua_ideck["amr"], "dtref_uniform",
      amr_deck.get< tag::dtref_uniform >(), false);

    // Mesh refinement frequency during time-stepping toggle
    storeIfSpecd< std::size_t >(lua_ideck["amr"], "dtfreq",
      amr_deck.get< tag::dtfreq >(), 3);

    // Maximum AMR levels
    storeIfSpecd< std::size_t >(lua_ideck["amr"], "maxlevels",
      amr_deck.get< tag::maxlevels >(), 2);

    // Initial AMR steps
    storeOptVecIfSpecd< inciter::ctr::AMRInitialType, inciter::ctr::AMRInitial >(
      lua_ideck["amr"], "initial", amr_deck.get< tag::initial >(), {});

    // Initial AMR coordinate based
    if (lua_ideck["amr"]["coords"].valid()) {
      auto rmax = std::numeric_limits< tk::real >::max() / 100;

      storeIfSpecd< tk::real >(lua_ideck["amr"]["coords"], "xminus",
        amr_deck.get< tag::coords, tag::xminus >(), rmax);
      storeIfSpecd< tk::real >(lua_ideck["amr"]["coords"], "xplus",
        amr_deck.get< tag::coords, tag::xplus >(), -rmax);

      storeIfSpecd< tk::real >(lua_ideck["amr"]["coords"], "yminus",
        amr_deck.get< tag::coords, tag::yminus >(), rmax);
      storeIfSpecd< tk::real >(lua_ideck["amr"]["coords"], "yplus",
        amr_deck.get< tag::coords, tag::yplus >(), -rmax);

      storeIfSpecd< tk::real >(lua_ideck["amr"]["coords"], "zminus",
        amr_deck.get< tag::coords, tag::zminus >(), rmax);
      storeIfSpecd< tk::real >(lua_ideck["amr"]["coords"], "zplus",
        amr_deck.get< tag::coords, tag::zplus >(), -rmax);
    }

    // Initial AMR edgelist based
    storeVecIfSpecd< std::size_t >(lua_ideck["amr"], "edgelist",
      amr_deck.get< tag::edgelist >(), {});
    if (amr_deck.get< tag::edgelist >().size() % 2 != 0)
      Throw("The number of edge-nodes, marking edges as pairs of nodes, used "
        "for explicit tagging of edges for initial mesh refineoment, is odd "
        "(it must be even).");

    // Error type for AMR
    storeOptIfSpecd< inciter::ctr::AMRErrorType, inciter::ctr::AMRError >(
      lua_ideck["amr"], "error", amr_deck.get< tag::error >(),
      inciter::ctr::AMRErrorType::JUMP);

    // Tolerances for refine/de-refine
    storeIfSpecd< tk::real >(lua_ideck["amr"], "tol_refine",
      amr_deck.get< tag::tol_refine >(), 0.2);
    storeIfSpecd< tk::real >(lua_ideck["amr"], "tol_derefine",
      amr_deck.get< tag::tol_derefine >(), 0.05);
  }

  // p-refinement block
  // ---------------------------------------------------------------------------
  gideck.get< tag::pref, tag::pref >() = false;
  if (lua_ideck["pref"].valid()) {
    auto& pref_deck = gideck.get< tag::pref >();
    pref_deck.get< tag::pref >() = true;

    // p-ref indicator type
    storeOptIfSpecd< inciter::ctr::PrefIndicatorType,
      inciter::ctr::PrefIndicator >(lua_ideck["pref"], "indicator",
      pref_deck.get< tag::indicator >(),
      inciter::ctr::PrefIndicatorType::SPECTRAL_DECAY);

    // p-ref max degrees-of-freedom per cell
    storeIfSpecd< std::size_t >(lua_ideck["pref"], "ndofmax",
      pref_deck.get< tag::ndofmax >(), 10);

    // p-ref tolerance
    storeIfSpecd< tk::real >(lua_ideck["pref"], "tolref",
      pref_deck.get< tag::tolref >(), 0.5);

    // error checking on the tolerance
    if (pref_deck.get< tag::tolref >() < 0.0 || pref_deck.get< tag::tolref >() > 1.0)
      Throw("The p-refinement tolerance must be a real number "
        "between 0.0 and 1.0, both inclusive.");
  }

  // Boundary conditions block
  // ---------------------------------------------------------------------------
  if (lua_ideck["bc"].valid()) {
    std::set< std::size_t > totalmesh;
    const sol::table& sol_bc = lua_ideck["bc"];
    auto& bc_deck = gideck.get< tag::bc >();
    bc_deck.resize(sol_bc.size());

    for (std::size_t i=0; i<bc_deck.size(); ++i) {

      checkBlock< inciter::ctr::bcList::Keys >(sol_bc[i+1], "bc");

      storeVecIfSpecd< std::size_t >(sol_bc[i+1], "mesh",
        bc_deck[i].get< tag::mesh >(), {1});
      // collect meshes for error checking
      totalmesh.insert(bc_deck[i].get< tag::mesh >().begin(),
        bc_deck[i].get< tag::mesh >().end());

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "dirichlet",
        bc_deck[i].get< tag::dirichlet >(), {});

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "symmetry",
        bc_deck[i].get< tag::symmetry >(), {});

      if (sol_bc[i+1]["inlet"].valid()) {
        const sol::table& sol_inbc = sol_bc[i+1]["inlet"];
        auto& inbc_deck = bc_deck[i].get< tag::inlet >();
        inbc_deck.resize(sol_inbc.size());

        for (std::size_t j=0; j<inbc_deck.size(); ++j) {
          storeVecIfSpecd< uint64_t >(sol_inbc[j+1], "sideset",
            inbc_deck[j].get< tag::sideset >(), {});

          storeVecIfSpecd< tk::real >(sol_inbc[j+1], "velocity",
            inbc_deck[j].get< tag::velocity >(), {0.0, 0.0, 0.0});
          if (inbc_deck[j].get< tag::velocity >().size() != 3)
            Throw("Inlet velocity requires 3 components.");

          storeIfSpecd< tk::real >(sol_inbc[j+1], "pressure",
            inbc_deck[j].get< tag::pressure >(), 0.0);

          storeIfSpecd< tk::real >(sol_inbc[j+1], "temperature",
            inbc_deck[j].get< tag::temperature >(), 0.0);

          storeIfSpecd< std::size_t >(sol_inbc[j+1], "materialid",
            inbc_deck[j].get< tag::materialid >(), 1);
        }
      }

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "outlet",
        bc_deck[i].get< tag::outlet >(), {});

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "farfield",
        bc_deck[i].get< tag::farfield >(), {});

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "extrapolate",
        bc_deck[i].get< tag::extrapolate >(), {});

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "noslipwall",
        bc_deck[i].get< tag::noslipwall >(), {});

      // Time-dependent BC
      if (sol_bc[i+1]["timedep"].valid()) {
        const sol::table& sol_tdbc = sol_bc[i+1]["timedep"];
        auto& tdbc_deck = bc_deck[i].get< tag::timedep >();
        tdbc_deck.resize(sol_tdbc.size());

        for (std::size_t j=0; j<tdbc_deck.size(); ++j) {
          storeVecIfSpecd< uint64_t >(sol_tdbc[j+1], "sideset",
            tdbc_deck[j].get< tag::sideset >(), {});
          storeVecIfSpecd< tk::real >(sol_tdbc[j+1], "fn",
            tdbc_deck[j].get< tag::fn >(), {});

          // error checking on user-def function
          if (tdbc_deck[j].get< tag::fn >().size() % 6 != 0)
            Throw("Incomplete user-defined function for time-dependent BC. An "
            "R->R^5 function is expected, the number of descrete entries must "
            "be divisible by 6: one 'column' for the abscissa, and 5 for the "
            "ordinate.");
        }
      }

      // Stagnation point
      storeVecIfSpecd< tk::real >(sol_bc[i+1], "stag_point",
        bc_deck[i].get< tag::stag_point >(), {});
      if (!bc_deck[i].get< tag::stag_point >().empty() &&
        bc_deck[i].get< tag::stag_point >().size() % 3 != 0)
        Throw("BC stagnation point requires 3 coordinate values for each "
          "point. Thus, this vector must be divisible by 3.");

      // Stagnation radius
      storeIfSpecd< tk::real >(sol_bc[i+1], "radius",
        bc_deck[i].get< tag::radius >(), 0.0);

      // Velocity for inlet/farfield
      storeVecIfSpecd< tk::real >(sol_bc[i+1], "velocity",
        bc_deck[i].get< tag::velocity >(), {0.0, 0.0, 0.0});
      if (bc_deck[i].get< tag::velocity >().size() != 3)
        Throw("BC velocity requires 3 components.");

      // Pressure for inlet/outlet/farfield
      storeIfSpecd< tk::real >(sol_bc[i+1], "pressure",
        bc_deck[i].get< tag::pressure >(), 0.0);

      // Density for inlet/outlet/farfield
      storeIfSpecd< tk::real >(sol_bc[i+1], "density",
        bc_deck[i].get< tag::density >(), 0.0);

      // Temperature for inlet/outlet/farfield
      storeIfSpecd< tk::real >(sol_bc[i+1], "temperature",
        bc_deck[i].get< tag::temperature >(), 0.0);

      // Mass fractions for inlet/farfield
      storeVecIfSpecd< tk::real >(sol_bc[i+1], "mass_fractions",
        bc_deck[i].get< tag::mass_fractions >(),
        std::vector< tk::real >(nspec, 1.0/static_cast<tk::real>(nspec)));
      if (bc_deck[i].get< tag::mass_fractions >().size() != nspec)
        Throw("BC mass fraction has incorrect number of species. "
          "Expected " + std::to_string(nspec));

      // Material-id for inlet/outlet/farfield
      storeIfSpecd< std::size_t >(sol_bc[i+1], "materialid",
        bc_deck[i].get< tag::materialid >(), 1);
    }

    // error checking on number of meshes
    if (totalmesh.size() != gideck.get< tag::mesh >().size())
      Throw("Total meshes (" + std::to_string(gideck.get< tag::mesh >().size()) +
        ") not equal to the meshes on which BC's are specified (" +
        std::to_string(totalmesh.size()));

    // error checking on mesh ids
    std::size_t ic(1);
    for (const auto& im : totalmesh) {
      if (im != ic) Throw("Non-contiguous mesh ids in BC-mesh");
      ++ic;
    }
  }
  else if (gideck.get< tag::scheme >() == inciter::ctr::SchemeType::ALECG ||
           gideck.get< tag::scheme >() == inciter::ctr::SchemeType::OversetFE)
    gideck.get< tag::bc >().resize(1);
  // error checking for unspecified BC's
  else
    Throw("No boundary conditions specified in input file.");

  // Initial condition block
  // ---------------------------------------------------------------------------
  if (lua_ideck["ic"].valid()) {
    auto& ic_deck = gideck.get< tag::ic >();

    checkBlock< inciter::ctr::icList::Keys >(lua_ideck["ic"], "ic");

    // background IC values
    storeIfSpecd< std::size_t >(lua_ideck["ic"], "materialid",
      ic_deck.get< tag::materialid >(), 1);

    storeIfSpecd< tk::real >(lua_ideck["ic"], "pressure",
      ic_deck.get< tag::pressure >(), 0.0);

    storeIfSpecd< tk::real >(lua_ideck["ic"], "temperature",
      ic_deck.get< tag::temperature >(), 0.0);

    storeVecIfSpecd< tk::real >(lua_ideck["ic"], "mass_fractions",
      ic_deck.get< tag::mass_fractions >(),
      std::vector< tk::real >(nspec, 1.0/static_cast<tk::real>(nspec)));
    if (ic_deck.get< tag::mass_fractions >().size() != nspec)
      Throw("IC mass fraction has incorrect number of species. "
        "Expected " + std::to_string(nspec));

    storeIfSpecd< tk::real >(lua_ideck["ic"], "density",
      ic_deck.get< tag::density >(), 0.0);

    storeIfSpecd< tk::real >(lua_ideck["ic"], "energy",
      ic_deck.get< tag::energy >(), 0.0);

    storeVecIfSpecd< tk::real >(lua_ideck["ic"], "velocity",
      ic_deck.get< tag::velocity >(), {0.0, 0.0, 0.0});
    if (ic_deck.get< tag::velocity >().size() != 3)
      Throw("Velocity in IC requires 3 components.");

    // IC box
    if (lua_ideck["ic"]["box"].valid()) {
      const sol::table& lua_box = lua_ideck["ic"]["box"];
      auto& box_deck = ic_deck.get< tag::box >();
      box_deck.resize(lua_box.size());

      for (std::size_t i=0; i<box_deck.size(); ++i) {

        checkBlock< inciter::ctr::boxList::Keys >(lua_box[i+1], "box");

        storeIfSpecd< std::size_t >(lua_box[i+1], "materialid",
          box_deck[i].get< tag::materialid >(), 1);

        storeIfSpecd< tk::real >(lua_box[i+1], "volume",
          box_deck[i].get< tag::volume >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "mass",
          box_deck[i].get< tag::mass >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "density",
          box_deck[i].get< tag::density >(), 0.0);

        storeVecIfSpecd< tk::real >(lua_box[i+1], "velocity",
          box_deck[i].get< tag::velocity >(), {0.0, 0.0, 0.0});
        if (box_deck[i].get< tag::velocity >().size() != 3)
          Throw("Velocity in IC box requires 3 components.");

        storeIfSpecd< tk::real >(lua_box[i+1], "pressure",
          box_deck[i].get< tag::pressure >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "energy",
          box_deck[i].get< tag::energy >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "energy_content",
          box_deck[i].get< tag::energy_content >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "temperature",
          box_deck[i].get< tag::temperature >(), 0.0);

        storeVecIfSpecd< tk::real >(lua_box[i+1], "mass_fractions",
          box_deck[i].get< tag::mass_fractions >(),
          std::vector< tk::real >(nspec, 1.0/static_cast<tk::real>(nspec)));
        if (box_deck[i].get< tag::mass_fractions >().size() != nspec)
          Throw("IC box mass fraction has incorrect number of species. "
            "Expected " + std::to_string(nspec));

        storeIfSpecd< tk::real >(lua_box[i+1], "xmin",
          box_deck[i].get< tag::xmin >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "xmax",
          box_deck[i].get< tag::xmax >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "ymin",
          box_deck[i].get< tag::ymin >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "ymax",
          box_deck[i].get< tag::ymax >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "zmin",
          box_deck[i].get< tag::zmin >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "zmax",
          box_deck[i].get< tag::zmax >(), 0.0);

        storeVecIfSpecd< tk::real >(lua_box[i+1], "orientation",
          box_deck[i].get< tag::orientation >(), {0.0, 0.0, 0.0});
        if (box_deck[i].get< tag::orientation >().size() != 3)
          Throw("Orientation in IC box requires 3 rotation angles.");

        storeOptIfSpecd< inciter::ctr::InitiateType, inciter::ctr::Initiate >(
          lua_box[i+1], "initiate", box_deck[i].get< tag::initiate >(),
          inciter::ctr::InitiateType::IMPULSE);

        storeVecIfSpecd< tk::real >(lua_box[i+1], "point",
          box_deck[i].get< tag::point >(), {0.0, 0.0, 0.0});
        if (box_deck[i].get< tag::point >().size() != 3)
          Throw("Point in IC box requires 3 coordinates.");

        storeIfSpecd< tk::real >(lua_box[i+1], "init_time",
          box_deck[i].get< tag::init_time >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "front_width",
          box_deck[i].get< tag::front_width >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "front_speed",
          box_deck[i].get< tag::front_speed >(), 0.0);
      }
    }

    // IC mesh-block
    if (lua_ideck["ic"]["meshblock"].valid()) {
      const sol::table& lua_meshblock = lua_ideck["ic"]["meshblock"];
      auto& mblk_deck = ic_deck.get< tag::meshblock >();
      mblk_deck.resize(lua_meshblock.size());

      for (std::size_t i=0; i<mblk_deck.size(); ++i) {

        checkBlock< inciter::ctr::meshblockList::Keys >(lua_meshblock[i+1],
          "meshblock");

        storeIfSpecd< std::size_t >(lua_meshblock[i+1], "blockid",
          mblk_deck[i].get< tag::blockid >(), 0);
        if (mblk_deck[i].get< tag::blockid >() == 0)
          Throw("Each IC mesh block must specify the mesh block id.");

        storeIfSpecd< std::size_t >(lua_meshblock[i+1], "materialid",
          mblk_deck[i].get< tag::materialid >(), 1);

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "energy_content",
          mblk_deck[i].get< tag::energy_content >(), 0.0);

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "volume",
          mblk_deck[i].get< tag::volume >(), 0.0);
        if (mblk_deck[i].get< tag::energy_content >() > 0.0 &&
          mblk_deck[i].get< tag::volume >() < 1e-12)
          Throw("Mesh block volume must be specified, if energy content is "
            "used to initialize block");

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "mass",
          mblk_deck[i].get< tag::mass >(), 0.0);

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "density",
          mblk_deck[i].get< tag::density >(), 0.0);

        storeVecIfSpecd< tk::real >(lua_meshblock[i+1], "velocity",
          mblk_deck[i].get< tag::velocity >(), {0.0, 0.0, 0.0});
        if (mblk_deck[i].get< tag::velocity >().size() != 3)
          Throw("Velocity in IC meshblock requires 3 components.");

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "pressure",
          mblk_deck[i].get< tag::pressure >(), 0.0);

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "energy",
          mblk_deck[i].get< tag::energy >(), 0.0);

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "temperature",
          mblk_deck[i].get< tag::temperature >(), 0.0);

        storeVecIfSpecd< tk::real >(lua_meshblock[i+1], "mass_fractions",
          mblk_deck[i].get< tag::mass_fractions >(),
          std::vector< tk::real >(nspec, 1.0/static_cast<tk::real>(nspec)));
        if (mblk_deck[i].get< tag::mass_fractions >().size() != nspec)
          Throw("IC meshblock mass fraction has incorrect number of species. "
            "Expected " + std::to_string(nspec));

        storeOptIfSpecd< inciter::ctr::InitiateType, inciter::ctr::Initiate >(
          lua_meshblock[i+1], "initiate", mblk_deck[i].get< tag::initiate >(),
          inciter::ctr::InitiateType::IMPULSE);

        storeVecIfSpecd< tk::real >(lua_meshblock[i+1], "point",
          mblk_deck[i].get< tag::point >(), {0.0, 0.0, 0.0});
        if (mblk_deck[i].get< tag::point >().size() != 3)
          Throw("Point in IC meshblock requires 3 coordinates.");

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "init_time",
          mblk_deck[i].get< tag::init_time >(), 0.0);

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "front_width",
          mblk_deck[i].get< tag::front_width >(), 0.0);

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "front_speed",
          mblk_deck[i].get< tag::front_speed >(), 0.0);
      }
    }
  }
  else {
    // TODO: remove double-specification of defaults
    auto& ic_deck = gideck.get< tag::ic >();
    ic_deck.get< tag::materialid >() = 1;
    ic_deck.get< tag::pressure >() = 0.0;
    ic_deck.get< tag::temperature >() = 1.0;
    ic_deck.get< tag::density >() = 0.0;
    ic_deck.get< tag::energy >() = 0.0;
    ic_deck.get< tag::velocity >() = {0.0, 0.0, 0.0};
    ic_deck.get< tag::mass_fractions >() =
      std::vector< tk::real >(nspec, 1.0/static_cast<tk::real>(nspec));
  }
}

void
LuaParser::checkStoreMatProp(
  const sol::table table,
  const std::string key,
  std::size_t vecsize,
  std::vector< tk::real >& storage )
// *****************************************************************************
//  Check and store material property into inpudeck storage
//! \param[in] table Sol-table which contains said property
//! \param[in] key Key for said property in Sol-table
//! \param[in] vecsize Number of said property in Sol-table (based on number of
//!   materials that are of the same eos type
//! \param[in,out] storage Storage space in inputdeck where said property is
//!   to be stored
// *****************************************************************************
{
  // check validity of table
  if (!table[key].valid())
    Throw("Material property '" + key + "' not specified");
  if (sol::table(table[key]).size() != vecsize)
    Throw("Incorrect number of '" + key + "'s specified. Expected " +
      std::to_string(vecsize));

  // store values from table to inputdeck
  storeVecIfSpecd< tk::real >(table, key, storage,
    std::vector< tk::real >(vecsize, 0.0));
}

void
LuaParser::checkStoreMatPropVec(
  const sol::table table,
  const std::string key,
  std::size_t nspec,
  std::size_t vecsize,
  std::vector<std::vector< tk::real >>& storage )
// *****************************************************************************
//  Check and store material property vector into inpudeck storage
//! \param[in] table Sol-table which contains said property
//! \param[in] key Key for said property in Sol-table
//! \param[in] nspec Number of species
//! \param[in] vecsize Number of said property in Sol-table (based on number of
//!   coefficients for the defined species)
//! \param[in,out] storage Storage space in inputdeck where said property is
//!   to be stored
// *****************************************************************************
{
  // check validity of table
  if (!table[key].valid())
    Throw("Material property '" + key + "' not specified");
  if (sol::table(table[key]).size() != nspec)
    Throw("Incorrect number of '" + key + "' vectors specified. Expected " +
      std::to_string(nspec) + " vectors");

  storage.resize(nspec);

  const auto& tableentry = table[key];
  for (std::size_t k=0; k < nspec; k++) {
    if (sol::table(tableentry[k+1]).size() != vecsize)
      Throw("Incorrect number of '" + key + "' entries in vector of species "
        + std::to_string(k+1) + " specified. Expected " +
        std::to_string(vecsize));

    // store values from table to inputdeck
    for (std::size_t i=0; i<vecsize; ++i)
      storage[k].push_back(tableentry[k+1][i+1]);
  }
}

void
LuaParser::checkStoreMatPropVecVec(
  const sol::table table,
  const std::string key,
  std::size_t nspec,
  std::size_t vecsize1,
  std::size_t vecsize2,
  std::vector<std::vector<std::vector< tk::real >>>& storage )
// *****************************************************************************
//  Check and store material property vector into inpudeck storage
//! \param[in] table Sol-table which contains said property
//! \param[in] key Key for said property in Sol-table
//! \param[in] nspec Number of species
//! \param[in] vecsize1 Outer number of said property in Sol-table (based on
//!   number of coefficients for the defined species)
//! \param[in] vecsize2 Inner number of said property in Sol-table (based on
//!   number of coefficients for the defined species)
//! \param[in,out] storage Storage space in inputdeck where said property is
//!   to be stored
// *****************************************************************************
{
  // check validity of table
  if (!table[key].valid())
    Throw("Material property '" + key + "' not specified");
  if (sol::table(table[key]).size() != nspec)
    Throw("Incorrect number of '" + key + "' vectors specified. Expected " +
      std::to_string(nspec) + " vectors");

  storage.resize(nspec);

  const auto& tableentry = table[key];
  for (std::size_t k=0; k < nspec; k++) {
    if (sol::table(tableentry[k+1]).size() != vecsize1)
      Throw("Incorrect outer number of '" + key + "' entries in vector of species "
        + std::to_string(k+1) + " specified. Expected " +
        std::to_string(vecsize1));

    // store values from table to inputdeck
    for (std::size_t i=0; i<vecsize1; ++i) {
      if (sol::table(tableentry[k+1][i+1]).size() != vecsize2)
        Throw("Incorrect inner number of '" + key + "' entries in vector of species "
          + std::to_string(k+1) + " specified. Expected " +
          std::to_string(vecsize2));
      std::vector< tk::real > temp_storage;
      for (std::size_t j=0; j<vecsize2; j++)
        temp_storage.push_back(tableentry[k+1][i+1][j+1]);
      storage[k].push_back(temp_storage);
    }
  }
}

void
LuaParser::addOutVar(
  const std::string& varname,
  const std::string& alias,
  std::vector< char >& depv,
  std::size_t nmat,
  std::size_t nspec,
  inciter::ctr::PDEType pde,
  tk::Centering c,
  std::vector< inciter::ctr::OutVar >& foutvar )
// *****************************************************************************
//  Check and store field output variables
//! \param[in] varname Name of variable requested
//! \param[in] alias User specified alias for output
//! \param[in] depv List of depvars
//! \param[in] nmat Number of materials configured
//! \param[in] nspec Number of species configured
//! \param[in] pde Type of PDE configured
//! \param[in] c Variable centering requested
//! \param[in,out] foutvar Input deck storage where output vars are stored
// *****************************************************************************
{
  // index-based quantity specification
  // ----------------------------------
  if (varname.length() == 2) {
    auto qty = varname.at(0);
    auto j = std::stoul(std::string{varname.at(1)}) - 1;

    if (pde == inciter::ctr::PDEType::MULTIMAT) {
    // multimat/matvar quantities
      if (qty == 'D') {  // density
        foutvar.emplace_back(
          inciter::ctr::OutVar(c, varname, alias, inciter::densityIdx(nmat,j)) );
      }
      else if (qty == 'F') {  // volume fraction
        foutvar.emplace_back(
          inciter::ctr::OutVar(c, varname, alias, inciter::volfracIdx(nmat,j)) );
      }
      else if (qty == 'M') {  // momentum
        foutvar.emplace_back(
          inciter::ctr::OutVar(c, varname, alias, inciter::momentumIdx(nmat,j)) );
      }
      else if (qty == 'E') {  // specific total energy
        foutvar.emplace_back(
          inciter::ctr::OutVar(c, varname, alias, inciter::energyIdx(nmat,j)) );
      }
      else if (qty == 'U') {  // velocity (primitive)
        foutvar.emplace_back(
          inciter::ctr::OutVar(c, varname, alias, inciter::velocityIdx(nmat,j)) );
      }
      else if (qty == 'P') {  // material pressure (primitive)
        foutvar.emplace_back(
          inciter::ctr::OutVar(c, varname, alias, inciter::pressureIdx(nmat,j)) );
      }
      else {
        // error out if incorrect matvar used
        Throw("field_output: matvar " + varname + " not found");
      }
    }
    else if (pde == inciter::ctr::PDEType::MULTISPECIES) {
    // multispecies/matvar quantities
      if (qty == 'D') {  // density
        foutvar.emplace_back(
          inciter::ctr::OutVar(c, varname, alias,
            inciter::multispecies::densityIdx(nspec,j)) );
      }
      else if (qty == 'M') {  // momentum
        foutvar.emplace_back(
          inciter::ctr::OutVar(c, varname, alias,
            inciter::multispecies::momentumIdx(nspec,j)) );
      }
      else if (qty == 'E') {  // specific total energy
        foutvar.emplace_back(
          inciter::ctr::OutVar(c, varname, alias,
            inciter::multispecies::energyIdx(nspec,j)) );
      }
      else {
        // error out if incorrect matvar used
        Throw("field_output: matvar " + varname + " not found");
      }
    }
    else {
    // quantities specified by depvar
      for (const auto& id : depv)
        if (std::tolower(qty) == id) {
          foutvar.emplace_back( inciter::ctr::OutVar(c, varname, alias, j) );
        }
    }
  }
  // analytic quantity specification
  // -------------------------------
  else if (varname.find("analytic") != std::string::npos) {
    foutvar.emplace_back( inciter::ctr::OutVar(c, varname, alias, 0) );
  }
  // name-based tensor quantity specification
  // ----------------------------------------
  else if (varname.find("_tensor") != std::string::npos) {
    std::string namet(varname);
    // remove the "_tensor" from the varname
    for (std::size_t i=0; i<7; ++i) namet.pop_back();

    for (std::size_t i=1; i<=3; ++i) {
      for (std::size_t j=1; j<=3; ++j) {
        std::string tij(namet + std::to_string(i) + std::to_string(j));
        foutvar.emplace_back( inciter::ctr::OutVar(c, tij, alias, 0, tij) );
      }
    }
  }
  // name-based quantity specification
  // ---------------------------------
  else {
    foutvar.emplace_back( inciter::ctr::OutVar(c, varname, alias, 0, varname) );
  }
}
