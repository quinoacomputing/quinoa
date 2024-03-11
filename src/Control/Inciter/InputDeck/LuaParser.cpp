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
#include "Inciter/Types.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "Inciter/InputDeck/NewInputDeck.hpp"
#include "Inciter/InputDeck/New2InputDeck.hpp"
#include "Inciter/InputDeck/LuaParser.hpp"
#include "Inciter/InputDeck/Grammar.hpp"

namespace tk {
namespace grm {

//TODO: remove this?
extern tk::Print g_print;

} // grm::
} // tk::

using inciter::LuaParser;

LuaParser::LuaParser( const tk::Print& print,
                      const ctr::CmdLine& cmdline,
                      ctr::New2InputDeck& inputdeck ) :
  FileParser( cmdline.get< tag::io, tag::control >() )
// *****************************************************************************
//  Constructor
//! \param[in] print Pretty printer
//! \param[in] cmdline Command line stack
//! \param[in,out] inputdeck Input deck stack where data is stored during
//!    parsing
// *****************************************************************************
{
  // Create InputDeck (a tagged tuple) to store parsed input
  ctr::New2InputDeck ideck( cmdline );

  //// Reset parser's output stream to that of print's. This is so that mild
  //// warnings emitted during parsing can be output using the pretty printer.
  //// Usually, errors and warnings are simply accumulated during parsing and
  //// printed during diagnostics after the parser has finished. Howver, in some
  //// special cases we can provide a more user-friendly message right during
  //// parsing since there is more information available to construct a more
  //// sensible message. This is done in e.g., tk::grm::store_option. Resetting
  //// the global g_print, to that of passed in as the constructor argument allows
  //// not to have to create a new pretty printer, but use the existing one.
  //tk::grm::g_print.reset( print.save() );

  //// Parse input file and populate the underlying tagged tuple
  //tao::pegtl::file_input<> in( m_filename );
  //tao::pegtl::parse< deck::read_file, tk::grm::action >( in, id );

  sol::state lua_deck;
  lua_deck.open_libraries(sol::lib::base);
  lua_deck.script_file(m_filename);

  // TODO: remove this cout
  //std::cout << " READ IN SUCCESSFUL " << std::endl;

  storeInputDeck(lua_deck["inciter"], ideck);

  // Echo errors and warnings accumulated during parsing
  diagnostics( print, std::vector< std::string >{}/*id.get< tag::error >()*/ );

  inputdeck = std::move( ideck );
}

void
LuaParser::storeInputDeck(
  const sol::table& lua_ideck,
  ctr::New2InputDeck& gideck )
// *****************************************************************************
//  Store lua inputdeck in custom struct
//! \param[in] lua_ideck Lua inputdeck parsed by sol2
//! \param[in,out] gideck Inciter's inputdeck storage
// *****************************************************************************
{
  // TODO: explore replacing storeIfSpecd() and storeOptIfSpecd with sol::get_or()

  storeIfSpecd< std::string >(
    lua_ideck, "title", gideck.get< newtag::title >(), "No title");

  // time stepping options
  // ---------------------------------------------------------------------------
  storeIfSpecd< uint64_t >(
    lua_ideck, "nstep", gideck.get< newtag::nstep >(),
    std::numeric_limits< uint64_t >::max());
  storeIfSpecd< tk::real >(
    lua_ideck, "term", gideck.get< newtag::term >(),
    std::numeric_limits< tk::real >::max());
  storeIfSpecd< tk::real >(
    lua_ideck, "t0", gideck.get< newtag::t0 >(), 0.0);
  storeIfSpecd< tk::real >(
    lua_ideck, "dt", gideck.get< newtag::dt >(), 0.0);
  storeIfSpecd< tk::real >(
    lua_ideck, "cfl", gideck.get< newtag::cfl >(), 0.0);
  storeIfSpecd< uint32_t >(
    lua_ideck, "ttyi", gideck.get< newtag::ttyi >(), 1);
  storeIfSpecd< bool >(
    lua_ideck, "steady_state", gideck.get< newtag::steady_state >(), false);
  storeIfSpecd< tk::real >(
    lua_ideck, "residual", gideck.get< newtag::residual >(), 1.0e-8);
  storeIfSpecd< uint32_t >(
    lua_ideck, "rescomp", gideck.get< newtag::rescomp >(), 1);

  // partitioning/reordering options
  // ---------------------------------------------------------------------------
  storeOptIfSpecd< tk::ctr::PartitioningAlgorithmType,
    tk::ctr::PartitioningAlgorithm >(
    lua_ideck, "partitioning", gideck.get< newtag::partitioning >(),
    tk::ctr::PartitioningAlgorithmType::RCB);
  storeIfSpecd< bool >(
    lua_ideck, "pelocal_reorder", gideck.get< newtag::pelocal_reorder >(),
    false);
  storeIfSpecd< bool >(
    lua_ideck, "operator_reorder", gideck.get< newtag::operator_reorder >(),
    false);

  // discretization scheme options
  // ---------------------------------------------------------------------------
  using inciter::ctr::SchemeType;
  storeOptIfSpecd< SchemeType, inciter::ctr::Scheme >(
    lua_ideck, "scheme", gideck.get< newtag::scheme >(), SchemeType::DiagCG);
  storeOptIfSpecd< inciter::ctr::LimiterType, inciter::ctr::Limiter >(
    lua_ideck, "limiter", gideck.get< newtag::limiter >(),
    inciter::ctr::LimiterType::NOLIMITER);
  storeIfSpecd< tk::real >(
    lua_ideck, "cweight", gideck.get< newtag::cweight >(), 1.0);
  storeIfSpecd< tk::real >(
    lua_ideck, "shock_detector_coeff",
    gideck.get< newtag::shock_detector_coeff >(), 1.0);
  storeIfSpecd< bool >(
    lua_ideck, "accuracy_test", gideck.get< newtag::accuracy_test >(), false);
  storeIfSpecd< bool >(
    lua_ideck, "limsol_projection", gideck.get< newtag::limsol_projection >(),
    true);
  storeIfSpecd< bool >(
    lua_ideck, "fct", gideck.get< newtag::fct >(), true);
  storeIfSpecd< bool >(
    lua_ideck, "fctclip", gideck.get< newtag::fctclip >(), false);
  storeIfSpecd< tk::real >(
    lua_ideck, "fcteps", gideck.get< newtag::fcteps >(),
    std::numeric_limits< tk::real >::epsilon());
  storeIfSpecd< tk::real >(
    lua_ideck, "ctau", gideck.get< newtag::ctau >(), 1.0);
  storeIfSpecd< bool >(
    lua_ideck, "sysfct", gideck.get< newtag::sysfct >(), false);
  storeVecIfSpecd< std::size_t >(
    lua_ideck, "sysfctvar", gideck.get< newtag::sysfctvar >(), {{0, 1, 2, 3, 4}});

  // configure solutions DOFs
  auto scheme = gideck.get< newtag::scheme >();
  auto& ndof = gideck.get< newtag::ndof >();
  auto& rdof = gideck.get< newtag::rdof >();
  ndof = rdof = 1;
  if (scheme == SchemeType::P0P1 || scheme == SchemeType::FV) {
    ndof = 1; rdof = 4;
  } else if (scheme == SchemeType::DGP1) {
    ndof = rdof = 4;
  } else if (scheme == SchemeType::DGP2) {
    ndof = rdof = 10;
  } else if (scheme == SchemeType::PDG) {
    ndof = rdof = 10;
    gideck.get< newtag::pref, newtag::pref >() = true;
  } else if (scheme != SchemeType::DG &&
      scheme != SchemeType::DiagCG &&
      scheme != SchemeType::ALECG &&
      scheme != SchemeType::OversetFE) {
    Throw("Scheme type not configured in configure_scheme");
  }

  // PDE options
  // ---------------------------------------------------------------------------

  // check transport
  if (lua_ideck["transport"].valid()) {
    gideck.get< newtag::pde >() = inciter::ctr::PDEType::TRANSPORT;
    storeIfSpecd< std::size_t >(
      lua_ideck["transport"], "ncomp",
      gideck.get< newtag::transport, newtag::ncomp >(), 1);
    storeOptIfSpecd< inciter::ctr::ProblemType, inciter::ctr::Problem >(
      lua_ideck["transport"], "problem",
      gideck.get< newtag::transport, newtag::problem >(),
      inciter::ctr::ProblemType::USER_DEFINED);
    storeIfSpecd< std::string >(
      lua_ideck, "depvar", gideck.get< newtag::depvar >(), "c");
    storeOptIfSpecd< inciter::ctr::FluxType, inciter::ctr::Flux >(
      lua_ideck, "flux", gideck.get< newtag::flux >(),
      inciter::ctr::FluxType::UPWIND);

    // store number of equations in PDE system
    gideck.get< newtag::ncomp >() =
      gideck.get< newtag::transport, newtag::ncomp >();
  }

  // check compflow
  if (lua_ideck["compflow"].valid()) {
    gideck.get< newtag::pde >() = inciter::ctr::PDEType::COMPFLOW;
    storeOptIfSpecd< inciter::ctr::ProblemType, inciter::ctr::Problem >(
      lua_ideck["compflow"], "problem",
      gideck.get< newtag::compflow, newtag::problem >(),
      inciter::ctr::ProblemType::USER_DEFINED);
    storeIfSpecd< std::string >(
      lua_ideck, "depvar", gideck.get< newtag::depvar >(), "a");
    storeOptIfSpecd< inciter::ctr::FluxType, inciter::ctr::Flux >(
      lua_ideck, "flux", gideck.get< newtag::flux >(),
      inciter::ctr::FluxType::HLLC);

    // store number of equations in PDE system
    gideck.get< newtag::ncomp >() = 5;
  }

  // check multimat
  if (lua_ideck["multimat"].valid()) {
    gideck.get< newtag::pde >() = inciter::ctr::PDEType::MULTIMAT;
    storeIfSpecd< std::size_t >(
      lua_ideck["multimat"], "nmat",
      gideck.get< newtag::multimat, newtag::nmat >(), 2);
    storeIfSpecd< uint64_t >(
      lua_ideck["multimat"], "prelax",
      gideck.get< newtag::multimat, newtag::prelax >(), 1);
    storeIfSpecd< tk::real >(
      lua_ideck["multimat"], "prelax_timescale",
      gideck.get< newtag::multimat, newtag::prelax_timescale >(), 0.25);
    storeIfSpecd< int >(
      lua_ideck["multimat"], "intsharp",
      gideck.get< newtag::multimat, newtag::intsharp >(), 0);
    storeIfSpecd< tk::real >(
      lua_ideck["multimat"], "intsharp_param",
      gideck.get< newtag::multimat, newtag::intsharp_param >(), 1.8);
    storeOptIfSpecd< inciter::ctr::ProblemType, inciter::ctr::Problem >(
      lua_ideck["multimat"], "problem",
      gideck.get< newtag::multimat, newtag::problem >(),
      inciter::ctr::ProblemType::USER_DEFINED);
    storeIfSpecd< std::string >(
      lua_ideck, "depvar", gideck.get< newtag::depvar >(), "a");
    storeOptIfSpecd< inciter::ctr::FluxType, inciter::ctr::Flux >(
      lua_ideck, "flux", gideck.get< newtag::flux >(),
      inciter::ctr::FluxType::AUSM);

    // number of equations in PDE system are determined based on materials
  }

  // physics
  // ---------------------------------------------------------------------------
  storeOptIfSpecd< inciter::ctr::PhysicsType, inciter::ctr::Physics >(
    lua_ideck, "physics", gideck.get< newtag::physics >(),
    inciter::ctr::PhysicsType::EULER);

  // Assemble material blocks
  // ---------------------------------------------------------------------------

  // solid counters
  std::size_t tmat(0), imatcntr(0), mtypei(0), isolcntr(0), isolidx(0);
  std::set< std::size_t > matidset;

  // material vector
  if (lua_ideck["material"].valid()) {

    // size material map vectors
    std::size_t nmat(1);
    if (gideck.get< newtag::pde >() == inciter::ctr::PDEType::MULTIMAT)
      nmat = gideck.get< newtag::multimat, newtag::nmat >();
    gideck.get< newtag::matidxmap, newtag::eosidx >().resize(nmat);
    gideck.get< newtag::matidxmap, newtag::matidx >().resize(nmat);
    gideck.get< newtag::matidxmap, newtag::solidx >().resize(nmat);

    // size material vector appropriately
    // size of the material vector is the number of distinct types of materials
    sol::table sol_mat = lua_ideck["material"];
    gideck.get< newtag::material >().resize(sol_mat.size());

    // store material properties
    for (std::size_t i=0; i<gideck.get< newtag::material >().size(); ++i) {
      auto& mati_deck = gideck.get< newtag::material >()[i];
      // eos
      storeOptIfSpecd< inciter::ctr::MaterialType, inciter::ctr::Material >(
        sol_mat[i+1], "eos", mati_deck.get< newtag::eos >(),
        inciter::ctr::MaterialType::STIFFENEDGAS);

      // material ids in this eos (default is for compflow i.e. single mat)
      storeVecIfSpecd< uint64_t >(
        sol_mat[i+1], "id", mati_deck.get< newtag::id >(),
        std::vector< uint64_t >(1,1));

      // Track total number of materials in multiple material blocks (eos's)
      tmat += mati_deck.get< newtag::id >().size();

      // Check for repeating user specified material ids
      for (auto midx : mati_deck.get< newtag::id >()) {
        if (!matidset.count(midx))
          matidset.insert(midx);
        else
          Throw("Repeating material id specified");
      }

      std::size_t ntype = mati_deck.get< newtag::id >().size();
      // cv
      if (!sol_mat[i+1]["cv"].valid())
        sol_mat[i+1]["cv"] = std::vector< tk::real >(ntype, 717.5);
      checkStoreMatProp(sol_mat[i+1], "cv", ntype, mati_deck.get< newtag::cv >());

      // reset solid-index
      isolidx = 0;

      // Stiffened-gas materials
      if (mati_deck.get< newtag::eos >() ==
        inciter::ctr::MaterialType::STIFFENEDGAS) {
        // gamma
        checkStoreMatProp(sol_mat[i+1], "gamma", ntype,
          mati_deck.get< newtag::gamma >());

        // pstiff
        if (!sol_mat[i+1]["pstiff"].valid())
          sol_mat[i+1]["pstiff"] = std::vector< tk::real >(ntype, 0.0);
        checkStoreMatProp(sol_mat[i+1], "pstiff", ntype,
          mati_deck.get< newtag::pstiff >());
      }
      // Small-shear solid materials
      else if (mati_deck.get< newtag::eos >() ==
        inciter::ctr::MaterialType::SMALLSHEARSOLID) {
        // gamma
        checkStoreMatProp(sol_mat[i+1], "gamma", ntype,
          mati_deck.get< newtag::gamma >());

        // pstiff
        if (!sol_mat[i+1]["pstiff"].valid())
          sol_mat[i+1]["pstiff"] = std::vector< tk::real >(ntype, 0.0);
        checkStoreMatProp(sol_mat[i+1], "pstiff", ntype,
          mati_deck.get< newtag::pstiff >());

        // mu
        checkStoreMatProp(sol_mat[i+1], "mu", ntype,
          mati_deck.get< newtag::mu >());

        // add to solid-counter
        ++isolcntr;
        // assign solid-counter value to solid-index
        isolidx = isolcntr;
      }
      // JWL materials
      else if (mati_deck.get< newtag::eos >() == inciter::ctr::MaterialType::JWL) {
        // w_gru
        checkStoreMatProp(sol_mat[i+1], "w_gru", ntype,
          mati_deck.get< newtag::w_gru >());

        // a_jwl
        checkStoreMatProp(sol_mat[i+1], "A_jwl", ntype,
          mati_deck.get< newtag::A_jwl >());

        // b_jwl
        checkStoreMatProp(sol_mat[i+1], "B_jwl", ntype,
          mati_deck.get< newtag::B_jwl >());

        // R1_jwl
        checkStoreMatProp(sol_mat[i+1], "R1_jwl", ntype,
          mati_deck.get< newtag::R1_jwl >());

        // R2_jwl
        checkStoreMatProp(sol_mat[i+1], "R2_jwl", ntype,
          mati_deck.get< newtag::R2_jwl >());

        // rho0_jwl
        checkStoreMatProp(sol_mat[i+1], "rho0_jwl", ntype,
          mati_deck.get< newtag::rho0_jwl >());

        // de_jwl
        checkStoreMatProp(sol_mat[i+1], "de_jwl", ntype,
          mati_deck.get< newtag::de_jwl >());

        // Pr_jwl
        checkStoreMatProp(sol_mat[i+1], "Pr_jwl", ntype,
          mati_deck.get< newtag::Pr_jwl >());

        // rhor_jwl
        if (sol_mat[i+1]["rhor_jwl"].valid()) {
          checkStoreMatProp(sol_mat[i+1], "rhor_jwl", ntype,
            mati_deck.get< newtag::rhor_jwl >());
        }
        // Tr_jwl
        else if (sol_mat[i+1]["Tr_jwl"].valid()) {
          checkStoreMatProp(sol_mat[i+1], "Tr_jwl", ntype,
            mati_deck.get< newtag::Tr_jwl >());
        }
        else
          Throw("Either reference density or reference temperature must be "
            "specified for JWL equation of state (EOS).");
      }

      // Generate mapping between material index and eos parameter index
      auto& eosmap = gideck.get< newtag::matidxmap, newtag::eosidx >();
      auto& idxmap = gideck.get< newtag::matidxmap, newtag::matidx >();
      auto& solidxmap = gideck.get< newtag::matidxmap, newtag::solidx >();
      for (auto midx : mati_deck.get< newtag::id >()) {
        midx -= 1;
        eosmap[midx] = mtypei;
        idxmap[midx] = imatcntr;
        solidxmap[midx] = isolidx;
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
    if (gideck.get< newtag::pde >() == inciter::ctr::PDEType::MULTIMAT) {
      auto ntot = nmat + nmat + 3 + nmat;
      // if solid EOS, add components
      const auto& solidx = gideck.get< newtag::matidxmap, newtag::solidx >();
      for (std::size_t i=0; i<solidx.size(); ++i) {
        if (solidx[i] > 0)
          ntot += 9;
      }
      gideck.get< newtag::ncomp >() = ntot;
    }
  }

  // Field output block
  // ---------------------------------------------------------------------------
  if (lua_ideck["field_output"].valid()) {
    auto& fo_deck = gideck.get< newtag::field_output >();

    // interval iteration
    storeIfSpecd< uint32_t >(
      lua_ideck["field_output"], "interval",
      fo_deck.get< newtag::iter_interval >(),
      std::numeric_limits< uint32_t >::max());

    // interval physical time
    storeIfSpecd< tk::real >(
      lua_ideck["field_output"], "time_interval",
      fo_deck.get< newtag::time_interval >(),
      std::numeric_limits< tk::real >::max());

    // interval time range
    storeVecIfSpecd< tk::real >(
      lua_ideck["field_output"], "time_range",
      fo_deck.get< newtag::time_range >(), {});

    // refined mesh field output
    storeIfSpecd< bool >(
      lua_ideck["field_output"], "refined", fo_deck.get< newtag::refined >(),
      false);

    // filetype
    storeOptIfSpecd< tk::ctr::FieldFileType, tk::ctr::FieldFile >(
      lua_ideck["field_output"], "filetype", fo_deck.get< newtag::filetype >(),
      tk::ctr::FieldFileType::EXODUSII);

    // sidesets for field output
    storeVecIfSpecd< uint64_t >(
      lua_ideck["field_output"], "sideset", fo_deck.get< newtag::sideset >(), {});

    // element variables
    storeVecIfSpecd< std::string >(
      lua_ideck["field_output"], "elemvar", fo_deck.get< newtag::elemvar >(), {});

    // node variables
    storeVecIfSpecd< std::string >(
      lua_ideck["field_output"], "nodevar", fo_deck.get< newtag::nodevar >(), {});
  }
  else {
    // TODO: remove double-specification of defaults
    auto& fo_deck = gideck.get< newtag::field_output >();
    fo_deck.get< newtag::iter_interval >() =
      std::numeric_limits< uint32_t >::max();
    fo_deck.get< newtag::time_interval >() =
      std::numeric_limits< tk::real >::max();
    fo_deck.get< newtag::time_range >() = {};
    fo_deck.get< newtag::refined >() = false;
    fo_deck.get< newtag::filetype >() = tk::ctr::FieldFileType::EXODUSII;
    fo_deck.get< newtag::sideset >() = {};
    fo_deck.get< newtag::elemvar >() = {};
    fo_deck.get< newtag::nodevar >() = {};
  }

  // Diagnostics output block
  // ---------------------------------------------------------------------------
  if (lua_ideck["diagnostics"].valid()) {
    auto& diag_deck = gideck.get< newtag::diagnostics >();

    // interval iteration
    storeIfSpecd< uint32_t >(
      lua_ideck["diagnostics"], "interval",
      diag_deck.get< newtag::iter_interval >(), 1);

    // error norm
    storeOptIfSpecd< tk::ctr::ErrorType, tk::ctr::Error >(
      lua_ideck["diagnostics"], "error", diag_deck.get< newtag::error >(),
      tk::ctr::ErrorType::L2);

    // float format
    storeOptIfSpecd< tk::ctr::TxtFloatFormatType, tk::ctr::TxtFloatFormat >(
      lua_ideck["diagnostics"], "format", diag_deck.get< newtag::format >(),
      tk::ctr::TxtFloatFormatType::DEFAULT);

    // precision
    storeIfSpecd< uint32_t >(
      lua_ideck["diagnostics"], "precision",
      diag_deck.get< newtag::precision >(), std::cout.precision());
  }
  else {
    // TODO: remove double-specification of defaults
    auto& diag_deck = gideck.get< newtag::diagnostics >();
    diag_deck.get< newtag::iter_interval >() = 1;
    diag_deck.get< newtag::error >() = tk::ctr::ErrorType::L2;
    diag_deck.get< newtag::format >() = tk::ctr::TxtFloatFormatType::DEFAULT;
    diag_deck.get< newtag::precision >() = std::cout.precision();
  }

  // History output block
  // ---------------------------------------------------------------------------
  if (lua_ideck["history_output"].valid()) {
    auto& hist_deck = gideck.get< newtag::history_output >();

    // interval iteration
    storeIfSpecd< uint32_t >(
      lua_ideck["history_output"], "interval",
      hist_deck.get< newtag::iter_interval >(),
      std::numeric_limits< uint32_t >::max());

    // interval time
    storeIfSpecd< tk::real >(
      lua_ideck["history_output"], "time_interval",
      hist_deck.get< newtag::time_interval >(),
      std::numeric_limits< tk::real >::max());

    // interval time range
    storeVecIfSpecd< tk::real >(
      lua_ideck["history_output"], "time_range",
      hist_deck.get< newtag::time_range >(), {});

    // point probes
    if (lua_ideck["history_output"]["point"].valid()) {
      const sol::table& sol_pt = lua_ideck["history_output"]["point"];
      hist_deck.get< newtag::point >().resize(sol_pt.size());

      for (std::size_t i=0; i<hist_deck.get< newtag::point >().size(); ++i) {
        auto& pti = hist_deck.get< newtag::point >()[i];
        storeIfSpecd< std::string >(
          sol_pt[i+1], "id", pti.get< newtag::id >(), "p");
        storeVecIfSpecd< tk::real >(
          sol_pt[i+1], "coord", pti.get< newtag::coord >(), {});
      }
    }

    // precision
    storeIfSpecd< uint32_t >(
      lua_ideck["history_output"], "precision",
      hist_deck.get< newtag::precision >(), std::cout.precision());

    // error check point
    for (std::size_t i=0; i<hist_deck.get< newtag::point >().size(); ++i) {
      if (hist_deck.get< newtag::point >()[i].get< newtag::coord >().size() != 3)
      Throw("Three reals required for point coordinates in history_output.");
    }
  }
  else {
    // TODO: remove double-specification of defaults
    auto& hist_deck = gideck.get< newtag::history_output >();
    hist_deck.get< newtag::iter_interval >() =
      std::numeric_limits< uint32_t >::max();
    hist_deck.get< newtag::time_interval >() =
      std::numeric_limits< tk::real >::max();
    hist_deck.get< newtag::time_range >() = {};
    hist_deck.get< newtag::precision >() = std::cout.precision();
    hist_deck.get< newtag::point >().resize(0);
  }

  // ALE block
  // ---------------------------------------------------------------------------
  gideck.get< newtag::ale, newtag::ale >() = false;
  if (lua_ideck["ale"].valid()) {
    auto& ale_deck = gideck.get< newtag::ale >();
    ale_deck.get< newtag::ale >() = true;

    // Mesh velocity smoother
    storeOptIfSpecd< inciter::ctr::MeshVelocitySmootherType,
      inciter::ctr::MeshVelocitySmoother >(lua_ideck["ale"], "smoother",
      ale_deck.get< newtag::smoother >(), inciter::ctr::MeshVelocitySmootherType::NONE);

    // Mesh velocity
    storeOptIfSpecd< inciter::ctr::MeshVelocityType,
      inciter::ctr::MeshVelocity >(lua_ideck["ale"], "mesh_velocity",
      ale_deck.get< newtag::mesh_velocity >(), inciter::ctr::MeshVelocityType::SINE);

    // Mesh motion direction
    storeVecIfSpecd< std::size_t >(lua_ideck["ale"], "mesh_motion",
      ale_deck.get< newtag::mesh_motion >(), { 0, 1, 2 });

    // Mesh force
    storeVecIfSpecd< tk::real >(lua_ideck["ale"], "meshforce",
      ale_deck.get< newtag::meshforce >(), { 0, 0, 0, 0 });

    // Move sidesets with user defined function
    if (lua_ideck["ale"]["move"].valid()) {
      const sol::table& sol_mv = lua_ideck["ale"]["move"];
      ale_deck.get< newtag::move >().resize(sol_mv.size());

      for (std::size_t i=0; i<ale_deck.get< newtag::move >().size(); ++i) {
        auto& mvi = ale_deck.get< newtag::move >()[i];
        storeOptIfSpecd< tk::ctr::UserTableType, tk::ctr::UserTable >(
          sol_mv[i+1], "fntype", mvi.get< newtag::fntype >(),
          tk::ctr::UserTableType::POSITION);
        storeIfSpecd< uint64_t >(
          sol_mv[i+1], "sideset", mvi.get< newtag::sideset >(), 0);
        storeVecIfSpecd< tk::real >(
          sol_mv[i+1], "fn", mvi.get< newtag::fn >(), {});

        // error checking on user-def function
        if (mvi.get< newtag::fn >().size() % 4 != 0)
          Throw("Incomplete user-defined function for ALE sideset movement. An "
          "R->R^3 function is expected, the number of descrete entries must be "
          "divisible by 4: one 'column' for the abscissa, and 3 for the "
          "ordinate.");
      }
    }

    // dv-CFL
    storeIfSpecd< tk::real >(lua_ideck["ale"], "dvcfl", ale_deck.get< newtag::dvcfl >(),
      0.01);

    // Vorticity multiplier
    storeIfSpecd< tk::real >(lua_ideck["ale"], "vortmult",
      ale_deck.get< newtag::vortmult >(), 0.0);

    // Mesh velocity max iterations
    storeIfSpecd< std::size_t >(lua_ideck["ale"], "maxit",
      ale_deck.get< newtag::maxit >(), 5);

    // Mesh velocity max iterations
    storeIfSpecd< tk::real >(lua_ideck["ale"], "tolerance",
      ale_deck.get< newtag::tolerance >(), 1e-2);
  }

  // AMR block
  // ---------------------------------------------------------------------------
  gideck.get< newtag::amr, newtag::amr >() = false;
  if (lua_ideck["amr"].valid()) {
    auto& amr_deck = gideck.get< newtag::amr >();
    amr_deck.get< newtag::amr >() = true;

    // Initial refinement toggle
    storeIfSpecd< bool >(lua_ideck["amr"], "t0ref", amr_deck.get< newtag::t0ref >(),
      false);

    // Mesh refinement during time-stepping toggle
    storeIfSpecd< bool >(lua_ideck["amr"], "dtref", amr_deck.get< newtag::dtref >(),
      false);

    // Uniform mesh refinement during time-stepping toggle
    storeIfSpecd< bool >(lua_ideck["amr"], "dtref_uniform",
      amr_deck.get< newtag::dtref_uniform >(), false);

    // Mesh refinement frequency during time-stepping toggle
    storeIfSpecd< std::size_t >(lua_ideck["amr"], "dtfreq",
      amr_deck.get< newtag::dtfreq >(), 3);

    // Maximum AMR levels
    storeIfSpecd< std::size_t >(lua_ideck["amr"], "maxlevels",
      amr_deck.get< newtag::maxlevels >(), 2);

    // Initial AMR steps
    storeOptVecIfSpecd< inciter::ctr::AMRInitialType, inciter::ctr::AMRInitial >(
      lua_ideck["amr"], "initial", amr_deck.get< newtag::initial >(), {});

    // Initial AMR coordinate based
    if (lua_ideck["amr"]["coords"].valid()) {
      auto rmax = std::numeric_limits< tk::real >::max() / 100;

      storeIfSpecd< tk::real >(lua_ideck["amr"]["coords"], "xminus",
        amr_deck.get< newtag::coords, newtag::xminus >(), rmax);
      storeIfSpecd< tk::real >(lua_ideck["amr"]["coords"], "xplus",
        amr_deck.get< newtag::coords, newtag::xplus >(), -rmax);

      storeIfSpecd< tk::real >(lua_ideck["amr"]["coords"], "yminus",
        amr_deck.get< newtag::coords, newtag::yminus >(), rmax);
      storeIfSpecd< tk::real >(lua_ideck["amr"]["coords"], "yplus",
        amr_deck.get< newtag::coords, newtag::yplus >(), -rmax);

      storeIfSpecd< tk::real >(lua_ideck["amr"]["coords"], "zminus",
        amr_deck.get< newtag::coords, newtag::zminus >(), rmax);
      storeIfSpecd< tk::real >(lua_ideck["amr"]["coords"], "zplus",
        amr_deck.get< newtag::coords, newtag::zplus >(), -rmax);
    }

    // Initial AMR edgelist based
    storeVecIfSpecd< std::size_t >(lua_ideck["amr"], "edgelist",
      amr_deck.get< newtag::edgelist >(), {});
    if (amr_deck.get< newtag::edgelist >().size() % 2 != 0)
      Throw("The number of edge-nodes, marking edges as pairs of nodes, used "
        "for explicit tagging of edges for initial mesh refineoment, is odd "
        "(it must be even).");

    // Error type for AMR
    storeOptIfSpecd< inciter::ctr::AMRErrorType, inciter::ctr::AMRError >(
      lua_ideck["amr"], "error", amr_deck.get< newtag::error >(),
      inciter::ctr::AMRErrorType::JUMP);

    // Tolerances for refine/de-refine
    storeIfSpecd< tk::real >(lua_ideck["amr"], "tol_refine",
      amr_deck.get< newtag::tol_refine >(), 0.2);
    storeIfSpecd< tk::real >(lua_ideck["amr"], "tol_derefine",
      amr_deck.get< newtag::tol_derefine >(), 0.05);
  }

  // p-refinement block
  // ---------------------------------------------------------------------------
  gideck.get< newtag::pref, newtag::pref >() = false;
  if (lua_ideck["pref"].valid()) {
    auto& pref_deck = gideck.get< newtag::pref >();
    pref_deck.get< newtag::pref >() = true;

    // p-ref indicator type
    storeOptIfSpecd< inciter::ctr::PrefIndicatorType,
      inciter::ctr::PrefIndicator >(lua_ideck["pref"], "indicator",
      pref_deck.get< newtag::indicator >(),
      inciter::ctr::PrefIndicatorType::SPECTRAL_DECAY);

    // p-ref max degrees-of-freedom per cell
    storeIfSpecd< std::size_t >(lua_ideck["pref"], "ndofmax",
      pref_deck.get< newtag::ndofmax >(), 10);

    // p-ref tolerance
    storeIfSpecd< tk::real >(lua_ideck["pref"], "tolref",
      pref_deck.get< newtag::tolref >(), 0.5);

    // error checking on the tolerance
    if (pref_deck.get< newtag::tolref >() < 0.0 || pref_deck.get< newtag::tolref >() > 1.0)
      Throw("The p-refinement tolerance must be a real number "
        "between 0.0 and 1.0, both inclusive.");
  }

  // Mesh specification block (for overset)
  // ---------------------------------------------------------------------------
  if (lua_ideck["mesh"].valid()) {
    const sol::table& lua_mesh = lua_ideck["mesh"];
    auto& mesh_deck = gideck.get< newtag::mesh >();
    mesh_deck.resize(lua_mesh.size());

    for (std::size_t i=0; i<mesh_deck.size(); ++i) {
      // filename
      storeIfSpecd< std::string >(lua_mesh[i+1], "filename",
        mesh_deck[i].get< newtag::filename >(), "");

      // location
      storeVecIfSpecd< tk::real >(lua_mesh[i+1], "location",
        mesh_deck[i].get< newtag::location >(), {0.0, 0.0, 0.0});
      if (mesh_deck[i].get< newtag::location >().size() != 3)
        Throw("Mesh location requires 3 coordinates.");

      // orientation
      storeVecIfSpecd< tk::real >(lua_mesh[i+1], "orientation",
        mesh_deck[i].get< newtag::orientation >(), {0.0, 0.0, 0.0});
      if (mesh_deck[i].get< newtag::orientation >().size() != 3)
        Throw("Mesh orientation requires 3 rotation angles.");

      // velocity
      storeVecIfSpecd< tk::real >(lua_mesh[i+1], "velocity",
        mesh_deck[i].get< newtag::velocity >(), {0.0, 0.0, 0.0});
      if (mesh_deck[i].get< newtag::velocity >().size() != 3)
        Throw("Mesh velocity requires 3 components.");
    }
  }
  else {
    // TODO: remove double-specification of defaults
    auto& mesh_deck = gideck.get< newtag::mesh >();
    mesh_deck.resize(1);
    mesh_deck[0].get< newtag::filename >() = "";
    mesh_deck[0].get< newtag::location >() = {0.0, 0.0, 0.0};
    mesh_deck[0].get< newtag::orientation >() = {0.0, 0.0, 0.0};
    mesh_deck[0].get< newtag::velocity >() = {0.0, 0.0, 0.0};
  }

  // Boundary conditions block
  // ---------------------------------------------------------------------------
  if (lua_ideck["bc"].valid()) {
    std::set< std::size_t > totalmesh;
    const sol::table& sol_bc = lua_ideck["bc"];
    auto& bc_deck = gideck.get< newtag::bc >();
    bc_deck.resize(sol_bc.size());

    for (std::size_t i=0; i<bc_deck.size(); ++i) {
      storeVecIfSpecd< std::size_t >(sol_bc[i+1], "mesh",
        bc_deck[i].get< newtag::mesh >(), {1});
      // collect meshes for error checking
      totalmesh.insert(bc_deck[i].get< newtag::mesh >().begin(),
        bc_deck[i].get< newtag::mesh >().end());

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "dirichlet",
        bc_deck[i].get< newtag::dirichlet >(), {});

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "symmetry",
        bc_deck[i].get< newtag::symmetry >(), {});

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "inlet",
        bc_deck[i].get< newtag::inlet >(), {});

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "outlet",
        bc_deck[i].get< newtag::outlet >(), {});

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "farfield",
        bc_deck[i].get< newtag::farfield >(), {});

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "extrapolate",
        bc_deck[i].get< newtag::extrapolate >(), {});

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "stag",
        bc_deck[i].get< newtag::stag >(), {});

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "skip",
        bc_deck[i].get< newtag::skip >(), {});

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "sponge",
        bc_deck[i].get< newtag::sponge >(), {});

      // Time-dependent BC
      if (sol_bc[i+1]["timedep"].valid()) {
        storeVecIfSpecd< uint64_t >(sol_bc[i+1]["timedep"], "sideset",
          bc_deck[i].get< newtag::timedep, newtag::sideset >(), {});
        storeVecIfSpecd< tk::real >(sol_bc[i+1]["timedep"], "fn",
          bc_deck[i].get< newtag::timedep, newtag::fn >(), {});

        // error checking on user-def function
        if (bc_deck[i].get< newtag::timedep, newtag::fn >().size() % 6 != 0)
          Throw("Incomplete user-defined function for time-dependent BC. An "
          "R->R^5 function is expected, the number of descrete entries must be "
          "divisible by 6: one 'column' for the abscissa, and 5 for the ordinate.");
      }

      // Stagnation point
      storeVecIfSpecd< tk::real >(sol_bc[i+1], "point",
        bc_deck[i].get< newtag::point >(), {0.0, 0.0, 0.0});
      if (bc_deck[i].get< newtag::point >().size() != 3)
        Throw("BC point requires 3 coordinates.");

      // Stagnation radius
      storeIfSpecd< tk::real >(sol_bc[i+1], "radius", bc_deck[i].get< newtag::radius >(),
        0.0);

      // Velocity for inlet/farfield
      storeVecIfSpecd< tk::real >(sol_bc[i+1], "velocity",
        bc_deck[i].get< newtag::velocity >(), {0.0, 0.0, 0.0});
      if (bc_deck[i].get< newtag::velocity >().size() != 3)
        Throw("BC velocity requires 3 components.");

      // Pressure for inlet/outlet/farfield
      storeIfSpecd< tk::real >(sol_bc[i+1], "pressure",
        bc_deck[i].get< newtag::pressure >(), 0.0);

      // Density for inlet/outlet/farfield
      storeIfSpecd< tk::real >(sol_bc[i+1], "density",
        bc_deck[i].get< newtag::density >(), 0.0);
    }

    // error checking on number of meshes
    if (totalmesh.size() != gideck.get< newtag::mesh >().size())
      Throw("Total meshes (" + std::to_string(gideck.get< newtag::mesh >().size()) +
        ") not equal to the meshes on which BC's are specified (" +
        std::to_string(totalmesh.size()));

    // error checking on mesh ids
    std::size_t ic(1);
    for (const auto& im : totalmesh) {
      if (im != ic) Throw("Non-contiguous mesh ids in BC-mesh");
      ++ic;
    }
  }
  // error checking for unspecified BC's
  else
    Throw("No boundary conditions specified in input file.");

  // Initial condition block
  // ---------------------------------------------------------------------------
  if (lua_ideck["ic"].valid()) {
    auto& ic_deck = gideck.get< newtag::ic >();

    // background IC values
    storeIfSpecd< std::size_t >(lua_ideck["ic"], "materialid",
      ic_deck.get< newtag::materialid >(), 1);

    storeIfSpecd< tk::real >(lua_ideck["ic"], "pressure",
      ic_deck.get< newtag::pressure >(), 0.0);

    storeIfSpecd< tk::real >(lua_ideck["ic"], "temperature",
      ic_deck.get< newtag::temperature >(), 0.0);

    storeVecIfSpecd< tk::real >(lua_ideck["ic"], "velocity",
      ic_deck.get< newtag::velocity >(), {0.0, 0.0, 0.0});
    if (ic_deck.get< newtag::velocity >().size() != 3)
      Throw("Velocity in IC requires 3 components.");

    // IC box
    if (lua_ideck["ic"]["box"].valid()) {
      const sol::table& lua_box = lua_ideck["ic"]["box"];
      auto& box_deck = ic_deck.get< newtag::box >();
      box_deck.resize(lua_box.size());

      for (std::size_t i=0; i<box_deck.size(); ++i) {
        storeIfSpecd< std::size_t >(lua_box[i+1], "materialid",
          box_deck[i].get< newtag::materialid >(), 1);

        storeIfSpecd< tk::real >(lua_box[i+1], "volume",
          box_deck[i].get< newtag::volume >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "mass",
          box_deck[i].get< newtag::mass >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "density",
          box_deck[i].get< newtag::density >(), 0.0);

        storeVecIfSpecd< tk::real >(lua_box[i+1], "velocity",
          box_deck[i].get< newtag::velocity >(), {0.0, 0.0, 0.0});
        if (box_deck[i].get< newtag::velocity >().size() != 3)
          Throw("Velocity in IC box requires 3 components.");

        storeIfSpecd< tk::real >(lua_box[i+1], "pressure",
          box_deck[i].get< newtag::pressure >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "energy",
          box_deck[i].get< newtag::energy >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "energy_content",
          box_deck[i].get< newtag::energy_content >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "temperature",
          box_deck[i].get< newtag::temperature >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "xmin",
          box_deck[i].get< newtag::xmin >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "xmax",
          box_deck[i].get< newtag::xmax >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "ymin",
          box_deck[i].get< newtag::ymin >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "ymax",
          box_deck[i].get< newtag::ymax >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "zmin",
          box_deck[i].get< newtag::zmin >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "zmax",
          box_deck[i].get< newtag::zmax >(), 0.0);

        storeVecIfSpecd< tk::real >(lua_box[i+1], "orientation",
          box_deck[i].get< newtag::orientation >(), {0.0, 0.0, 0.0});
        if (box_deck[i].get< newtag::orientation >().size() != 3)
          Throw("Orientation in IC box requires 3 rotation angles.");

        storeOptIfSpecd< inciter::ctr::InitiateType, inciter::ctr::Initiate >(
          lua_box[i+1], "initiate", box_deck[i].get< newtag::initiate >(),
          inciter::ctr::InitiateType::IMPULSE);

        storeVecIfSpecd< tk::real >(lua_box[i+1], "point",
          box_deck[i].get< newtag::point >(), {0.0, 0.0, 0.0});
        if (box_deck[i].get< newtag::point >().size() != 3)
          Throw("Point in IC box requires 3 coordinates.");

        storeIfSpecd< tk::real >(lua_box[i+1], "init_time",
          box_deck[i].get< newtag::init_time >(), 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "front_width",
          box_deck[i].get< newtag::front_width >(), 0.0);
      }
    }

    // IC mesh-block
    if (lua_ideck["ic"]["meshblock"].valid()) {
      const sol::table& lua_meshblock = lua_ideck["ic"]["meshblock"];
      auto& mblk_deck = ic_deck.get< newtag::meshblock >();
      mblk_deck.resize(lua_meshblock.size());

      for (std::size_t i=0; i<mblk_deck.size(); ++i) {
        storeIfSpecd< std::size_t >(lua_meshblock[i+1], "blockid",
          mblk_deck[i].get< newtag::blockid >(), 0);
        if (mblk_deck[i].get< newtag::blockid >() == 0)
          Throw("Each IC mesh block must specify the mesh block id.");

        storeIfSpecd< std::size_t >(lua_meshblock[i+1], "materialid",
          mblk_deck[i].get< newtag::materialid >(), 1);

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "energy_content",
          mblk_deck[i].get< newtag::energy_content >(), 0.0);

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "volume",
          mblk_deck[i].get< newtag::volume >(), 0.0);
        if (mblk_deck[i].get< newtag::energy_content >() > 0.0 &&
          mblk_deck[i].get< newtag::volume >() < 1e-12)
          Throw("Mesh block volume must be specified, if energy content is "
            "used to initialize block");

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "mass",
          mblk_deck[i].get< newtag::mass >(), 0.0);

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "density",
          mblk_deck[i].get< newtag::density >(), 0.0);

        storeVecIfSpecd< tk::real >(lua_meshblock[i+1], "velocity",
          mblk_deck[i].get< newtag::velocity >(), {0.0, 0.0, 0.0});
        if (mblk_deck[i].get< newtag::velocity >().size() != 3)
          Throw("Velocity in IC meshblock requires 3 components.");

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "pressure",
          mblk_deck[i].get< newtag::pressure >(), 0.0);

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "energy",
          mblk_deck[i].get< newtag::energy >(), 0.0);

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "temperature",
          mblk_deck[i].get< newtag::temperature >(), 0.0);

        storeOptIfSpecd< inciter::ctr::InitiateType, inciter::ctr::Initiate >(
          lua_meshblock[i+1], "initiate", mblk_deck[i].get< newtag::initiate >(),
          inciter::ctr::InitiateType::IMPULSE);

        storeVecIfSpecd< tk::real >(lua_meshblock[i+1], "point",
          mblk_deck[i].get< newtag::point >(), {0.0, 0.0, 0.0});
        if (mblk_deck[i].get< newtag::point >().size() != 3)
          Throw("Point in IC meshblock requires 3 coordinates.");

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "init_time",
          mblk_deck[i].get< newtag::init_time >(), 0.0);

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "front_width",
          mblk_deck[i].get< newtag::front_width >(), 0.0);
      }
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // Testing
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------

  //std::cout << gideck.get< newtag::title >() << std::endl;
  //std::cout << "cfl = " << gideck.get< newtag::cfl >() << std::endl;
  //std::cout << "scheme = " << static_cast< std::size_t >(gideck.get< newtag::scheme >()) << std::endl;
  //auto nmat = gideck.get< newtag::multimat, newtag::nmat >();
  //std::cout << "nmat " << nmat << std::endl;
  //std::cout << "ncomp " << gideck.get< newtag::ncomp >() << std::endl;

  //for (std::size_t k=0; k<nmat; ++k)
  //  std::cout << gideck.get< newtag::matidxmap, newtag::eosidx >()[k] << ", ";
  //std::cout << std::endl;
  //for (std::size_t k=0; k<nmat; ++k)
  //  std::cout << gideck.get< newtag::matidxmap, newtag::matidx >()[k] << ", ";
  //std::cout << std::endl;
  //for (std::size_t k=0; k<nmat; ++k)
  //  std::cout << gideck.get< newtag::matidxmap, newtag::solidx >()[k] << ", ";
  //std::cout << std::endl;

  //inciter::ctr::Material matclass;
  //for (std::size_t i=0; i<gideck.get< newtag::material >().size(); ++i) {
  //  const auto& mati_deck = gideck.get< newtag::material >()[i];
  //  std::cout << "MatType "
  //    << matclass.name(mati_deck.get< newtag::eos >()) << std::endl;
  //  for (std::size_t k=0; k<mati_deck.get< newtag::id >().size(); ++k) {
  //    std::cout << mati_deck.get< newtag::id >()[k] << std::endl;

  //    //std::cout << ", cv = " << mati_deck.get< newtag::cv >()[k];

  //    //if (mati_deck.get< newtag::eos >() == inciter::ctr::MaterialType::STIFFENEDGAS) {
  //    //  std::cout << ", gamma = " << mati_deck.get< newtag::gamma >()[k] << ", "
  //    //    << ", pstiff = " << mati_deck.get< newtag::pstiff >()[k] << ", ";
  //    //}

  //    //if (mati_deck.get< newtag::eos >() == inciter::ctr::MaterialType::SMALLSHEARSOLID) {
  //    //  std::cout << ", gamma = " << mati_deck.get< newtag::gamma >()[k] << ", "
  //    //    << ", pstiff = " << mati_deck.get< newtag::pstiff >()[k] << ", "
  //    //    << ", mu = " << mati_deck.get< newtag::mu >()[k] << ", ";
  //    //}

  //    //if (mati_deck.get< newtag::eos >() == inciter::ctr::MaterialType::JWL) {
  //    //  std::cout << ", w_gru = " << mati_deck.get< newtag::w_gru >()[k] << ", "
  //    //    << ", A_jwl = " << mati_deck.get< newtag::A_jwl >()[k] << ", "
  //    //    << ", B_jwl = " << mati_deck.get< newtag::B_jwl >()[k] << ", "
  //    //    << ", R1_jwl = " << mati_deck.get< newtag::R1_jwl >()[k] << ", "
  //    //    << ", R2_jwl = " << mati_deck.get< newtag::R2_jwl >()[k] << ", "
  //    //    << ", rho0_jwl = " << mati_deck.get< newtag::rho0_jwl >()[k] << ", "
  //    //    << ", de_jwl = " << mati_deck.get< newtag::de_jwl >()[k] << ", "
  //    //    << ", Tr_jwl = " << mati_deck.get< newtag::Tr_jwl >()[k] << ", "
  //    //    << ", Pr_jwl = " << mati_deck.get< newtag::Pr_jwl >()[k] << ", ";
  //    //}

  //    //std::cout << std::endl;
  //  }
  //}

  //std::cout << " diag iter: " << gideck.get< newtag::diagnostics, newtag::iter_interval >() << std::endl;
  //std::cout << " diag norm: " << static_cast< std::size_t >(gideck.get< newtag::diagnostics, newtag::error >()) << std::endl;
  //std::cout << " diag format: " << static_cast< std::size_t >(gideck.get< newtag::diagnostics, newtag::format >()) << std::endl;
  //std::cout << " diag prec: " << gideck.get< newtag::diagnostics, newtag::precision >() << std::endl;

  //std::cout << " F-O/P ref: "<< gideck.get< newtag::field_output, newtag::refined >() << std::endl;
  //std::cout << " F-O/P time: "<< gideck.get< newtag::field_output, newtag::time_interval >() << std::endl;
  //std::cout << " F-O/P iter: "<< gideck.get< newtag::field_output, newtag::iter_interval >() << std::endl;
  //std::cout << " F-O/P filetype: "<< static_cast< std::size_t >(gideck.get< newtag::field_output, newtag::filetype >()) << std::endl;
  //std::cout << " F-O/P sidesets: ";
  //for (std::size_t i=0; i< gideck.get< newtag::field_output, newtag::sideset >().size(); ++i)
  //  std::cout << gideck.get< newtag::field_output, newtag::sideset >()[i] << ", ";
  //std::cout << std::endl;
  //std::cout << " F-O/P elemvars: ";
  //for (std::size_t i=0; i< gideck.get< newtag::field_output, newtag::elemvar >().size(); ++i)
  //  std::cout << gideck.get< newtag::field_output, newtag::elemvar >()[i] << ", ";
  //std::cout << std::endl;
  //std::cout << " F-O/P nodevars: ";
  //for (std::size_t i=0; i< gideck.get< newtag::field_output, newtag::nodevar >().size(); ++i)
  //  std::cout << gideck.get< newtag::field_output, newtag::nodevar >()[i] << ", ";
  //std::cout << std::endl;

  //std::cout << " hist time: " << gideck.get< newtag::history_output, newtag::time_interval >() << std::endl;
  //std::cout << " hist iter: " << gideck.get< newtag::history_output, newtag::iter_interval >() << std::endl;
  //std::cout << " hist range: ";
  //for (std::size_t i=0; i< gideck.get< newtag::history_output, newtag::time_range >().size(); ++i)
  //  std::cout << gideck.get< newtag::history_output, newtag::time_range >()[i] << ", ";
  //std::cout << std::endl;
  //std::cout << " hist point: " << std::endl;
  //for (std::size_t i=0; i< gideck.get< newtag::history_output, newtag::point >().size(); ++i) {
  //  std::cout << gideck.get< newtag::history_output, newtag::point >()[i].get< newtag::id >() << ": ";
  //  for (std::size_t j=0; j< gideck.get< newtag::history_output, newtag::point >()[i].get< newtag::coord >().size(); ++j)
  //    std::cout << gideck.get< newtag::history_output, newtag::point >()[i].get< newtag::coord >()[j] << ", ";
  //  std::cout << std::endl;
  //}
  //std::cout << std::endl;

  //std::cout << " ALE " << gideck.get< newtag::ale >().get< newtag::ale >() << std::endl;
  //for (std::size_t i=0; i< gideck.get< newtag::ale >().get< newtag::move >().size(); ++i) {
  //  const auto& mvi = gideck.get< newtag::ale >().get< newtag::move >()[i];
  //  std::cout << static_cast<std::size_t>(mvi.get< newtag::fntype >()) << ": "
  //    << mvi.get< newtag::sideset >() << std::endl;
  //  for (std::size_t j=0; j< mvi.get< newtag::fn >().size(); ++j)
  //    std::cout << mvi.get< newtag::fn >()[j] << ", ";
  //  std::cout << std::endl;
  //}

  //std::cout << " AMR " << gideck.get< newtag::amr, newtag::amr >() << std::endl;
  //std::cout << " t0ref " << gideck.get< newtag::amr, newtag::t0ref >() << std::endl;
  //std::cout << " dtref " << gideck.get< newtag::amr, newtag::dtref >() << std::endl;
  //std::cout << " dtref_uniform " << gideck.get< newtag::amr, newtag::dtref_uniform >() << std::endl;
  //std::cout << " dtfreq " << gideck.get< newtag::amr, newtag::dtfreq >() << std::endl;
  //std::cout << " maxlevels " << gideck.get< newtag::amr, newtag::maxlevels >() << std::endl;
  //std::cout << " initial ";
  //for (std::size_t i=0; i<gideck.get< newtag::amr, newtag::initial >().size(); ++i)
  //  std::cout << static_cast< std::size_t >(gideck.get< newtag::amr, newtag::initial >()[i]) << ", ";
  //std::cout << std::endl;
  //std::cout << " error type " << static_cast< std::size_t >(gideck.get< newtag::amr, newtag::error >())
  //  << std::endl;

  //std::cout << " pref " << gideck.get< newtag::pref, newtag::pref >() << std::endl;
  //std::cout << " indicator type " <<
  //  static_cast< std::size_t >(gideck.get< newtag::pref, newtag::indicator >()) << std::endl;
  //std::cout << " ndofmax " << gideck.get< newtag::pref, newtag::ndofmax >() << std::endl;
  //std::cout << " tolref " << gideck.get< newtag::pref, newtag::tolref >() << std::endl;

  //std::cout << " Meshes: " << std::endl;
  //for (std::size_t j=0; j< gideck.get< newtag::mesh >().size(); ++j) {
  //  std::cout << " " << gideck.get< newtag::mesh >()[j].get< newtag::filename >() << std::endl;
  //  std::cout << " location: ";
  //  for (std::size_t i=0; i<gideck.get< newtag::mesh >()[j].get< newtag::location >().size(); ++i)
  //    std::cout << gideck.get< newtag::mesh >()[j].get< newtag::location >()[i] << ", ";
  //  std::cout << std::endl;
  //  std::cout << " orientation: ";
  //  for (std::size_t i=0; i<gideck.get< newtag::mesh >()[j].get< newtag::orientation >().size(); ++i)
  //    std::cout << gideck.get< newtag::mesh >()[j].get< newtag::orientation >()[i] << ", ";
  //  std::cout << std::endl;
  //  std::cout << " velocity: ";
  //  for (std::size_t i=0; i<gideck.get< newtag::mesh >()[j].get< newtag::velocity >().size(); ++i)
  //    std::cout << gideck.get< newtag::mesh >()[j].get< newtag::velocity >()[i] << ", ";
  //  std::cout << std::endl;
  //}

  //std::cout << " BCs: " << std::endl;
  //for (std::size_t j=0; j< gideck.get< newtag::bc >().size(); ++j) {
  //  const auto& bcj = gideck.get< newtag::bc >()[j];
  //  std::cout << " bc-meshes: ";
  //  for (std::size_t i=0; i<bcj.get< newtag::mesh >().size(); ++i)
  //    std::cout << bcj.get< newtag::mesh >()[i] << ", ";
  //  std::cout << std::endl;
  //  std::cout << " dirichlet: ";
  //  for (std::size_t i=0; i<bcj.get< newtag::dirichlet >().size(); ++i)
  //    std::cout << bcj.get< newtag::dirichlet >()[i] << ", ";
  //  std::cout << std::endl;
  //  std::cout << " symmetry: ";
  //  for (std::size_t i=0; i<bcj.get< newtag::symmetry >().size(); ++i)
  //    std::cout << bcj.get< newtag::symmetry >()[i] << ", ";
  //  std::cout << std::endl;
  //  std::cout << " farfield: ";
  //  for (std::size_t i=0; i<bcj.get< newtag::farfield >().size(); ++i)
  //    std::cout << bcj.get< newtag::farfield >()[i] << ", ";
  //  std::cout << std::endl;
  //  std::cout << " inlet: ";
  //  for (std::size_t i=0; i<bcj.get< newtag::inlet >().size(); ++i)
  //    std::cout << bcj.get< newtag::inlet >()[i] << ", ";
  //  std::cout << std::endl;
  //  std::cout << " pressure " << bcj.get< newtag::pressure >() << std::endl;
  //  std::cout << " density " << bcj.get< newtag::density >() << std::endl;
  //  std::cout << " velocity " << bcj.get< newtag::velocity >()[0] << ", "
  //    << bcj.get< newtag::velocity >()[1] << ", "
  //    << bcj.get< newtag::velocity >()[2] << std::endl;
  //  std::cout << " timedep: ";
  //  for (std::size_t i=0; i<bcj.get< newtag::timedep, newtag::sideset >().size(); ++i)
  //    std::cout << bcj.get< newtag::timedep, newtag::sideset >()[i] << ", ";
  //  std::cout << std::endl;
  //  for (std::size_t i=0; i<bcj.get< newtag::timedep, newtag::fn >().size(); ++i)
  //    std::cout << bcj.get< newtag::timedep, newtag::fn >()[i] << ", ";
  //  std::cout << std::endl;
  //}

  //std::cout << " ICs: " << std::endl;
  //std::cout << " matid " << gideck.get< newtag::ic, newtag::materialid >() << std::endl;
  //std::cout << " pressure " << gideck.get< newtag::ic, newtag::pressure >() << std::endl;
  //std::cout << " temperature " << gideck.get< newtag::ic, newtag::temperature >() << std::endl;
  //std::cout << " velocity " << gideck.get< newtag::ic, newtag::velocity >()[0] << ", "
  //  << gideck.get< newtag::ic, newtag::velocity >()[1] << ", "
  //  << gideck.get< newtag::ic, newtag::velocity >()[2] << std::endl;
  //std::cout << " IC box: " << std::endl;
  //for (std::size_t i=0; i<gideck.get< newtag::ic, newtag::box >().size(); ++i) {
  //  const auto& bxi = gideck.get< newtag::ic, newtag::box >()[i];
  //  std::cout << "  box " << i << std::endl;
  //  std::cout << " energy_content " <<
  //    bxi.get< newtag::energy_content >() << std::endl;
  //  std::cout << " xmin " <<
  //    bxi.get< newtag::xmin >() << std::endl;
  //  std::cout << " xmax " <<
  //    bxi.get< newtag::xmax >() << std::endl;
  //  std::cout << " ymin " <<
  //    bxi.get< newtag::ymin >() << std::endl;
  //  std::cout << " ymax " <<
  //    bxi.get< newtag::ymax >() << std::endl;
  //  std::cout << " zmin " <<
  //    bxi.get< newtag::zmin >() << std::endl;
  //  std::cout << " zmax " <<
  //    bxi.get< newtag::zmax >() << std::endl;
  //  std::cout << " orientation "
  //    << bxi.get< newtag::orientation >()[0] << ", "
  //    << bxi.get< newtag::orientation >()[1] << ", "
  //    << bxi.get< newtag::orientation >()[2] << ", "
  //    << std::endl;
  //  std::cout << " initiate " << static_cast< std::size_t >(bxi.get< newtag::initiate >()) << std::endl;
  //  std::cout << " point "
  //    << bxi.get< newtag::point >()[0] << ", "
  //    << bxi.get< newtag::point >()[1] << ", "
  //    << bxi.get< newtag::point >()[2] << ", "
  //    << std::endl;
  //  std::cout << " init_time " <<
  //    bxi.get< newtag::init_time >() << std::endl;
  //  std::cout << " front_width " <<
  //    bxi.get< newtag::front_width >() << std::endl;
  //}
  //std::cout << " IC block: " << std::endl;
  //for (std::size_t i=0; i<gideck.get< newtag::ic, newtag::meshblock >().size(); ++i) {
  //  const auto& mblki = gideck.get< newtag::ic, newtag::meshblock >()[i];
  //  std::cout << "  meshblock " << i << std::endl;
  //  std::cout << " blockid " <<
  //    mblki.get< newtag::blockid >() << std::endl;
  //  std::cout << " materid " <<
  //    mblki.get< newtag::materialid >() << std::endl;
  //  std::cout << " energy_content " <<
  //    mblki.get< newtag::energy_content >() << std::endl;
  //  std::cout << " initiate " << static_cast< std::size_t >(mblki.get< newtag::initiate >()) << std::endl;
  //  std::cout << " point "
  //    << mblki.get< newtag::point >()[0] << ", "
  //    << mblki.get< newtag::point >()[1] << ", "
  //    << mblki.get< newtag::point >()[2] << ", "
  //    << std::endl;
  //  std::cout << " init_time " <<
  //    mblki.get< newtag::init_time >() << std::endl;
  //  std::cout << " front_width " <<
  //    mblki.get< newtag::front_width >() << std::endl;
  //}
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
    Throw("Material property " + key + " not specified");
  if (sol::table(table[key]).size() != vecsize)
    Throw("Incorrect number of " + key + "'s specified. Expected " +
      std::to_string(vecsize));

  // store values from table to inputdeck
  storeVecIfSpecd< tk::real >(table, key, storage,
    std::vector< tk::real >(vecsize, 0.0));
}
