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
                      ctr::NewInputDeck& inputdeck ) :
  FileParser( cmdline.get< tag::io, tag::control >() )
// *****************************************************************************
//  Constructor
//! \param[in] print Pretty printer
//! \param[in] cmdline Command line stack
//! \param[in,out] inputdeck Input deck stack where data is stored during
//!    parsing
// *****************************************************************************
{
  //// Create InputDeck (a tagged tuple) to store parsed input
  //ctr::InputDeck id( cmdline );

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

  inputdeck = storeInputDeck(lua_deck["inciter"]);

  // Echo errors and warnings accumulated during parsing
  diagnostics( print, std::vector< std::string >{}/*id.get< tag::error >()*/ );
}

inciter::ctr::NewInputDeck
LuaParser::storeInputDeck(
  const sol::table& lua_ideck )
// *****************************************************************************
//  Store lua inputdeck in custom struct
//! \param[in] lua_ideck Lua inputdeck parsed by sol2
// *****************************************************************************
{
  inciter::ctr::NewInputDeck gideck;

  // TODO: explore replacing storeIfSpecd() and storeOptIfSpecd with sol::get_or()

  storeIfSpecd< std::string >(
    lua_ideck, "title", gideck.title.data, "No title");

  // time stepping options
  // ---------------------------------------------------------------------------
  storeIfSpecd< uint64_t >(
    lua_ideck, "nstep", gideck.nstep.data,
    std::numeric_limits< uint64_t >::max());
  storeIfSpecd< tk::real >(
    lua_ideck, "term", gideck.term.data,
    std::numeric_limits< tk::real >::max());
  storeIfSpecd< tk::real >(
    lua_ideck, "t0", gideck.t0.data, 0.0);
  storeIfSpecd< tk::real >(
    lua_ideck, "dt", gideck.dt.data, 0.0);
  storeIfSpecd< tk::real >(
    lua_ideck, "cfl", gideck.cfl.data, 0.0);
  storeIfSpecd< uint32_t >(
    lua_ideck, "ttyi", gideck.ttyi.data, 1);
  storeIfSpecd< bool >(
    lua_ideck, "steady_state", gideck.steady_state.data, false);
  storeIfSpecd< tk::real >(
    lua_ideck, "residual", gideck.residual.data, 1.0e-8);
  storeIfSpecd< uint32_t >(
    lua_ideck, "rescomp", gideck.rescomp.data, 1);

  // partitioning/reordering options
  // ---------------------------------------------------------------------------
  storeOptIfSpecd< tk::ctr::PartitioningAlgorithmType,
    tk::ctr::PartitioningAlgorithm >(
    lua_ideck, "partitioning", gideck.partitioning.data,
    tk::ctr::PartitioningAlgorithmType::RCB);
  storeIfSpecd< bool >(
    lua_ideck, "pelocal_reorder", gideck.pelocal_reorder.data, false);
  storeIfSpecd< bool >(
    lua_ideck, "operator_reorder", gideck.operator_reorder.data, false);

  // discretization scheme options
  // ---------------------------------------------------------------------------
  using inciter::ctr::SchemeType;
  storeOptIfSpecd< SchemeType, inciter::ctr::Scheme >(
    lua_ideck, "scheme", gideck.scheme.data, SchemeType::DiagCG);
  storeOptIfSpecd< inciter::ctr::LimiterType, inciter::ctr::Limiter >(
    lua_ideck, "limiter", gideck.limiter.data,
    inciter::ctr::LimiterType::NOLIMITER);
  storeIfSpecd< tk::real >(
    lua_ideck, "cweight", gideck.cweight.data, 1.0);
  storeIfSpecd< tk::real >(
    lua_ideck, "shock_detector_coeff", gideck.shock_detector_coeff.data, 1.0);
  storeIfSpecd< bool >(
    lua_ideck, "accuracy_test", gideck.accuracy_test.data, false);
  storeIfSpecd< bool >(
    lua_ideck, "limsol_projection", gideck.limsol_projection.data, true);
  storeIfSpecd< bool >(
    lua_ideck, "fct", gideck.fct.data, true);
  storeIfSpecd< bool >(
    lua_ideck, "fctclip", gideck.fctclip.data, false);
  storeIfSpecd< tk::real >(
    lua_ideck, "fcteps", gideck.fcteps.data,
    std::numeric_limits< tk::real >::epsilon());
  storeIfSpecd< tk::real >(
    lua_ideck, "ctau", gideck.ctau.data, 1.0);
  storeIfSpecd< bool >(
    lua_ideck, "sysfct", gideck.sysfct.data, false);
  storeVecIfSpecd< std::size_t >(
    lua_ideck, "sysfctvar", gideck.sysfctvar.data, {{0, 1, 2, 3, 4}});

  // configure solutions DOFs
  auto scheme = gideck.scheme.data;
  auto& ndof = gideck.ndof.data;
  auto& rdof = gideck.rdof.data;
  ndof = rdof = 1;
  if (scheme == SchemeType::P0P1 || scheme == SchemeType::FV) {
    ndof = 1; rdof = 4;
  } else if (scheme == SchemeType::DGP1) {
    ndof = rdof = 4;
  } else if (scheme == SchemeType::DGP2) {
    ndof = rdof = 10;
  } else if (scheme == SchemeType::PDG) {
    ndof = rdof = 10;
    gideck.pref.pref.data = true;
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
    gideck.pde.data = inciter::ctr::PDEType::TRANSPORT;
    storeIfSpecd< std::size_t >(
      lua_ideck["transport"], "ncomp", gideck.transport.ncomp.data, 1);
    storeIfSpecd< std::string >(
      lua_ideck["transport"], "depvar", gideck.transport.depvar.data, "c");
    storeOptIfSpecd< inciter::ctr::ProblemType, inciter::ctr::Problem >(
      lua_ideck["transport"], "problem", gideck.transport.problem.data,
      inciter::ctr::ProblemType::USER_DEFINED);
  storeOptIfSpecd< inciter::ctr::FluxType, inciter::ctr::Flux >(
    lua_ideck, "flux", gideck.flux.data, inciter::ctr::FluxType::UPWIND);

    // store number of equations in PDE system
    gideck.ncomp.data = gideck.transport.ncomp.data;
  }

  // check compflow
  if (lua_ideck["compflow"].valid()) {
    gideck.pde.data = inciter::ctr::PDEType::COMPFLOW;
    storeOptIfSpecd< inciter::ctr::ProblemType, inciter::ctr::Problem >(
      lua_ideck["compflow"], "problem", gideck.compflow.problem.data,
      inciter::ctr::ProblemType::USER_DEFINED);
    storeOptIfSpecd< inciter::ctr::FluxType, inciter::ctr::Flux >(
      lua_ideck, "flux", gideck.flux.data, inciter::ctr::FluxType::HLLC);

    // store number of equations in PDE system
    gideck.ncomp.data = 5;
  }

  // check multimat
  if (lua_ideck["multimat"].valid()) {
    gideck.pde.data = inciter::ctr::PDEType::MULTIMAT;
    storeIfSpecd< std::size_t >(
      lua_ideck["multimat"], "nmat", gideck.multimat.nmat.data, 2);
    storeIfSpecd< uint64_t >(
      lua_ideck["multimat"], "prelax", gideck.multimat.prelax.data, 1);
    storeIfSpecd< tk::real >(
      lua_ideck["multimat"], "prelax_timescale",
      gideck.multimat.prelax_timescale.data, 0.25);
    storeIfSpecd< int >(
      lua_ideck["multimat"], "intsharp", gideck.multimat.intsharp.data, 0);
    storeIfSpecd< tk::real >(
      lua_ideck["multimat"], "intsharp_param",
      gideck.multimat.intsharp_param.data, 1.8);
    storeOptIfSpecd< inciter::ctr::ProblemType, inciter::ctr::Problem >(
      lua_ideck["multimat"], "problem", gideck.multimat.problem.data,
      inciter::ctr::ProblemType::USER_DEFINED);
    storeOptIfSpecd< inciter::ctr::FluxType, inciter::ctr::Flux >(
      lua_ideck, "flux", gideck.flux.data, inciter::ctr::FluxType::AUSM);

    // number of equations in PDE system are determined based on materials
  }

  // physics
  // ---------------------------------------------------------------------------
  storeOptIfSpecd< inciter::ctr::PhysicsType, inciter::ctr::Physics >(
    lua_ideck, "physics", gideck.physics.data, inciter::ctr::PhysicsType::EULER);

  // Assemble material blocks
  // ---------------------------------------------------------------------------

  // solid counters
  std::size_t tmat(0), imatcntr(0), mtypei(0), isolcntr(0), isolidx(0);
  std::set< std::size_t > matidset;

  // material vector
  if (lua_ideck["material"].valid()) {

    // size material map vectors
    std::size_t nmat(1);
    if (gideck.pde.data == inciter::ctr::PDEType::MULTIMAT)
      nmat = gideck.multimat.nmat.data;
    gideck.matidxmap.eosidx.resize(nmat);
    gideck.matidxmap.matidx.resize(nmat);
    gideck.matidxmap.solidx.resize(nmat);

    // size material vector appropriately
    // size of the material vector is the number of distinct types of materials
    sol::table sol_mat = lua_ideck["material"];
    gideck.material.resize(sol_mat.size());

    // store material properties
    for (std::size_t i=0; i<gideck.material.size(); ++i) {
      // eos
      storeOptIfSpecd< inciter::ctr::MaterialType, inciter::ctr::Material >(
        sol_mat[i+1], "eos", gideck.material[i].eos.data,
        inciter::ctr::MaterialType::STIFFENEDGAS);

      // material ids in this eos (default is for compflow i.e. single mat)
      storeVecIfSpecd< uint64_t >(
        sol_mat[i+1], "id", gideck.material[i].id.data,
        std::vector< uint64_t >(1,1));

      // Track total number of materials in multiple material blocks (eos's)
      tmat += gideck.material[i].id.data.size();

      // Check for repeating user specified material ids
      for (auto midx : gideck.material[i].id.data) {
        if (!matidset.count(midx))
          matidset.insert(midx);
        else
          Throw("Repeating material id specified");
      }

      std::size_t ntype = gideck.material[i].id.data.size();
      // cv
      if (!sol_mat[i+1]["cv"].valid())
        sol_mat[i+1]["cv"] = std::vector< tk::real >(ntype, 717.5);
      checkStoreMatProp(sol_mat[i+1], "cv", ntype, gideck.material[i].cv.data);

      // reset solid-index
      isolidx = 0;

      // Stiffened-gas materials
      if (gideck.material[i].eos.data ==
        inciter::ctr::MaterialType::STIFFENEDGAS) {
        // gamma
        checkStoreMatProp(sol_mat[i+1], "gamma", ntype,
          gideck.material[i].gamma.data);

        // pstiff
        if (!sol_mat[i+1]["pstiff"].valid())
          sol_mat[i+1]["pstiff"] = std::vector< tk::real >(ntype, 0.0);
        checkStoreMatProp(sol_mat[i+1], "pstiff", ntype,
          gideck.material[i].pstiff.data);
      }
      // Small-shear solid materials
      else if (gideck.material[i].eos.data ==
        inciter::ctr::MaterialType::SMALLSHEARSOLID) {
        // gamma
        checkStoreMatProp(sol_mat[i+1], "gamma", ntype,
          gideck.material[i].gamma.data);

        // pstiff
        if (!sol_mat[i+1]["pstiff"].valid())
          sol_mat[i+1]["pstiff"] = std::vector< tk::real >(ntype, 0.0);
        checkStoreMatProp(sol_mat[i+1], "pstiff", ntype,
          gideck.material[i].pstiff.data);

        // mu
        checkStoreMatProp(sol_mat[i+1], "mu", ntype,
          gideck.material[i].mu.data);

        // add to solid-counter
        ++isolcntr;
        // assign solid-counter value to solid-index
        isolidx = isolcntr;
      }
      // JWL materials
      else if (gideck.material[i].eos.data == inciter::ctr::MaterialType::JWL) {
        // w_gru
        checkStoreMatProp(sol_mat[i+1], "w_gru", ntype,
          gideck.material[i].w_gru.data);

        // a_jwl
        checkStoreMatProp(sol_mat[i+1], "A_jwl", ntype,
          gideck.material[i].A_jwl.data);

        // b_jwl
        checkStoreMatProp(sol_mat[i+1], "B_jwl", ntype,
          gideck.material[i].B_jwl.data);

        // R1_jwl
        checkStoreMatProp(sol_mat[i+1], "R1_jwl", ntype,
          gideck.material[i].R1_jwl.data);

        // R2_jwl
        checkStoreMatProp(sol_mat[i+1], "R2_jwl", ntype,
          gideck.material[i].R2_jwl.data);

        // rho0_jwl
        checkStoreMatProp(sol_mat[i+1], "rho0_jwl", ntype,
          gideck.material[i].rho0_jwl.data);

        // de_jwl
        checkStoreMatProp(sol_mat[i+1], "de_jwl", ntype,
          gideck.material[i].de_jwl.data);

        // Pr_jwl
        checkStoreMatProp(sol_mat[i+1], "Pr_jwl", ntype,
          gideck.material[i].Pr_jwl.data);

        // rhor_jwl
        if (sol_mat[i+1]["rhor_jwl"].valid()) {
          checkStoreMatProp(sol_mat[i+1], "rhor_jwl", ntype,
            gideck.material[i].rhor_jwl.data);
        }
        // Tr_jwl
        else if (sol_mat[i+1]["Tr_jwl"].valid()) {
          checkStoreMatProp(sol_mat[i+1], "Tr_jwl", ntype,
            gideck.material[i].Tr_jwl.data);
        }
        else
          Throw("Either reference density or reference temperature must be "
            "specified for JWL equation of state (EOS).");
      }

      // Generate mapping between material index and eos parameter index
      auto& eosmap = gideck.matidxmap.eosidx;
      auto& idxmap = gideck.matidxmap.matidx;
      auto& solidxmap = gideck.matidxmap.solidx;
      for (auto midx : gideck.material[i].id.data) {
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
    if (gideck.pde.data == inciter::ctr::PDEType::MULTIMAT) {
      auto ntot = nmat + nmat + 3 + nmat;
      // if solid EOS, add components
      for (std::size_t i=0; i<gideck.matidxmap.solidx.size(); ++i) {
        if (gideck.matidxmap.solidx[i] > 0)
          ntot += 9;
      }
      gideck.ncomp.data = ntot;
    }
  }

  // Field output block
  // ---------------------------------------------------------------------------
  if (lua_ideck["field_output"].valid()) {

    // interval iteration
    storeIfSpecd< uint32_t >(
      lua_ideck["field_output"], "interval",
      gideck.field_output.iter_interval.data,
      std::numeric_limits< uint32_t >::max());

    // interval physical time
    storeIfSpecd< tk::real >(
      lua_ideck["field_output"], "time_interval",
      gideck.field_output.time_interval.data,
      std::numeric_limits< tk::real >::max());

    // interval time range
    storeVecIfSpecd< tk::real >(
      lua_ideck["field_output"], "time_range",
      gideck.field_output.time_range.data, {});

    // refined mesh field output
    storeIfSpecd< bool >(
      lua_ideck["field_output"], "refined", gideck.field_output.refined.data,
      false);

    // filetype
    storeOptIfSpecd< tk::ctr::FieldFileType, tk::ctr::FieldFile >(
      lua_ideck["field_output"], "filetype", gideck.field_output.filetype.data,
      tk::ctr::FieldFileType::EXODUSII);

    // sidesets for field output
    storeVecIfSpecd< uint64_t >(
      lua_ideck["field_output"], "sideset", gideck.field_output.sideset.data, {});

    // element variables
    storeVecIfSpecd< std::string >(
      lua_ideck["field_output"], "elemvar", gideck.field_output.elemvar.data, {});

    // node variables
    storeVecIfSpecd< std::string >(
      lua_ideck["field_output"], "nodevar", gideck.field_output.nodevar.data, {});
  }
  else {
    // TODO: remove double-specification of defaults
    gideck.field_output.iter_interval.data =
      std::numeric_limits< uint32_t >::max();
    gideck.field_output.time_interval.data =
      std::numeric_limits< tk::real >::max();
    gideck.field_output.time_range.data = {};
    gideck.field_output.refined.data = false;
    gideck.field_output.filetype.data = tk::ctr::FieldFileType::EXODUSII;
    gideck.field_output.sideset.data = {};
    gideck.field_output.elemvar.data = {};
    gideck.field_output.nodevar.data = {};
  }

  // Diagnostics output block
  // ---------------------------------------------------------------------------
  if (lua_ideck["diagnostics"].valid()) {

    // interval iteration
    storeIfSpecd< uint32_t >(
      lua_ideck["diagnostics"], "interval",
      gideck.diagnostics.iter_interval.data, 1);

    // error norm
    storeOptIfSpecd< tk::ctr::ErrorType, tk::ctr::Error >(
      lua_ideck["diagnostics"], "error",
      gideck.diagnostics.error.data, tk::ctr::ErrorType::L2);

    // float format
    storeOptIfSpecd< tk::ctr::TxtFloatFormatType, tk::ctr::TxtFloatFormat >(
      lua_ideck["diagnostics"], "format",
      gideck.diagnostics.format.data, tk::ctr::TxtFloatFormatType::DEFAULT);

    // precision
    storeIfSpecd< uint32_t >(
      lua_ideck["diagnostics"], "precision",
      gideck.diagnostics.precision.data, std::cout.precision());
  }
  else {
    // TODO: remove double-specification of defaults
    gideck.diagnostics.iter_interval.data = 1;
    gideck.diagnostics.error.data = tk::ctr::ErrorType::L2;
    gideck.diagnostics.format.data = tk::ctr::TxtFloatFormatType::DEFAULT;
    gideck.diagnostics.precision.data = std::cout.precision();
  }

  // History output block
  // ---------------------------------------------------------------------------
  if (lua_ideck["history_output"].valid()) {

    // interval iteration
    storeIfSpecd< uint32_t >(
      lua_ideck["history_output"], "interval",
      gideck.history_output.iter_interval.data,
      std::numeric_limits< uint32_t >::max());

    // interval time
    storeIfSpecd< tk::real >(
      lua_ideck["history_output"], "time_interval",
      gideck.history_output.time_interval.data,
      std::numeric_limits< tk::real >::max());

    // interval time range
    storeVecIfSpecd< tk::real >(
      lua_ideck["history_output"], "time_range",
      gideck.history_output.time_range.data, {});

    // point probes
    if (lua_ideck["history_output"]["point"].valid()) {
      const sol::table& sol_pt = lua_ideck["history_output"]["point"];
      gideck.history_output.point.resize(sol_pt.size());

      for (std::size_t i=0; i<gideck.history_output.point.size(); ++i) {
        storeIfSpecd< std::string >(
          sol_pt[i+1], "id", gideck.history_output.point[i].id.data, "p");
        storeVecIfSpecd< tk::real >(
          sol_pt[i+1], "coord", gideck.history_output.point[i].coord.data, {});
      }
    }

    // precision
    storeIfSpecd< uint32_t >(
      lua_ideck["history_output"], "precision",
      gideck.history_output.precision.data, std::cout.precision());

    // error check point
    for (std::size_t i=0; i<gideck.history_output.point.size(); ++i) {
      if (gideck.history_output.point[i].coord.data.size() != 3)
      Throw("Three reals required for point coordinates in history_output.");
    }
  }
  else {
    // TODO: remove double-specification of defaults
    gideck.history_output.iter_interval.data =
      std::numeric_limits< uint32_t >::max();
    gideck.history_output.time_interval.data =
      std::numeric_limits< tk::real >::max();
    gideck.history_output.time_range.data = {};
    gideck.history_output.precision.data = std::cout.precision();
    gideck.history_output.point.resize(0);
  }

  // ALE block
  // ---------------------------------------------------------------------------
  gideck.ale.ale.data = false;
  if (lua_ideck["ale"].valid()) {
    gideck.ale.ale.data = true;

    // Mesh velocity smoother
    storeOptIfSpecd< inciter::ctr::MeshVelocitySmootherType,
      inciter::ctr::MeshVelocitySmoother >(lua_ideck["ale"], "smoother",
      gideck.ale.smoother.data, inciter::ctr::MeshVelocitySmootherType::NONE);

    // Mesh velocity
    storeOptIfSpecd< inciter::ctr::MeshVelocityType,
      inciter::ctr::MeshVelocity >(lua_ideck["ale"], "mesh_velocity",
      gideck.ale.mesh_velocity.data, inciter::ctr::MeshVelocityType::SINE);

    // Mesh motion direction
    storeVecIfSpecd< std::size_t >(lua_ideck["ale"], "mesh_motion",
      gideck.ale.mesh_motion.data, { 0, 1, 2 });

    // Mesh force
    storeVecIfSpecd< tk::real >(lua_ideck["ale"], "meshforce",
      gideck.ale.meshforce.data, { 0, 0, 0, 0 });

    // Move sidesets with user defined function
    if (lua_ideck["ale"]["move"].valid()) {
      const sol::table& sol_mv = lua_ideck["ale"]["move"];
      gideck.ale.move.resize(sol_mv.size());

      for (std::size_t i=0; i<gideck.ale.move.size(); ++i) {
        storeOptIfSpecd< tk::ctr::UserTableType, tk::ctr::UserTable >(
          sol_mv[i+1], "fntype", gideck.ale.move[i].fntype.data,
          tk::ctr::UserTableType::POSITION);
        storeIfSpecd< uint64_t >(
          sol_mv[i+1], "sideset", gideck.ale.move[i].sideset.data, 0);
        storeVecIfSpecd< tk::real >(
          sol_mv[i+1], "fn", gideck.ale.move[i].fn.data, {});

        // error checking on user-def function
        if (gideck.ale.move[i].fn.data.size() % 4 != 0)
          Throw("Incomplete user-defined function for ALE sideset movement. An "
          "R->R^3 function is expected, the number of descrete entries must be "
          "divisible by 4: one 'column' for the abscissa, and 3 for the "
          "ordinate.");
      }
    }

    // dv-CFL
    storeIfSpecd< tk::real >(lua_ideck["ale"], "dvcfl", gideck.ale.dvcfl.data,
      0.01);

    // Vorticity multiplier
    storeIfSpecd< tk::real >(lua_ideck["ale"], "vortmult",
      gideck.ale.vortmult.data, 0.0);

    // Mesh velocity max iterations
    storeIfSpecd< std::size_t >(lua_ideck["ale"], "maxit",
      gideck.ale.maxit.data, 5);

    // Mesh velocity max iterations
    storeIfSpecd< tk::real >(lua_ideck["ale"], "tolerance",
      gideck.ale.tolerance.data, 1e-2);
  }

  // AMR block
  // ---------------------------------------------------------------------------
  gideck.amr.amr.data = false;
  if (lua_ideck["amr"].valid()) {
    gideck.amr.amr.data = true;

    // Initial refinement toggle
    storeIfSpecd< bool >(lua_ideck["amr"], "t0ref", gideck.amr.t0ref.data,
      false);

    // Mesh refinement during time-stepping toggle
    storeIfSpecd< bool >(lua_ideck["amr"], "dtref", gideck.amr.dtref.data,
      false);

    // Uniform mesh refinement during time-stepping toggle
    storeIfSpecd< bool >(lua_ideck["amr"], "dtref_uniform",
      gideck.amr.dtref_uniform.data, false);

    // Mesh refinement frequency during time-stepping toggle
    storeIfSpecd< std::size_t >(lua_ideck["amr"], "dtfreq",
      gideck.amr.dtfreq.data, 3);

    // Maximum AMR levels
    storeIfSpecd< std::size_t >(lua_ideck["amr"], "maxlevels",
      gideck.amr.maxlevels.data, 2);

    // Initial AMR steps
    storeOptVecIfSpecd< inciter::ctr::AMRInitialType, inciter::ctr::AMRInitial >(
      lua_ideck["amr"], "initial", gideck.amr.initial.data, {});

    // Initial AMR coordinate based
    if (lua_ideck["amr"]["coords"].valid()) {
      auto rmax = std::numeric_limits< tk::real >::max() / 100;

      storeIfSpecd< tk::real >(lua_ideck["amr"]["coords"], "xminus",
        gideck.amr.coords.xminus.data, rmax);
      storeIfSpecd< tk::real >(lua_ideck["amr"]["coords"], "xplus",
        gideck.amr.coords.xplus.data, -rmax);

      storeIfSpecd< tk::real >(lua_ideck["amr"]["coords"], "yminus",
        gideck.amr.coords.yminus.data, rmax);
      storeIfSpecd< tk::real >(lua_ideck["amr"]["coords"], "yplus",
        gideck.amr.coords.yplus.data, -rmax);

      storeIfSpecd< tk::real >(lua_ideck["amr"]["coords"], "zminus",
        gideck.amr.coords.zminus.data, rmax);
      storeIfSpecd< tk::real >(lua_ideck["amr"]["coords"], "zplus",
        gideck.amr.coords.zplus.data, -rmax);
    }

    // Initial AMR edgelist based
    storeVecIfSpecd< std::size_t >(lua_ideck["amr"], "edgelist",
      gideck.amr.edgelist.data, {});

    // Error type for AMR
    storeOptIfSpecd< inciter::ctr::AMRErrorType, inciter::ctr::AMRError >(
      lua_ideck["amr"], "error", gideck.amr.error.data,
      inciter::ctr::AMRErrorType::JUMP);

    // Tolerances for refine/de-refine
    storeIfSpecd< tk::real >(lua_ideck["amr"], "tol_refine",
      gideck.amr.tol_refine.data, 0.2);
    storeIfSpecd< tk::real >(lua_ideck["amr"], "tol_derefine",
      gideck.amr.tol_derefine.data, 0.05);
  }

  // p-refinement block
  // ---------------------------------------------------------------------------
  gideck.pref.pref.data = false;
  if (lua_ideck["pref"].valid()) {
    gideck.pref.pref.data = true;

    // p-ref indicator type
    storeOptIfSpecd< inciter::ctr::PrefIndicatorType,
      inciter::ctr::PrefIndicator >(lua_ideck["pref"], "indicator",
      gideck.pref.indicator.data,
      inciter::ctr::PrefIndicatorType::SPECTRAL_DECAY);

    // p-ref max degrees-of-freedom per cell
    storeIfSpecd< std::size_t >(lua_ideck["pref"], "ndofmax",
      gideck.pref.ndofmax.data, 10);

    // p-ref tolerance
    storeIfSpecd< tk::real >(lua_ideck["pref"], "tolref",
      gideck.pref.tolref.data, 0.5);

    // error checking on the tolerance
    if (gideck.pref.tolref.data < 0.0 || gideck.pref.tolref.data > 1.0)
      Throw("The p-refinement tolerance must be a real number "
        "between 0.0 and 1.0, both inclusive.");
  }

  // Mesh specification block (for overset)
  // ---------------------------------------------------------------------------
  if (lua_ideck["mesh"].valid()) {
    const sol::table& lua_mesh = lua_ideck["mesh"];
    gideck.mesh.resize(lua_mesh.size());

    for (std::size_t i=0; i<gideck.mesh.size(); ++i) {
      // filename
      storeIfSpecd< std::string >(lua_mesh[i+1], "filename",
        gideck.mesh[i].filename.data, "");

      // location
      storeVecIfSpecd< tk::real >(lua_mesh[i+1], "location",
        gideck.mesh[i].location.data, {0.0, 0.0, 0.0});
      if (gideck.mesh[i].location.data.size() != 3)
        Throw("Mesh location requires 3 coordinates.");

      // orientation
      storeVecIfSpecd< tk::real >(lua_mesh[i+1], "orientation",
        gideck.mesh[i].orientation.data, {0.0, 0.0, 0.0});
      if (gideck.mesh[i].orientation.data.size() != 3)
        Throw("Mesh orientation requires 3 rotation angles.");

      // velocity
      storeVecIfSpecd< tk::real >(lua_mesh[i+1], "velocity",
        gideck.mesh[i].velocity.data, {0.0, 0.0, 0.0});
      if (gideck.mesh[i].velocity.data.size() != 3)
        Throw("Mesh velocity requires 3 components.");
    }
  }
  else {
    // TODO: remove double-specification of defaults
    gideck.mesh.resize(1);
    gideck.mesh[0].filename.data = "";
    gideck.mesh[0].location.data = {0.0, 0.0, 0.0};
    gideck.mesh[0].orientation.data = {0.0, 0.0, 0.0};
    gideck.mesh[0].velocity.data = {0.0, 0.0, 0.0};
  }

  // Boundary conditions block
  // ---------------------------------------------------------------------------
  if (lua_ideck["bc"].valid()) {
    std::set< std::size_t > totalmesh;
    const sol::table& sol_bc = lua_ideck["bc"];
    gideck.bc.resize(sol_bc.size());

    for (std::size_t i=0; i<gideck.bc.size(); ++i) {
      storeVecIfSpecd< std::size_t >(sol_bc[i+1], "mesh",
        gideck.bc[i].mesh.data, {1});
      // collect meshes for error checking
      totalmesh.insert(gideck.bc[i].mesh.data.begin(),
        gideck.bc[i].mesh.data.end());

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "dirichlet",
        gideck.bc[i].dirichlet.data, {});

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "symmetry",
        gideck.bc[i].symmetry.data, {});

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "inlet",
        gideck.bc[i].inlet.data, {});

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "outlet",
        gideck.bc[i].outlet.data, {});

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "farfield",
        gideck.bc[i].farfield.data, {});

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "extrapolate",
        gideck.bc[i].extrapolate.data, {});

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "stag",
        gideck.bc[i].stag.data, {});

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "skip",
        gideck.bc[i].skip.data, {});

      storeVecIfSpecd< uint64_t >(sol_bc[i+1], "sponge",
        gideck.bc[i].sponge.data, {});

      // Time-dependent BC
      if (sol_bc[i+1]["timedep"].valid()) {
        storeVecIfSpecd< uint64_t >(sol_bc[i+1]["timedep"], "sideset",
          gideck.bc[i].timedep.sideset.data, {});
        storeVecIfSpecd< tk::real >(sol_bc[i+1]["timedep"], "fn",
          gideck.bc[i].timedep.fn.data, {});

        // error checking on user-def function
        if (gideck.bc[i].timedep.fn.data.size() % 6 != 0)
          Throw("Incomplete user-defined function for time-dependent BC. An "
          "R->R^5 function is expected, the number of descrete entries must be "
          "divisible by 6: one 'column' for the abscissa, and 5 for the ordinate.");
      }

      // Stagnation point
      storeVecIfSpecd< tk::real >(sol_bc[i+1], "point", gideck.bc[i].point.data,
        {0.0, 0.0, 0.0});
      if (gideck.bc[i].point.data.size() != 3)
        Throw("BC point requires 3 coordinates.");

      // Stagnation radius
      storeIfSpecd< tk::real >(sol_bc[i+1], "radius", gideck.bc[i].radius.data,
        0.0);

      // Velocity for inlet/farfield
      storeVecIfSpecd< tk::real >(sol_bc[i+1], "velocity",
        gideck.bc[i].velocity.data, {0.0, 0.0, 0.0});
      if (gideck.bc[i].velocity.data.size() != 3)
        Throw("BC velocity requires 3 components.");

      // Pressure for inlet/outlet/farfield
      storeIfSpecd< tk::real >(sol_bc[i+1], "pressure",
        gideck.bc[i].pressure.data, 0.0);

      // Density for inlet/outlet/farfield
      storeIfSpecd< tk::real >(sol_bc[i+1], "density",
        gideck.bc[i].density.data, 0.0);
    }

    // error checking on number of meshes
    if (totalmesh.size() != gideck.mesh.size())
      Throw("Total meshes (" + std::to_string(gideck.mesh.size()) +
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

    // background IC values
    storeIfSpecd< std::size_t >(lua_ideck["ic"], "materialid",
      gideck.ic.materialid.data, 1);

    storeIfSpecd< tk::real >(lua_ideck["ic"], "pressure",
      gideck.ic.pressure.data, 0.0);

    storeIfSpecd< tk::real >(lua_ideck["ic"], "temperature",
      gideck.ic.temperature.data, 0.0);

    storeVecIfSpecd< tk::real >(lua_ideck["ic"], "velocity",
      gideck.ic.velocity.data, {0.0, 0.0, 0.0});
    if (gideck.ic.velocity.data.size() != 3)
      Throw("Velocity in IC requires 3 components.");

    // IC box
    if (lua_ideck["ic"]["box"].valid()) {
      const sol::table& lua_box = lua_ideck["ic"]["box"];
      gideck.ic.box.resize(lua_box.size());

      for (std::size_t i=0; i<gideck.ic.box.size(); ++i) {
        storeIfSpecd< std::size_t >(lua_box[i+1], "materialid",
          gideck.ic.box[i].materialid.data, 1);

        storeIfSpecd< tk::real >(lua_box[i+1], "volume",
          gideck.ic.box[i].volume.data, 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "mass",
          gideck.ic.box[i].mass.data, 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "density",
          gideck.ic.box[i].density.data, 0.0);

        storeVecIfSpecd< tk::real >(lua_box[i+1], "velocity",
          gideck.ic.box[i].velocity.data, {0.0, 0.0, 0.0});
        if (gideck.ic.box[i].velocity.data.size() != 3)
          Throw("Velocity in IC box requires 3 components.");

        storeIfSpecd< tk::real >(lua_box[i+1], "pressure",
          gideck.ic.box[i].pressure.data, 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "energy",
          gideck.ic.box[i].energy.data, 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "energy_content",
          gideck.ic.box[i].energy_content.data, 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "temperature",
          gideck.ic.box[i].temperature.data, 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "xmin",
          gideck.ic.box[i].xmin.data, 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "xmax",
          gideck.ic.box[i].xmax.data, 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "ymin",
          gideck.ic.box[i].ymin.data, 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "ymax",
          gideck.ic.box[i].ymax.data, 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "zmin",
          gideck.ic.box[i].zmin.data, 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "zmax",
          gideck.ic.box[i].zmax.data, 0.0);

        storeVecIfSpecd< tk::real >(lua_box[i+1], "orientation",
          gideck.ic.box[i].orientation.data, {0.0, 0.0, 0.0});
        if (gideck.ic.box[i].orientation.data.size() != 3)
          Throw("Orientation in IC box requires 3 rotation angles.");

        storeOptIfSpecd< inciter::ctr::InitiateType, inciter::ctr::Initiate >(
          lua_box[i+1], "initiate", gideck.ic.box[i].initiate.data,
          inciter::ctr::InitiateType::IMPULSE);

        storeVecIfSpecd< tk::real >(lua_box[i+1], "point",
          gideck.ic.box[i].point.data, {0.0, 0.0, 0.0});
        if (gideck.ic.box[i].point.data.size() != 3)
          Throw("Point in IC box requires 3 coordinates.");

        storeIfSpecd< tk::real >(lua_box[i+1], "init_time",
          gideck.ic.box[i].init_time.data, 0.0);

        storeIfSpecd< tk::real >(lua_box[i+1], "front_width",
          gideck.ic.box[i].front_width.data, 0.0);
      }
    }

    // IC mesh-block
    if (lua_ideck["ic"]["meshblock"].valid()) {
      const sol::table& lua_meshblock = lua_ideck["ic"]["meshblock"];
      gideck.ic.meshblock.resize(lua_meshblock.size());

      for (std::size_t i=0; i<gideck.ic.meshblock.size(); ++i) {
        storeIfSpecd< std::size_t >(lua_meshblock[i+1], "blockid",
          gideck.ic.meshblock[i].blockid.data, 0);
        if (gideck.ic.meshblock[i].blockid.data == 0)
          Throw("Each IC mesh block must specify the mesh block id.");

        storeIfSpecd< std::size_t >(lua_meshblock[i+1], "materialid",
          gideck.ic.meshblock[i].materialid.data, 1);

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "energy_content",
          gideck.ic.meshblock[i].energy_content.data, 0.0);

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "volume",
          gideck.ic.meshblock[i].volume.data, 0.0);
        if (gideck.ic.meshblock[i].energy_content.data > 0.0 &&
          gideck.ic.meshblock[i].volume.data < 1e-12)
          Throw("Mesh block volume must be specified, if energy content is "
            "used to initialize block");

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "mass",
          gideck.ic.meshblock[i].mass.data, 0.0);

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "density",
          gideck.ic.meshblock[i].density.data, 0.0);

        storeVecIfSpecd< tk::real >(lua_meshblock[i+1], "velocity",
          gideck.ic.meshblock[i].velocity.data, {0.0, 0.0, 0.0});
        if (gideck.ic.meshblock[i].velocity.data.size() != 3)
          Throw("Velocity in IC meshblock requires 3 components.");

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "pressure",
          gideck.ic.meshblock[i].pressure.data, 0.0);

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "energy",
          gideck.ic.meshblock[i].energy.data, 0.0);

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "temperature",
          gideck.ic.meshblock[i].temperature.data, 0.0);

        storeOptIfSpecd< inciter::ctr::InitiateType, inciter::ctr::Initiate >(
          lua_meshblock[i+1], "initiate", gideck.ic.meshblock[i].initiate.data,
          inciter::ctr::InitiateType::IMPULSE);

        storeVecIfSpecd< tk::real >(lua_meshblock[i+1], "point",
          gideck.ic.meshblock[i].point.data, {0.0, 0.0, 0.0});
        if (gideck.ic.meshblock[i].point.data.size() != 3)
          Throw("Point in IC meshblock requires 3 coordinates.");

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "init_time",
          gideck.ic.meshblock[i].init_time.data, 0.0);

        storeIfSpecd< tk::real >(lua_meshblock[i+1], "front_width",
          gideck.ic.meshblock[i].front_width.data, 0.0);
      }
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // Testing
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------

  //std::cout << gideck.title.data << std::endl;
  //std::cout << "cfl = " << gideck.cfl.data << std::endl;
  //std::cout << "scheme = " << static_cast< std::size_t >(gideck.scheme.data) << std::endl;
  //auto nmat = gideck.multimat.nmat.data;
  //std::cout << "nmat " << nmat << std::endl;
  //std::cout << "ncomp " << gideck.ncomp.data << std::endl;

  //for (std::size_t k=0; k<nmat; ++k)
  //  std::cout << gideck.matidxmap.eosidx[k] << ", ";
  //std::cout << std::endl;
  //for (std::size_t k=0; k<nmat; ++k)
  //  std::cout << gideck.matidxmap.matidx[k] << ", ";
  //std::cout << std::endl;
  //for (std::size_t k=0; k<nmat; ++k)
  //  std::cout << gideck.matidxmap.solidx[k] << ", ";
  //std::cout << std::endl;

  //inciter::ctr::Material matclass;
  //for (std::size_t i=0; i<gideck.material.size(); ++i) {
  //  std::cout << "MatType "
  //    << matclass.name(gideck.material[i].eos.data) << std::endl;
  //  for (std::size_t k=0; k<gideck.material[i].id.data.size(); ++k) {
  //    std::cout << gideck.material[i].id.data[k] << std::endl;

  //    //std::cout << ", cv = " << gideck.material[i].cv.data[k];

  //    //if (gideck.material[i].eos.data == "stiffenedgas") {
  //    //  std::cout << ", gamma = " << gideck.material[i].gamma.data[k] << ", "
  //    //    << ", pstiff = " << gideck.material[i].pstiff.data[k] << ", ";
  //    //}

  //    //if (gideck.material[i].eos.data == "smallshearsolid") {
  //    //  std::cout << ", gamma = " << gideck.material[i].gamma.data[k] << ", "
  //    //    << ", pstiff = " << gideck.material[i].pstiff.data[k] << ", "
  //    //    << ", mu = " << gideck.material[i].mu.data[k] << ", ";
  //    //}

  //    //if (gideck.material[i].eos.data == "jwl") {
  //    //  std::cout << ", w_gru = " << gideck.material[i].w_gru.data[k] << ", "
  //    //    << ", A_jwl = " << gideck.material[i].A_jwl.data[k] << ", "
  //    //    << ", B_jwl = " << gideck.material[i].B_jwl.data[k] << ", "
  //    //    << ", R1_jwl = " << gideck.material[i].R1_jwl.data[k] << ", "
  //    //    << ", R2_jwl = " << gideck.material[i].R2_jwl.data[k] << ", "
  //    //    << ", rho0_jwl = " << gideck.material[i].rho0_jwl.data[k] << ", "
  //    //    << ", de_jwl = " << gideck.material[i].de_jwl.data[k] << ", "
  //    //    << ", Tr_jwl = " << gideck.material[i].Tr_jwl.data[k] << ", "
  //    //    << ", Pr_jwl = " << gideck.material[i].Pr_jwl.data[k] << ", ";
  //    //}

  //    //std::cout << std::endl;
  //  }
  //}

  //std::cout << " diag iter: " << gideck.diagnostics.iter_interval.data << std::endl;
  //std::cout << " diag norm: " << static_cast< std::size_t >(gideck.diagnostics.error.data) << std::endl;
  //std::cout << " diag format: " << static_cast< std::size_t >(gideck.diagnostics.format.data) << std::endl;
  //std::cout << " diag prec: " << gideck.diagnostics.precision.data << std::endl;

  //std::cout << " F-O/P ref: "<< gideck.field_output.refined.data << std::endl;
  //std::cout << " F-O/P time: "<< gideck.field_output.time_interval.data << std::endl;
  //std::cout << " F-O/P iter: "<< gideck.field_output.iter_interval.data << std::endl;
  //std::cout << " F-O/P filetype: "<< static_cast< std::size_t >(gideck.field_output.filetype.data) << std::endl;
  //std::cout << " F-O/P sidesets: ";
  //for (std::size_t i=0; i< gideck.field_output.sideset.data.size(); ++i)
  //  std::cout << gideck.field_output.sideset.data[i] << ", ";
  //std::cout << std::endl;
  //std::cout << " F-O/P elemvars: ";
  //for (std::size_t i=0; i< gideck.field_output.elemvar.data.size(); ++i)
  //  std::cout << gideck.field_output.elemvar.data[i] << ", ";
  //std::cout << std::endl;
  //std::cout << " F-O/P nodevars: ";
  //for (std::size_t i=0; i< gideck.field_output.nodevar.data.size(); ++i)
  //  std::cout << gideck.field_output.nodevar.data[i] << ", ";
  //std::cout << std::endl;

  //std::cout << " hist time: " << gideck.history_output.time_interval.data << std::endl;
  //std::cout << " hist iter: " << gideck.history_output.iter_interval.data << std::endl;
  //std::cout << " hist range: ";
  //for (std::size_t i=0; i< gideck.history_output.time_range.data.size(); ++i)
  //  std::cout << gideck.history_output.time_range.data[i] << ", ";
  //std::cout << std::endl;
  //std::cout << " hist point: " << std::endl;
  //for (std::size_t i=0; i< gideck.history_output.point.size(); ++i) {
  //  std::cout << gideck.history_output.point[i].id.data << ": ";
  //  for (std::size_t j=0; j< gideck.history_output.point[i].coord.data.size(); ++j)
  //    std::cout << gideck.history_output.point[i].coord.data[j] << ", ";
  //  std::cout << std::endl;
  //}
  //std::cout << std::endl;

  //std::cout << " ALE " << gideck.ale.ale.data << std::endl;
  //for (std::size_t i=0; i< gideck.ale.move.size(); ++i) {
  //  std::cout << static_cast<std::size_t>(gideck.ale.move[i].fntype.data) << ": "
  //    << gideck.ale.move[i].sideset.data << std::endl;
  //  for (std::size_t j=0; j< gideck.ale.move[i].fn.data.size(); ++j)
  //    std::cout << gideck.ale.move[i].fn.data[j] << ", ";
  //  std::cout << std::endl;
  //}

  //std::cout << " AMR " << gideck.amr.amr.data << std::endl;
  //std::cout << " t0ref " << gideck.amr.t0ref.data << std::endl;
  //std::cout << " dtref " << gideck.amr.dtref.data << std::endl;
  //std::cout << " dtref_uniform " << gideck.amr.dtref_uniform.data << std::endl;
  //std::cout << " dtfreq " << gideck.amr.dtfreq.data << std::endl;
  //std::cout << " maxlevels " << gideck.amr.maxlevels.data << std::endl;
  //std::cout << " initial ";
  //for (std::size_t i=0; i<gideck.amr.initial.data.size(); ++i)
  //  std::cout << static_cast< std::size_t >(gideck.amr.initial.data[i]) << ", ";
  //std::cout << std::endl;
  //std::cout << " error type " << static_cast< std::size_t >(gideck.amr.error.data)
  //  << std::endl;

  //std::cout << " pref " << gideck.pref.pref.data << std::endl;
  //std::cout << " indicator type " <<
  //  static_cast< std::size_t >(gideck.pref.indicator.data) << std::endl;
  //std::cout << " ndofmax " << gideck.pref.ndofmax.data << std::endl;
  //std::cout << " tolref " << gideck.pref.tolref.data << std::endl;

  //std::cout << " Meshes: " << std::endl;
  //for (std::size_t j=0; j< gideck.mesh.size(); ++j) {
  //  std::cout << " " << gideck.mesh[j].filename.data << std::endl;
  //  std::cout << " location: ";
  //  for (std::size_t i=0; i<gideck.mesh[j].location.data.size(); ++i)
  //    std::cout << gideck.mesh[j].location.data[i] << ", ";
  //  std::cout << std::endl;
  //  std::cout << " orientation: ";
  //  for (std::size_t i=0; i<gideck.mesh[j].orientation.data.size(); ++i)
  //    std::cout << gideck.mesh[j].orientation.data[i] << ", ";
  //  std::cout << std::endl;
  //  std::cout << " velocity: ";
  //  for (std::size_t i=0; i<gideck.mesh[j].velocity.data.size(); ++i)
  //    std::cout << gideck.mesh[j].velocity.data[i] << ", ";
  //  std::cout << std::endl;
  //}

  //std::cout << " BCs: " << std::endl;
  //for (std::size_t j=0; j< gideck.bc.size(); ++j) {
  //  std::cout << " bc-meshes: ";
  //  for (std::size_t i=0; i<gideck.bc[j].mesh.data.size(); ++i)
  //    std::cout << gideck.bc[j].mesh.data[i] << ", ";
  //  std::cout << std::endl;
  //  std::cout << " dirichlet: ";
  //  for (std::size_t i=0; i<gideck.bc[j].dirichlet.data.size(); ++i)
  //    std::cout << gideck.bc[j].dirichlet.data[i] << ", ";
  //  std::cout << std::endl;
  //  std::cout << " symmetry: ";
  //  for (std::size_t i=0; i<gideck.bc[j].symmetry.data.size(); ++i)
  //    std::cout << gideck.bc[j].symmetry.data[i] << ", ";
  //  std::cout << std::endl;
  //  std::cout << " farfield: ";
  //  for (std::size_t i=0; i<gideck.bc[j].farfield.data.size(); ++i)
  //    std::cout << gideck.bc[j].farfield.data[i] << ", ";
  //  std::cout << std::endl;
  //  std::cout << " inlet: ";
  //  for (std::size_t i=0; i<gideck.bc[j].inlet.data.size(); ++i)
  //    std::cout << gideck.bc[j].inlet.data[i] << ", ";
  //  std::cout << std::endl;
  //  std::cout << " pressure " << gideck.bc[j].pressure.data << std::endl;
  //  std::cout << " density " << gideck.bc[j].density.data << std::endl;
  //  std::cout << " velocity " << gideck.bc[j].velocity.data[0] << ", "
  //    << gideck.bc[j].velocity.data[1] << ", "
  //    << gideck.bc[j].velocity.data[2] << std::endl;
  //  std::cout << " timedep: ";
  //  for (std::size_t i=0; i<gideck.bc[j].timedep.sideset.data.size(); ++i)
  //    std::cout << gideck.bc[j].timedep.sideset.data[i] << ", ";
  //  std::cout << std::endl;
  //  for (std::size_t i=0; i<gideck.bc[j].timedep.fn.data.size(); ++i)
  //    std::cout << gideck.bc[j].timedep.fn.data[i] << ", ";
  //  std::cout << std::endl;
  //}

  //std::cout << " ICs: " << std::endl;
  //std::cout << " matid " << gideck.ic.materialid.data << std::endl;
  //std::cout << " pressure " << gideck.ic.pressure.data << std::endl;
  //std::cout << " temperature " << gideck.ic.temperature.data << std::endl;
  //std::cout << " velocity " << gideck.ic.velocity.data[0] << ", "
  //  << gideck.ic.velocity.data[1] << ", "
  //  << gideck.ic.velocity.data[2] << std::endl;
  //std::cout << " IC box: " << std::endl;
  //for (std::size_t i=0; i<gideck.ic.box.size(); ++i) {
  //  std::cout << "  box " << i << std::endl;
  //  std::cout << " energy_content " <<
  //    gideck.ic.box[i].energy_content.data << std::endl;
  //  std::cout << " xmin " <<
  //    gideck.ic.box[i].xmin.data << std::endl;
  //  std::cout << " xmax " <<
  //    gideck.ic.box[i].xmax.data << std::endl;
  //  std::cout << " ymin " <<
  //    gideck.ic.box[i].ymin.data << std::endl;
  //  std::cout << " ymax " <<
  //    gideck.ic.box[i].ymax.data << std::endl;
  //  std::cout << " zmin " <<
  //    gideck.ic.box[i].zmin.data << std::endl;
  //  std::cout << " zmax " <<
  //    gideck.ic.box[i].zmax.data << std::endl;
  //  std::cout << " orientation "
  //    << gideck.ic.box[i].orientation.data[0] << ", "
  //    << gideck.ic.box[i].orientation.data[1] << ", "
  //    << gideck.ic.box[i].orientation.data[2] << ", "
  //    << std::endl;
  //  std::cout << " initiate " << static_cast< std::size_t >(gideck.ic.box[i].initiate.data) << std::endl;
  //  std::cout << " point "
  //    << gideck.ic.box[i].point.data[0] << ", "
  //    << gideck.ic.box[i].point.data[1] << ", "
  //    << gideck.ic.box[i].point.data[2] << ", "
  //    << std::endl;
  //  std::cout << " init_time " <<
  //    gideck.ic.box[i].init_time.data << std::endl;
  //  std::cout << " front_width " <<
  //    gideck.ic.box[i].front_width.data << std::endl;
  //}
  //std::cout << " IC block: " << std::endl;
  //for (std::size_t i=0; i<gideck.ic.meshblock.size(); ++i) {
  //  std::cout << "  meshblock " << i << std::endl;
  //  std::cout << " blockid " <<
  //    gideck.ic.meshblock[i].blockid.data << std::endl;
  //  std::cout << " materid " <<
  //    gideck.ic.meshblock[i].materialid.data << std::endl;
  //  std::cout << " energy_content " <<
  //    gideck.ic.meshblock[i].energy_content.data << std::endl;
  //  std::cout << " initiate " << static_cast< std::size_t >(gideck.ic.meshblock[i].initiate.data) << std::endl;
  //  std::cout << " point "
  //    << gideck.ic.meshblock[i].point.data[0] << ", "
  //    << gideck.ic.meshblock[i].point.data[1] << ", "
  //    << gideck.ic.meshblock[i].point.data[2] << ", "
  //    << std::endl;
  //  std::cout << " init_time " <<
  //    gideck.ic.meshblock[i].init_time.data << std::endl;
  //  std::cout << " front_width " <<
  //    gideck.ic.meshblock[i].front_width.data << std::endl;
  //}

  return gideck;
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
