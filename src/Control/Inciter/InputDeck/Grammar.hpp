// *****************************************************************************
/*!
  \file      src/Control/Inciter/InputDeck/Grammar.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter's input deck grammar definition
  \details   Inciter's input deck grammar definition. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Word of advice: read from the bottom up.
*/
// *****************************************************************************
#ifndef InciterInputDeckGrammar_h
#define InciterInputDeckGrammar_h

#include <limits>
#include <cmath>

#include "CommonGrammar.hpp"
#include "CartesianProduct.hpp"
#include "Keywords.hpp"
#include "ContainerUtil.hpp"
#include "Centering.hpp"
#include "PDE/MultiMat/MultiMatIndexing.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck_defaults;

//! Inciter input deck facilitating user input for computing shock hydrodynamics
namespace deck {

  //! \brief Specialization of tk::grm::use for Inciter's input deck parser
  template< typename keyword >
  using use = tk::grm::use< keyword, ctr::InputDeck::keywords >;

  // Inciter's InputDeck state

  //! \brief Number of registered equations
  //! \details Counts the number of parsed equation blocks during parsing.
  static tk::TaggedTuple< brigand::list<
             tag::transport,   std::size_t
           , tag::compflow,    std::size_t
           , tag::multimat,    std::size_t
         > > neq;

  //! \brief Parser-lifetime storage for point names
  //! \details Used to track the point names registered so that parsing new ones
  //!    can be required to be unique.
  static std::set< std::string > pointnames;

  //! Parser-lifetime storage of elem or node centering
  static tk::Centering centering = tk::Centering::NODE;

  //! Assign depvars starting from this value
  static char depvar_cnt = 'a';
  //! Assign solver coupling starting from this value
  static std::size_t transfer_cnt = 0;

  //! Accepted multimat output variable labels and associated index functions
  //! \details The keys are a list of characters accepted as labels for
  //! denoting (matvar-style) output variables used for multi-material variable
  //! output. We use a case- insesitive comparitor, since when this set is used
  //! we only care about whether the variable is selected or not and not whether
  //! it denotes a full variable (upper case) or a fluctuation (lower case).
  //! This is true when matching these labels for output variables with
  //! instantaenous variables as well terms of products in parsing requested
  //! statistics (for turbulence). The values are associated indexing functions
  //! used to index into the state of the multimaterial system, all must follow
  //! the same signature.
  static std::map< char, tk::MultiMatIdxFn,
                   tk::ctr::CaseInsensitiveCharLess >
    multimatvars{
      { 'd', densityIdx }       // density
    , { 'f', volfracIdx }       // volume fraction
    , { 'm', momentumIdx }      // momentum
    , { 'e', energyIdx }        // specific total energy
    , { 'u', velocityIdx }      // velocity (primitive)
    , { 'p', pressureIdx }      // material pressure (primitive)
  };

} // ::deck
} // ::inciter

namespace tk {
namespace grm {

  using namespace tao;

  // Note that PEGTL action specializations must be in the same namespace as the
  // template being specialized. See http://stackoverflow.com/a/3052604.

  // Inciter's InputDeck actions

  //! Rule used to trigger action
  template< class eq > struct register_inciter_eq : pegtl::success {};
  //! \brief Register differential equation after parsing its block
  //! \details This is used by the error checking functors (check_*) during
  //!    parsing to identify the recently-parsed block.
  template< class eq >
  struct action< register_inciter_eq< eq > > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& ) {
      using inciter::deck::neq;
      ++neq.get< eq >();
    }
  };

  //! Rule used to trigger action
  template< class eq > struct check_transport : pegtl::success {};
  //! \brief Set defaults and do error checking on the transport equation block
  //! \details This is error checking that only the transport equation block
  //!   must satisfy. Besides error checking we also set defaults here as
  //!   this block is called when parsing of a transport...end block has
  //!   just finished.
  template< class eq >
  struct action< check_transport< eq > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using inciter::deck::neq;
      using tag::param;

      // Set depvar
      auto& depvar = stack.template get< param, eq, tag::depvar >();
      if (depvar.empty() || depvar.size() != neq.get< eq >()) {
        depvar.push_back('c');
        // add depvar to deck::depvars so it can be selected as outvar later
        tk::grm::depvars.insert( 'c' );
      }

      // If no number of components has been selected, default to 1
      auto& ncomp = stack.template get< tag::component, eq >();
      if (ncomp.empty() || ncomp.size() != neq.get< eq >())
        ncomp.push_back( 1 );

      // If physics type is not given, default to 'advection'
      auto& physics = stack.template get< param, eq, tag::physics >();
      if (physics.empty() || physics.size() != neq.get< eq >())
        physics.push_back( inciter::ctr::PhysicsType::ADVECTION );

      // If problem type is not given, error out
      auto& problem = stack.template get< param, eq, tag::problem >();
      if (problem.empty() || problem.size() != neq.get< eq >())
        Message< Stack, ERROR, MsgKey::NOPROBLEM >( stack, in );
      // Error check Dirichlet boundary condition block for all transport eq
      // configurations
      const auto& bc = stack.template get< param, eq, tag::bc, tag::bcdir >();
      for (const auto& s : bc)
        if (s.empty()) Message< Stack, ERROR, MsgKey::BC_EMPTY >( stack, in );

      // If interface compression is not specified, default to 'false'
      auto& intsharp = stack.template get< param, eq, tag::intsharp >();
      if (intsharp.empty() || intsharp.size() != neq.get< eq >())
        intsharp.push_back( 0 );

      // If interface compression parameter is not specified, default to 1.8
      auto& intsharp_p = stack.template get< param, eq,
                                            tag::intsharp_param >();
      if (intsharp_p.empty() || intsharp_p.size() != neq.get< eq >())
        intsharp_p.push_back( 1.8 );
    }
  };

  //! Rule used to trigger action
  template< class eq > struct check_compflow : pegtl::success {};
  //! \brief Set defaults and do error checking on the compressible flow
  //!   equation block
  //! \details This is error checking that only the compressible flow equation
  //!   block must satisfy. Besides error checking we also set defaults here as
  //!   this block is called when parsing of a compflow...end block has
  //!   just finished.
  template< class eq >
  struct action< check_compflow< eq > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using inciter::deck::neq;
      using tag::param;

      // Set default flux to HLLC if not specified
      auto& flux = stack.template get< tag::param, eq, tag::flux >();
      if (flux.empty() || flux.size() != neq.get< eq >())
        flux.push_back( inciter::ctr::FluxType::HLLC );

      auto sysfct = stack.template get< param, eq, tag::sysfct >();
      auto& sysfctvar = stack.template get< param, eq, tag::sysfctvar >();
      // If sysfctvar is not specified, use all variables for system FCT
      if (sysfct && sysfctvar.empty()) sysfctvar = { 0,1,2,3,4 };

      // Set number of components to 5 (mass, 3 x mom, energy) if not coupled
      stack.template get< tag::component, eq >().push_back( 5 );

      if (stack.template get< tag::couple, tag::transfer >().empty()) {
        auto& problem = stack.template get< param, eq, tag::problem >();
        if (problem.empty()) {
          problem.push_back( inciter::ctr::ProblemType::USER_DEFINED );
        }
        auto& depvar = stack.template get< param, eq, tag::depvar >();
        if (depvar.empty() || depvar.size() != neq.get< eq >()) {
          depvar.push_back('a');
          // add depvar to deck::depvars so it can be selected as outvar later
          tk::grm::depvars.insert( 'a' );
        }
      }

      // Verify correct number of material properties configured
      auto& matprop = stack.template get< param, eq, tag::material >().back();
      auto& matidxmap = stack.template get< param, eq, tag::matidxmap >();
      matidxmap.template get< tag::eosidx >().resize(1);
      matidxmap.template get< tag::matidx >().resize(1);
      auto& meos = matprop.template get< tag::eos >();
      auto& mat_id = matprop.template get< tag::id >();

      if (mat_id.empty())
        mat_id.push_back(0);
      else if (mat_id.size() != 1)
        Message< Stack, ERROR, MsgKey::NUMMAT >( stack, in );
      else
        mat_id[0] = 0;

      if (meos == inciter::ctr::MaterialType::STIFFENEDGAS) {
        const auto& gamma = matprop.template get< tag::gamma >();
        // If gamma vector is wrong size, error out
        if (gamma.empty() || gamma.size() != 1)
          Message< Stack, ERROR, MsgKey::EOSGAMMA >( stack, in );

        auto& cv = matprop.template get< tag::cv >();
        // As a default, the specific heat of air (717.5 J/Kg-K) is used
        if (cv.empty())
          cv.push_back(717.5);
        // If specific heat vector is wrong size, error out
        if (cv.size() != 1)
          Message< Stack, ERROR, MsgKey::EOSCV >( stack, in );

        auto& pstiff = matprop.template get< tag::pstiff >();
        // As a default, a stiffness coefficient of 0.0 is used
        if (pstiff.empty())
          pstiff.push_back(0.0);
        // If stiffness coefficient vector is wrong size, error out
        if (pstiff.size() != 1)
          Message< Stack, ERROR, MsgKey::EOSPSTIFF >( stack, in );
      }
      else {
        Message< Stack, ERROR, MsgKey::NOEOS >( stack, in );
      }

      // Generate mapping between material index and eos parameter index
      auto& eosmap = matidxmap.template get< tag::eosidx >();
      auto& idxmap = matidxmap.template get< tag::matidx >();
      eosmap[mat_id[0]] = static_cast< std::size_t >(matprop.template get<
        tag::eos >());
      idxmap[mat_id[0]] = 0;

      // If problem type is not given, default to 'user_defined'
      auto& problem = stack.template get< param, eq, tag::problem >();
      if (problem.empty() || problem.size() != neq.get< eq >()) {
        problem.push_back( inciter::ctr::ProblemType::USER_DEFINED );
      }

      // Error check on user-defined problem type
      auto& ic = stack.template get< param, eq, tag::ic >();
      //auto& bgdensityic = ic.template get< tag::density >();
      //auto& bgvelocityic = ic.template get< tag::velocity >();
      //auto& bgpressureic = ic.template get< tag::pressure >();
      //auto& bgenergyic = ic.template get< tag::energy >();
      //auto& bgtemperatureic = ic.template get< tag::temperature >();
      if (problem.back() == inciter::ctr::ProblemType::USER_DEFINED) {
        // must have defined background ICs for user-defined ICs
        //auto n = neq.get< eq >();
        //if ( bgdensityic.size() != n || bgvelocityic.size() != n ||
        //     ( bgpressureic.size() != n && bgenergyic.size() != n &&
        //       bgtemperatureic.size() != n ) )
        //{
        //  Message< Stack, ERROR, MsgKey::BGICMISSING >( stack, in );
        //}

        // Error check for ic box
        auto& box = ic.template get< tag::box >();
        if (!box.empty()) {
          for (auto& b : box.back()) {   // for all boxes
            auto& boxorient = b.template get< tag::orientation >();
            if (boxorient.size() == 0)
              boxorient.resize(3, 0.0);
            else if (boxorient.size() != 3)
              Message< Stack, ERROR, MsgKey::BOXORIENTWRONG >(stack, in);
          }
        }

        // Error check for ic mesh block
        const auto& mblock = ic.template get< tag::meshblock >();
        if (!mblock.empty()) {
          for (const auto& b : mblock.back()) {   // for all blocks
            if (stack.template get< tag::discr, tag::scheme >() ==
              inciter::ctr::SchemeType::ALECG) {
              Message< Stack, ERROR, MsgKey::MESHBLOCKSUPPORT >(stack, in);
            }
            else {
              const auto& blkid = b.template get< tag::blockid >();
              if (blkid == 0)
                Message< Stack, ERROR, MsgKey::MESHBLOCKIDMISSING >(stack, in);

              // if energy content is used to initialize block, then volume must
              // be specified
              const auto& blkenc = b.template get< tag::energy_content >();
              const auto& blkvol = b.template get< tag::volume >();
              if (blkenc > 0.0 && blkvol < 1e-12)
                Message< Stack, ERROR, MsgKey::MESHBLOCKVOL >(stack, in);
            }
          }
        }

        // Error check Dirichlet boundary condition block for all compflow
        // configurations
        const auto& bc = stack.template get< param, eq, tag::bc, tag::bcdir >();
        for (const auto& s : bc)
          if (s.empty()) Message< Stack, ERROR, MsgKey::BC_EMPTY >( stack, in );

        // Error check stagnation BC block
        const auto& stag = stack.template get<tag::param, eq, tag::stag>();
        const auto& spoint = stag.template get< tag::point >();
        const auto& sradius = stag.template get< tag::radius >();
        if ( (!spoint.empty() && !spoint.back().empty() &&
              !sradius.empty() && !sradius.back().empty() &&
              spoint.back().size() != 3*sradius.back().size())
         || (!sradius.empty() && !sradius.back().empty() &&
              !spoint.empty() && !spoint.back().empty() &&
              spoint.back().size() != 3*sradius.back().size())
         || (!spoint.empty() && !spoint.back().empty() &&
              (sradius.empty() || (!sradius.empty() && sradius.back().empty())))
         || (!sradius.empty() && !sradius.back().empty() &&
              (spoint.empty() || (!spoint.empty() && spoint.back().empty()))) )
        {
          Message< Stack, ERROR, MsgKey::STAGBCWRONG >( stack, in );
        }

        // Error check skip BC block
        const auto& skip = stack.template get<tag::param, eq, tag::skip>();
        const auto& kpoint = skip.template get< tag::point >();
        const auto& kradius = skip.template get< tag::radius >();
        if ( (!kpoint.empty() && !kpoint.back().empty() &&
              !kradius.empty() && !kradius.back().empty() &&
              kpoint.back().size() != 3*kradius.back().size())
          || (!kradius.empty() && !kradius.back().empty() &&
              !kpoint.empty() && !kpoint.back().empty() &&
              kpoint.back().size() != 3*kradius.back().size())
          || (!kpoint.empty() && !kpoint.back().empty() &&
              (kradius.empty() || (!kradius.empty() && kradius.back().empty())))
          || (!kradius.empty() && !kradius.back().empty() &&
              (kpoint.empty() || (!kpoint.empty() && kpoint.back().empty()))) )
        {
          Message< Stack, ERROR, MsgKey::SKIPBCWRONG >( stack, in );
        }

        // Error check sponge BC parameter vectors for symmetry BC block
        const auto& sponge =
          stack.template get< tag::param, eq, tag::sponge >();
        const auto& ss = sponge.template get< tag::sideset >();

        const auto& spvel = sponge.template get< tag::velocity >();
        if ( !spvel.empty()) {
          if (spvel.size() != ss.size()) {
            Message< Stack, ERROR, MsgKey::SPONGEBCWRONG >( stack, in );
          }
          for (const auto& s : spvel) {
            if ( s < 0.0 || s > 1.0 ) {
              Message< Stack, ERROR, MsgKey::SPONGEBCWRONG >( stack, in );
            }
          }
        }

        const auto& sppre = sponge.template get< tag::velocity >();
        if ( !sppre.empty()) {
          if (sppre.size() != ss.size()) {
            Message< Stack, ERROR, MsgKey::SPONGEBCWRONG >( stack, in );
          }
          for (const auto& s : sppre) {
            if ( s < 0.0 || s > 1.0 ) {
              Message< Stack, ERROR, MsgKey::SPONGEBCWRONG >( stack, in );
            }
          }
        }

        // Error check user defined time dependent BC for this system
        const auto& tdepbc =
          stack.template get< tag::param, eq, tag::bctimedep >().back();
        // multiple time dependent BCs can be specified on different side sets
        for (const auto& bndry : tdepbc) {
          const auto& s = bndry.template get< tag::sideset >();
          if (s.empty()) Message< Stack, ERROR, MsgKey::BC_EMPTY >( stack, in );
          const auto& f = bndry.template get< tag::fn >();
          if (f.empty() or f.size() % 6 != 0)
            Message< Stack, ERROR, MsgKey::INCOMPLETEUSERFN>( stack, in );
        }
      }
    }
  };

  //! Rule used to trigger action
  template< class eq > struct check_multimat : pegtl::success {};
  //! \brief Set defaults and do error checking on the multimaterial
  //!    compressible flow equation block
  //! \details This is error checking that only the multimaterial compressible
  //!   flow equation block must satisfy. Besides error checking we also set
  //!   defaults here as this block is called when parsing of a
  //!   multimat...end block has just finished.
  template< class eq >
  struct action< check_multimat< eq > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using inciter::deck::neq;
      using tag::param;

      // Set depvar
      auto& depvar = stack.template get< param, eq, tag::depvar >();
      if (depvar.empty() || depvar.size() != neq.get< eq >()) {
        depvar.push_back('a');
        // add depvar to deck::depvars so it can be selected as outvar later
        tk::grm::depvars.insert( 'a' );
      }

      // If physics type is not given, default to 'euler'
      auto& physics = stack.template get< param, eq, tag::physics >();
      if (physics.empty() || physics.size() != neq.get< eq >())
        physics.push_back( inciter::ctr::PhysicsType::EULER );

      // Set default flux to AUSM if not specified
      auto& flux = stack.template get< tag::param, eq, tag::flux >();
      if (flux.empty() || flux.size() != neq.get< eq >())
        flux.push_back( inciter::ctr::FluxType::AUSM );

      // Set number of scalar components based on number of materials
      auto& nmat = stack.template get< param, eq, tag::nmat >();
      auto& ncomp = stack.template get< tag::component, eq >();
      if (physics.back() == inciter::ctr::PhysicsType::EULER ||
        physics.back() == inciter::ctr::PhysicsType::ENERGYPILL) {
        // physics = euler/energy pill: m-material compressible flow
        // scalar components: volfrac:m + mass:m + momentum:3 + energy:m
        // if nmat is unspecified, configure it be 2
        if (nmat == inciter::g_inputdeck_defaults.get<param, eq, tag::nmat>()) {
          Message< Stack, WARNING, MsgKey::NONMAT >( stack, in );
          nmat = 2;
        }

        // set ncomp based on nmat
        // if solid EOS, add components
        auto ntot = nmat + nmat + 3 + nmat;
        const auto& matprop = stack.template get< param, eq, tag::material >();
        for (const auto& mtype : matprop) {
          if (mtype.template get< tag::eos >() ==
            inciter::ctr::MaterialType::SMALLSHEARSOLID) {
            ntot += 9;
          }
        }

        ncomp.push_back( ntot );
      }

      // Verify correct number of multi-material properties (gamma, cv, pstiff)
      // have been configured
      auto& matprop = stack.template get< param, eq, tag::material >();
      auto& matidxmap = stack.template get< param, eq, tag::matidxmap >();
      matidxmap.template get< tag::eosidx >().resize(nmat);
      matidxmap.template get< tag::matidx >().resize(nmat);
      matidxmap.template get< tag::solidx >().resize(nmat);
      std::size_t tmat(0), i(0), mtypei(0), isolcntr(0), isolidx(0);
      std::set< std::size_t > matidset;

      for (auto& mtype : matprop) {
        const auto& meos = mtype.template get< tag::eos >();
        const auto& mat_id = mtype.template get< tag::id >();

        if (meos == inciter::ctr::MaterialType::STIFFENEDGAS ||
          meos == inciter::ctr::MaterialType::SMALLSHEARSOLID) {
          const auto& gamma = mtype.template get< tag::gamma >();
          // If gamma vector is wrong size, error out
          if (gamma.empty() || gamma.size() != mat_id.size())
            Message< Stack, ERROR, MsgKey::EOSGAMMA >( stack, in );

          auto& cv = mtype.template get< tag::cv >();
          // As a default, the specific heat of air (717.5 J/Kg-K) is used
          if (cv.empty()) {
            for (std::size_t k=0; k<mat_id.size(); ++k) {
              cv.push_back(717.5);
            }
          }
          // If specific heat vector is wrong size, error out
          if (cv.size() != mat_id.size())
            Message< Stack, ERROR, MsgKey::EOSCV >( stack, in );

          auto& pstiff = mtype.template get< tag::pstiff >();
          // As a default, a stiffness coefficient of 0.0 is used
          if (pstiff.empty()) {
            for (std::size_t k=0; k<mat_id.size(); ++k) {
              pstiff.push_back(0.0);
            }
          }
          // If stiffness coefficient vector is wrong size, error out
          if (pstiff.size() != mat_id.size())
            Message< Stack, ERROR, MsgKey::EOSPSTIFF >( stack, in );

          // Check shear modulus vector size
          if (meos == inciter::ctr::MaterialType::SMALLSHEARSOLID) {
            const auto& mu = mtype.template get< tag::mu >();
            if (mu.size() != mat_id.size())
              Message< Stack, ERROR, MsgKey::EOSMU >( stack, in );

            // add to solid-counter
            ++isolcntr;
            // assign solid-counter value to solid-index
            isolidx = isolcntr;
          }
          else {
            // since not solid, assign 0 to solid-index
            isolidx = 0;
          }
        }
        else if (meos == inciter::ctr::MaterialType::JWL) {
          auto& cv = mtype.template get< tag::cv >();
          // As a default, the specific heat of air (717.5 J/Kg-K) is used
          if (cv.empty()) {
            for (std::size_t k=0; k<mat_id.size(); ++k) {
              cv.push_back(717.5);
            }
          }
          // If specific heat vector is wrong size, error out
          if (cv.size() != mat_id.size())
            Message< Stack, ERROR, MsgKey::EOSCV >( stack, in );

          // If JWL parameter vectors are wrong size, error out
          const auto& w_gru = mtype.template get< tag::w_gru >();
          const auto& a_jwl = mtype.template get< tag::A_jwl >();
          const auto& b_jwl = mtype.template get< tag::B_jwl >();
          const auto& r1_jwl = mtype.template get< tag::R1_jwl >();
          const auto& r2_jwl = mtype.template get< tag::R2_jwl >();
          const auto& rho0_jwl = mtype.template get< tag::rho0_jwl >();
          const auto& de_jwl = mtype.template get< tag::de_jwl >();
          const auto& rhor_jwl = mtype.template get< tag::rhor_jwl >();
          const auto& Tr_jwl = mtype.template get< tag::Tr_jwl >();
          const auto& pr_jwl = mtype.template get< tag::Pr_jwl >();
          if (w_gru.size() != mat_id.size() || a_jwl.size() != mat_id.size() ||
            b_jwl.size() != mat_id.size() || r1_jwl.size() != mat_id.size() ||
            r2_jwl.size() != mat_id.size() || rho0_jwl.size() != mat_id.size()
            || de_jwl.size() != mat_id.size() || pr_jwl.size() != mat_id.size())
            Message< Stack, ERROR, MsgKey::EOSJWLPARAM >( stack, in );

          if (rhor_jwl.size() != mat_id.size() && Tr_jwl.size() != mat_id.size())
            Message< Stack, ERROR, MsgKey::EOSJWLREFSTATE >( stack, in );
        }

        // Track total number of materials in multiple material blocks
        tmat += mat_id.size();

        // Check for repeating user specified material ids
        for (auto midx : mat_id) {
          if (!matidset.count(midx))
            matidset.insert(midx);
          else
            Message< Stack, ERROR, MsgKey::REPMATID >( stack, in );
        }

        // Generate mapping between material index and eos parameter index
        auto& eosmap = matidxmap.template get< tag::eosidx >();
        auto& idxmap = matidxmap.template get< tag::matidx >();
        auto& solidxmap = matidxmap.template get< tag::solidx >();
        for (auto midx : mat_id) {
          midx -= 1;
          //eosmap[midx] = static_cast< std::size_t >(mtype.template get<
          //  tag::eos >());
          eosmap[midx] = mtypei;
          idxmap[midx] = i;
          solidxmap[midx] = isolidx;
          ++i;
        }
        // end of materials for this eos, thus reset index counter
        i = 0;
        // increment material-type/eos-type index counter
        ++mtypei;
      }

      // If total number of materials is incorrect, error out
      if (tmat != nmat)
        Message< Stack, ERROR, MsgKey::NUMMAT >( stack, in );

      // Check if material ids are contiguous and 1-based
      if (!matidset.count(1))
        Message< Stack, ERROR, MsgKey::ONEMATID >( stack, in );
      std::size_t icount(1);
      for (auto midx : matidset) {
        if (midx != icount)
          Message< Stack, ERROR, MsgKey::GAPMATID >( stack, in );
        ++icount;
      }

      // If pressure relaxation is not specified, default to 'true'
      auto& prelax = stack.template get< param, eq, tag::prelax >();
      if (prelax.empty() || prelax.size() != neq.get< eq >())
        prelax.push_back( 1 );

      // If pressure relaxation time-scale is not specified, default to 0.25
      auto& prelax_ts = stack.template get< param, eq,
                                            tag::prelax_timescale >();
      if (prelax_ts.empty() || prelax_ts.size() != neq.get< eq >())
        prelax_ts.push_back( 0.25 );

      // If interface compression is not specified, default to 'false'
      auto& intsharp = stack.template get< param, eq, tag::intsharp >();
      if (intsharp.empty() || intsharp.size() != neq.get< eq >())
        intsharp.push_back( 0 );

      // If interface compression parameter is not specified, default to 1.8
      auto& intsharp_p = stack.template get< param, eq,
                                            tag::intsharp_param >();
      if (intsharp_p.empty() || intsharp_p.size() != neq.get< eq >())
        intsharp_p.push_back( 1.8 );

      // If problem type is not given, default to 'user_defined'
      auto& problem = stack.template get< param, eq, tag::problem >();
      if (problem.empty() || problem.size() != neq.get< eq >()) {
        problem.push_back( inciter::ctr::ProblemType::USER_DEFINED );
      }

      // Error check on user-defined problem type
      auto& ic = stack.template get< param, eq, tag::ic >();
      auto& bgmatid = ic.template get< tag::materialid >();
      auto& bgdensityic = ic.template get< tag::density >();
      auto& bgvelocityic = ic.template get< tag::velocity >();
      auto& bgpressureic = ic.template get< tag::pressure >();
      auto& bgenergyic = ic.template get< tag::energy >();
      auto& bgtemperatureic = ic.template get< tag::temperature >();
      if (problem.back() == inciter::ctr::ProblemType::USER_DEFINED) {
        // must have defined background ICs for user-defined ICs
        auto n = neq.get< eq >();
        if (bgmatid.size() != n) {
          Message< Stack, ERROR, MsgKey::BGMATIDMISSING >( stack, in );
        }

        if ( bgdensityic.size() != n || bgvelocityic.size() != n ||
             ( bgpressureic.size() != n && bgenergyic.size() != n &&
               bgtemperatureic.size() != n ) )
        {
          Message< Stack, ERROR, MsgKey::BGICMISSING >( stack, in );
        }

        // each IC box should have material id specified, and it should be
        // within nmat
        auto& icbox = ic.template get< tag::box >();

        if (!icbox.empty()) {
          for (auto& b : icbox.back()) {   // for all boxes
            auto boxmatid = b.template get< tag::materialid >();
            if (boxmatid == 0) {
              Message< Stack, ERROR, MsgKey::BOXMATIDMISSING >( stack, in );
            }
            else if (boxmatid > nmat) {
              Message< Stack, ERROR, MsgKey::BOXMATIDWRONG >( stack, in );
            }
            auto& boxorient = b.template get< tag::orientation >();
            if (boxorient.size() == 0)
              boxorient.resize(3, 0.0);
            else if (boxorient.size() != 3)
              Message< Stack, ERROR, MsgKey::BOXORIENTWRONG >(stack, in);
          }
        }

        // each IC mesh block should have block id and material id specified,
        // and the material id should be within nmat
        const auto& mblock = ic.template get< tag::meshblock >();
        if (!mblock.empty()) {
          for (const auto& b : mblock.back()) {   // for all blocks
            const auto& blkid = b.template get< tag::blockid >();
            if (blkid == 0)
              Message< Stack, ERROR, MsgKey::MESHBLOCKIDMISSING >(stack, in);
            const auto& blkmatid = b.template get< tag::materialid >();
            if (blkmatid == 0) {
              Message< Stack, ERROR, MsgKey::BOXMATIDMISSING >( stack, in );
            }
            else if (blkmatid > nmat) {
              Message< Stack, ERROR, MsgKey::BOXMATIDWRONG >( stack, in );
            }

            // if energy content is used to initialize block, then volume must
            // be specified
            const auto& blkenc = b.template get< tag::energy_content >();
            const auto& blkvol = b.template get< tag::volume >();
            if (blkenc > 0.0 && blkvol < 1e-12)
              Message< Stack, ERROR, MsgKey::MESHBLOCKVOL >(stack, in);
          }
        }
      }

      // Error check Dirichlet boundary condition block for all multimat
      // configurations
      const auto& bc = stack.template get< param, eq, tag::bc, tag::bcdir >();
      for (const auto& s : bc)
        if (s.empty()) Message< Stack, ERROR, MsgKey::BC_EMPTY >( stack, in );
    }
  };

  //! Rule used to trigger action
  template< class Option, typename...tags >
  struct store_inciter_option : pegtl::success {};
  //! \brief Put option in state at position given by tags
  //! \details This is simply a wrapper around tk::grm::store_option passing the
  //!    stack defaults for inciter.
  template< class Option, typename... tags >
  struct action< store_inciter_option< Option, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      store_option< Stack, inciter::deck::use, Option, inciter::ctr::InputDeck,
                    Input, tags... >
                  ( stack, in, inciter::g_inputdeck_defaults );
    }
  };

  //! Function object to ensure disjoint side sets for all boundary conditions
  //! \details This is instantiated using a Cartesian product of all PDE types
  //!    and all BC types at compile time. It goes through all side sets
  //!    configured by the user and triggers an error if a side set is assigned
  //!    a BC more than once (within a solver).
  template< typename Input, typename Stack >
  struct ensure_disjoint {
    const Input& m_input;
    Stack& m_stack;
    explicit ensure_disjoint( const Input& in, Stack& stack ) :
      m_input( in ), m_stack( stack ) {}
    template< typename U > void operator()( brigand::type_<U> ) {
      using Eq = typename brigand::front< U >;
      using BC = typename brigand::back< U >;
      const auto& bc = m_stack.template get< tag::param, Eq, tag::bc, BC >();
      std::unordered_set< int > bcset;
      for (const auto& s : bc) {
        auto id = std::stoi(s);
        if (bcset.find(id) != end(bcset))
          Message< Stack, ERROR, MsgKey::NONDISJOINTBC >( m_stack, m_input );
        else
          bcset.insert( id );
      }
    }
  };

  //! Rule used to trigger action
  struct configure_scheme : pegtl::success {};
  //! Configure scheme selected by user
  //! \details This grammar action configures the number of degrees of freedom
  //! (NDOF) used for DG methods. For finite volume (or DGP0), the DOF are the
  //! cell-averages. This implies ndof=1 for DGP0. Similarly, ndof=4 and 10 for
  //! DGP1 and DGP2 respectively, since they evolve higher (>1) order solution
  //! information (e.g. gradients) as well. "rdof" includes degrees of freedom
  //! that are both, evolved and reconstructed. For rDGPnPm methods (e.g. P0P1
  //! and P1P2), "n" denotes the evolved solution-order and "m" denotes the
  //! reconstructed solution-order; i.e. P0P1 has ndof=1 and rdof=4, whereas
  //! P1P2 has ndof=4 and rdof=10. For a pure DG method without reconstruction
  //! (DGP0, DGP1, DGP2), rdof=ndof. For more information about rDGPnPm methods,
  //! ref. Luo, H. et al. (2013).  A reconstructed discontinuous Galerkin method
  //! based on a hierarchical WENO reconstruction for compressible flows on
  //! tetrahedral grids. Journal of Computational Physics, 236, 477-492.
  template<> struct action< configure_scheme > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& stack ) {
      using inciter::ctr::SchemeType;
      auto& discr = stack.template get< tag::discr >();
      auto& ndof = discr.template get< tag::ndof >();
      auto& rdof = discr.template get< tag::rdof >();
      auto scheme = discr.template get< tag::scheme >();
      if (scheme == SchemeType::P0P1 || scheme == SchemeType::FV) {
        ndof = 1; rdof = 4;
      } else if (scheme == SchemeType::DGP1) {
        ndof = rdof = 4;
      } else if (scheme == SchemeType::DGP2) {
        ndof = rdof = 10;
      } else if (scheme == SchemeType::PDG) {
        ndof = rdof = 10;
        stack.template get< tag::pref, tag::pref >() = true;
      } else if (scheme != SchemeType::DG &&
          scheme != SchemeType::DiagCG &&
          scheme != SchemeType::ALECG) {
        Throw("Scheme type not configured in configure_scheme");
      }
    }
  };

  //! Function object to do error checking on output time ranges
  template< typename Stack, typename Input >
  struct range_errchk {
    Stack& stack;
    const Input& input;
    explicit range_errchk( Stack& s, const Input& in ) : stack(s), input(in) {}
    template< typename U > void operator()( brigand::type_<U> ) {
      for (const auto& r : stack.template get< tag::output, tag::range, U >())
        if ( r.size() != 3 or r[0] > r[1] or r[2] < 0.0 or r[2] > r[1]-r[0] )
          Message< Stack, ERROR, MsgKey::BADRANGE >( stack, input );
    }
  };

  //! Rule used to trigger action
  struct check_inciter : pegtl::success {};
  //! \brief Do error checking on the inciter block
  template<> struct action< check_inciter > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using inciter::deck::neq;
      using inciter::g_inputdeck_defaults;

      // Error out if no dt policy has been selected
      const auto& dt = stack.template get< tag::discr, tag::dt >();
      const auto& cfl = stack.template get< tag::discr, tag::cfl >();
      if ( std::abs(dt - g_inputdeck_defaults.get< tag::discr, tag::dt >()) <
            std::numeric_limits< tk::real >::epsilon() &&
          std::abs(cfl - g_inputdeck_defaults.get< tag::discr, tag::cfl >()) <
            std::numeric_limits< tk::real >::epsilon() )
        Message< Stack, ERROR, MsgKey::NODT >( stack, in );

      // If both dt and cfl are given, warn that dt wins over cfl
      if ( std::abs(dt - g_inputdeck_defaults.get< tag::discr, tag::dt >()) >
            std::numeric_limits< tk::real >::epsilon() &&
          std::abs(cfl - g_inputdeck_defaults.get< tag::discr, tag::cfl >()) >
            std::numeric_limits< tk::real >::epsilon() )
        Message< Stack, WARNING, MsgKey::MULDT >( stack, in );

      // Do error checking on time history points
      const auto& hist = stack.template get< tag::history, tag::point >();
      if (std::any_of( begin(hist), end(hist),
           [](const auto& p){ return p.size() != 3; } ) )
      {
        Message< Stack, ERROR, MsgKey::WRONGSIZE >( stack, in );
      }

      // Ensure no different BC types are assigned to the same side set
      using PDETypes = inciter::ctr::parameters::Keys;
      using BCTypes = inciter::ctr::bc::Keys;
      brigand::for_each< tk::cartesian_product< PDETypes, BCTypes > >(
        ensure_disjoint< Input, Stack >( in, stack ) );

      // Do error checking on output time range configuration parameters: they
      // all must be a 3 reals: mintime, maxtime, and dt with maxtime >
      // mintime, and dt<maxtime-mintime.
      brigand::for_each< inciter::ctr::time_range::Keys >
                       ( range_errchk< Stack, Input >( stack, in ) );

      // Do error checking on time history point names (this is a programmer
      // error if triggers, hence assert)
      Assert(
        (stack.template get< tag::history, tag::id >().size() == hist.size()),
        "Number of history points and ids must equal" );
    }
  };

  //! Rule used to trigger action
  struct check_ale : pegtl::success {};
  //! \brief Do error checking on the inciter block
  template<> struct action< check_ale > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using inciter::g_inputdeck_defaults;

     // Trigger error if steady state + ALE are both enabled
      auto steady = stack.template get< tag::discr, tag::steady_state >();
      auto ale = stack.template get< tag::ale, tag::ale >();
      if (steady && ale) {
        Message< Stack, ERROR, MsgKey::STEADYALE >( stack, in );
      }

      // Set a sensible default for dvCFL if ALE is enabled and if dvcfl not set
      auto& dvcfl = stack.template get< tag::ale, tag::dvcfl >();
      auto dvcfl_default = g_inputdeck_defaults.get< tag::ale, tag::dvcfl >();
      auto eps = std::numeric_limits< tk::real >::epsilon();
      if (ale && std::abs(dvcfl - dvcfl_default) < eps) dvcfl = 0.01;

      // Set a default of zeros for the mesh force ALE parameters
      auto& meshforce = stack.template get< tag::ale, tag::meshforce >();
      if (ale && meshforce.size() != 4) meshforce = { 0, 0, 0, 0 };

      // Set a default for the ALE mesh motion dimensions
      auto& mesh_motion = stack.template get< tag::ale, tag::mesh_motion >();
      if (ale && mesh_motion.empty()) mesh_motion = { 0, 1, 2 };

      // Error out if mesh motion dimensions are wrong
      if (ale && (mesh_motion.size() > 3 ||
                  std::any_of( begin(mesh_motion), end(mesh_motion),
                     [](auto d){return d > 2;} )) )
      {
        Message< Stack, ERROR, MsgKey::WRONGMESHMOTION >( stack, in );
      }

      // Error checking on user-defined function for ALE's moving sides
      const auto& move = stack.template get< tag::ale, tag::move >();
      for (const auto& s : move) {
        const auto& f = s.template get< tag::fn >();
        if (f.empty() or f.size() % 4 != 0)
          Message< Stack, ERROR, MsgKey::INCOMPLETEUSERFN>( stack, in );
      }
    }
  };

  //! Rule used to trigger action
  template< typename Feature >
  struct enable : pegtl::success {};
  //! Enable adaptive mesh refinement (AMR)
  template< typename Feature >
  struct action< enable< Feature > > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& stack ) {
      stack.template get< Feature, Feature >() = true;
    }
  };

  //! Rule used to trigger action
  struct compute_refvar_idx : pegtl::success {};
  //! Compute indices of refinement variables
  //! \details This functor computes the indices in the unknown vector for all
  //!   refinement variables in the system of systems of dependent variables
  //!   after the refvar...end block has been parsed in the amr...end block.
  //!   After basic error checking, the vector at stack.get<tag::amr,tag::id>()
  //!   is filled.
  template<>
  struct action< compute_refvar_idx > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // reference variables just parsed by refvar...end block
      const auto& refvar = stack.template get< tag::amr, tag::refvar >();
      // get ncomponents object from this input deck
      const auto& ncomps = stack.template get< tag::component >();
      // compute number of components associated to dependent variabels
      auto ncompmap = ncomps.ncompmap( stack );
      // reference variable index vector to fill
      auto& refidx = stack.template get< tag::amr, tag::id >();
      // Compute indices for all refvars
      for (const auto& v : refvar) {    // for all reference variables parsed
        // depvar is the first char of a refvar
        auto depvar = v[0];
        // the field ID is optional and is the rest of the depvar string
        std::size_t f = (v.size()>1 ? std::stoul(v.substr(1)) : 1) - 1;
        // field ID must be less than or equal to the number of scalar
        // components configured for the eq system for this dependent variable
        if (f >= tk::cref_find( ncompmap, depvar ))
          Message< Stack, ERROR, MsgKey::NOSUCHCOMPONENT >( stack, in );
        // the index is the eq field ID
        auto idx = f;
        // save refvar index in system of all systems
        refidx.push_back( idx );
      }
    }
  };

  //! Rule used to trigger action
  struct check_amr_errors : pegtl::success {};
  //! Do error checking for the amr...end block
  //! \details This is error checking that only the amr...end block
  //!   must satisfy. Besides error checking this can also set defaults
  //!   as this block is called when parsing of a amr...end block has
  //!   just finished.
  template<>
  struct action< check_amr_errors > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // Error out if refvar size does not equal refidx size (programmer error)
      Assert( (stack.template get< tag::amr, tag::refvar >().size() ==
               stack.template get< tag::amr, tag::id >().size()),
              "The size of refvar and refidx vectors must equal" );
      const auto& initref = stack.template get< tag::amr, tag::init >();
      const auto& refvar = stack.template get< tag::amr, tag::refvar >();
      const auto& edgelist = stack.template get< tag::amr, tag::edge >();
      // Error out if initref edge list is not divisible by 2 (user error)
      if (edgelist.size() % 2 == 1)
        Message< Stack, ERROR, MsgKey::T0REFODD >( stack, in );
      // Warn if initial AMR will be a no-op
      if ( stack.template get< tag::amr, tag::t0ref >() && initref.empty() )
        Message< Stack, WARNING, MsgKey::T0REFNOOP >( stack, in );
      // Error out if timestepping AMR will be a no-op (user error)
      if ( stack.template get< tag::amr, tag::dtref >() && refvar.empty() )
        Message< Stack, ERROR, MsgKey::DTREFNOOP >( stack, in );
      // Error out if mesh refinement frequency is zero (programmer error)
      Assert( (stack.template get< tag::amr, tag::dtfreq >() > 0),
              "Mesh refinement frequency must be positive" );
    }
  };

  //! Rule used to trigger action
  struct check_pref_errors : pegtl::success {};
  //! Do error checking for the pref...end block
  template<>
  struct action< check_pref_errors > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      auto& tolref = stack.template get< tag::pref, tag::tolref >();
      if (tolref < 0.0 || tolref > 1.0)
        Message< Stack, ERROR, MsgKey::PREFTOL >( stack, in );
    }
  };

  //! Rule used to trigger action
  struct match_pointname : pegtl::success {};
  //! \brief Match PDF name to the registered ones
  //! \details This is used to check the set of PDF names dependent previously
  //!    registered to make sure all are unique.
  template<>
  struct action< match_pointname > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using inciter::deck::pointnames;
      // find matched name in set of registered ones
      if (pointnames.find( in.string() ) == pointnames.end()) {
        pointnames.insert( in.string() );
        stack.template get< tag::history, tag::id >().push_back( in.string() );
      }
      else  // error out if name matched var is already registered
        Message< Stack, ERROR, MsgKey::POINTEXISTS >( stack, in );
    }
  };

  //! Rule used to trigger action
  struct push_depvar : pegtl::success {};
  //! Add matched outvar based on depvar into vector of outvars
  //! \details Push outvar based on depvar: use first char of matched token as
  //! OutVar::var, OutVar::name = "" by default. OutVar::name being empty will
  //! be used to differentiate a depvar-based outvar from a human-readable
  //! outvar. Depvar-based outvars can directly access solution arrays using
  //! their field. Human-readable outvars need a mechanism (a function) to read
  //! and compute their variables from solution arrays. The 'getvar' function,
  //! used to compute a physics variable from the numerical solution is assigned
  //! after initial migration and thus not assigned here (during parsing).
  template<>
  struct action< push_depvar > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using inciter::deck::centering;
      using inciter::ctr::OutVar;
      auto& vars = stack.template get< tag::cmd, tag::io, tag::outvar >();
      vars.emplace_back(OutVar(in.string()[0], field, centering));
      field = 0;        // reset field
    }
  };

  //! Rule used to trigger action
  struct push_matvar : pegtl::success {};
  //! Add matched outvar based on matvar into vector of outvars
  //! \details Push outvar based on matvar: use depvar char as OutVar::var,
  //! OutVar::name = "" by default. Matvar-based outvars are similar to
  //! depvar-base outvars, in that the OutVar has empty name and the OutVar::var
  //! is a depvar, but instead of having the user try to guess the field id, the
  //! grammar accepts a physics label (accepted multimatvars) and a material
  //! index, which are then converted to a depvar + field index.
  template<>
  struct action< push_matvar > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using tag::param;
      using tag::multimat;
      using inciter::deck::centering;
      using inciter::deck::multimatvars;
      using inciter::ctr::OutVar;
      auto& vars = stack.template get< tag::cmd, tag::io, tag::outvar >();
      auto nmat = stack.template get< param, multimat, tag::nmat >();
      auto depvar = stack.template get< param, multimat, tag::depvar >().back();
      // first char of matched token: accepted multimatvar label char
      char v = static_cast<char>( in.string()[0] );
      // Since multimat outvars are configured based on an acceptable character
      // label (see inciter::deck::multimatvars, denoting a physics variable),
      // and a material index (instead of a depvar + a component index),
      // multimat material bounds are checked here. Note that for momentum and
      // velocity, the field id is the spatial direction not the material id.
      // Also note that field (in grammar's state) starts from 0.
      if ( ((v=='u'||v=='U'||v=='m'||v=='M') && field>2) ||
           ((v!='u'&&v!='U'&&v!='m'&&v!='M') && field>=nmat) )
        Message< Stack, ERROR, MsgKey::NOSUCHCOMPONENT >( stack, in );
      // field contains material id, compute multiat component index
      auto comp = tk::cref_find( multimatvars, v )( nmat, field );
      // save depvar + component index based on physics label + material id,
      // also save physics label + material id as matvar
      vars.emplace_back(
        OutVar(depvar, comp, centering, {}, {}, v+std::to_string(field+1)) );
      field = 0;        // reset field
    }
  };

  //! Function object for adding a human-readable output variable
  //! \details Since human-readable outvars do not necessarily have any
  //! reference to the depvar of their system they refer to, nor which system
  //! they refer to, we configure them for all of the systems they are preceded
  //! by. If there is only a single system of the type the outvar is configured,
  //! we simply look up the depvar and use that as OutVar::var. If there are
  //! multiple systems configured upstream to which the outvar could refer to,
  //! we configure an outvar for all systems configured, and postfix the
  //! human-readable OutVar::name with '_' + depvar. Hence this function object
  //! so the code below can be invoked for all equation types.
  template< typename Stack >
  struct AddOutVarHuman {
    Stack& stack;
    const std::string& in_string;
    explicit AddOutVarHuman( Stack& s, const std::string& ins )
      : stack(s), in_string(ins) {}
    template< typename Eq > void operator()( brigand::type_<Eq> ) {
      using inciter::deck::centering;
      using inciter::ctr::OutVar;
      const auto& depvar = stack.template get< tag::param, Eq, tag::depvar >();
      auto& vars = stack.template get< tag::cmd, tag::io, tag::outvar >();
      if (depvar.size() == 1)
        vars.emplace_back( OutVar( depvar[0], 0, centering, in_string ) );
      else
        for (auto d : depvar)
          vars.emplace_back( OutVar( d, 0, centering, in_string + '_' + d ) );
    }
  };

  //! Rule used to trigger action
  struct push_humanvar : pegtl::success {};
  //! Add matched outvar based on depvar into vector of vector of outvars
  //! \details Push outvar based on human readable string for which
  //! OutVar::name = matched token. OutVar::name being not empty will be used to
  //! differentiate a depvar-based outvar from a human-readable outvar.
  //! Depvar-based outvars can directly access solution arrays using their
  //! field. Human-readable outvars need a mechanism (a function) to read and
  //! compute their variables from solution arrays. The 'getvar' function, used
  //! to compute a physics variable from the numerical solution is assigned
  //! after initial migration and thus not assigned here (during parsing).
  template<>
  struct action< push_humanvar > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      brigand::for_each< inciter::ctr::parameters::Keys >
                       ( AddOutVarHuman< Stack >( stack, in.string() ) );
    }
  };

  //! Rule used to trigger action
  struct set_outvar_alias : pegtl::success {};
  //! Set alias of last pushed output variable
  template<>
  struct action< set_outvar_alias > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // Set alias of last pushed outvar:
      auto& vars = stack.template get< tag::cmd, tag::io, tag::outvar >();
      if (!vars.empty()) vars.back().alias = in.string();
    }
  };

  //! Function object for error checking outvar bounds for each equation type
  template< typename Stack >
  struct OutVarBounds {
    const Stack& stack;
    bool& inbounds;
    explicit OutVarBounds( const Stack& s, bool& i )
      : stack(s), inbounds(i) { inbounds = false; }
    template< typename U > void operator()( brigand::type_<U> ) {
      if (std::is_same_v< U, tag::multimat >) inbounds = true;  // skip multimat
      const auto& depvar = stack.template get< tag::param, U, tag::depvar >();
      const auto& ncomp = stack.template get< tag::component, U >();
      Assert( depvar.size() == ncomp.size(), "Size mismatch" );
      // called after matching each outvar, so only check the last one
      auto& vars = stack.template get< tag::cmd, tag::io, tag::outvar >();
      const auto& last_outvar = vars.back();
      const auto& v = static_cast<char>( std::tolower(last_outvar.var) );
      for (std::size_t e=0; e<depvar.size(); ++e)
        if (v == depvar[e] && last_outvar.field < ncomp[e]) inbounds = true;
    }
  };

  //! Rule used to trigger action
  struct set_centering : pegtl::success {};
  //! Set variable centering in parser's state
  template<>
  struct action< set_centering > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& ) {
      inciter::deck::centering =
        (in.string() == "node") ? tk::Centering::NODE : tk::Centering::ELEM;
    }
  };

  //! Rule used to trigger action
  struct match_outvar : pegtl::success {};
  //! Match output variable based on depvar
  template<>
  struct action< match_outvar > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using inciter::deck::neq;
      using inciter::deck::multimatvars;
      // convert matched string to char
      auto var = stack.template convert< char >( in.string() );
      if (neq.get< tag::multimat >() == 0) {    // if not multimat
        // find matched variable in set of selected ones
        if (depvars.find(var) != end(depvars))
          action< push_depvar >::apply( in, stack );
        else  // error out if matched var is not selected
          Message< Stack, ERROR, MsgKey::NOSUCHOUTVAR >( stack, in );
      } else {    // if multimat
        // find matched variable in set accepted for multimat
        if (multimatvars.find(var) != end(multimatvars))
          action< push_matvar >::apply( in, stack );
        else
          Message< Stack, ERROR, MsgKey::NOSUCHMULTIMATVAR >( stack, in );
      }
    }
  };

  //! Rule used to trigger action
  struct couple_mesh : pegtl::success {};
  //! Couple solvers to use mesh2mesh transfer
  template<> struct action< couple_mesh > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& stack ) {
      using namespace inciter::ctr;
      if (inciter::deck::depvar_cnt != 'a') {   // dst only
        // set number of scalar components
        stack.template get< tag::component, tag::compflow >().push_back( 5 );
        // configure physics policy
        stack.template get< tag::param, tag::compflow, tag::physics >().
          push_back( PhysicsType::EULER );
        // add new PDE to instantiate
        stack.template get< tag::selected, tag::pde >().
          push_back( PDEType::COMPFLOW );
        // configure its problem policy
         stack.template get< tag::param, tag::compflow, tag::problem >().
           push_back( ProblemType::USER_DEFINED );
        // setup coupling: src:0 -> dst:1...N
        auto& t = stack.template get< tag::couple, tag::transfer >();
        t.emplace_back( 0, ++inciter::deck::transfer_cnt );
      }
      // add depvar to deck::depvars so it can be selected as outvar later
      tk::grm::depvars.insert( inciter::deck::depvar_cnt );
      // assign depvar
      stack.template get< tag::param, tag::compflow, tag::depvar >().
        push_back( inciter::deck::depvar_cnt++ );
    }
  };

} // ::grm
} // ::tk

namespace inciter {

//! Inciter input deck facilitating user input for computing shock hydrodynamics
namespace deck {

  using namespace tao;

  // Inciter's InputDeck grammar

  //! scan and store_back equation keyword and option
  template< typename keyword, class eq >
  struct scan_eq :
         tk::grm::scan< typename keyword::pegtl_string,
                        tk::grm::store_back_option< use,
                                                    ctr::PDE,
                                                    tag::selected,
                                                    tag::pde > > {};

  //! Error checks after an equation...end block has been parsed
  template< class eq, template< class > class eqchecker >
  struct check_errors :
         pegtl::seq<
           // register differential equation block
           tk::grm::register_inciter_eq< eq >,
           // do error checking on this block
           eqchecker< eq > > {};

  //! Match discretization option
  template< template< class > class use, class keyword, class Option,
            class Tag >
  struct discroption :
         tk::grm::process< use< keyword >,
                           tk::grm::store_inciter_option<
                             Option, tag::discr, Tag >,
                           pegtl::alpha > {};

  //! Discretization parameters
  struct discretization :
         pegtl::sor<
           tk::grm::discrparam< use, kw::nstep, tag::nstep >,
           tk::grm::discrparam< use, kw::term, tag::term >,
           tk::grm::discrparam< use, kw::t0, tag::t0 >,
           tk::grm::discrparam< use, kw::dt, tag::dt >,
           tk::grm::discrparam< use, kw::cfl, tag::cfl >,
           tk::grm::discrparam< use, kw::residual, tag::residual >,
           tk::grm::discrparam< use, kw::rescomp, tag::rescomp >,
           tk::grm::process< use< kw::fcteps >,
                             tk::grm::Store< tag::discr, tag::fcteps > >,
           tk::grm::process< use< kw::fctclip >,
                             tk::grm::Store< tag::discr, tag::fctclip >,
                             pegtl::alpha >,
           tk::grm::process< use< kw::fct >,
                             tk::grm::Store< tag::discr, tag::fct >,
                             pegtl::alpha >,
           tk::grm::process< use< kw::ctau >,
                             tk::grm::Store< tag::discr, tag::ctau > >,
           tk::grm::process< use< kw::pelocal_reorder >,
                             tk::grm::Store< tag::discr, tag::pelocal_reorder >,
                             pegtl::alpha >,
           tk::grm::process< use< kw::operator_reorder >,
                             tk::grm::Store< tag::discr, tag::operator_reorder >,
                             pegtl::alpha >,
           tk::grm::process< use< kw::steady_state >,
                             tk::grm::Store< tag::discr, tag::steady_state >,
                             pegtl::alpha >,
           tk::grm::interval_iter< use< kw::ttyi >,
                                   tag::output, tag::iter, tag::tty >,
           tk::grm::process_alpha< use< kw::scheme >,
                                   tk::grm::store_inciter_option<
                                     inciter::ctr::Scheme,
                                     tag::discr,
                                     tag::scheme >,
                                   tk::grm::configure_scheme >,
           discroption< use, kw::limiter, inciter::ctr::Limiter, tag::limiter >,
           tk::grm::discrparam< use, kw::cweight, tag::cweight >,
           tk::grm::process< use< kw::accuracy_test >,
                             tk::grm::Store< tag::discr, tag::accuracy_test >,
                             pegtl::alpha >,
           tk::grm::process< use< kw::limsol_projection >,
             tk::grm::Store< tag::discr, tag::limsol_projection >,
             pegtl::alpha >,
           tk::grm::discrparam< use, kw::shock_detector_coeff,
             tag::shock_detector_coeff >
         > {};

  //! PDE parameter vector
  template< class keyword, class eq, class param, class... xparams >
  struct pde_parameter_vector :
         tk::grm::parameter_vector< use,
                                    use< keyword >,
                                    tk::grm::Store_back_back,
                                    tk::grm::start_vector,
                                    tk::grm::check_vector,
                                    eq, param, xparams... > {};

   //! Match box parameter
  template< class eq, typename keyword, typename target, typename icunit >
  struct box_parameter :
           pegtl::if_must<
             tk::grm::readkw< typename use<keyword>::pegtl_string >,
             tk::grm::scan<
               pegtl::sor< tk::grm::number,
                           tk::grm::msg< tk::grm::ERROR,
                                         tk::grm::MsgKey::MISSING > >,
               tk::grm::Back_back_store< target,
                 tag::param, eq, tag::ic, icunit > > > {};

   //! Match box parameter and store deep
  template< class eq, typename keyword, typename target, typename subtarget,
    typename icunit >
  struct box_deep_parameter :
           pegtl::if_must<
             tk::grm::readkw< typename use<keyword>::pegtl_string >,
             tk::grm::scan<
               pegtl::sor< tk::grm::number,
                           tk::grm::msg< tk::grm::ERROR,
                                         tk::grm::MsgKey::MISSING > >,
               tk::grm::Back_back_deep_store< target, subtarget,
                 tag::param, eq, tag::ic, icunit > > > {};

   //! Match box parameter vector
  template< class eq, typename keyword, typename target, typename icunit >
  struct box_vector :
         tk::grm::vector< use< keyword >,
                          tk::grm::Back_back_store_back< target,
                            tag::param, eq, tag::ic, icunit >,
                          use< kw::end > > {};

   //! Match box parameter vector and store deep
  template< class eq, typename keyword, typename target, typename subtarget,
    typename icunit >
  struct box_deep_vector :
         tk::grm::vector< use< keyword >,
                          tk::grm::Back_back_deep_store_back< target, subtarget,
                            tag::param, eq, tag::ic, icunit >,
                          use< kw::end > > {};


   //! Match box option
  template< class eq, typename Option, typename keyword, typename target,
            typename subtarget, typename icunit >
  struct box_option :
         tk::grm::process<
           use< keyword >,
           tk::grm::back_back_deep_store_option< target, subtarget, use,
             Option, tag::param, eq, tag::ic, icunit >,
           pegtl::alpha > {};

   //! Match material option
  template< class eq, typename Option, typename keyword, typename target >
  struct material_option :
         tk::grm::process<
           use< keyword >,
           tk::grm::back_store_option< target, use, Option,
             tag::param, eq, tag::material >,
           pegtl::alpha > {};

  //! Match material parameter vector
  template< class eq, typename keyword, typename target >
  struct material_vector :
         tk::grm::vector< use< keyword >,
                          tk::grm::Back_store_back< target,
                            tag::param, eq, tag::material >,
                          use< kw::end > > {};

  //! put in PDE parameter for equation matching keyword
  template< typename eq, typename keyword, typename param,
            class kw_type = tk::grm::number >
  struct parameter :
         tk::grm::process< use< keyword >,
                           tk::grm::Store_back< tag::param, eq, param >,
                           kw_type > {};

  //! put in PDE parameter for equation matching keyword
  template< typename eq, typename keyword, typename param,
            class kw_type = tk::grm::number >
  struct store_parameter :
         tk::grm::process< use< keyword >,
                           tk::grm::Store< tag::param, eq, param >,
                           kw_type > {};

  //! put in PDE bool parameter for equation matching keyword into vector< int >
  template< typename eq, typename keyword, typename p >
  struct parameter_bool :
         tk::grm::process< use< keyword >,
                           tk::grm::Store_bool< tag::param, eq, p >,
                           pegtl::alpha > {};

  //! Boundary conditions block
  template< class keyword, class eq, class param >
  struct bc :
         pegtl::if_must<
           tk::grm::readkw< typename use< keyword >::pegtl_string >,
           tk::grm::block<
             use< kw::end >,
             tk::grm::parameter_vector< use,
                                        use< kw::sideset >,
                                        tk::grm::Store_back,
                                        tk::grm::noop,
                                        tk::grm::check_vector,
                                        eq, tag::bc, param > > > {};

  //! Match user-defined function as a discrete list of real numbers
  template< class target, template< class... > class insert, class tag,
    class... tags >
  struct user_fn :
         pegtl::if_must<
           tk::grm::readkw< use< kw::fn >::pegtl_string >,
           tk::grm::block< use< kw::end >,
             tk::grm::scan< tk::grm::number,
               insert< target, tag, tags... > > > > {};

  //! User defined time dependent BC bc_timedep...end block
  template< class eq >
  struct timedep_bc :
         pegtl::if_must<
           tk::grm::readkw< use< kw::bc_timedep >::pegtl_string >,
           tk::grm::start_vector_back< tag::param, eq, tag::bctimedep >,
           tk::grm::block< use< kw::end >,
             user_fn< tag::fn, tk::grm::Back_back_store_back, tag::param, eq,
               tag::bctimedep >,
             pegtl::if_must< tk::grm::vector< use< kw::sideset >,
               tk::grm::Back_back_store_back< tag::sideset, tag::param, eq,
                 tag::bctimedep >,
               use< kw::end > > > > > {};

  //! Stagnation boundary conditions block
  template< class eq, class bc, class kwbc >
  struct bc_spec :
         pegtl::if_must<
           tk::grm::readkw< typename kwbc::pegtl_string >,
           tk::grm::block<
             use< kw::end >,
             tk::grm::parameter_vector< use,
                                        use< kw::radius >,
                                        tk::grm::Store_back_back,
                                        tk::grm::start_vector,
                                        tk::grm::check_vector,
                                        eq, bc, tag::radius >,
             tk::grm::parameter_vector< use,
                                        use< kw::point >,
                                        tk::grm::Store_back_back,
                                        tk::grm::start_vector,
                                        tk::grm::check_vector,
                                        eq, bc, tag::point > > > {};

  //! Boundary conditions block
  template< class eq >
  struct sponge :
         pegtl::if_must<
           tk::grm::readkw< typename use< kw::sponge >::pegtl_string >,
           tk::grm::block<
             use< kw::end >,
             tk::grm::parameter_vector< use,
                                        use< kw::velocity >,
                                        tk::grm::Store_back,
                                        tk::grm::noop,
                                        tk::grm::check_vector,
                                        eq, tag::sponge, tag::velocity >,
             tk::grm::parameter_vector< use,
                                        use< kw::pressure >,
                                        tk::grm::Store_back,
                                        tk::grm::noop,
                                        tk::grm::check_vector,
                                        eq, tag::sponge, tag::pressure >,
             tk::grm::parameter_vector< use,
                                        use< kw::sideset >,
                                        tk::grm::Store_back,
                                        tk::grm::noop,
                                        tk::grm::check_vector,
                                        eq, tag::sponge, tag::sideset > > > {};

  //! Farfield boundary conditions block
  template< class keyword, class eq, class param >
  struct farfield_bc :
         pegtl::if_must<
           tk::grm::readkw< typename use< keyword >::pegtl_string >,
           tk::grm::block<
             use< kw::end >,
             store_parameter< eq, kw::pressure, tag::farfield_pressure >,
             store_parameter< eq, kw::density, tag::farfield_density >,
             tk::grm::parameter_vector< use,
                                        use< kw::velocity >,
                                        tk::grm::Store_back,
                                        tk::grm::noop,
                                        tk::grm::check_vector,
                                        eq, tag::farfield_velocity >,
             tk::grm::parameter_vector< use,
                                        use< kw::sideset >,
                                        tk::grm::Store_back,
                                        tk::grm::noop,
                                        tk::grm::check_vector,
                                        eq, tag::bc, param > > > {};

  //! edgelist ... end block
  struct edgelist :
         tk::grm::vector< use< kw::amr_edgelist >,
                          tk::grm::Store_back< tag::amr, tag::edge >,
                          use< kw::end >,
                          tk::grm::check_vector< tag::amr, tag::edge > > {};

  //! xminus configuring coordinate-based edge tagging for mesh refinement
  template< typename keyword, typename Tag >
  struct half_world :
         tk::grm::control< use< keyword >, pegtl::digit, tk::grm::Store,
                           tag::amr, Tag > {};

  //! coords ... end block
  struct coords :
           pegtl::if_must<
             tk::grm::readkw< use< kw::amr_coords >::pegtl_string >,
             tk::grm::block< use< kw::end >,
                             half_world< kw::amr_xminus, tag::xminus >,
                             half_world< kw::amr_xplus, tag::xplus >,
                             half_world< kw::amr_yminus, tag::yminus >,
                             half_world< kw::amr_yplus, tag::yplus >,
                             half_world< kw::amr_zminus, tag::zminus >,
                             half_world< kw::amr_zplus, tag::zplus > > > {};

  //! initial conditins box block
  template< class eq >
  struct box :
         pegtl::if_must<
           tk::grm::readkw< use< kw::box >::pegtl_string >,
           tk::grm::start_vector_back< tag::param, eq, tag::ic, tag::box >,
           tk::grm::block< use< kw::end >
             , box_parameter< eq, kw::xmin, tag::xmin, tag::box >
             , box_parameter< eq, kw::xmax, tag::xmax, tag::box >
             , box_parameter< eq, kw::ymin, tag::ymin, tag::box >
             , box_parameter< eq, kw::ymax, tag::ymax, tag::box >
             , box_parameter< eq, kw::zmin, tag::zmin, tag::box >
             , box_parameter< eq, kw::zmax, tag::zmax, tag::box >
             , box_vector< eq, kw::orientation, tag::orientation, tag::box >
             , box_parameter< eq, kw::materialid, tag::materialid, tag::box >
             , box_parameter< eq, kw::density, tag::density, tag::box >
             , box_parameter< eq, kw::pressure, tag::pressure, tag::box >
             , box_parameter< eq, kw::temperature, tag::temperature, tag::box >
             , box_parameter< eq, kw::energy_content, tag::energy_content, tag::box >
             , box_parameter< eq, kw::energy, tag::energy, tag::box >
             , box_parameter< eq, kw::mass, tag::mass, tag::box >
             , box_vector< eq, kw::velocity, tag::velocity, tag::box >
             , box_option< eq, ctr::Initiate, kw::initiate, tag::initiate,
                           tag::init, tag::box >
             , pegtl::if_must<
                 tk::grm::readkw< use< kw::linear >::pegtl_string >,
                 tk::grm::block< use< kw::end >
                   , box_deep_vector< eq, kw::point, tag::initiate, tag::point,
                      tag::box >
                   , box_deep_parameter< eq, kw::init_time, tag::initiate,
                                         tag::init_time, tag::box >
                   , box_deep_parameter< eq, kw::front_width, tag::initiate,
                                         tag::front_width, tag::box >
                   , box_deep_parameter< eq, kw::velocity, tag::initiate,
                                         tag::velocity, tag::box > > >
             > > {};

  //! initial conditins meshblock block
  template< class eq >
  struct meshblock :
         pegtl::if_must<
           tk::grm::readkw< use< kw::meshblock >::pegtl_string >,
           tk::grm::start_vector_back< tag::param, eq, tag::ic, tag::meshblock >,
           tk::grm::block< use< kw::end >
             , box_parameter< eq, kw::blockid, tag::blockid, tag::meshblock >
             , box_parameter< eq, kw::materialid, tag::materialid,
                tag::meshblock >
             , box_parameter< eq, kw::volume, tag::volume, tag::meshblock >
             , box_parameter< eq, kw::density, tag::density, tag::meshblock >
             , box_parameter< eq, kw::pressure, tag::pressure, tag::meshblock >
             , box_parameter< eq, kw::temperature, tag::temperature,
                tag::meshblock >
             , box_parameter< eq, kw::energy_content, tag::energy_content,
                tag::meshblock >
             , box_parameter< eq, kw::energy, tag::energy, tag::meshblock >
             , box_parameter< eq, kw::mass, tag::mass, tag::meshblock >
             , box_vector< eq, kw::velocity, tag::velocity, tag::meshblock >
             , box_option< eq, ctr::Initiate, kw::initiate, tag::initiate,
                           tag::init, tag::meshblock >
             , pegtl::if_must<
                 tk::grm::readkw< use< kw::linear >::pegtl_string >,
                 tk::grm::block< use< kw::end >
                   , box_deep_vector< eq, kw::point, tag::initiate, tag::point,
                      tag::meshblock >
                   , box_deep_parameter< eq, kw::init_time, tag::initiate,
                                         tag::init_time, tag::meshblock >
                   , box_deep_parameter< eq, kw::front_width, tag::initiate,
                                         tag::front_width, tag::meshblock >
                   , box_deep_parameter< eq, kw::velocity, tag::initiate,
                                         tag::velocity, tag::meshblock > > >
             > > {};

  //! initial conditions block for compressible flow
  template< class eq >
  struct ic :
         pegtl::if_must<
           tk::grm::readkw< use< kw::ic >::pegtl_string >,
           tk::grm::block< use< kw::end >,
             pegtl::sor<
               pegtl::if_must<
                 tk::grm::vector< kw::density,
                   tk::grm::Store_back_back< tag::param, eq, tag::ic,
                                             tag::density >,
                   use< kw::end > > >,
               pegtl::if_must<
                 tk::grm::vector< kw::materialid,
                   tk::grm::Store_back_back< tag::param, eq, tag::ic,
                                             tag::materialid >,
                   use< kw::end > > >,
               pegtl::if_must<
                 tk::grm::vector< kw::velocity,
                   tk::grm::Store_back_back< tag::param, eq, tag::ic,
                                             tag::velocity >,
                   use< kw::end > > >,
               pegtl::if_must<
                 tk::grm::vector< kw::pressure,
                   tk::grm::Store_back_back< tag::param, eq, tag::ic,
                                             tag::pressure >,
                   use< kw::end > > >,
               pegtl::if_must<
                 tk::grm::vector< kw::temperature,
                   tk::grm::Store_back_back< tag::param, eq, tag::ic,
                                             tag::temperature >,
                   use< kw::end > > >,
               pegtl::if_must<
                 tk::grm::vector< kw::energy,
                   tk::grm::Store_back_back< tag::param, eq, tag::ic,
                                             tag::energy >,
                   use< kw::end > > >,
               pegtl::seq< box< eq > >,
               pegtl::seq< meshblock< eq > >
           > > > {};

  //! put in material property for equation matching keyword
  template< typename eq, typename keyword, typename property >
  struct material_property :
         pde_parameter_vector< keyword, eq, property > {};

  //! Material properties block for compressible flow
  template< class eq >
  struct material_properties :
         pegtl::seq<
          pegtl::if_must<
            tk::grm::readkw< use< kw::material >::pegtl_string >,
            tk::grm::block< use< kw::end >,
                material_vector< eq, kw::id, tag::id >
              , material_vector< eq, kw::mat_gamma, tag::gamma >
              , material_vector< eq, kw::mat_mu, tag::mu >
              , material_vector< eq, kw::mat_pstiff, tag::pstiff >
              , material_vector< eq, kw::w_gru, tag::w_gru >
              , material_vector< eq, kw::A_jwl, tag::A_jwl >
              , material_vector< eq, kw::B_jwl, tag::B_jwl >
              , material_vector< eq, kw::C_jwl, tag::C_jwl >
              , material_vector< eq, kw::R1_jwl, tag::R1_jwl >
              , material_vector< eq, kw::R2_jwl, tag::R2_jwl >
              , material_vector< eq, kw::rho0_jwl, tag::rho0_jwl >
              , material_vector< eq, kw::de_jwl, tag::de_jwl >
              , material_vector< eq, kw::rhor_jwl, tag::rhor_jwl >
              , material_vector< eq, kw::Tr_jwl, tag::Tr_jwl >
              , material_vector< eq, kw::Pr_jwl, tag::Pr_jwl >
              , material_vector< eq, kw::mat_cv, tag::cv >
              , material_vector< eq, kw::mat_k, tag::k >
              , material_option< eq, ctr::Material, kw::eos, tag::eos >
            > > > {};

  //! mesh ... end block
  struct mesh :
         pegtl::if_must<
           tk::grm::readkw< use< kw::mesh >::pegtl_string >,
           tk::grm::block< use< kw::end >,
             pegtl::if_must<
               tk::grm::filename< use, tag::param, tag::compflow, tag::mesh,
                                  tag::filename >,
               tk::grm::couple_mesh > > > {};

  //! transport equation for scalars
  struct transport :
         pegtl::if_must<
           scan_eq< use< kw::transport >, tag::transport >,
           tk::grm::block< use< kw::end >,
                           tk::grm::policy< use,
                                            use< kw::physics >,
                                            ctr::Physics,
                                            tag::transport,
                                            tag::physics >,
                           tk::grm::policy< use,
                                            use< kw::problem >,
                                            ctr::Problem,
                                            tag::transport,
                                            tag::problem >,
                           tk::grm::depvar< use,
                                            tag::transport,
                                            tag::depvar >,
                           tk::grm::component< use< kw::ncomp >,
                                               tag::transport >,
                           tk::grm::parameter_vector< use,
                                                      use<kw::pde_diffusivity>,
                                                      tk::grm::Store_back,
                                                      tk::grm::noop,
                                                      tk::grm::check_vector,
                                                      tag::transport,
                                                      tag::diffusivity >,
                           tk::grm::parameter_vector< use,
                                                      use< kw::pde_lambda >,
                                                      tk::grm::Store_back,
                                                      tk::grm::noop,
                                                      tk::grm::check_vector,
                                                      tag::transport,
                                                      tag::lambda >,
                           tk::grm::parameter_vector< use,
                                                      use< kw::pde_u0 >,
                                                      tk::grm::Store_back,
                                                      tk::grm::noop,
                                                      tk::grm::check_vector,
                                                      tag::transport,
                                                      tag::u0 >,
                           bc< kw::bc_dirichlet, tag::transport, tag::bcdir >,
                           bc< kw::bc_sym, tag::transport, tag::bcsym >,
                           bc< kw::bc_inlet, tag::transport, tag::bcinlet >,
                           bc< kw::bc_outlet, tag::transport, tag::bcoutlet >,
                           bc< kw::bc_extrapolate, tag::transport,
                               tag::bcextrapolate >,
                           parameter< tag::transport,
                                      kw::intsharp_param,
                                      tag::intsharp_param >,
                           parameter< tag::transport,
                                      kw::intsharp,
                                      tag::intsharp > >,
           check_errors< tag::transport, tk::grm::check_transport > > {};

  //! compressible flow
  struct compflow :
         pegtl::if_must<
           scan_eq< use< kw::compflow >, tag::compflow >,
           tk::grm::start_vector< tag::param, tag::compflow, tag::ic,
                                  tag::density >,
           tk::grm::start_vector< tag::param, tag::compflow, tag::ic,
                                  tag::velocity >,
           tk::grm::start_vector< tag::param, tag::compflow, tag::ic,
                                  tag::pressure >,
           tk::grm::start_vector< tag::param, tag::compflow, tag::ic,
                                  tag::temperature >,
           tk::grm::start_vector< tag::param, tag::compflow, tag::ic,
                                  tag::energy >,
           tk::grm::start_vector< tag::param, tag::compflow, tag::ic,
                                  tag::box >,
           tk::grm::start_vector< tag::param, tag::compflow, tag::ic,
                                  tag::meshblock >,
           tk::grm::start_vector< tag::param, tag::compflow, tag::material >,
           tk::grm::start_vector< tag::param, tag::compflow, tag::bctimedep >,
           tk::grm::block< use< kw::end >,
                           tk::grm::policy< use,
                                            use< kw::physics >,
                                            ctr::Physics,
                                            tag::compflow,
                                            tag::physics >,
                           tk::grm::policy< use,
                                            use< kw::problem >,
                                            ctr::Problem,
                                            tag::compflow,
                                            tag::problem >,
                           mesh,
                           tk::grm::process<
                             use< kw::flux >,
                               tk::grm::store_back_option< use,
                                                           ctr::Flux,
                                                           tag::param,
                                                           tag::compflow,
                                                           tag::flux >,
                             pegtl::alpha >,
                           ic< tag::compflow >,
                           tk::grm::lua< use, tag::param, tag::compflow >,
                           material_properties< tag::compflow >,
                           parameter< tag::compflow, kw::sysfctvar,
                                      tag::sysfctvar >,
                           parameter_bool< tag::compflow, kw::sysfct,
                                           tag::sysfct >,
                           store_parameter< tag::compflow, kw::pde_alpha,
                                            tag::alpha >,
                           store_parameter< tag::compflow, kw::pde_p0,
                                            tag::p0 >,
                           store_parameter< tag::compflow, kw::pde_betax,
                                            tag::betax >,
                           store_parameter< tag::compflow, kw::pde_betay,
                                            tag::betay >,
                           store_parameter< tag::compflow, kw::pde_betaz,
                                            tag::betaz >,
                           store_parameter< tag::compflow, kw::pde_beta,
                                            tag::beta >,
                           store_parameter< tag::compflow, kw::pde_r0,
                                            tag::r0 >,
                           store_parameter< tag::compflow, kw::pde_ce,
                                            tag::ce >,
                           store_parameter< tag::compflow, kw::pde_kappa,
                                            tag::kappa >,
                           bc< kw::bc_dirichlet, tag::compflow, tag::bcdir >,
                           bc< kw::bc_sym, tag::compflow, tag::bcsym >,
                           bc_spec< tag::compflow, tag::stag, kw::bc_stag >,
                           bc_spec< tag::compflow, tag::skip, kw::bc_skip >,
                           bc< kw::bc_inlet, tag::compflow, tag::bcinlet >,
                           sponge< tag::compflow >,
                           farfield_bc< kw::bc_farfield,
                                        tag::compflow,
                                        tag::bcfarfield >,
                           bc< kw::bc_extrapolate, tag::compflow,
                               tag::bcextrapolate >,
                           timedep_bc< tag::compflow >
                           >,
           check_errors< tag::compflow, tk::grm::check_compflow > > {};

  //! compressible multi-material flow
  struct multimat :
         pegtl::if_must<
           scan_eq< use< kw::multimat >, tag::multimat >,
           tk::grm::start_vector< tag::param, tag::multimat, tag::ic,
                                  tag::density >,
           tk::grm::start_vector< tag::param, tag::multimat, tag::ic,
                                  tag::materialid >,
           tk::grm::start_vector< tag::param, tag::multimat, tag::ic,
                                  tag::velocity >,
           tk::grm::start_vector< tag::param, tag::multimat, tag::ic,
                                  tag::pressure >,
           tk::grm::start_vector< tag::param, tag::multimat, tag::ic,
                                  tag::temperature >,
           tk::grm::start_vector< tag::param, tag::multimat, tag::ic,
                                  tag::energy >,
           tk::grm::start_vector< tag::param, tag::multimat, tag::ic,
                                  tag::box >,
           tk::grm::start_vector< tag::param, tag::multimat, tag::ic,
                                  tag::meshblock >,
           tk::grm::start_vector< tag::param, tag::multimat, tag::material >,
           tk::grm::block< use< kw::end >,
                           tk::grm::policy< use,
                                            use< kw::physics >,
                                            ctr::Physics,
                                            tag::multimat,
                                            tag::physics >,
                           tk::grm::policy< use,
                                            use< kw::problem >,
                                            ctr::Problem,
                                            tag::multimat,
                                            tag::problem >,
                           store_parameter<tag::multimat, kw::nmat, tag::nmat>,
                           tk::grm::process<
                             use< kw::flux >,
                               tk::grm::store_back_option< use,
                                                           ctr::Flux,
                                                           tag::param,
                                                           tag::multimat,
                                                           tag::flux >,
                             pegtl::alpha >,
                           ic< tag::multimat >,
                           material_properties< tag::multimat >,
                           store_parameter< tag::multimat,
                                            kw::pde_alpha,
                                            tag::alpha >,
                           store_parameter< tag::multimat,
                                            kw::pde_p0,
                                            tag::p0 >,
                           store_parameter< tag::multimat,
                                            kw::pde_beta,
                                            tag::beta >,
                           bc< kw::bc_dirichlet,
                               tag::multimat,
                               tag::bcdir >,
                           bc< kw::bc_sym,
                               tag::multimat,
                               tag::bcsym >,
                           bc< kw::bc_inlet,
                               tag::multimat,
                               tag::bcinlet >,
                           bc< kw::bc_outlet,
                               tag::multimat,
                               tag::bcoutlet >,
                           bc< kw::bc_extrapolate,
                               tag::multimat,
                               tag::bcextrapolate >,
                           farfield_bc< kw::bc_farfield,
                                        tag::multimat,
                                        tag::bcfarfield >,
                           parameter< tag::multimat,
                                      kw::prelax_timescale,
                                      tag::prelax_timescale >,
                           parameter< tag::multimat,
                                      kw::prelax,
                                      tag::prelax >,
                           parameter< tag::multimat,
                                      kw::intsharp_param,
                                      tag::intsharp_param >,
                           parameter< tag::multimat,
                                      kw::intsharp,
                                      tag::intsharp > >,
           check_errors< tag::multimat, tk::grm::check_multimat > > {};

  //! partitioning ... end block
  struct partitioning :
         pegtl::if_must<
           tk::grm::readkw< use< kw::partitioning >::pegtl_string >,
           tk::grm::block< use< kw::end >,
                           tk::grm::process<
                             use< kw::algorithm >,
                             tk::grm::store_inciter_option<
                               tk::ctr::PartitioningAlgorithm,
                               tag::selected,
                               tag::partitioner >,
                             pegtl::alpha > > > {};

  //! equation types
  struct equations :
         pegtl::sor< transport, compflow, multimat > {};

  //! refinement variable(s) (refvar) ... end block
  struct refvars :
         pegtl::if_must<
           tk::grm::vector< use< kw::amr_refvar >,
                            tk::grm::match_depvar<
                              tk::grm::Store_back< tag::amr, tag::refvar > >,
                            use< kw::end >,
                            tk::grm::check_vector< tag::amr, tag::refvar >,
                            tk::grm::fieldvar< pegtl::alpha > >,
           tk::grm::compute_refvar_idx > {};

  //! adaptive mesh refinement (AMR) amr...end block
  struct amr :
         pegtl::if_must<
           tk::grm::readkw< use< kw::amr >::pegtl_string >,
           // enable AMR if amr...end block encountered
           tk::grm::enable< tag::amr >,
           tk::grm::block< use< kw::end >,
                           refvars,
                           edgelist,
                           coords,
                           tk::grm::process<
                             use< kw::amr_initial >,
                             tk::grm::store_back_option< use,
                                                         ctr::AMRInitial,
                                                         tag::amr,
                                                         tag::init >,
                             pegtl::alpha >,
                           tk::grm::process<
                             use< kw::amr_error >,
                             tk::grm::store_inciter_option<
                               ctr::AMRError,
                               tag::amr, tag::error >,
                             pegtl::alpha >,
                           tk::grm::control< use< kw::amr_tolref >,
                                             pegtl::digit,
                                             tk::grm::Store,
                                             tag::amr,
                                             tag::tolref >,
                           tk::grm::control< use< kw::amr_tolderef >,
                                             pegtl::digit,
                                             tk::grm::Store,
                                             tag::amr,
                                             tag::tolderef >,
                           tk::grm::process< use< kw::amr_t0ref >,
                             tk::grm::Store< tag::amr, tag::t0ref >,
                             pegtl::alpha >,
                           tk::grm::process< use< kw::amr_dtref_uniform >,
                             tk::grm::Store< tag::amr, tag::dtref_uniform >,
                             pegtl::alpha >,
                           tk::grm::process< use< kw::amr_dtref >,
                             tk::grm::Store< tag::amr, tag::dtref >,
                             pegtl::alpha >,
                           tk::grm::process< use< kw::amr_dtfreq >,
                             tk::grm::Store< tag::amr, tag::dtfreq >,
                             pegtl::digit >,
                           tk::grm::process< use< kw::amr_maxlevels >,
                             tk::grm::Store< tag::amr, tag::maxlevels >,
                             pegtl::digit >
                         >,
           tk::grm::check_amr_errors > {};


  //! Arbitrary-Lagrangian-Eulerian (ALE) move...end block
  struct moving_sides :
         pegtl::if_must<
           tk::grm::readkw< use< kw::move >::pegtl_string >,
           tk::grm::start_vector< tag::ale, tag::move >,
           tk::grm::block< use< kw::end >,
             tk::grm::process<
                use< kw::fntype >,
                tk::grm::back_store_option< tag::fntype,
                                            use,
                                            tk::ctr::UserTable,
                                            tag::ale, tag::move >,
                pegtl::alpha >,
             user_fn< tag::fn, tk::grm::Back_store_back, tag::ale, tag::move >,
             pegtl::if_must< tk::grm::vector< use< kw::sideset >,
               tk::grm::Back_store_back< tag::sideset, tag::ale, tag::move >,
               use< kw::end > > > > > {};

  //! Arbitrary-Lagrangian-Eulerian (ALE) ale...end block
  struct ale :
         pegtl::if_must<
           tk::grm::readkw< use< kw::ale >::pegtl_string >,
           // enable ALE if ale ...end block encountered
           tk::grm::enable< tag::ale >,
           tk::grm::block< use< kw::end >,
              tk::grm::control< use< kw::dvcfl >,
                                pegtl::digit,
                                tk::grm::Store,
                                tag::ale, tag::dvcfl >,
              tk::grm::control< use< kw::vortmult >,
                                pegtl::digit,
                                tk::grm::Store,
                                tag::ale, tag::vortmult >,
              tk::grm::control< use< kw::meshvel_maxit >,
                                pegtl::digit,
                                tk::grm::Store,
                                tag::ale, tag::maxit >,
              tk::grm::control< use< kw::meshvel_tolerance >,
                                pegtl::digit,
                                tk::grm::Store,
                                tag::ale, tag::tolerance >,
              moving_sides,
              tk::grm::process<
                use< kw::meshvelocity >,
                tk::grm::store_inciter_option< ctr::MeshVelocity,
                                               tag::ale, tag::meshvelocity >,
                pegtl::alpha >,
              tk::grm::process<
                use< kw::smoother >,
                tk::grm::store_inciter_option< ctr::MeshVelocitySmoother,
                                               tag::ale, tag::smoother >,
                pegtl::alpha >,
              pegtl::if_must< tk::grm::dimensions< use< kw::mesh_motion >,
                              tk::grm::Store_back< tag::ale, tag::mesh_motion >,
                              use< kw::end > > >,
              pegtl::if_must< tk::grm::vector< use< kw::meshforce >,
                              tk::grm::Store_back< tag::ale, tag::meshforce >,
                              use< kw::end > > >,
              pegtl::if_must<
                tk::grm::readkw< use< kw::bc_dirichlet >::pegtl_string >,
                tk::grm::block< use< kw::end >,
                  pegtl::if_must< tk::grm::vector< use< kw::sideset >,
                                  tk::grm::Store_back< tag::ale, tag::bcdir >,
                                  use< kw::end > > > > >,
              pegtl::if_must<
                tk::grm::readkw< use< kw::bc_sym >::pegtl_string >,
                tk::grm::block< use< kw::end >,
                  pegtl::if_must< tk::grm::vector< use< kw::sideset >,
                                  tk::grm::Store_back< tag::ale, tag::bcsym >,
                                  use< kw::end > > > > > >,
              tk::grm::check_ale > {};

  //! p-adaptive refinement (pref) ...end block
  struct pref :
         pegtl::if_must<
           tk::grm::readkw< use< kw::pref >::pegtl_string >,
           tk::grm::block< use< kw::end >,
                           tk::grm::control< use< kw::pref_tolref >,
                                             pegtl::digit,
                                             tk::grm::Store,
                                             tag::pref,
                                             tag::tolref  >,
                           tk::grm::control< use< kw::pref_ndofmax >,
                                             pegtl::digit,
                                             tk::grm::Store,
                                             tag::pref,
                                             tag::ndofmax >,
                           tk::grm::process<
                             use< kw::pref_indicator >,
                             tk::grm::store_inciter_option<
                               ctr::PrefIndicator,
                               tag::pref, tag::indicator >,
                             pegtl::alpha >
                         >,
           tk::grm::check_pref_errors > {};

  //! Match output variable alias
  struct outvar_alias :
         tk::grm::quoted< tk::grm::set_outvar_alias > {};

  //! Match an output variable in a human readable form: var must be a keyword
  template< class var >
  struct outvar_human :
         tk::grm::exact_scan< use< var >, tk::grm::push_humanvar > {};

  //! Match an output variable based on depvar defined upstream of input file
  struct outvar_depvar :
           tk::grm::scan< tk::grm::fieldvar< pegtl::upper >,
             tk::grm::match_outvar > {};

  //! Parse a centering token and if matches, set centering in parser's state
  struct outvar_centering :
         pegtl::sor<
           tk::grm::exact_scan< use< kw::node >, tk::grm::set_centering >,
           tk::grm::exact_scan< use< kw::elem >, tk::grm::set_centering > > {};

  //! outvar ... end block
  struct outvar_block :
         pegtl::if_must<
           tk::grm::readkw< use< kw::outvar >::pegtl_string >,
           tk::grm::block<
             use< kw::end >
           , outvar_centering
           , outvar_depvar
           , outvar_alias
           , outvar_human< kw::outvar_density >
           , outvar_human< kw::outvar_xmomentum >
           , outvar_human< kw::outvar_ymomentum >
           , outvar_human< kw::outvar_zmomentum >
           , outvar_human< kw::outvar_specific_total_energy >
           , outvar_human< kw::outvar_volumetric_total_energy >
           , outvar_human< kw::outvar_xvelocity >
           , outvar_human< kw::outvar_yvelocity >
           , outvar_human< kw::outvar_zvelocity >
           , outvar_human< kw::outvar_pressure >
           , outvar_human< kw::outvar_material_indicator >
           , outvar_human< kw::outvar_analytic >
           > > {};

  //! field_output ... end block
  struct field_output :
         pegtl::if_must<
           tk::grm::readkw< use< kw::field_output >::pegtl_string >,
           tk::grm::block<
             use< kw::end >,
             outvar_block,
             tk::grm::process< use< kw::filetype >,
                               tk::grm::store_inciter_option<
                                 tk::ctr::FieldFile,
                                 tag::selected,
                                 tag::filetype >,
                               pegtl::alpha >,
             tk::grm::interval_iter< use< kw::interval_iter >,
                                     tag::output, tag::iter, tag::field >,
             tk::grm::interval_time< use< kw::interval_time >,
                                     tag::output, tag::time, tag::field >,
             tk::grm::time_range< use, kw::time_range,
                                  tag::output, tag::range, tag::field >,
             tk::grm::process<
               use< kw::refined >,
               tk::grm::Store< tag::cmd, tag::io, tag::refined >,
               pegtl::alpha >,
             pegtl::if_must<
               tk::grm::vector<
                 use< kw::sideset >,
                 tk::grm::Store_back< tag::cmd, tag::io, tag::surface >,
                 use< kw::end > > > > > {};

  //! history_output ... end block
  struct history_output :
         pegtl::if_must<
           tk::grm::readkw< use< kw::history_output >::pegtl_string >,
           tk::grm::block<
             use< kw::end >,
             outvar_block,
             tk::grm::interval_iter< use< kw::interval_iter >,
               tag::output, tag::iter, tag::history >,
             tk::grm::interval_time< use< kw::interval_time >,
               tag::output, tag::time, tag::history >,
             tk::grm::time_range< use, kw::time_range,
                                  tag::output, tag::range, tag::history >,
             tk::grm::precision< use, tag::history >,
             tk::grm::process<
               use< kw::txt_float_format >,
               tk::grm::store_inciter_option< tk::ctr::TxtFloatFormat,
                                              tag::flformat,
                                              tag::history >,
               pegtl::alpha >,
             pegtl::if_must<
               tk::grm::readkw< use< kw::point >::pegtl_string >,
               tk::grm::act< pegtl::identifier, tk::grm::match_pointname >,
               pegtl::seq<
                 tk::grm::start_vector< tag::history, tag::point >,
                 tk::grm::block<
                   use< kw::end >,
                   tk::grm::scan< tk::grm::number,
                     tk::grm::Store_back_back< tag::history, tag::point > > >
               > > > > {};

  //! 'inciter' block
  struct inciter :
         pegtl::if_must<
           tk::grm::readkw< use< kw::inciter >::pegtl_string >,
           pegtl::sor<
             pegtl::seq< tk::grm::block<
                           use< kw::end >,
                           discretization,
                           equations,
                           amr,
                           ale,
                           pref,
                           partitioning,
                           field_output,
                           history_output,
                           tk::grm::diagnostics<
                             use,
                             tk::grm::store_inciter_option > >,
                         tk::grm::check_inciter >,
            tk::grm::msg< tk::grm::MsgType::ERROR,
                          tk::grm::MsgKey::UNFINISHED > > > {};

  //! \brief All keywords
  struct keywords :
         pegtl::sor< tk::grm::title< use >, inciter > {};

  //! \brief Grammar entry point: parse keywords and ignores until eof
  struct read_file :
         tk::grm::read_file< keywords, tk::grm::ignore > {};

} // deck::
} // inciter::

#endif // InciterInputDeckGrammar_h
