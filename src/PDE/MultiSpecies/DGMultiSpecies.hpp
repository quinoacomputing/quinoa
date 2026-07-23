// *****************************************************************************
/*!
  \file      src/PDE/MultiSpecies/DGMultiSpecies.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Compressible multi-species flow using discontinuous Galerkin
    finite elements
  \details   This file implements calls to the physics operators governing
    compressible multi-species flow using discontinuous Galerkin discretization.
*/
// *****************************************************************************
#ifndef DGMultiSpecies_h
#define DGMultiSpecies_h

#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <map>
#include <array>
#include <iostream>

#include "Macro.hpp"
#include "Exception.hpp"
#include "Vector.hpp"
#include "ContainerUtil.hpp"
#include "UnsMesh.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "Integrate/Basis.hpp"
#include "Integrate/Quadrature.hpp"
#include "Integrate/Initialize.hpp"
#include "Integrate/Surface.hpp"
#include "Integrate/Boundary.hpp"
#include "Integrate/Volume.hpp"
#include "RiemannChoice.hpp"
#include "MultiSpecies/MultiSpeciesIndexing.hpp"
#include "Reconstruction.hpp"
#include "Limiter.hpp"
#include "Problem/FieldOutput.hpp"
#include "Problem/BoxInitialization.hpp"
#include "PrefIndicator.hpp"
#include "MultiSpecies/BCFunctions.hpp"
#include "MultiSpecies/MiscMultiSpeciesFns.hpp"
#include "EoS/GetMatProp.hpp"
#include "Mixture/Mixture.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

namespace dg {

//! \brief MultiSpecies used polymorphically with tk::DGPDE
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Physics - physics configuration, see PDE/MultiSpecies/Physics.h
//!   - Problem - problem configuration, see PDE/MultiSpecies/Problem.h
//! \note The default physics is Euler, which is set in
//!   inciter::LuaParser::storeInputDeck()
template< class Physics, class Problem >
class MultiSpecies {

  private:
    using eq = tag::multispecies;

  public:
    //! Constructor
    explicit MultiSpecies() :
      m_physics(),
      m_ncomp( g_inputdeck.get< tag::ncomp >() ),
      m_nprim(nprim()),
      m_riemann( multispeciesRiemannSolver( g_inputdeck.get< tag::flux >() ) ),
      m_riemannjac( multispeciesRiemannJac( g_inputdeck.get< tag::flux >() ) )
    {
      // associate boundary condition configurations with state functions
      brigand::for_each< ctr::bclist::Keys >( ConfigBC( m_bc,
        // BC State functions
        { dirichlet
        , symmetry
        , invalidBC         // Outlet BC not implemented
        , farfield
        , extrapolate
        , noslipwall
        , symmetry    // Slip-wall equivalent to symmetry without mesh motion
        , isothermal_wall },
        // BC Gradient functions
        { noOpGrad
        , symmetryGrad
        , noOpGrad
        , noOpGrad
        , noOpGrad
        , noOpGrad
        , symmetryGrad
        , noOpGrad }
        ) );

      // EoS initialization
      initializeSpeciesEoS( m_mat_blk );
    }

    //! Find the number of primitive quantities required for this PDE system
    //! \return The number of primitive quantities required to be stored for
    //!   this PDE system
    std::size_t nprim() const
    {
      // mixture temperature
      std::size_t np(1);

      return np;
    }

    //! Find the number of materials set up for this PDE system
    //! \return The number of materials set up for this PDE system
    std::size_t nmat() const { return 1; }

    //! Assign number of DOFs per equation in the PDE system
    //! \param[in,out] numEqDof Array storing number of Dofs for each PDE
    //!   equation
    void numEquationDofs(std::vector< std::size_t >& numEqDof) const
    {
      // all equation-dofs initialized to ndofs
      for (std::size_t i=0; i<m_ncomp; ++i) {
        numEqDof.push_back(g_inputdeck.get< tag::ndof >());
      }
    }

    //! Determine elements that lie inside the user-defined IC box
    //! \param[in] geoElem Element geometry array
    //! \param[in] nielem Number of internal elements
    //! \param[in,out] inbox List of nodes at which box user ICs are set for
    //!    each IC box
    void IcBoxElems( const tk::Fields& geoElem,
      std::size_t nielem,
      std::vector< std::unordered_set< std::size_t > >& inbox ) const
    {
      tk::BoxElems< eq >(geoElem, nielem, inbox);
    }

    //! Find how many 'stiff equations' in this PDE system
    //! \return number of stiff equations. Zero for now, but will need to change
    //!   as chemical non-equilibrium is added
    std::size_t nstiffeq() const
    {
      return 0;
    }

    //! Find how many 'non-stiff equations' in this PDE system
    //! \return number of non-stiff equations
    std::size_t nnonstiffeq() const
    {
      return m_ncomp-nstiffeq();
    }

    //! Locate the stiff equations.
    //! \param[out] stiffEqIdx list with pointers to stiff equations. Empty
    //!   for now but will have to index to chemical non-equilibrium when added
    void setStiffEqIdx( std::vector< std::size_t >& stiffEqIdx ) const
    {
      stiffEqIdx.resize(0);
    }

    //! Locate the nonstiff equations.
    //! \param[out] nonStiffEqIdx list with pointers to nonstiff equations
    void setNonStiffEqIdx( std::vector< std::size_t >& nonStiffEqIdx ) const
    {
      nonStiffEqIdx.resize(0);
    }

    //! Initialize the compressible flow equations, prepare for time integration
    //! \param[in] geoElem Element geometry array
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] inbox List of elements at which box user ICs are set for
    //!   each IC box
    //! \param[in] elemblkid Element ids associated with mesh block ids where
    //!   user ICs are set
    //! \param[in,out] unk Array of unknowns
    //! \param[in] t Physical time
    //! \param[in] nielem Number of internal elements
    void initialize( const tk::Fields& geoElem,
      const std::vector< std::size_t >& inpoel,
      const tk::UnsMesh::Coords& coord,
      const std::vector< std::unordered_set< std::size_t > >& inbox,
      const std::unordered_map< std::size_t, std::set< std::size_t > >&
        elemblkid,
      tk::Fields& unk,
      tk::real t,
      const std::size_t nielem ) const
    {
      tk::initialize( m_ncomp, m_mat_blk, geoElem, inpoel, coord,
                      Problem::initialize, unk, t, nielem );

      const auto rdof = g_inputdeck.get< tag::rdof >();
      const auto& ic = g_inputdeck.get< tag::ic >();
      const auto& icbox = ic.get< tag::box >();
      const auto& icmbk = ic.get< tag::meshblock >();

      // Set initial conditions inside user-defined IC boxes and mesh blocks
      std::vector< tk::real > s(m_ncomp, 0.0);
      for (std::size_t e=0; e<nielem; ++e) {
        // inside user-defined box
        if (!icbox.empty()) {
          std::size_t bcnt = 0;
          for (const auto& b : icbox) {   // for all boxes
            if (inbox.size() > bcnt && inbox[bcnt].find(e) != inbox[bcnt].end())
            {
              std::vector< tk::real > box
                { b.template get< tag::xmin >(), b.template get< tag::xmax >(),
                  b.template get< tag::ymin >(), b.template get< tag::ymax >(),
                  b.template get< tag::zmin >(), b.template get< tag::zmax >() };
              auto V_ex = (box[1]-box[0]) * (box[3]-box[2]) * (box[5]-box[4]);
              for (std::size_t c=0; c<m_ncomp; ++c) {
                auto mark = c*rdof;
                s[c] = unk(e,mark);
                // set high-order DOFs to zero
                for (std::size_t i=1; i<rdof; ++i)
                  unk(e,mark+i) = 0.0;
              }
              initializeBox<ctr::boxList>( m_mat_blk, V_ex, t, b, s );
              // store box-initialization in solution vector
              for (std::size_t c=0; c<m_ncomp; ++c) {
                auto mark = c*rdof;
                unk(e,mark) = s[c];
              }
            }
            ++bcnt;
          }
        }

        // inside user-specified mesh blocks
        if (!icmbk.empty()) {
          for (const auto& b : icmbk) { // for all blocks
            auto blid = b.get< tag::blockid >();
            auto V_ex = b.get< tag::volume >();
            if (elemblkid.find(blid) != elemblkid.end()) {
              const auto& elset = tk::cref_find(elemblkid, blid);
              if (elset.find(e) != elset.end()) {
                initializeBox<ctr::meshblockList>( m_mat_blk, V_ex, t, b, s );
                // store initialization in solution vector
                for (std::size_t c=0; c<m_ncomp; ++c) {
                  auto mark = c*rdof;
                  unk(e,mark) = s[c];
                }
              }
            }
          }
        }
      }
    }

    //! Compute average plastic deformation on each element
    // //! \param[in] nelem Number of elements
    // //! \param[in] unk Array of unknowns
    // //! \param[in] pri Array of primitives
    //! \param[out] plasticDeformation Frobenius norm of Lp matrix
    void computePlasticDeformation( std::size_t /*nelem*/,
                                    tk::Fields& /*unk*/,
                                    tk::Fields& /*pri*/,
                                    std::vector< tk::real >& plasticDeformation) const
    {
      plasticDeformation.resize(0);
    }

    //! Update the interface cells to first order dofs. No-op.
    // //! \param[in] unk Array of unknowns
    // //! \param[in] nielem Number of internal elements
    // //! \param[in,out] ndofel Array of dofs
    // //! \param[in,out] interface Vector of interface marker
    void updateInterfaceCells( tk::Fields& /*unk*/,
      std::size_t /*nielem*/,
      std::vector< std::size_t >& /*ndofel*/,
      std::vector< std::size_t >& /*interface*/ ) const {}

    //! Update the primitives for this PDE system.
    //! \param[in] unk Array of unknowns
    //! \param[in] geoElem Element geometry array
    //! \param[in,out] prim Array of primitives
    //! \param[in] nielem Number of internal elements
    //! \param[in] ndofel Array of dofs
    //! \details This function computes and stores the dofs for primitive
    //!   quantities, which are required for obtaining reconstructed states used
    //!   in the Riemann solver. See for eg. /PDE/Riemann/AUSMMultiSpecies.hpp,
    //!   where temperature is independently reconstructed.
    void updatePrimitives( const tk::Fields& unk,
                           const tk::Fields& geoElem,
                           tk::Fields& prim,
                           std::size_t nielem,
                           const std::vector< std::size_t >& ndofel ) const
    {
      const auto rdof = g_inputdeck.get< tag::rdof >();
      const auto ndof = g_inputdeck.get< tag::ndof >();
      auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();

      Assert( unk.nunk() == prim.nunk(), "Number of unknowns in solution "
              "vector and primitive vector at recent time step incorrect" );
      Assert( unk.nprop() == rdof*m_ncomp, "Number of components in solution "
              "vector must equal "+ std::to_string(rdof*m_ncomp) );
      Assert( prim.nprop() == rdof*m_nprim, "Number of components in vector of "
              "primitive quantities must equal "+ std::to_string(rdof*m_nprim) );

      auto mass_m = tk::massMatrixDubiner();
      std::vector< tk::real > state(m_ncomp, 0.0);

      for (std::size_t e=0; e<nielem; ++e)
      {
        std::vector< tk::real > R(m_nprim*ndof, 0.0);

        auto ng = tk::NGvol(ndof);

        // arrays for quadrature points
        std::array< std::vector< tk::real >, 3 > coordgp;
        std::vector< tk::real > wgp;

        coordgp[0].resize( ng );
        coordgp[1].resize( ng );
        coordgp[2].resize( ng );
        wgp.resize( ng );

        tk::GaussQuadratureTet( ng, coordgp, wgp );

        // Local degree of freedom
        auto dof_el = ndofel[e];

        auto vole = geoElem(e, 0);

        std::vector< tk::real > B(dof_el);

        // Loop over quadrature points in element e
        for (std::size_t igp=0; igp<ng; ++igp)
        {
          // Compute the basis function
          tk::eval_basis( dof_el, coordgp[0][igp], coordgp[1][igp],
            coordgp[2][igp], B );

          auto w = wgp[igp] * vole;

          tk::eval_state( m_ncomp, rdof, dof_el, e, unk, B, state.data() );

          // Mixture state at quadrature point
          Mixture mixgp(nspec, state, m_mat_blk);

          // Mixture density at quadrature point
          tk::real rhob = mixgp.get_mix_density();

          // velocity vector at quadrature point
          std::array< tk::real, 3 >
            vel{ state[multispecies::momentumIdx(nspec, 0)]/rhob,
                 state[multispecies::momentumIdx(nspec, 1)]/rhob,
                 state[multispecies::momentumIdx(nspec, 2)]/rhob };

          std::vector< tk::real > pri(m_nprim, 0.0);

          // Evaluate mixture temperature at quadrature point
          auto rhoE0 = state[multispecies::energyIdx(nspec, 0)];
          int iconv(0);
          pri[multispecies::temperatureIdx(nspec,0)] =
            mixgp.temperature(rhob, vel[0], vel[1], vel[2], rhoE0, m_mat_blk,
              iconv, prim(e,multispecies::temperatureDofIdx(nspec, 0, rdof, 0)));

          // Check for unphysical state and error out with information
          if (!iconv) {
            std::cout << "Element centroid:  " << geoElem(e,1) << ", "
              << geoElem(e,2) << ", " << geoElem(e,3) << std::endl;
            std::cout << "Mass fractions:    ";
            for (std::size_t k=0; k<nspec; ++k) {
              std::cout << state[multispecies::densityIdx(nspec,k)]/rhob << ", ";
            }
            std::cout << std::endl;
            std::cout << "Mixture density:   " << rhob << std::endl;
            std::cout << "Total energy:      " << rhoE0 << std::endl;
            std::cout << "Velocity:          " <<
              vel[0] << ", " << vel[1] << ", " << vel[2] << std::endl;

            Throw( "Newton iteration for temperature failed to converge "
              " with temperature at final iteration = " +
              std::to_string(pri[multispecies::temperatureIdx(nspec,0)]) );
          }

          // Clip temperature
          pri[multispecies::temperatureIdx(nspec,0)] = constrain_temperature(
            pri[multispecies::temperatureIdx(nspec,0)]);

          for(std::size_t k = 0; k < m_nprim; k++)
          {
            auto mark = k * ndof;
            for(std::size_t idof = 0; idof < dof_el; idof++)
              R[mark+idof] += w * pri[k] * B[idof];
          }
        }

        // Update the DG solution of primitive variables
        for(std::size_t k = 0; k < m_nprim; k++)
        {
          auto mark = k * ndof;
          auto rmark = k * rdof;
          for(std::size_t idof = 0; idof < dof_el; idof++)
          {
            prim(e, rmark+idof) = R[mark+idof] / (mass_m[idof]*vole);
          }
        }
      }
    }

    //! Clean up the state of trace materials for this PDE system. No-op.
    // //! \param[in] t Physical time
    // //! \param[in] geoElem Element geometry array
    // //! \param[in,out] unk Array of unknowns
    // //! \param[in,out] prim Array of primitives
    // //! \param[in] nielem Number of internal elements
    void cleanTraceMaterial( tk::real /*t*/,
                             const tk::Fields& /*geoElem*/,
                             tk::Fields& /*unk*/,
                             tk::Fields& /*prim*/,
                             std::size_t /*nielem*/ ) const {}

    //! Reconstruct second-order solution from first-order
    //! \param[in] geoElem Element geometry array
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] esup Elements-surrounding-nodes connectivity
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in,out] U Solution vector at recent time step
    //! \param[in,out] P Vector of primitives at recent time step
    //! \param[in] pref Indicator for p-adaptive algorithm
    //! \param[in] ndofel Vector of local number of degrees of freedome
    void reconstruct( tk::real,
                      const tk::Fields&,
                      const tk::Fields& geoElem,
                      const inciter::FaceData& fd,
                      const std::map< std::size_t, std::vector< std::size_t > >&
                        esup,
                      const std::vector< std::size_t >& inpoel,
                      const tk::UnsMesh::Coords& coord,
                      tk::Fields& U,
                      tk::Fields& P,
                      const bool pref,
                      const std::vector< std::size_t >& ndofel ) const
    {
      const auto rdof = g_inputdeck.get< tag::rdof >();
      const auto ndof = g_inputdeck.get< tag::ndof >();

      bool is_p0p1(false);
      if (rdof == 4 && ndof == 1)
        is_p0p1 = true;

      const auto nelem = fd.Esuel().size()/4;

      Assert( U.nprop() == rdof*m_ncomp, "Number of components in solution "
              "vector must equal "+ std::to_string(rdof*m_ncomp) );
      Assert( P.nprop() == rdof*m_nprim, "Number of components in primitive "
              "vector must equal "+ std::to_string(rdof*m_nprim) );

      //----- reconstruction of conserved quantities -----
      //--------------------------------------------------

      for (std::size_t e=0; e<nelem; ++e)
      {
        std::vector< std::size_t > vars;
        // check if element is marked as p0p1
        if ( (pref && ndofel[e] == 1) || is_p0p1 ) {
          // 1. specify how many variables need to be reconstructed
          for (std::size_t c=0; c<m_ncomp; ++c) vars.push_back(c);

          // 2. solve 3x3 least-squares system
          // Reconstruct second-order dofs in Taylor space using nodal-stencils
          tk::recoLeastSqExtStencil( rdof, e, esup, inpoel, geoElem, U, vars );

          // 3. transform reconstructed derivatives to Dubiner dofs
          tk::transform_P0P1( rdof, e, inpoel, coord, U, vars );
        }
      }

      //----- reconstruction of primitive quantities -----
      //--------------------------------------------------
      // For multispecies, conserved and primitive quantities are reconstructed
      // separately.

      for (std::size_t e=0; e<nelem; ++e)
      {
        // There are two conditions that requires the reconstruction of the
        // primitive variables:
        //   1. p-adaptive is triggered and P0P1 scheme is applied to specific
        //      elements
        //   2. p-adaptive is not triggered and P0P1 scheme is applied to the
        //      whole computation domain
        if ((pref && ndofel[e] == 1) || (!pref && is_p0p1)) {
          std::vector< std::size_t > vars;
          for (std::size_t c=0; c<m_nprim; ++c) vars.push_back(c);

          // 1.
          // Reconstruct second-order dofs in Taylor space using nodal-stencils
          tk::recoLeastSqExtStencil( rdof, e, esup, inpoel, geoElem, P, vars );

          // 2.
          tk::transform_P0P1(rdof, e, inpoel, coord, P, vars);
        }
      }
    }

    //! Limit second-order solution, and primitive quantities separately
    // //! \param[in] pref Indicator for p-adaptive algorithm
    //! \param[in] geoFace Face geometry array
    //! \param[in] geoElem Element geometry array
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] esup Elements-surrounding-nodes connectivity
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] ndofel Vector of local number of degrees of freedome
    // //! \param[in] gid Local->global node id map
    // //! \param[in] bid Local chare-boundary node ids (value) associated to
    // //!   global node ids (key)
    // //! \param[in] mtInv Inverse of Taylor mass matrix
    //! \param[in,out] U Solution vector at recent time step
    //! \param[in,out] P Vector of primitives at recent time step
    //! \param[in,out] shockmarker Vector of shock-marker values
    void limit( [[maybe_unused]] tk::real,
                const bool /*pref*/,
                const tk::Fields& geoFace,
                const tk::Fields& geoElem,
                const inciter::FaceData& fd,
                const std::map< std::size_t, std::vector< std::size_t > >& esup,
                const std::vector< std::size_t >& inpoel,
                const tk::UnsMesh::Coords& coord,
                const std::vector< std::size_t >& ndofel,
                const std::vector< std::size_t >& /*gid*/,
                const std::unordered_map< std::size_t, std::size_t >& /*bid*/,
                const std::vector< std::vector<tk::real> >& /*mtInv*/,
                tk::Fields& U,
                tk::Fields& P,
                std::vector< std::size_t >& shockmarker ) const
    {
      const auto limiter = g_inputdeck.get< tag::limiter >();
      auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();
      const auto rdof = g_inputdeck.get< tag::rdof >();
      const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();

      // limit vectors of conserved and primitive quantities
      if (limiter == ctr::LimiterType::VERTEXBASEDP1 && rdof == 4)
      {
        VertexBasedMultiSpecies_P1( esup, inpoel, ndofel, fd.Esuel().size()/4,
          m_mat_blk, fd, geoFace, geoElem, coord, flux, solidx, U, P, nspec,
          shockmarker );
      }
      else if (limiter == ctr::LimiterType::VERTEXBASEDP1 && rdof == 10)
      {
        VertexBasedMultiSpecies_P2( esup, inpoel, ndofel, fd.Esuel().size()/4,
          m_mat_blk, fd, geoFace, geoElem, coord, flux, solidx, U, P, nspec,
          shockmarker );
      }
      else if (limiter != ctr::LimiterType::NOLIMITER)
      {
        Throw("Limiter type not configured for multispecies.");
      }
    }

    //! Apply CPL to the conservative variable solution for this PDE system
    //! \param[in] prim Array of primitive variables
    //! \param[in] geoElem Element geometry array
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in,out] unk Array of conservative variables
    //! \param[in] nielem Number of internal elements
    //! \details This function applies CPL to obtain consistent dofs for
    //!   conservative quantities based on the limited primitive quantities.
    //!   See appendix of paper: Pandare et al. (2023). On the Design of Stable,
    //!   Consistent, and Conservative High-Order Methods for Multi-Material
    //!   Hydrodynamics. J Comp Phys, 112313.
    void CPL( const tk::Fields& prim,
      const tk::Fields& geoElem,
      const std::vector< std::size_t >& inpoel,
      const tk::UnsMesh::Coords& coord,
      tk::Fields& unk,
      std::size_t nielem ) const
    {
      [[maybe_unused]] const auto rdof = g_inputdeck.get< tag::rdof >();
      auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();

      Assert( unk.nunk() == prim.nunk(), "Number of unknowns in solution "
              "vector and primitive vector at recent time step incorrect" );
      Assert( unk.nprop() == rdof*m_ncomp, "Number of components in solution "
              "vector must equal "+ std::to_string(rdof*m_ncomp) );
      Assert( prim.nprop() == rdof*m_nprim, "Number of components in vector of "
              "primitive quantities must equal "+ std::to_string(rdof*m_nprim) );

      correctLimConservMultiSpecies(nielem, m_mat_blk, nspec, inpoel,
        coord, geoElem, prim, unk);
    }

    //! Return cell-average deformation gradient tensor. No-op.
    std::array< std::vector< tk::real >, 9 > cellAvgDeformGrad(
      const tk::Fields&,
      std::size_t ) const
    { return {}; }

    //! Reset the high order solution for p-adaptive scheme
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in,out] unk Solution vector at recent time step
    //! \param[in,out] prim Primitive vector at recent time step
    //! \param[in] ndofel Vector of local number of degrees of freedome
    //! \details This function reset the high order coefficient for p-adaptive
    //!   solution polynomials. Unlike compflow class, the high order of fv
    //!   solution will not be reset since p0p1 is the base scheme for
    //!   multi-species p-adaptive DG method.
    void resetAdapSol( const inciter::FaceData& fd,
                       tk::Fields& unk,
                       tk::Fields& prim,
                       const std::vector< std::size_t >& ndofel ) const
    {
      const auto rdof = g_inputdeck.get< tag::rdof >();
      const auto ncomp = unk.nprop() / rdof;
      const auto nprim = prim.nprop() / rdof;

      for(std::size_t e = 0; e < fd.Esuel().size()/4; e++)
      {
        if(ndofel[e] < 10)
        {
          for (std::size_t c=0; c<ncomp; ++c)
          {
            auto mark = c*rdof;
            unk(e, mark+4) = 0.0;
            unk(e, mark+5) = 0.0;
            unk(e, mark+6) = 0.0;
            unk(e, mark+7) = 0.0;
            unk(e, mark+8) = 0.0;
            unk(e, mark+9) = 0.0;
          }
          for (std::size_t c=0; c<nprim; ++c)
          {
            auto mark = c*rdof;
            prim(e, mark+4) = 0.0;
            prim(e, mark+5) = 0.0;
            prim(e, mark+6) = 0.0;
            prim(e, mark+7) = 0.0;
            prim(e, mark+8) = 0.0;
            prim(e, mark+9) = 0.0;
          }
        }
      }
    }

    //! Compute right hand side
    //! \param[in] t Physical time
    //! \param[in] pref Indicator for p-adaptive algorithm
    //! \param[in] geoFace Face geometry array
    //! \param[in] geoElem Element geometry array
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] U Solution vector at recent time step
    //! \param[in] P Primitive vector at recent time step
    //! \param[in] ndofel Vector of local number of degrees of freedom
    //! \param[in] dt Delta time
    //! \param[in,out] R Right-hand side vector computed
    void rhs( tk::real t,
              const bool pref,
              const tk::Fields& geoFace,
              const tk::Fields& geoElem,
              const inciter::FaceData& fd,
              const std::vector< std::size_t >& inpoel,
              const std::vector< std::unordered_set< std::size_t > >&,
              const tk::UnsMesh::Coords& coord,
              const tk::Fields& U,
              const tk::Fields& P,
              const std::vector< std::size_t >& ndofel,
              const tk::real dt,
              tk::Fields& R ) const
    {
      const auto ndof = g_inputdeck.get< tag::ndof >();
      const auto rdof = g_inputdeck.get< tag::rdof >();
      const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();
      auto viscous = g_inputdeck.get< tag::multispecies, tag::viscous >();
      auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();

      const auto nelem = fd.Esuel().size()/4;

      Assert( U.nunk() == P.nunk(), "Number of unknowns in solution "
              "vector and primitive vector at recent time step incorrect" );
      Assert( U.nunk() == R.nunk(), "Number of unknowns in solution "
              "vector and right-hand side at recent time step incorrect" );
      Assert( U.nprop() == rdof*m_ncomp, "Number of components in solution "
              "vector must equal "+ std::to_string(rdof*m_ncomp) );
      Assert( P.nprop() == rdof*m_nprim, "Number of components in primitive "
              "vector must equal "+ std::to_string(rdof*m_nprim) );
      Assert( R.nprop() == ndof*m_ncomp, "Number of components in right-hand "
              "side vector must equal "+ std::to_string(ndof*m_ncomp) );
      Assert( fd.Inpofa().size()/3 == fd.Esuf().size()/2,
              "Mismatch in inpofa size" );

      // set rhs to zero
      R.fill(0.0);

      // empty vector for non-conservative terms. This vector is unused for
      // multi-species flow since, there are no non-conservative terms
      // in the system of PDEs.
      std::vector< std::vector< tk::real > > riemannDeriv;

      std::vector< std::vector< tk::real > > vriem;
      std::vector< std::vector< tk::real > > riemannLoc;

      // configure a no-op lambda for prescribed velocity
      auto velfn = []( ncomp_t, tk::real, tk::real, tk::real, tk::real ){
        return tk::VelFn::result_type(); };

      // p-adaptive DG
      if (!pref) {
        // compute internal surface flux integrals
        tk::surfInt_constP( 1, m_mat_blk, t, ndof, rdof, inpoel, solidx,
          coord, fd, geoFace, geoElem, m_riemann, velfn, U, P,
          dt, R, riemannDeriv );

        // compute boundary surface flux integrals
        for (const auto& b : m_bc)
          tk::bndSurfInt_constP( 1, m_mat_blk, ndof, rdof, std::get<0>(b), fd,
            geoFace, geoElem, inpoel, coord, t, m_riemann, velfn,
            std::get<1>(b), U, P, R, riemannDeriv );

        // compute volume integrals
        tk::volInt_constP( 1, t, m_mat_blk, ndof, rdof, nelem, inpoel, coord,
          geoElem, flux, velfn, Problem::src, U, P, R );
      }
      else {
        // compute internal surface flux integrals
        tk::surfInt( pref, 1, m_mat_blk, t, ndof, rdof, inpoel, solidx,
                     coord, fd, geoFace, geoElem, m_riemann, velfn, U, P, ndofel,
                     dt, R, riemannDeriv );

        // compute boundary surface flux integrals
        for (const auto& b : m_bc)
          tk::bndSurfInt( pref, 1, m_mat_blk, ndof, rdof, std::get<0>(b), fd,
                          geoFace, geoElem, inpoel, coord, t, m_riemann, velfn,
                          std::get<1>(b), U, P, ndofel, R, riemannDeriv );

        // compute volume integrals
        tk::volInt( 1, t, m_mat_blk, ndof, rdof, nelem, inpoel, coord, geoElem,
          flux, velfn, Problem::src, U, P, ndofel, R );
      }

      if (viscous) {
        tk::surfIntViscousMultiSpecies( nspec, m_mat_blk, ndof, rdof, inpoel,
          coord, fd, geoFace, geoElem, U, P, R );
        for (const auto& b : m_bc)
          tk::bndSurfIntViscousMultiSpecies( nspec, m_mat_blk, ndof, rdof,
            std::get<0>(b), fd, geoFace, geoElem, inpoel, coord, t,
            std::get<1>(b), std::get<2>(b), U, P, R );
      }

      // compute external (energy) sources
      //m_physics.physSrc(nspec, t, geoElem, {}, R, {});
    }

    //  Assemble element-local analytic residual Jacobian for point implicit
    //! \param[in] t Physical time
    //! \param[in] geoFace Face geometry array
    //! \param[in] geoElem Element geometry array
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] inpoel Element connectivity
    //! \param[in] coord Mesh coordinates
    //! \param[in] U Conservative solution
    //! \param[in] P Primitive solution
    //! \param[in] ndofel Number of DOFs per element
    //! \return Dense matrix dR_e/du_e for all owned elements
    //! \details First P0/source-free MultiSpecies implementation. Surface
    //!   contributions are assembled analytically using the existing Riemann
    //!   Jacobian. Only symmetry boundary terms are included in this first 
    //!   version.
    std::vector< std::vector< std::vector< tk::real > > >
    point_implicit_jacobian_analytic(
      tk::real t,
      const tk::Fields& geoFace,
      [[maybe_unused]] const tk::Fields& geoElem,
      const inciter::FaceData& fd,
      [[maybe_unused]] const std::vector< std::size_t >& inpoel,
      const tk::UnsMesh::Coords& coord,
      const tk::Fields& U,
      const tk::Fields& P,
      [[maybe_unused]] const std::vector< std::size_t >& ndofel ) const
    {
      const auto rdof = g_inputdeck.get< tag::rdof >();
      const auto ndof = g_inputdeck.get< tag::ndof >();

      const auto nelem = fd.Esuel().size()/4;
      const auto ncomp = static_cast< std::size_t >( m_ncomp );
      const auto nprim = static_cast< std::size_t >( m_nprim );
      const auto nunk = ncomp*ndof;

      const auto& inpofa = fd.Inpofa();

      const auto& cx = coord[0];
      const auto& cy = coord[1];
      const auto& cz = coord[2];

      // Current exceptions for this method
      if (ndof != 1) {
        Throw(
          "MultiSpecies analytic point-implicit Jacobian currently only "
          "implemented for P0 (ndof=1)");
      }

      if (rdof != 1) {
        Throw(
          "MultiSpecies analytic point-implicit Jacobian currently requires "
          "rdof=1");
      }

      // Full RHS - conserved Jacobian
      std::vector< std::vector< std::vector< tk::real > > >
        dRdu( nelem,
          std::vector< std::vector< tk::real > >(
            nunk, std::vector< tk::real >( nunk, 0.0 ) ) );

      // Calculation of analytic jacobian requires a vector of both conserved
      // variables and primitives.
      // This is required a few times (L and R state, BCs), so put in a lambda
      auto load_state =
        [&]( std::size_t e, std::vector< tk::real >& state )
        {
          state.assign( ncomp+nprim, 0.0 );

          for (std::size_t c=0; c<ncomp; ++c) {
            state[c] = U(e,c*rdof);
          }

          for (std::size_t p=0; p<nprim; ++p) {
            state[ncomp+p] = P(e,p*rdof);
          }
        };

      const auto& esuf = fd.Esuf();

      // Keep the quadrature structure consistent with surfInt_constP().
      const auto ng = tk::NGfa( ndof );

      std::array< std::vector< tk::real >, 2 > coordgp;
      std::vector< tk::real > wgp;

      coordgp[0].resize( ng );
      coordgp[1].resize( ng );
      wgp.resize( ng );

      // Gauss weights required for full analytic Jacobian
      tk::GaussQuadratureTri( ng, coordgp, wgp );

      std::array< std::vector< tk::real >, 2 > state;
      state[0].resize( ncomp+nprim );
      state[1].resize( ncomp+nprim );

      // Boundary faces: Only symmetry BC is implemented for now.
      // Other BCs are skipped in this first version.
      const auto& bface = fd.Bface();

      for (const auto& b : m_bc) {
        const auto& bcconfig = std::get< 0 >(b);
        const auto& bstatefn = std::get< 1 >(b);

        // Check if this BC is symmetry by comparing function pointers.
        // Both symmetry and slip-wall are configured with the symmetry function.
        const auto* target = bstatefn.target< decltype(&symmetry) >();
        if (target == nullptr || *target != &symmetry) {
          continue;
        }

        // Loop over side-sets for this BC
        for (const auto sideset : bcconfig) {
          const auto it = bface.find( static_cast< int >( sideset ) );
          if (it == bface.end()) continue;

          // Loop over actual boundary faces in this side-set
          for (const auto f : it->second) {
            // Face coordinates
            std::array< std::array< tk::real, 3 >, 3 > coordfa{{
              {{
                cx[ inpofa[3*f] ], cy[ inpofa[3*f] ], cz[ inpofa[3*f] ]
              }},
              {{
                cx[ inpofa[3*f+1] ], cy[ inpofa[3*f+1] ], cz[ inpofa[3*f+1] ]
              }},
              {{
                cx[ inpofa[3*f+2] ], cy[ inpofa[3*f+2] ], cz[ inpofa[3*f+2] ]
              }}
            }};

            if( esuf[2*f] < 0 ) Throw("Inside boundary element is invalid");
            if( esuf[2*f+1] != -1 ) 
              Throw("Outside boundary element is not -1");

            const auto el = static_cast< std::size_t >( esuf[2*f] );

            if( el >= nelem ) Throw("Boundary element index out of bounds");

            std::array< tk::real, 3 >
              fn{{ geoFace(f,1), geoFace(f,2), geoFace(f,3) }};

            // Accumulate Jacobian contributions
            // Apply chain rule:
            //   dF/dU_internal = dF/dU_left + dF/dU_right * dU_right/dU_left
            for (std::size_t igp=0; igp<ng; ++igp) {
              const auto r = coordgp[0][igp];
              const auto s = coordgp[1][igp];

              // Map reference-triangle point (r,s) to the physical face.
              std::array< tk::real, 3 > gp;

              // This is overkill if we are restricted to P0 since we only have
              // one Gauss point at the centroid, but hopefully it's forward-
              // looking for higher-order solvers when we need the GP's
              // at not just the centroid.
              for (std::size_t i=0; i<3; ++i) {
                gp[i] =
                  (1.0 - r - s) * coordfa[0][i]
                  + r * coordfa[1][i]
                  + s * coordfa[2][i];
              }

              const auto wt = wgp[igp] * geoFace(f,0);

              // Load internal (left) state
              load_state( el, state[0] );

              // Compute both left and right states using BC
              state = bstatefn( ncomp, m_mat_blk, state[0],
                                gp[0], gp[1], gp[2], t, fn );

              auto dUdP = conservedPrimitiveJac(m_mat_blk, {}, state, {});

              // Compute flux Jacobian
              const auto dFdU = m_riemannjac( m_mat_blk, fn, dUdP, state, {} );

              // Compute ghost state Jacobian: dU_ghost/dU_internal
              const auto dUgdUi = symmetryJacobian( ncomp, fn );

              for (std::size_t row=0; row<ncomp; ++row) {
                const auto rmark = row*ndof;
                for (std::size_t col=0; col<ncomp; ++col) {
                  const auto cmark = col*ndof;
                  // Compute total flux derivative: J_left + J_right * dUgdUi
                  auto dflux = dFdU[0][row][col];
                  for (std::size_t k=0; k<ncomp; ++k) {
                    dflux += dFdU[1][row][k] * dUgdUi[k][col];
                  }

                  dRdu[el][rmark][cmark] -= wt * dflux;
                }
              }
            }
          }
        }
      }

      // Interior faces
      for (auto f=fd.Nbfac(); f<esuf.size()/2; ++f) {
        if( esuf[2*f] < 0 || esuf[2*f+1] < 0) {
          Throw("Interior element detected as -1 in "
            "analytic point-implicit Jacobian");
        }

        const auto el = static_cast< std::size_t >( esuf[2*f] );
        const auto er = static_cast< std::size_t >( esuf[2*f+1] );

        if( el >= U.nunk()) Throw("Left element index out of bounds" );
        if( er >= U.nunk()) Throw("Right element index out of bounds" );

        std::array< tk::real, 3 >
          fn{{ geoFace(f,1), geoFace(f,2), geoFace(f,3) }};

        // Need state, containing conserved and primitive variables
        load_state( el, state[0] );
        load_state( er, state[1] );

        auto dUdP = conservedPrimitiveJac(m_mat_blk, {}, state, {});

        const auto dFdU = m_riemannjac( m_mat_blk, fn, dUdP, state, {} );

        for (std::size_t igp=0; igp<ng; ++igp) {
          const auto wt = wgp[igp] * geoFace(f,0);

          for (std::size_t row=0; row<ncomp; ++row) {
            const auto rmark = row*ndof;

            for (std::size_t col=0; col<ncomp; ++col) {
              const auto cmark = col*ndof;

              // Left residual receives -F.
              if (el < nelem) {
                dRdu[el][rmark][cmark] -= wt * dFdU[0][row][col];
              }

              // Right residual receives +F.
              if (er < nelem) {
                dRdu[er][rmark][cmark] += wt * dFdU[1][row][col];
              }
            }
          }
        }
      }

      return dRdu;
    }

    //! Evaluate the adaptive indicator and mark the ndof for each element
    //! \param[in] nunk Number of unknowns
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] unk Array of unknowns
    // //! \param[in] prim Array of primitive quantities
    //! \param[in] indicator p-refinement indicator type
    //! \param[in] ndof Number of degrees of freedom in the solution
    //! \param[in] ndofmax Max number of degrees of freedom for p-refinement
    //! \param[in] tolref Tolerance for p-refinement
    //! \param[in,out] ndofel Vector of local number of degrees of freedome
    void eval_ndof( std::size_t nunk,
                    [[maybe_unused]] const tk::UnsMesh::Coords& coord,
                    [[maybe_unused]] const std::vector< std::size_t >& inpoel,
                    const inciter::FaceData& fd,
                    const tk::Fields& unk,
                    const tk::Fields& /*prim*/,
                    inciter::ctr::PrefIndicatorType indicator,
                    std::size_t ndof,
                    std::size_t ndofmax,
                    tk::real tolref,
                    std::vector< std::size_t >& ndofel ) const
    {
      const auto& esuel = fd.Esuel();

      if(indicator == inciter::ctr::PrefIndicatorType::SPECTRAL_DECAY)
        spectral_decay(1, nunk, esuel, unk, ndof, ndofmax, tolref, ndofel);
      else
        Throw( "No such adaptive indicator type" );
    }

    //! Compute the minimum time step size
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] geoFace Face geometry array
    //! \param[in] geoElem Element geometry array
//    //! \param[in] ndofel Vector of local number of degrees of freedom
    //! \param[in] U Solution vector at recent time step
    //! \param[in] P Vector of primitive quantities at recent time step
    //! \param[in] nielem Number of internal elements
    //! \param[in,out] local_dte Time step size for each element (for local
    //!   time stepping)
    //! \return Minimum time step size
    //! \details The allowable dt is calculated by looking at the maximum
    //!   wave-speed in elements surrounding each face, times the area of that
    //!   face. Once the maximum of this quantity over the mesh is determined,
    //!   the volume of each cell is divided by this quantity. A minimum of this
    //!   ratio is found over the entire mesh, which gives the allowable dt.
    tk::real dt( const std::array< std::vector< tk::real >, 3 >&,
                 const std::vector< std::size_t >&,
                 const inciter::FaceData& fd,
                 const tk::Fields& geoFace,
                 const tk::Fields& geoElem,
                 const std::vector< std::size_t >& /*ndofel*/,
                 const tk::Fields& U,
                 const tk::Fields& P,
                 const std::size_t nielem,
                 std::vector< tk::real >& local_dte ) const
    {
      const auto ndof = g_inputdeck.get< tag::ndof >();
      auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();

      auto mindt = timeStepSizeMultiSpecies( m_mat_blk, fd.Esuf(), geoFace,
        geoElem, nielem, nspec, U, P, local_dte);

      //mindt = std::min(mindt, m_physics.dtRestriction(geoElem, nielem, {}));

      tk::real dgp = 0.0;
      if (ndof == 4)
      {
        dgp = 1.0;
      }
      else if (ndof == 10)
      {
        dgp = 2.0;
      }

      // Scale smallest dt with CFL coefficient and the CFL is scaled by (2*p+1)
      // where p is the order of the DG polynomial by linear stability theory.
      mindt /= (2.0*dgp + 1.0);
      for (std::size_t e=0; e<nielem; ++e)
        local_dte[e] /= (2.0*dgp + 1.0);
      return mindt;
    }

    //! Balances elastic energy after plastic update. Not implemented here.
    // //! \param[in] e Element number
    // //! \param[in] x_star Stiff variables before implicit update
    // //! \param[in] x Stiff variables after implicit update
    // //! \param[in] U Field of conserved variables
    void balance_plastic_energy( std::size_t /*e*/,
                                 std::vector< tk::real > /*x_star*/,
                                 std::vector< tk::real > /*x*/,
                                 tk::Fields& /*U*/ ) const {}

    //! Compute stiff terms for a single element. No-op until chem sources added
    // //! \param[in] e Element number
    // //! \param[in] geoElem Element geometry array
    // //! \param[in] U Solution vector at recent time step
    // //! \param[in] ndofel Vector of local number of degrees of freedom
    // //! \param[in,out] R Right-hand side vector computed
    void stiff_rhs( std::size_t /*e*/,
                    const tk::Fields& /*geoElem*/,
                    const tk::Fields& /*U*/,
                    const std::vector< std::size_t >& /*ndofel*/,
                    tk::Fields& /*R*/ ) const {}

    //! Extract the velocity field at cell nodes. Currently unused.
    // //! \param[in] U Solution vector at recent time step
    // //! \param[in] N Element node indices
    //! \return Array of the four values of the velocity field
    std::array< std::array< tk::real, 4 >, 3 >
    velocity( const tk::Fields& /*U*/,
              const std::array< std::vector< tk::real >, 3 >&,
              const std::array< std::size_t, 4 >& /*N*/ ) const
    {
      std::array< std::array< tk::real, 4 >, 3 > v;
      return v;
    }

    //! Return a map that associates user-specified strings to functions
    //! \return Map that associates user-specified strings to functions that
    //!   compute relevant quantities to be output to file
    std::map< std::string, tk::GetVarFn > OutVarFn() const
    { return MultiSpeciesOutVarFn(); }

    //! Return analytic field names to be output to file
    //! \return Vector of strings labelling analytic fields output in file
    std::vector< std::string > analyticFieldNames() const {
      auto nspec = g_inputdeck.get< eq, tag::nspec >();

      return MultiSpeciesFieldNames(nspec);
    }

    //! Return time history field names to be output to file
    //! \return Vector of strings labelling time history fields output in file
    std::vector< std::string > histNames() const {
      return MultiSpeciesHistNames();
    }

    //! Return surface field names to be output to file
    //! \return Vector of strings labelling surface fields output in file
    std::vector< std::string > surfNames() const
    {
      return MultiSpeciesSurfNames();
    }

    //! Return surface field output going to file
    std::vector< std::vector< tk::real > >
    surfOutput( const inciter::FaceData& fd,
      const tk::Fields& geoFace,
      const std::vector< std::size_t >& inpoel,
      const tk::UnsMesh::Coords& coord,
      const tk::Fields& U,
      const tk::Fields& P ) const
    {
      const auto rdof = g_inputdeck.get< tag::rdof >();
      const auto nspec = g_inputdeck.get< eq, tag::nspec >();

      return MultiSpeciesSurfOutput( nspec, rdof, fd, geoFace, inpoel, coord,
        U, P );
    }

    //! Return time history field output evaluated at time history points
    //! \param[in] h History point data
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] U Array of unknowns
    //! \param[in] P Array of primitive quantities
    //! \return Vector of time history output of bulk flow quantities (density,
    //!   velocity, total energy, and pressure) evaluated at time history points
    std::vector< std::vector< tk::real > >
    histOutput( const std::vector< HistData >& h,
                const std::vector< std::size_t >& inpoel,
                const tk::UnsMesh::Coords& coord,
                const tk::Fields& U,
                const tk::Fields& P ) const
    {
      const auto rdof = g_inputdeck.get< tag::rdof >();
      auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();

      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      std::vector< std::vector< tk::real > > Up(h.size());
      std::vector< tk::real > B(rdof), uhp(m_ncomp, 0.0), php(m_nprim, 0.0);

      std::size_t j = 0;
      for (const auto& p : h) {
        auto e = p.get< tag::elem >();
        auto chp = p.get< tag::coord >();

        // Evaluate inverse Jacobian
        std::array< std::array< tk::real, 3>, 4 > cp{{
          {{ x[inpoel[4*e  ]], y[inpoel[4*e  ]], z[inpoel[4*e  ]] }},
          {{ x[inpoel[4*e+1]], y[inpoel[4*e+1]], z[inpoel[4*e+1]] }},
          {{ x[inpoel[4*e+2]], y[inpoel[4*e+2]], z[inpoel[4*e+2]] }},
          {{ x[inpoel[4*e+3]], y[inpoel[4*e+3]], z[inpoel[4*e+3]] }} }};
        auto J = tk::inverseJacobian( cp[0], cp[1], cp[2], cp[3] );

        // evaluate solution at history-point
        std::array< tk::real, 3 > dc{{chp[0]-cp[0][0], chp[1]-cp[0][1],
          chp[2]-cp[0][2]}};
        tk::eval_basis(rdof, tk::dot(J[0],dc), tk::dot(J[1],dc),
          tk::dot(J[2],dc), B);
        eval_state(m_ncomp, rdof, rdof, e, U, B, uhp.data());
        eval_state(m_nprim, rdof, rdof, e, P, B, php.data());

        // Mixture calculations, initialized
        Mixture mix(nspec, uhp, m_mat_blk);

        // store solution in history output vector
        Up[j].resize(6+nspec, 0.0);
        Up[j][0] = mix.get_mix_density();
        Up[j][1] = uhp[multispecies::momentumIdx(nspec,0)]/Up[j][0];
        Up[j][2] = uhp[multispecies::momentumIdx(nspec,1)]/Up[j][0];
        Up[j][3] = uhp[multispecies::momentumIdx(nspec,2)]/Up[j][0];
        Up[j][4] = uhp[multispecies::energyIdx(nspec,0)];
        Up[j][5] = mix.pressure( Up[j][0],
          php[multispecies::temperatureIdx(nspec,0)] );
        for (std::size_t k=0; k<nspec; ++k) {
          Up[j][6+k] = uhp[multispecies::densityIdx(nspec,k)]/Up[j][0];
        }
        ++j;
      }

      return Up;
    }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    std::vector< std::string > names() const
    {
      auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();
      return MultiSpeciesDiagNames(nspec);
    }

    //! Return analytic solution (if defined by Problem) at xi, yi, zi, t
    //! \param[in] xi X-coordinate at which to evaluate the analytic solution
    //! \param[in] yi Y-coordinate at which to evaluate the analytic solution
    //! \param[in] zi Z-coordinate at which to evaluate the analytic solution
    //! \param[in] t Physical time at which to evaluate the analytic solution
    //! \return Vector of analytic solution at given location and time
    std::vector< tk::real >
    analyticSolution( tk::real xi, tk::real yi, tk::real zi, tk::real t ) const
    { return Problem::analyticSolution( m_ncomp, m_mat_blk, xi, yi, zi, t ); }

    //! Return analytic solution for conserved variables
    //! \param[in] xi X-coordinate at which to evaluate the analytic solution
    //! \param[in] yi Y-coordinate at which to evaluate the analytic solution
    //! \param[in] zi Z-coordinate at which to evaluate the analytic solution
    //! \param[in] t Physical time at which to evaluate the analytic solution
    //! \return Vector of analytic solution at given location and time
    std::vector< tk::real >
    solution( tk::real xi, tk::real yi, tk::real zi, tk::real t ) const
    { return Problem::initialize( m_ncomp, m_mat_blk, xi, yi, zi, t ); }

    //! Return cell-averaged specific total energy for an element
    //! \param[in] e Element id for which total energy is required
    //! \param[in] unk Vector of conserved quantities
    //! \return Cell-averaged specific total energy for given element
    tk::real sp_totalenergy(std::size_t e, const tk::Fields& unk) const
    {
      const auto rdof = g_inputdeck.get< tag::rdof >();
      auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();

      return unk(e, multispecies::energyDofIdx(nspec,0,rdof,0));
    }

  private:
    //! Physics policy
    const Physics m_physics;
    //! Number of components in this PDE system
    const ncomp_t m_ncomp;
    //! Number of primitive quantities stored in this PDE system
    const ncomp_t m_nprim;
    //! Riemann solver
    tk::RiemannFluxFn m_riemann;
    //! Riemann solver Jacboian
    tk::RiemannFluxJacFn m_riemannjac;
    //! BC configuration
    BCStateFn m_bc;
    //! EOS material block
    std::vector< EOS > m_mat_blk;

    //! Evaluate conservative part of physical flux function for this PDE system
    //! \param[in] ncomp Number of scalar components in this PDE system
    //! \param[in] ugp Numerical solution at the Gauss point at which to
    //!   evaluate the flux
    //! \return Flux vectors for all components in this PDE system
    //! \note The function signature must follow tk::FluxFn
    static tk::FluxFn::result_type
    flux( [[maybe_unused]] ncomp_t ncomp,
          const std::vector< EOS >& mat_blk,
          const std::vector< tk::real >& ugp,
          const std::vector< std::array< tk::real, 3 > >& )
    {
      auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();

      std::vector< std::array< tk::real, 3 > > fl( ugp.size() );

      Mixture mix(nspec, ugp, mat_blk);
      auto rhob = mix.get_mix_density();

      std::array< tk::real, 3 > u{{
        ugp[multispecies::momentumIdx(nspec,0)] / rhob,
        ugp[multispecies::momentumIdx(nspec,1)] / rhob,
        ugp[multispecies::momentumIdx(nspec,2)] / rhob }};
      auto p = mix.pressure(rhob,
        ugp[ncomp+multispecies::temperatureIdx(nspec,0)]);

      // density flux
      for (std::size_t k=0; k<nspec; ++k) {
        auto idx = multispecies::densityIdx(nspec, k);
        for (std::size_t j=0; j<3; ++j) {
          fl[idx][j] = ugp[idx] * u[j];
        }
      }

      // momentum flux
      for (std::size_t i=0; i<3; ++i) {
        auto idx = multispecies::momentumIdx(nspec,i);
        for (std::size_t j=0; j<3; ++j) {
          fl[idx][j] = ugp[idx] * u[j];
          if (i == j) fl[idx][j] += p;
        }
      }

      // energy flux
      auto idx = multispecies::energyIdx(nspec,0);
      for (std::size_t j=0; j<3; ++j) {
        fl[idx][j] = u[j] * (ugp[idx] + p);
      }

      return fl;
    }

  //! Calculates the Jacobian of the converved variables with respect to the
  //! primitive variables
  //! \param[in] u Left and right unknown/state vector
  //! \return Derivatives of conserved variables with respect to primitive
  //! variables for the left and right states
  //! \note The function signature must follow tk::RiemannFluxJacFn
  static tk::RiemannFluxJacFn::result_type
  conservedPrimitiveJac(
    const std::vector< EOS >& mat_blk,
    const std::array< tk::real, 3 >&,
    const std::array< std::vector< tk::real >, 2 >& u,
    const std::vector< std::array< tk::real, 3 > >& = {} )
  {
    auto ncomp = u[0].size()-1;
    auto nspec = g_inputdeck.get< tag::multispecies, tag::nspec >();
    std::array dUdP{ std::vector(ncomp, std::vector< tk::real >(ncomp)),
                     std::vector(ncomp, std::vector< tk::real >(ncomp)) };
    tk::real rhol(0.0), rhor(0.0), Tl(0.0), Tr(0.0);
    std::size_t Tid = ncomp - 1;
    const auto uid = multispecies::momentumIdx( nspec, 0 );
    const auto vid = multispecies::momentumIdx( nspec, 1 );
    const auto wid = multispecies::momentumIdx( nspec, 2 );
    const auto eidx = multispecies::energyIdx( nspec, 0 );

    // Initialize mixtures
    Mixture mixl(nspec, u[0], mat_blk);
    Mixture mixr(nspec, u[1], mat_blk);

    Tl = u[0][ncomp+multispecies::temperatureIdx(nspec, 0)];
    Tr = u[1][ncomp+multispecies::temperatureIdx(nspec, 0)];
    rhol = mixl.get_mix_density();
    rhor = mixr.get_mix_density();
    if (!(rhol > 0.0) || !(rhor > 0.0)) {
      Throw(
        "Non-positive mixture density in MultiSpecies "
        "conserved-to-primitive Jacobian");
    }

    // Velocities
    auto ul = u[0][multispecies::momentumIdx(nspec, 0)]/rhol;
    auto vl = u[0][multispecies::momentumIdx(nspec, 1)]/rhol;
    auto wl = u[0][multispecies::momentumIdx(nspec, 2)]/rhol;
    auto ur = u[1][multispecies::momentumIdx(nspec, 0)]/rhor;
    auto vr = u[1][multispecies::momentumIdx(nspec, 1)]/rhor;
    auto wr = u[1][multispecies::momentumIdx(nspec, 2)]/rhor;

    // Partials of species density w.r.t. species density
    for (std::size_t k=0; k<nspec; ++k)
    {
      dUdP[0][k][k] = 1.0;
      dUdP[1][k][k] = 1.0;
    }

    // Partials of momentum...
    // ... w.r.t. species density
    for (std::size_t k=0; k<nspec; ++k)
    {
      dUdP[0][uid][k] = ul;
      dUdP[0][vid][k] = vl;
      dUdP[0][wid][k] = wl;
      dUdP[1][uid][k] = ur;
      dUdP[1][vid][k] = vr;
      dUdP[1][wid][k] = wr;
    }

    // ... w.r.t. velocity
    dUdP[0][uid][uid] = rhol;
    dUdP[1][uid][uid] = rhor;
    dUdP[0][vid][vid] = rhol;
    dUdP[1][vid][vid] = rhor;
    dUdP[0][wid][wid] = rhol;
    dUdP[1][wid][wid] = rhor;

    // Partials of energy...
    // ... w.r.t. species density
    for (std::size_t k=0; k<nspec; ++k)
    {
      // This assumes de_s/drho_s = 0
      dUdP[0][eidx][k] = mat_blk[k].compute< EOS::internalenergy >(Tl)
                      + 0.5 * (ul*ul + vl*vl + wl*wl);
      dUdP[1][eidx][k] = mat_blk[k].compute< EOS::internalenergy >(Tr)
                      + 0.5 * (ur*ur + vr*vr + wr*wr);
    }
    // ... w.r.t. velocity
    dUdP[0][eidx][uid] = rhol * ul;
    dUdP[0][eidx][vid] = rhol * vl;
    dUdP[0][eidx][wid] = rhol * wl;
    dUdP[1][eidx][uid] = rhor * ur;
    dUdP[1][eidx][vid] = rhor * vr;
    dUdP[1][eidx][wid] = rhor * wr;

    // ... w.r.t. temperature
    dUdP[0][eidx][Tid] = rhol * mixl.mix_Cv(Tl, mat_blk);
    dUdP[1][eidx][Tid] = rhor * mixr.mix_Cv(Tr, mat_blk);

    return dUdP;
  }

    //! \brief Boundary state function providing the left and right state of a
    //!   face at Dirichlet boundaries
    //! \param[in] ncomp Number of scalar components in this PDE system
    //! \param[in] mat_blk EOS material block
    //! \param[in] ul Left (domain-internal) state
    //! \param[in] x X-coordinate at which to compute the states
    //! \param[in] y Y-coordinate at which to compute the states
    //! \param[in] z Z-coordinate at which to compute the states
    //! \param[in] t Physical time
    //! \return Left and right states for all scalar components in this PDE
    //!   system
    //! \note The function signature must follow tk::StateFn.
    static tk::StateFn::result_type
    dirichlet( ncomp_t ncomp,
               const std::vector< EOS >& mat_blk,
               const std::vector< tk::real >& ul, tk::real x, tk::real y,
               tk::real z, tk::real t, const std::array< tk::real, 3 >& )
    {
      return {{ ul, Problem::initialize( ncomp, mat_blk, x, y, z, t ) }};
    }

    // Other boundary condition types that do not depend on "Problem" should be
    // added in BCFunctions.hpp
};

} // dg::

} // inciter::

#endif // DGMultiSpecies_h
