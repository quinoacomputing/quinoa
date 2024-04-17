// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/DGMultiMat.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Compressible multi-material flow using discontinuous Galerkin
    finite elements
  \details   This file implements calls to the physics operators governing
    compressible multi-material flow (with velocity equilibrium) using
    discontinuous Galerkin discretizations.
*/
// *****************************************************************************
#ifndef DGMultiMat_h
#define DGMultiMat_h

#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <map>
#include <array>

#include "Macro.hpp"
#include "Exception.hpp"
#include "Vector.hpp"
#include "ContainerUtil.hpp"
#include "UnsMesh.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "Integrate/Basis.hpp"
#include "Integrate/Quadrature.hpp"
#include "Integrate/Initialize.hpp"
#include "Integrate/Mass.hpp"
#include "Integrate/Surface.hpp"
#include "Integrate/Boundary.hpp"
#include "Integrate/Volume.hpp"
#include "Integrate/MultiMatTerms.hpp"
#include "Integrate/Source.hpp"
#include "Integrate/SolidTerms.hpp"
#include "RiemannChoice.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "Reconstruction.hpp"
#include "Limiter.hpp"
#include "Problem/FieldOutput.hpp"
#include "Problem/BoxInitialization.hpp"
#include "PrefIndicator.hpp"
#include "MultiMat/BCFunctions.hpp"
#include "MultiMat/MiscMultiMatFns.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

namespace dg {

//! \brief MultiMat used polymorphically with tk::DGPDE
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Physics - physics configuration, see PDE/MultiMat/Physics.h
//!   - Problem - problem configuration, see PDE/MultiMat/Problem.h
//! \note The default physics is Euler, set in inciter::deck::check_multimat()
template< class Physics, class Problem >
class MultiMat {

  private:
    using eq = tag::multimat;

  public:
    //! Constructor
    explicit MultiMat() :
      m_ncomp( g_inputdeck.get< tag::ncomp >() ),
      m_nprim(nprim()),
      m_riemann( multimatRiemannSolver(
        g_inputdeck.get< tag::flux >() ) )
    {
      // associate boundary condition configurations with state functions
      brigand::for_each< ctr::bclist::Keys >( ConfigBC( m_bc,
        { dirichlet
        , symmetry
        , invalidBC         // Inlet BC not implemented
        , invalidBC         // Outlet BC not implemented
        , farfieldOutlet
        , extrapolate } ) );

      // EoS initialization
      initializeMaterialEoS( m_mat_blk );
    }

    //! Find the number of primitive quantities required for this PDE system
    //! \return The number of primitive quantities required to be stored for
    //!   this PDE system
    std::size_t nprim() const
    {
      auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
      const auto& solidx = inciter::g_inputdeck.get<
        tag::matidxmap, tag::solidx >();

      // individual material pressures and three velocity components
      std::size_t np(nmat+3);

      for (std::size_t k=0; k<nmat; ++k) {
        if (solidx[k] > 0) {
          // individual material Cauchy stress tensor components
          np += 6;
        }
      }

      return np;
    }

    //! Find the number of materials set up for this PDE system
    //! \return The number of materials set up for this PDE system
    std::size_t nmat() const
    {
      auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
      return nmat;
    }

    //! Assign number of DOFs per equation in the PDE system
    //! \param[in,out] numEqDof Array storing number of Dofs for each PDE
    //!   equation
    void numEquationDofs(std::vector< std::size_t >& numEqDof) const
    {
      // all equation-dofs initialized to ndofs first
      for (std::size_t i=0; i<m_ncomp; ++i) {
        numEqDof.push_back(g_inputdeck.get< tag::ndof >());
      }

      // volume fractions are P0Pm (ndof = 1) for multi-material simulations
      auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
      if(nmat > 1)
        for (std::size_t k=0; k<nmat; ++k)
          numEqDof[volfracIdx(nmat, k)] = 1;
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

    //! Initalize the compressible flow equations, prepare for time integration
    //! \param[in] L Block diagonal mass matrix
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] inbox List of elements at which box user ICs are set for
    //!   each IC box
    //! \param[in] elemblkid Element ids associated with mesh block ids where
    //!   user ICs are set
    //! \param[in,out] unk Array of unknowns
    //! \param[in] t Physical time
    //! \param[in] nielem Number of internal elements
    void initialize( const tk::Fields& L,
      const std::vector< std::size_t >& inpoel,
      const tk::UnsMesh::Coords& coord,
      const std::vector< std::unordered_set< std::size_t > >& inbox,
      const std::unordered_map< std::size_t, std::set< std::size_t > >&
        elemblkid,
      tk::Fields& unk,
      tk::real t,
      const std::size_t nielem ) const
    {
      tk::initialize( m_ncomp, m_mat_blk, L, inpoel, coord,
                      Problem::initialize, unk, t, nielem );

      const auto rdof = g_inputdeck.get< tag::rdof >();
      const auto& ic = g_inputdeck.get< tag::ic >();
      const auto& icbox = ic.get< tag::box >();
      const auto& icmbk = ic.get< tag::meshblock >();

      const auto& bgpre = ic.get< tag::pressure >();
      const auto& bgtemp = ic.get< tag::temperature >();

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
              initializeBox<ctr::boxList>( m_mat_blk, V_ex, t, b, bgpre,
                bgtemp, s );
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
                initializeBox<ctr::meshblockList>( m_mat_blk, V_ex, t, b,
                  bgpre, bgtemp, s );
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

    //! Compute the left hand side block-diagonal mass matrix
    //! \param[in] geoElem Element geometry array
    //! \param[in,out] l Block diagonal mass matrix
    void lhs( const tk::Fields& geoElem, tk::Fields& l ) const {
      const auto ndof = g_inputdeck.get< tag::ndof >();
      // Unlike Compflow and Transport, there is a weak reconstruction about
      // conservative variable after limiting function which will require the
      // size of left hand side vector to be rdof
      tk::mass( m_ncomp, ndof, geoElem, l );
    }

    //! Update the interface cells to first order dofs
    //! \param[in] unk Array of unknowns
    //! \param[in] nielem Number of internal elements
//    //! \param[in,out] ndofel Array of dofs
    //! \details This function resets the high-order terms in interface cells.
    void updateInterfaceCells( tk::Fields& unk,
      std::size_t nielem,
      std::vector< std::size_t >& /*ndofel*/ ) const
    {
      auto intsharp =
        g_inputdeck.get< tag::multimat, tag::intsharp >();
      // If this cell is not material interface, return this function
      if(not intsharp)  return;

      auto rdof = g_inputdeck.get< tag::rdof >();
      auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
      const auto& solidx = g_inputdeck.get<
        tag::matidxmap, tag::solidx >();

      for (std::size_t e=0; e<nielem; ++e) {
        std::vector< std::size_t > matInt(nmat, 0);
        std::vector< tk::real > alAvg(nmat, 0.0);
        for (std::size_t k=0; k<nmat; ++k)
          alAvg[k] = unk(e, volfracDofIdx(nmat,k,rdof,0));
        auto intInd = interfaceIndicator(nmat, alAvg, matInt);

        // interface cells cannot be high-order
        if (intInd) {
          //ndofel[e] = 1;
          for (std::size_t k=0; k<nmat; ++k) {
            if (matInt[k]) {
              for (std::size_t i=1; i<rdof; ++i) {
                unk(e, densityDofIdx(nmat,k,rdof,i)) = 0.0;
                unk(e, energyDofIdx(nmat,k,rdof,i)) = 0.0;
              }
              if (solidx[k] > 0) {
                for (std::size_t i=0; i<3; ++i)
                  for (std::size_t j=0; j<3; ++j)
                    for (std::size_t idof=1; idof<rdof; ++idof) {
                      unk(e, deformDofIdx(nmat,solidx[k],i,j,rdof,idof)) = 0.0;
                    }
              }
            }
          }
          for (std::size_t idir=0; idir<3; ++idir) {
            for (std::size_t i=1; i<rdof; ++i) {
              unk(e, momentumDofIdx(nmat,idir,rdof,i)) = 0.0;
            }
          }
        }
      }
    }

    //! Update the primitives for this PDE system
    //! \param[in] unk Array of unknowns
    //! \param[in] L The left hand side block-diagonal mass matrix
    //! \param[in] geoElem Element geometry array
    //! \param[in,out] prim Array of primitives
    //! \param[in] nielem Number of internal elements
    //! \details This function computes and stores the dofs for primitive
    //!   quantities, which are required for obtaining reconstructed states used
    //!   in the Riemann solver. See /PDE/Riemann/AUSM.hpp, where the
    //!   normal velocity for advection is calculated from independently
    //!   reconstructed velocities.
    void updatePrimitives( const tk::Fields& unk,
                           const tk::Fields& L,
                           const tk::Fields& geoElem,
                           tk::Fields& prim,
                           std::size_t nielem ) const
    {
      const auto rdof = g_inputdeck.get< tag::rdof >();
      const auto ndof = g_inputdeck.get< tag::ndof >();
      auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
      const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();

      Assert( unk.nunk() == prim.nunk(), "Number of unknowns in solution "
              "vector and primitive vector at recent time step incorrect" );
      Assert( unk.nprop() == rdof*m_ncomp, "Number of components in solution "
              "vector must equal "+ std::to_string(rdof*m_ncomp) );
      Assert( prim.nprop() == rdof*m_nprim, "Number of components in vector of "
              "primitive quantities must equal "+ std::to_string(rdof*m_nprim) );

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

        // Loop over quadrature points in element e
        for (std::size_t igp=0; igp<ng; ++igp)
        {
          // Compute the basis function
          auto B =
            tk::eval_basis( ndof, coordgp[0][igp], coordgp[1][igp], coordgp[2][igp] );

          auto w = wgp[igp] * geoElem(e, 0);

          auto state = tk::eval_state( m_ncomp, rdof, ndof, e, unk, B );

          // bulk density at quadrature point
          tk::real rhob(0.0);
          for (std::size_t k=0; k<nmat; ++k)
            rhob += state[densityIdx(nmat, k)];

          // velocity vector at quadrature point
          std::array< tk::real, 3 >
            vel{ state[momentumIdx(nmat, 0)]/rhob,
                 state[momentumIdx(nmat, 1)]/rhob,
                 state[momentumIdx(nmat, 2)]/rhob };

          std::vector< tk::real > pri(m_nprim, 0.0);

          // Evaluate material pressure at quadrature point
          for(std::size_t imat = 0; imat < nmat; imat++)
          {
            auto alphamat = state[volfracIdx(nmat, imat)];
            auto arhomat = state[densityIdx(nmat, imat)];
            auto arhoemat = state[energyIdx(nmat, imat)];
            auto gmat = getDeformGrad(nmat, imat, state);
            pri[pressureIdx(nmat,imat)] = m_mat_blk[imat].compute<
              EOS::pressure >( arhomat, vel[0], vel[1], vel[2], arhoemat,
              alphamat, imat, gmat );

            pri[pressureIdx(nmat,imat)] = constrain_pressure( m_mat_blk,
              pri[pressureIdx(nmat,imat)], arhomat, alphamat, imat);

            if (solidx[imat] > 0) {
              auto asigmat = m_mat_blk[imat].computeTensor< EOS::CauchyStress >(
              arhomat, vel[0], vel[1], vel[2], arhoemat,
              alphamat, imat, gmat );

              pri[stressIdx(nmat,solidx[imat],0)] = asigmat[0][0];
              pri[stressIdx(nmat,solidx[imat],1)] = asigmat[1][1];
              pri[stressIdx(nmat,solidx[imat],2)] = asigmat[2][2];
              pri[stressIdx(nmat,solidx[imat],3)] = asigmat[0][1];
              pri[stressIdx(nmat,solidx[imat],4)] = asigmat[0][2];
              pri[stressIdx(nmat,solidx[imat],5)] = asigmat[1][2];
            }
          }

          // Evaluate bulk velocity at quadrature point
          for (std::size_t idir=0; idir<3; ++idir) {
            pri[velocityIdx(nmat,idir)] = vel[idir];
          }

          for(std::size_t k = 0; k < m_nprim; k++)
          {
            auto mark = k * ndof;
            for(std::size_t idof = 0; idof < ndof; idof++)
              R[mark+idof] += w * pri[k] * B[idof];
          }
        }

        // Update the DG solution of primitive variables
        for(std::size_t k = 0; k < m_nprim; k++)
        {
          auto mark = k * ndof;
          auto rmark = k * rdof;
          for(std::size_t idof = 0; idof < ndof; idof++)
          {
            prim(e, rmark+idof) = R[mark+idof] / L(e, mark+idof);
            if(fabs(prim(e, rmark+idof)) < 1e-16)
              prim(e, rmark+idof) = 0;
          }
        }
      }
    }

    //! Clean up the state of trace materials for this PDE system
    //! \param[in] t Physical time
    //! \param[in] geoElem Element geometry array
    //! \param[in,out] unk Array of unknowns
    //! \param[in,out] prim Array of primitives
    //! \param[in] nielem Number of internal elements
    //! \details This function cleans up the state of materials present in trace
    //!   quantities in each cell. Specifically, the state of materials with
    //!   very low volume-fractions in a cell is replaced by the state of the
    //!   material which is present in the largest quantity in that cell. This
    //!   becomes necessary when shocks pass through cells which contain a very
    //!   small amount of material. The state of that tiny material might
    //!   become unphysical and cause solution to diverge; thus requiring such
    //!   a "reset".
    void cleanTraceMaterial( tk::real t,
                             const tk::Fields& geoElem,
                             tk::Fields& unk,
                             tk::Fields& prim,
                             std::size_t nielem ) const
    {
      [[maybe_unused]] const auto rdof = g_inputdeck.get< tag::rdof >();
      auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();

      Assert( unk.nunk() == prim.nunk(), "Number of unknowns in solution "
              "vector and primitive vector at recent time step incorrect" );
      Assert( unk.nprop() == rdof*m_ncomp, "Number of components in solution "
              "vector must equal "+ std::to_string(rdof*m_ncomp) );
      Assert( prim.nprop() == rdof*m_nprim, "Number of components in vector of "
              "primitive quantities must equal "+ std::to_string(rdof*m_nprim) );

      auto neg_density = cleanTraceMultiMat(t, nielem, m_mat_blk, geoElem, nmat,
        unk, prim);

      if (neg_density) Throw("Negative partial density.");
    }

    //! Reconstruct second-order solution from first-order
    //! \param[in] geoElem Element geometry array
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] esup Elements-surrounding-nodes connectivity
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in,out] U Solution vector at recent time step
    //! \param[in,out] P Vector of primitives at recent time step
    void reconstruct( tk::real,
                      const tk::Fields&,
                      const tk::Fields& geoElem,
                      const inciter::FaceData& fd,
                      const std::map< std::size_t, std::vector< std::size_t > >&
                        esup,
                      const std::vector< std::size_t >& inpoel,
                      const tk::UnsMesh::Coords& coord,
                      tk::Fields& U,
                      tk::Fields& P ) const
    {
      const auto rdof = g_inputdeck.get< tag::rdof >();
      const auto ndof = g_inputdeck.get< tag::ndof >();

      bool is_p0p1(false);
      if (rdof == 4 && ndof == 1)
        is_p0p1 = true;

      const auto nelem = fd.Esuel().size()/4;
      auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();

      Assert( U.nprop() == rdof*m_ncomp, "Number of components in solution "
              "vector must equal "+ std::to_string(rdof*m_ncomp) );

      //----- reconstruction of conserved quantities -----
      //--------------------------------------------------
      // specify how many variables need to be reconstructed
      std::vector< std::size_t > vars;
      for (std::size_t c=0; c<m_ncomp; ++c) vars.push_back(c);
      // If DG is applied, reconstruct only volume fractions
      if (!is_p0p1 && ndof > 1)
      {
        vars.clear();
        for (std::size_t k=0; k<nmat; ++k) vars.push_back(volfracIdx(nmat, k));
      }

      // 1. solve 3x3 least-squares system
      for (std::size_t e=0; e<nelem; ++e)
      {
        // Reconstruct second-order dofs of volume-fractions in Taylor space
        // using nodal-stencils, for a good interface-normal estimate
        tk::recoLeastSqExtStencil( rdof, e, esup, inpoel, geoElem, U, vars );
      }

      // 2. transform reconstructed derivatives to Dubiner dofs
      tk::transform_P0P1(rdof, nelem, inpoel, coord, U, vars);

      //----- reconstruction of primitive quantities -----
      //--------------------------------------------------
      // For multimat, conserved and primitive quantities are reconstructed
      // separately.
      if (is_p0p1) {
        vars.clear();
        for (std::size_t c=0; c<m_nprim; ++c) vars.push_back(c);

        // 1.
        for (std::size_t e=0; e<nelem; ++e)
        {
          // Reconstruct second-order dofs of volume-fractions in Taylor space
          // using nodal-stencils, for a good interface-normal estimate
          tk::recoLeastSqExtStencil( rdof, e, esup, inpoel, geoElem, P, vars );
        }

        // 2.
        tk::transform_P0P1(rdof, nelem, inpoel, coord, P, vars);
      }
    }

    //! Limit second-order solution, and primitive quantities separately
    //! \param[in] t Physical time
    //! \param[in] geoFace Face geometry array
    //! \param[in] geoElem Element geometry array
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] esup Elements-surrounding-nodes connectivity
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] ndofel Vector of local number of degrees of freedome
    //! \param[in] gid Local->global node id map
    //! \param[in] bid Local chare-boundary node ids (value) associated to
    //!   global node ids (key)
    //! \param[in] uNodalExtrm Chare-boundary nodal extrema for conservative
    //!   variables
    //! \param[in] pNodalExtrm Chare-boundary nodal extrema for primitive
    //!   variables
    //! \param[in] mtInv Inverse of Taylor mass matrix
    //! \param[in,out] U Solution vector at recent time step
    //! \param[in,out] P Vector of primitives at recent time step
    //! \param[in,out] shockmarker Vector of shock-marker values
    void limit( [[maybe_unused]] tk::real t,
                const tk::Fields& geoFace,
                const tk::Fields& geoElem,
                const inciter::FaceData& fd,
                const std::map< std::size_t, std::vector< std::size_t > >& esup,
                const std::vector< std::size_t >& inpoel,
                const tk::UnsMesh::Coords& coord,
                const std::vector< std::size_t >& ndofel,
                const std::vector< std::size_t >& gid,
                const std::unordered_map< std::size_t, std::size_t >& bid,
                const std::vector< std::vector<tk::real> >& uNodalExtrm,
                const std::vector< std::vector<tk::real> >& pNodalExtrm,
                const std::vector< std::vector<tk::real> >& mtInv,
                tk::Fields& U,
                tk::Fields& P,
                std::vector< std::size_t >& shockmarker ) const
    {
      Assert( U.nunk() == P.nunk(), "Number of unknowns in solution "
              "vector and primitive vector at recent time step incorrect" );

      const auto limiter = g_inputdeck.get< tag::limiter >();
      auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
      const auto rdof = g_inputdeck.get< tag::rdof >();
      const auto& solidx = g_inputdeck.get<
        tag::matidxmap, tag::solidx >();

      // limit vectors of conserved and primitive quantities
      if (limiter == ctr::LimiterType::SUPERBEEP1)
      {
        SuperbeeMultiMat_P1( fd.Esuel(), inpoel, ndofel,
          coord, solidx, U, P, nmat );
      }
      else if (limiter == ctr::LimiterType::VERTEXBASEDP1 && rdof == 4)
      {
        VertexBasedMultiMat_P1( esup, inpoel, ndofel, fd.Esuel().size()/4,
          m_mat_blk, fd, geoFace, geoElem, coord, flux, solidx, U, P,
          nmat, shockmarker );
      }
      else if (limiter == ctr::LimiterType::VERTEXBASEDP1 && rdof == 10)
      {
        VertexBasedMultiMat_P2( esup, inpoel, ndofel, fd.Esuel().size()/4,
          m_mat_blk, fd, geoFace, geoElem, coord, gid, bid,
          uNodalExtrm, pNodalExtrm, mtInv, flux, solidx, U, P, nmat,
          shockmarker );
      }
      else if (limiter != ctr::LimiterType::NOLIMITER)
      {
        Throw("Limiter type not configured for multimat.");
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
    //!   See Pandare et al. (2023). On the Design of Stable,
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
      auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();

      Assert( unk.nunk() == prim.nunk(), "Number of unknowns in solution "
              "vector and primitive vector at recent time step incorrect" );
      Assert( unk.nprop() == rdof*m_ncomp, "Number of components in solution "
              "vector must equal "+ std::to_string(rdof*m_ncomp) );
      Assert( prim.nprop() == rdof*m_nprim, "Number of components in vector of "
              "primitive quantities must equal "+ std::to_string(rdof*m_nprim) );

      correctLimConservMultiMat(nielem, m_mat_blk, nmat, inpoel,
        coord, geoElem, prim, unk);
    }

    //! Return cell-average deformation gradient tensor
    //! \param[in] unk Solution vector at recent time step
    //! \param[in] nielem Number of internal elements
    //! \details This function returns the bulk cell-average inverse
    //!   deformation gradient tensor
    std::array< std::vector< tk::real >, 9 > cellAvgDeformGrad(
      const tk::Fields& unk,
      std::size_t nielem ) const
    {
      const auto rdof = g_inputdeck.get< tag::rdof >();
      auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
      const auto& solidx = g_inputdeck.get< tag::matidxmap, tag::solidx >();

      std::array< std::vector< tk::real >, 9 > gb;
      if (inciter::haveSolid(nmat, solidx)) {
        for (auto& gij : gb)
          gij.resize(nielem, 0.0);
        for (std::size_t e=0; e<nielem; ++e) {
          for (std::size_t k=0; k<nmat; ++k) {
            if (solidx[k] > 0) {
              for (std::size_t i=0; i<3; ++i)
                for (std::size_t j=0; j<3; ++j)
                  gb[3*i+j][e] += unk(e, volfracDofIdx(nmat,k,rdof,0)) *
                    unk(e,deformDofIdx(nmat,solidx[k],i,j,rdof,0));
            }
          }
        }
      }

      return gb;
    }


    //! Compute right hand side
    //! \param[in] t Physical time
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
      auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
      const auto intsharp =
        g_inputdeck.get< tag::multimat, tag::intsharp >();
      const auto& solidx = inciter::g_inputdeck.get<
        tag::matidxmap, tag::solidx >();
      auto nsld = numSolids(nmat, solidx);

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

      // Allocate space for Riemann derivatives used in non-conservative terms.
      // The following Riemann derivatives are stored, in order:
      // 1) 3*nmat terms: derivatives of partial pressure of each material,
      //    for the energy equations.
      // 2) ndof terms: derivatives of Riemann velocity times the basis
      //    function, for the volume fraction equations.
      // 3) nmat*3*3*9 terms: 3 derivatives of u_l*g_ij for each material, for
      //    the deformation gradient equations.
      // 4) 3*nsld terms: 3 derivatives of \alpha \sigma_ij for each solid
      //    material, for the energy equations.
      std::vector< std::vector< tk::real > >
        riemannDeriv(3*nmat+ndof+3*nsld, std::vector<tk::real>(U.nunk(),0.0));

      // vectors to store the data of riemann velocity used for reconstruction
      // in volume fraction equation
      std::vector< std::vector< tk::real > > vriem( U.nunk() );
      std::vector< std::vector< tk::real > > riemannLoc( U.nunk() );

      // configure a no-op lambda for prescribed velocity
      auto velfn = []( ncomp_t, tk::real, tk::real, tk::real, tk::real ){
        return tk::VelFn::result_type(); };

      // compute internal surface flux integrals
      tk::surfInt( nmat, m_mat_blk, t, ndof, rdof, inpoel, solidx,
                   coord, fd, geoFace, geoElem, m_riemann, velfn, U, P, ndofel,
                   dt, R, vriem, riemannLoc, riemannDeriv, intsharp );

      // compute optional source term
      tk::srcInt( m_mat_blk, t, ndof, fd.Esuel().size()/4, inpoel,
                  coord, geoElem, Problem::src, ndofel, R, nmat );

      if(ndof > 1)
        // compute volume integrals
        tk::volInt( nmat, t, m_mat_blk, ndof, rdof, nelem,
                    inpoel, coord, geoElem, flux, velfn, U, P, ndofel, R,
                    intsharp );

      // compute boundary surface flux integrals
      for (const auto& b : m_bc)
        tk::bndSurfInt( nmat, m_mat_blk, ndof, rdof,
                        b.first, fd, geoFace, geoElem, inpoel, coord, t,
                        m_riemann, velfn, b.second, U, P, ndofel, R, vriem,
                        riemannLoc, riemannDeriv, intsharp );

      Assert( riemannDeriv.size() == 3*nmat+ndof+3*nsld, "Size of "
              "Riemann derivative vector incorrect" );

      // get derivatives from riemannDeriv
      for (std::size_t k=0; k<riemannDeriv.size(); ++k)
      {
        Assert( riemannDeriv[k].size() == U.nunk(), "Riemann derivative vector "
                "for non-conservative terms has incorrect size" );
        for (std::size_t e=0; e<U.nunk(); ++e)
          riemannDeriv[k][e] /= geoElem(e, 0);
      }

      // compute volume integrals of non-conservative terms
      tk::nonConservativeInt( nmat, m_mat_blk, ndof, rdof, nelem,
                              inpoel, coord, geoElem, U, P, riemannDeriv,
                              ndofel, R, intsharp );

      // Code below commented until details about the form of these terms in the
      // \alpha_k g_k equations are sorted out.
      // // Compute integrals for inverse deformation in solid materials
      // if (inciter::haveSolid(nmat, solidx))
      //   tk::solidTermsVolInt( nmat, m_mat_blk, ndof, rdof, nelem,
      //                         inpoel, coord, geoElem, U, P, ndofel, dt, R);

      // compute finite pressure relaxation terms
      if (g_inputdeck.get< tag::multimat, tag::prelax >())
      {
        const auto ct = g_inputdeck.get< tag::multimat,
                                         tag::prelax_timescale >();
        tk::pressureRelaxationInt( nmat, m_mat_blk, ndof,
                                   rdof, nelem, inpoel, coord, geoElem, U, P,
                                   ndofel, ct, R, intsharp );
      }
    }

    //! Evaluate the adaptive indicator and mark the ndof for each element
    //! \param[in] nunk Number of unknowns
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] unk Array of unknowns
    //! \param[in] prim Array of primitive quantities
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
                    const tk::Fields& prim,
                    inciter::ctr::PrefIndicatorType indicator,
                    std::size_t ndof,
                    std::size_t ndofmax,
                    tk::real tolref,
                    std::vector< std::size_t >& ndofel ) const
    {
      const auto& esuel = fd.Esuel();
      auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();

      if(indicator == inciter::ctr::PrefIndicatorType::SPECTRAL_DECAY)
        spectral_decay(nmat, nunk, esuel, unk, prim, ndof, ndofmax, tolref,
          ndofel);
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
                 const std::size_t nielem ) const
    {
      const auto ndof = g_inputdeck.get< tag::ndof >();
      auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();

      auto mindt = timeStepSizeMultiMat( m_mat_blk, fd.Esuf(), geoFace, geoElem,
        nielem, nmat, U, P);

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
      return mindt;
    }

    //! Extract the velocity field at cell nodes. Currently unused.
    //! \param[in] U Solution vector at recent time step
    //! \param[in] N Element node indices
    //! \return Array of the four values of the velocity field
    std::array< std::array< tk::real, 4 >, 3 >
    velocity( const tk::Fields& U,
              const std::array< std::vector< tk::real >, 3 >&,
              const std::array< std::size_t, 4 >& N ) const
    {
      const auto rdof = g_inputdeck.get< tag::rdof >();
      auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();

      std::array< std::array< tk::real, 4 >, 3 > v;
      v[0] = U.extract( momentumDofIdx(nmat, 0, rdof, 0), N );
      v[1] = U.extract( momentumDofIdx(nmat, 1, rdof, 0), N );
      v[2] = U.extract( momentumDofIdx(nmat, 2, rdof, 0), N );

      std::vector< std::array< tk::real, 4 > > ar;
      ar.resize(nmat);
      for (std::size_t k=0; k<nmat; ++k)
        ar[k] = U.extract( densityDofIdx(nmat, k, rdof, 0), N );

      std::array< tk::real, 4 > r{{ 0.0, 0.0, 0.0, 0.0 }};
      for (std::size_t i=0; i<r.size(); ++i) {
        for (std::size_t k=0; k<nmat; ++k)
          r[i] += ar[k][i];
      }

      std::transform( r.begin(), r.end(), v[0].begin(), v[0].begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      std::transform( r.begin(), r.end(), v[1].begin(), v[1].begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      std::transform( r.begin(), r.end(), v[2].begin(), v[2].begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      return v;
    }

    //! Return a map that associates user-specified strings to functions
    //! \return Map that associates user-specified strings to functions that
    //!   compute relevant quantities to be output to file
    std::map< std::string, tk::GetVarFn > OutVarFn() const
    { return MultiMatOutVarFn(); }

    //! Return analytic field names to be output to file
    //! \return Vector of strings labelling analytic fields output in file
    std::vector< std::string > analyticFieldNames() const {
      auto nmat = g_inputdeck.get< eq, tag::nmat >();

      return MultiMatFieldNames(nmat);
    }

    //! Return time history field names to be output to file
    //! \return Vector of strings labelling time history fields output in file
    std::vector< std::string > histNames() const {
      return MultiMatHistNames();
    }

    //! Return surface field output going to file
    std::vector< std::vector< tk::real > >
    surfOutput( const std::map< int, std::vector< std::size_t > >&,
                tk::Fields& ) const
    {
      std::vector< std::vector< tk::real > > s; // punt for now
      return s;
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
      auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();

      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      std::vector< std::vector< tk::real > > Up(h.size());

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
        auto B = tk::eval_basis(rdof, tk::dot(J[0],dc), tk::dot(J[1],dc),
          tk::dot(J[2],dc));
        auto uhp = eval_state(m_ncomp, rdof, rdof, e, U, B);
        auto php = eval_state(m_nprim, rdof, rdof, e, P, B);

        // store solution in history output vector
        Up[j].resize(6, 0.0);
        for (std::size_t k=0; k<nmat; ++k) {
          Up[j][0] += uhp[densityIdx(nmat,k)];
          Up[j][4] += uhp[energyIdx(nmat,k)];
          Up[j][5] += php[pressureIdx(nmat,k)];
        }
        Up[j][1] = php[velocityIdx(nmat,0)];
        Up[j][2] = php[velocityIdx(nmat,1)];
        Up[j][3] = php[velocityIdx(nmat,2)];
        ++j;
      }

      return Up;
    }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    std::vector< std::string > names() const
    {
      auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
      return MultiMatDiagNames(nmat);
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
      auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();

      tk::real sp_te(0.0);
      // sum each material total energy
      for (std::size_t k=0; k<nmat; ++k) {
        sp_te += unk(e, energyDofIdx(nmat,k,rdof,0));
      }
      return sp_te;
    }

  private:
    //! Number of components in this PDE system
    const ncomp_t m_ncomp;
    //! Number of primitive quantities stored in this PDE system
    const ncomp_t m_nprim;
    //! Riemann solver
    tk::RiemannFluxFn m_riemann;
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
    flux( ncomp_t ncomp,
          const std::vector< EOS >& mat_blk,
          const std::vector< tk::real >& ugp,
          const std::vector< std::array< tk::real, 3 > >& )
    {
      auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();

      return tk::fluxTerms(ncomp, nmat, mat_blk, ugp);
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
    //! \note The function signature must follow tk::StateFn. For multimat, the
    //!   left or right state is the vector of conserved quantities, followed by
    //!   the vector of primitive quantities appended to it.
    static tk::StateFn::result_type
    dirichlet( ncomp_t ncomp,
               const std::vector< EOS >& mat_blk,
               const std::vector< tk::real >& ul, tk::real x, tk::real y,
               tk::real z, tk::real t, const std::array< tk::real, 3 >& )
    {
      auto nmat = g_inputdeck.get< tag::multimat, tag::nmat >();
      const auto& solidx = g_inputdeck.get<
        tag::matidxmap, tag::solidx >();

      [[maybe_unused]] auto nsld = numSolids(nmat, solidx);

      auto ur = Problem::initialize( ncomp, mat_blk, x, y, z, t );
      Assert( ur.size() == ncomp, "Incorrect size for boundary state vector" );

      ur.resize(ul.size());

      tk::real rho(0.0);
      for (std::size_t k=0; k<nmat; ++k)
        rho += ur[densityIdx(nmat, k)];

      // get primitives in boundary state

      // velocity
      ur[ncomp+velocityIdx(nmat, 0)] = ur[momentumIdx(nmat, 0)] / rho;
      ur[ncomp+velocityIdx(nmat, 1)] = ur[momentumIdx(nmat, 1)] / rho;
      ur[ncomp+velocityIdx(nmat, 2)] = ur[momentumIdx(nmat, 2)] / rho;

      // material pressures
      for (std::size_t k=0; k<nmat; ++k)
      {
        auto gk = getDeformGrad(nmat, k, ur);
        ur[ncomp+pressureIdx(nmat, k)] = mat_blk[k].compute< EOS::pressure >(
          ur[densityIdx(nmat, k)], ur[ncomp+velocityIdx(nmat, 0)],
          ur[ncomp+velocityIdx(nmat, 1)], ur[ncomp+velocityIdx(nmat, 2)],
          ur[energyIdx(nmat, k)], ur[volfracIdx(nmat, k)], k, gk );
      }

      Assert( ur.size() == ncomp+nmat+3+nsld*6, "Incorrect size for appended "
              "boundary state vector" );

      return {{ std::move(ul), std::move(ur) }};
    }

    // Other boundary condition types that do not depend on "Problem" should be
    // added in BCFunctions.hpp
};

} // dg::

} // inciter::

#endif // DGMultiMat_h
