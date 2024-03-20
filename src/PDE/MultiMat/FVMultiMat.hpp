// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/FVMultiMat.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Compressible multi-material flow using finite volumes
  \details   This file implements calls to the physics operators governing
    compressible multi-material flow (with velocity equilibrium) using finite
    volume discretizations.
*/
// *****************************************************************************
#ifndef FVMultiMat_h
#define FVMultiMat_h

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
#include "Inciter/InputDeck/New2InputDeck.hpp"
#include "Integrate/Basis.hpp"
#include "Integrate/Quadrature.hpp"
#include "Integrate/Initialize.hpp"
#include "Integrate/Mass.hpp"
#include "Integrate/Surface.hpp"
#include "Integrate/Boundary.hpp"
#include "Integrate/Volume.hpp"
#include "Integrate/MultiMatTerms.hpp"
#include "Integrate/Source.hpp"
#include "RiemannChoice.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "Reconstruction.hpp"
#include "Limiter.hpp"
#include "Problem/FieldOutput.hpp"
#include "Problem/BoxInitialization.hpp"
#include "MultiMat/BCFunctions.hpp"
#include "MultiMat/MiscMultiMatFns.hpp"

namespace inciter {

extern ctr::New2InputDeck g_inputdeck;

namespace fv {

//! \brief MultiMat used polymorphically with tk::FVPDE
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Physics - physics configuration, see PDE/MultiMat/Physics.h
//!   - Problem - problem configuration, see PDE/MultiMat/Problem.h
//! \note The default physics is Euler, set in inciter::deck::check_multimat()
template< class Physics, class Problem >
class MultiMat {

  private:
    using eq = newtag::multimat;

  public:
    //! Constructor
    explicit MultiMat() :
      m_physics(),
      m_ncomp( g_inputdeck.get< newtag::ncomp >() ),
      m_riemann( multimatRiemannSolver(
        g_inputdeck.get< newtag::flux >() ) )
    {
      // associate boundary condition configurations with state functions
      brigand::for_each< newtag::bclist::Keys >( ConfigBC( m_bc,
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
      auto nmat = g_inputdeck.get< newtag::multimat, newtag::nmat >();
      // multimat needs individual material pressures and velocities currently
      return (nmat+3);
    }

    //! Find the number of materials set up for this PDE system
    //! \return The number of materials set up for this PDE system
    std::size_t nmat() const
    {
      auto nmat = g_inputdeck.get< newtag::multimat, newtag::nmat >();
      return nmat;
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

      const auto rdof = g_inputdeck.get< newtag::rdof >();
      const auto& ic = g_inputdeck.get< newtag::ic >();
      const auto& icbox = ic.get< newtag::box >();
      const auto& icmbk = ic.get< newtag::meshblock >();

      const auto& bgpre = ic.get< newtag::pressure >();
      const auto& bgtemp = ic.get< newtag::temperature >();

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
                { b.template get< newtag::xmin >(), b.template get< newtag::xmax >(),
                  b.template get< newtag::ymin >(), b.template get< newtag::ymax >(),
                  b.template get< newtag::zmin >(), b.template get< newtag::zmax >() };
              auto V_ex = (box[1]-box[0]) * (box[3]-box[2]) * (box[5]-box[4]);
              for (std::size_t c=0; c<m_ncomp; ++c) {
                auto mark = c*rdof;
                s[c] = unk(e,mark);
                // set high-order DOFs to zero
                for (std::size_t i=1; i<rdof; ++i)
                  unk(e,mark+i) = 0.0;
              }
              initializeBox<newtag::newbox>( m_mat_blk, V_ex, t, b, bgpre,
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
        for (const auto& b : icmbk) { // for all blocks
          auto blid = b.get< newtag::blockid >();
          auto V_ex = b.get< newtag::volume >();
          if (elemblkid.find(blid) != elemblkid.end()) {
            const auto& elset = tk::cref_find(elemblkid, blid);
            if (elset.find(e) != elset.end()) {
              initializeBox<newtag::newmeshblock>( m_mat_blk, V_ex, t, b,
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

    //! Compute the left hand side block-diagonal mass matrix
    //! \param[in] geoElem Element geometry array
    //! \param[in,out] l Block diagonal mass matrix
    void lhs( const tk::Fields& geoElem, tk::Fields& l ) const {
      const auto nelem = geoElem.nunk();
      for (std::size_t e=0; e<nelem; ++e)
        for (ncomp_t c=0; c<m_ncomp; ++c)
          l(e, c) = geoElem(e,0);
    }

    //! Update the primitives for this PDE system
    //! \param[in] unk Array of unknowns
    //! \param[in,out] prim Array of primitives
    //! \param[in] nielem Number of internal elements
    //! \details This function computes and stores the dofs for primitive
    //!   quantities, which are required for obtaining reconstructed states used
    //!   in the Riemann solver. See /PDE/Riemann/AUSM.hpp, where the
    //!   normal velocity for advection is calculated from independently
    //!   reconstructed velocities.
    void updatePrimitives( const tk::Fields& unk,
                           tk::Fields& prim,
                           std::size_t nielem ) const
    {
      const auto rdof = g_inputdeck.get< newtag::rdof >();
      auto nmat = g_inputdeck.get< newtag::multimat, newtag::nmat >();

      Assert( unk.nunk() == prim.nunk(), "Number of unknowns in solution "
              "vector and primitive vector at recent time step incorrect" );
      Assert( unk.nprop() == rdof*m_ncomp, "Number of components in solution "
              "vector must equal "+ std::to_string(rdof*m_ncomp) );
      Assert( prim.nprop() == rdof*nprim(), "Number of components in vector of "
              "primitive quantities must equal "+ std::to_string(rdof*nprim()) );

      for (std::size_t e=0; e<nielem; ++e)
      {
        // cell-average bulk density
        tk::real rhob(0.0);
        for (std::size_t k=0; k<nmat; ++k)
        {
          rhob += unk(e, densityDofIdx(nmat, k, rdof, 0));
        }

        // cell-average velocity
        std::array< tk::real, 3 >
          vel{{ unk(e, momentumDofIdx(nmat, 0, rdof, 0))/rhob,
                unk(e, momentumDofIdx(nmat, 1, rdof, 0))/rhob,
                unk(e, momentumDofIdx(nmat, 2, rdof, 0))/rhob }};

        for (std::size_t idir=0; idir<3; ++idir)
        {
          prim(e, velocityDofIdx(nmat, idir, rdof, 0)) = vel[idir];
          for (std::size_t idof=1; idof<rdof; ++idof)
            prim(e, velocityDofIdx(nmat, idir, rdof, idof)) = 0.0;
        }

        // cell-average material pressure
        for (std::size_t k=0; k<nmat; ++k)
        {
          tk::real arhomat = unk(e, densityDofIdx(nmat, k, rdof, 0));
          tk::real arhoemat = unk(e, energyDofIdx(nmat, k, rdof, 0));
          tk::real alphamat = unk(e, volfracDofIdx(nmat, k, rdof, 0));
          auto agmat = getDeformGrad(nmat, k, unk.extract(e));
          prim(e, pressureDofIdx(nmat, k, rdof, 0)) =
            m_mat_blk[k].compute< EOS::pressure >( arhomat, vel[0], vel[1],
            vel[2], arhoemat, alphamat, k, agmat );
          prim(e, pressureDofIdx(nmat, k, rdof, 0)) =
            constrain_pressure( m_mat_blk,
            prim(e, pressureDofIdx(nmat, k, rdof, 0)), arhomat, alphamat, k);
          for (std::size_t idof=1; idof<rdof; ++idof)
            prim(e, pressureDofIdx(nmat, k, rdof, idof)) = 0.0;
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
      [[maybe_unused]] const auto rdof = g_inputdeck.get< newtag::rdof >();
      auto nmat = g_inputdeck.get< newtag::multimat, newtag::nmat >();

      Assert( unk.nunk() == prim.nunk(), "Number of unknowns in solution "
              "vector and primitive vector at recent time step incorrect" );
      Assert( unk.nprop() == rdof*m_ncomp, "Number of components in solution "
              "vector must equal "+ std::to_string(rdof*m_ncomp) );
      Assert( prim.nprop() == rdof*nprim(), "Number of components in vector of "
              "primitive quantities must equal "+ std::to_string(rdof*nprim()) );
      Assert( (g_inputdeck.get< newtag::ndof >()) <= 4, "High-order "
              "discretizations not set up for multimat cleanTraceMaterial()" );

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
    void reconstruct( const tk::Fields& geoElem,
                      const inciter::FaceData& fd,
                      const std::map< std::size_t, std::vector< std::size_t > >&
                        esup,
                      const std::vector< std::size_t >& inpoel,
                      const tk::UnsMesh::Coords& coord,
                      tk::Fields& U,
                      tk::Fields& P ) const
    {
      const auto rdof = g_inputdeck.get< newtag::rdof >();
      const auto nelem = fd.Esuel().size()/4;
      auto nmat = g_inputdeck.get< newtag::multimat, newtag::nmat >();

      Assert( U.nprop() == rdof*m_ncomp, "Number of components in solution "
              "vector must equal "+ std::to_string(rdof*m_ncomp) );

      //----- reconstruction of conserved quantities -----
      //--------------------------------------------------
      // specify how many variables need to be reconstructed
      std::vector< std::size_t > vars;
      for (std::size_t k=0; k<nmat; ++k) {
        vars.push_back(volfracIdx(nmat,k));
        vars.push_back(densityIdx(nmat,k));
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
      vars.clear();
      for (std::size_t c=0; c<nprim(); ++c) vars.push_back(c);
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

    //! Limit second-order solution, and primitive quantities separately
    //! \param[in] geoFace Face geometry array
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] esup Elements-surrounding-nodes connectivity
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] srcFlag Whether the energy source was added
    //! \param[in,out] U Solution vector at recent time step
    //! \param[in,out] P Vector of primitives at recent time step
    void limit( const tk::Fields& geoFace,
                const inciter::FaceData& fd,
                const std::map< std::size_t, std::vector< std::size_t > >& esup,
                const std::vector< std::size_t >& inpoel,
                const tk::UnsMesh::Coords& coord,
                const std::vector< int >& srcFlag,
                tk::Fields& U,
                tk::Fields& P ) const
    {
      Assert( U.nunk() == P.nunk(), "Number of unknowns in solution "
              "vector and primitive vector at recent time step incorrect" );

      const auto limiter = g_inputdeck.get< newtag::limiter >();
      auto nmat = g_inputdeck.get< newtag::multimat, newtag::nmat >();
      const auto& solidx = g_inputdeck.get<
        newtag::matidxmap, newtag::solidx >();

      // limit vectors of conserved and primitive quantities
      if (limiter == ctr::LimiterType::VERTEXBASEDP1)
      {
        VertexBasedMultiMat_FV( esup, inpoel, fd.Esuel().size()/4,
          coord, srcFlag, solidx, U, P, nmat );
        PositivityPreservingMultiMat_FV( inpoel, fd.Esuel().size()/4, nmat,
          m_mat_blk, coord, geoFace, U, P );
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
      [[maybe_unused]] const auto rdof = g_inputdeck.get< newtag::rdof >();
      auto nmat = g_inputdeck.get< newtag::multimat, newtag::nmat >();

      Assert( unk.nunk() == prim.nunk(), "Number of unknowns in solution "
              "vector and primitive vector at recent time step incorrect" );
      Assert( unk.nprop() == rdof*m_ncomp, "Number of components in solution "
              "vector must equal "+ std::to_string(rdof*m_ncomp) );
      Assert( prim.nprop() == rdof*nprim(), "Number of components in vector of "
              "primitive quantities must equal "+ std::to_string(rdof*nprim()) );

      correctLimConservMultiMat(nielem, m_mat_blk, nmat, inpoel,
        coord, geoElem, prim, unk);
    }

    //! Compute right hand side
    //! \param[in] t Physical time
    //! \param[in] geoFace Face geometry array
    //! \param[in] geoElem Element geometry array
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] elemblkid Element ids associated with mesh block ids where
    //!   user ICs are set
    //! \param[in] U Solution vector at recent time step
    //! \param[in] P Primitive vector at recent time step
    //! \param[in,out] R Right-hand side vector computed
    //! \param[in,out] srcFlag Whether the energy source was added
    void rhs( tk::real t,
              const tk::Fields& geoFace,
              const tk::Fields& geoElem,
              const inciter::FaceData& fd,
              const std::vector< std::size_t >& inpoel,
              const tk::UnsMesh::Coords& coord,
              const std::unordered_map< std::size_t, std::set< std::size_t > >&
                elemblkid,
              const tk::Fields& U,
              const tk::Fields& P,
              tk::Fields& R,
              std::vector< int >& srcFlag ) const
    {
      const auto rdof = g_inputdeck.get< newtag::rdof >();
      auto nmat = g_inputdeck.get< newtag::multimat, newtag::nmat >();
      const auto intsharp =
        g_inputdeck.get< newtag::multimat, newtag::intsharp >();

      const auto nelem = fd.Esuel().size()/4;

      Assert( U.nunk() == P.nunk(), "Number of unknowns in solution "
              "vector and primitive vector at recent time step incorrect" );
      Assert( U.nunk() == R.nunk(), "Number of unknowns in solution "
              "vector and right-hand side at recent time step incorrect" );
      Assert( U.nprop() == rdof*m_ncomp, "Number of components in solution "
              "vector must equal "+ std::to_string(rdof*m_ncomp) );
      Assert( P.nprop() == rdof*nprim(), "Number of components in primitive "
              "vector must equal "+ std::to_string(rdof*nprim()) );
      Assert( fd.Inpofa().size()/3 == fd.Esuf().size()/2,
              "Mismatch in inpofa size" );

      // set rhs to zero
      R.fill(0.0);

      // configure a no-op lambda for prescribed velocity
      auto velfn = []( ncomp_t, tk::real, tk::real, tk::real, tk::real ){
        return tk::VelFn::result_type(); };

      // compute internal surface flux (including non-conservative) integrals
      tk::surfIntFV( nmat, m_mat_blk, t, rdof, inpoel,
                     coord, fd, geoFace, geoElem, m_riemann, velfn, U, P,
                     srcFlag, R, intsharp );

      // compute boundary surface flux (including non-conservative) integrals
      for (const auto& b : m_bc)
        tk::bndSurfIntFV( nmat, m_mat_blk, rdof, b.first,
                          fd, geoFace, geoElem, inpoel, coord, t, m_riemann,
                          velfn, b.second, U, P, srcFlag, R, intsharp );

      // compute optional source term
      tk::srcIntFV( m_mat_blk, t, fd.Esuel().size()/4,
                    geoElem, Problem::src, R, nmat );

      // compute finite pressure relaxation terms
      if (g_inputdeck.get< newtag::multimat, newtag::prelax >())
      {
        const auto ct = g_inputdeck.get< newtag::multimat,
                                         newtag::prelax_timescale >();
        tk::pressureRelaxationIntFV( nmat, m_mat_blk, rdof,
                                     nelem, inpoel, coord, geoElem, U, P, ct,
                                     R );
      }

      // compute external (energy) sources
      m_physics.physSrc(nmat, t, geoElem, elemblkid, R, srcFlag);
    }

    //! Compute the minimum time step size
//    //! \param[in] fd Face connectivity and boundary conditions object
//    //! \param[in] geoFace Face geometry array
    //! \param[in] geoElem Element geometry array
    //! \param[in] U Solution vector at recent time step
    //! \param[in] P Vector of primitive quantities at recent time step
    //! \param[in] nielem Number of internal elements
    //! \param[in] srcFlag Whether the energy source was added
    //! \param[in,out] local_dte Time step size for each element (for local
    //!   time stepping)
    //! \return Minimum time step size
    //! \details The allowable dt is calculated by looking at the maximum
    //!   wave-speed in elements surrounding each face, times the area of that
    //!   face. Once the maximum of this quantity over the mesh is determined,
    //!   the volume of each cell is divided by this quantity. A minimum of this
    //!   ratio is found over the entire mesh, which gives the allowable dt.
    tk::real dt( const inciter::FaceData& /*fd*/,
                 const tk::Fields& /*geoFace*/,
                 const tk::Fields& geoElem,
                 const tk::Fields& U,
                 const tk::Fields& P,
                 const std::size_t nielem,
                 const std::vector< int >& srcFlag,
                 std::vector< tk::real >& local_dte ) const
    {
      auto nmat = g_inputdeck.get< newtag::multimat, newtag::nmat >();

      // obtain dt restrictions from all physics
      auto dt_e = timeStepSizeMultiMatFV(m_mat_blk, geoElem, nielem, nmat, U,
        P, local_dte);
      auto dt_p = m_physics.dtRestriction(geoElem, nielem, srcFlag);

      return std::min(dt_e, dt_p);
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
      const auto rdof = g_inputdeck.get< newtag::rdof >();
      auto nmat = g_inputdeck.get< newtag::multimat, newtag::nmat >();

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

    //! Return analytic field names to be output to file
    //! \return Vector of strings labelling analytic fields output in file
    std::vector< std::string > analyticFieldNames() const {
      auto nmat = g_inputdeck.get< eq, newtag::nmat >();

      return MultiMatFieldNames(nmat);
    }

    //! Return surface field names to be output to file
    //! \return Vector of strings labelling surface fields output in file
    std::vector< std::string > surfNames() const
    { return MultiMatSurfNames(); }

    //! Return time history field names to be output to file
    //! \return Vector of strings labelling time history fields output in file
    std::vector< std::string > histNames() const {
      return MultiMatHistNames();
    }

    //! Return surface field output going to file
    std::vector< std::vector< tk::real > >
    surfOutput( const inciter::FaceData& fd,
      const tk::Fields& U,
      const tk::Fields& P ) const
    {
      const auto rdof = g_inputdeck.get< newtag::rdof >();
      auto nmat = g_inputdeck.get< newtag::multimat, newtag::nmat >();

      return MultiMatSurfOutput( nmat, rdof, fd, U, P );
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
      const auto rdof = g_inputdeck.get< newtag::rdof >();
      auto nmat = g_inputdeck.get< newtag::multimat, newtag::nmat >();

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
        auto php = eval_state(nprim(), rdof, rdof, e, P, B);

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
      auto nmat = g_inputdeck.get< newtag::multimat, newtag::nmat >();
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
      const auto rdof = g_inputdeck.get< newtag::rdof >();
      auto nmat = g_inputdeck.get< newtag::multimat, newtag::nmat >();

      tk::real sp_te(0.0);
      // sum each material total energy
      for (std::size_t k=0; k<nmat; ++k) {
        sp_te += unk(e, energyDofIdx(nmat,k,rdof,0));
      }
      return sp_te;
    }

    //! Compute relevant sound speed for output
    //! \param[in] nielem Number of internal elements
    //! \param[in] U Solution vector at recent time step
    //! \param[in] P Primitive vector at recent time step
    //! \param[in,out] ss Sound speed vector
    void soundspeed(
      std::size_t nielem,
      const tk::Fields& U,
      const tk::Fields& P,
      std::vector< tk::real >& ss) const
    {
      Assert( ss.size() == nielem, "Size of sound speed vector incorrect " );

      const auto ndof = g_inputdeck.get< newtag::ndof >();
      const auto rdof = g_inputdeck.get< newtag::rdof >();
      auto nmat = g_inputdeck.get< newtag::multimat, newtag::nmat >();
      std::size_t ncomp = U.nprop()/rdof;
      std::size_t nprim = P.nprop()/rdof;

      std::vector< tk::real > ugp(ncomp, 0.0), pgp(nprim, 0.0);

      for (std::size_t e=0; e<nielem; ++e) {
        // basis function at centroid
        std::vector< tk::real > B(rdof, 0.0);
        B[0] = 1.0;

        // get conserved quantities
        ugp = eval_state(ncomp, rdof, ndof, e, U, B);
        // get primitive quantities
        pgp = eval_state(nprim, rdof, ndof, e, P, B);

        // acoustic speed (this should be consistent with time-step calculation)
        ss[e] = 0.0;
        for (std::size_t k=0; k<nmat; ++k)
        {
          if (ugp[volfracIdx(nmat, k)] > 1.0e-04) {
            ss[e] = std::max( ss[e], m_mat_blk[k].compute< EOS::soundspeed >(
              ugp[densityIdx(nmat, k)], pgp[pressureIdx(nmat, k)],
              ugp[volfracIdx(nmat, k)], k ) );
          }
        }
      }
    }

  private:
    //! Physics policy
    const Physics m_physics;
    //! Number of components in this PDE system
    const ncomp_t m_ncomp;
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
      auto nmat = g_inputdeck.get< newtag::multimat, newtag::nmat >();

      return tk::fluxTerms(ncomp, nmat, mat_blk, ugp);
    }

    //! \brief Boundary state function providing the left and right state of a
    //!   face at Dirichlet boundaries
    //! \param[in] ncomp Number of scalar components in this PDE system
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
      auto nmat = g_inputdeck.get< newtag::multimat, newtag::nmat >();

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
        auto agk = getDeformGrad(nmat, k, ur);
        ur[ncomp+pressureIdx(nmat, k)] = mat_blk[k].compute< EOS::pressure >(
          ur[densityIdx(nmat, k)], ur[ncomp+velocityIdx(nmat, 0)],
          ur[ncomp+velocityIdx(nmat, 1)], ur[ncomp+velocityIdx(nmat, 2)],
          ur[energyIdx(nmat, k)], ur[volfracIdx(nmat, k)], k, agk );
      }

      Assert( ur.size() == ncomp+nmat+3, "Incorrect size for appended "
              "boundary state vector" );

      return {{ std::move(ul), std::move(ur) }};
    }

    // Other boundary condition types that do not depend on "Problem" should be
    // added in BCFunctions.hpp
};

} // fv::

} // inciter::

#endif // FVMultiMat_h
