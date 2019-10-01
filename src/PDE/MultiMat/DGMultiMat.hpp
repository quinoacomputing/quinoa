// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/DGMultiMat.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Compressible multi-material flow using discontinuous Galerkin
    finite elements
  \details   This file implements calls to the physics operators governing
    compressible multi-material flow using discontinuous Galerkin
    discretizations.
*/
// *****************************************************************************
#ifndef MultiMatDG_h
#define MultiMatDG_h

#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <map>

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
#include "Integrate/Riemann/AUSM.hpp"
#include "EoS/EoS.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "Reconstruction.hpp"
#include "Limiter.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

namespace dg {

//! \brief MultiMat used polymorphically with tk::DGPDE
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Physics - physics configuration, see PDE/MultiMat/Physics.h
//!   - Problem - problem configuration, see PDE/MultiMat/Problem.h
//! \note The default physics is velocity equilibrium (veleq), set in
//!   inciter::deck::check_multimat()
template< class Physics, class Problem >
class MultiMat {

  private:
    using ncomp_t = kw::ncomp::info::expect::type;
    using bcconf_t = kw::sideset::info::expect::type;
    using eq = tag::multimat;

    //! Extract BC configuration ignoring if BC not specified
    //! \param[in] c Equation system index (among multiple systems configured)
    //! \return Vector of BC config of type bcconf_t used to apply BCs for all
    //!   scalar components this Euler eq system is configured for
    //! \note A more preferable way of catching errors such as this function
    //!   hides is during parsing, so that we don't even get here if BCs are not
    //!   correctly specified. For now we simply ignore if BCs are not
    //!   specified by allowing empty BC vectors from the user input.
    template< typename bctag >
    std::vector< bcconf_t >
    config( ncomp_t c ) {
      std::vector< bcconf_t > bc;
      const auto& v = g_inputdeck.get< tag::param, eq, bctag >();
      if (v.size() > c) bc = v[c];
      return bc;
    }

  public:
    //! Constructor
    //! \param[in] c Equation system index (among multiple systems configured)
    explicit MultiMat( ncomp_t c ) :
      m_system( c ),
      m_ncomp( g_inputdeck.get< tag::component, eq >().at(c) ),
      m_offset( g_inputdeck.get< tag::component >().offset< eq >(c) ),
      //m_riemann( tk::cref_find( RiemannSolvers(),
      //             g_inputdeck.get< tag::discr, tag::flux >() ) ),
      m_bcdir( config< tag::bcdir >( c ) ),
      m_bcsym( config< tag::bcsym >( c ) ),
      m_bcextrapolate( config< tag::bcextrapolate >( c ) )
    {}

    //! Find the number of primitive quantities required for this PDE system
    //! \return The number of primitive quantities required to be stored for
    //!   this PDE system
    std::size_t nprim() const
    {
      const auto nmat =
        g_inputdeck.get< tag::param, tag::multimat, tag::nmat >()[m_system];
      // multimat needs individual material pressures and velocities currently
      return (nmat+3);
    }

    //! Initalize the compressible flow equations, prepare for time integration
    //! \param[in] L Block diagonal mass matrix
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in,out] unk Array of unknowns
    //! \param[in] t Physical time
    //! \param[in] nielem Number of internal elements
    void initialize( const tk::Fields& L,
                     const std::vector< std::size_t >& inpoel,
                     const tk::UnsMesh::Coords& coord,
                     tk::Fields& unk,
                     tk::real t,
                     const std::size_t nielem ) const
    {
      tk::initialize( m_system, m_ncomp, m_offset, L, inpoel, coord,
                      Problem::solution, unk, t, nielem );
    }

    //! Compute the left hand side block-diagonal mass matrix
    //! \param[in] geoElem Element geometry array
    //! \param[in,out] l Block diagonal mass matrix
    void lhs( const tk::Fields& geoElem, tk::Fields& l ) const {
      tk::mass( m_ncomp, m_offset, geoElem, l );
    }

    //! Update the primitives for this PDE system
    //! \param[in] unk Array of unknowns
    //! \param[in,out] prim Array of primitives
    //! \param[in] nielem Number of internal elements
    //! \details This function computes and stores the dofs for primitive
    //!   quantities, which are required for obtaining reconstructed states used
    //!   in the Riemann solver. See /PDE/Integrate/Riemann/AUSM.hpp, where the
    //!   normal velocity for advection is calculated from independently
    //!   reconstructed velocities.
    void updatePrimitives( const tk::Fields& unk,
                           tk::Fields& prim,
                           std::size_t nielem ) const
    {
      const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();
      const auto nmat =
        g_inputdeck.get< tag::param, tag::multimat, tag::nmat >()[m_system];

      Assert( unk.nunk() == prim.nunk(), "Number of unknowns in solution "
              "vector and primitive vector at recent time step incorrect" );
      Assert( unk.nprop() == rdof*m_ncomp, "Number of components in solution "
              "vector must equal "+ std::to_string(rdof*m_ncomp) );
      Assert( prim.nprop() == rdof*nprim(), "Number of components in vector of "
              "primitive quantities must equal "+ std::to_string(rdof*nprim()) );
      Assert( (g_inputdeck.get< tag::discr, tag::ndof >()) == 1, "High-order "
              "discretizations not set up for multimat updatePrimitives()" );

      for (std::size_t e=0; e<nielem; ++e)
      {
        // cell-average bulk density
        tk::real rhob(0.0);
        for (std::size_t k=0; k<nmat; ++k)
        {
          rhob += unk(e, densityIdx(nmat, k)*rdof, m_offset);
        }

        // cell-average velocity
        std::array< tk::real, 3 >
          vel{{ unk(e, momentumIdx(nmat, 0)*rdof, m_offset)/rhob,
                unk(e, momentumIdx(nmat, 1)*rdof, m_offset)/rhob,
                unk(e, momentumIdx(nmat, 2)*rdof, m_offset)/rhob }};

        for (std::size_t idir=0; idir<3; ++idir)
        {
          prim(e, velocityIdx(nmat, idir)*rdof, m_offset) = vel[idir];
        }

        // cell-average material pressure
        for (std::size_t k=0; k<nmat; ++k)
        {
          tk::real rhomat = unk(e, densityIdx(nmat, k)*rdof, m_offset)
            / unk(e, volfracIdx(nmat, k)*rdof, m_offset);
          tk::real rhoemat = unk(e, energyIdx(nmat, k)*rdof, m_offset)
            / unk(e, volfracIdx(nmat, k)*rdof, m_offset);
          prim(e, pressureIdx(nmat, k)*rdof, m_offset) =
            eos_pressure< tag::multimat >( m_system, rhomat, vel[0], vel[1],
              vel[2], rhoemat, k );
        }
      }
    }

    //! Reconstruct second-order solution from first-order
    //! \param[in] t Physical time
    //! \param[in] geoFace Face geometry array
    //! \param[in] geoElem Element geometry array
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in,out] U Solution vector at recent time step
    //! \param[in,out] P Vector of primitives at recent time step
    void reconstruct( tk::real t,
                      const tk::Fields& geoFace,
                      const tk::Fields& geoElem,
                      const inciter::FaceData& fd,
                      const std::vector< std::size_t >& inpoel,
                      const tk::UnsMesh::Coords& coord,
                      tk::Fields& U,
                      tk::Fields& P ) const
    {
      const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();

      Assert( U.nprop() == rdof*m_ncomp, "Number of components in solution "
              "vector must equal "+ std::to_string(rdof*m_ncomp) );
      Assert( inpoel.size()/4 == U.nunk(), "Connectivity inpoel has incorrect "
              "size" );
      Assert( fd.Inpofa().size()/3 == fd.Esuf().size()/2,
              "Mismatch in inpofa size" );

      // supported boundary condition types and associated state functions
      std::vector< std::pair< std::vector< bcconf_t >, tk::StateFn > >
        bctypes{{
          { m_bcdir, Dirichlet },
          { m_bcsym, Symmetry },
          { m_bcextrapolate, Extrapolate } }};

      // allocate and initialize matrix and vector for reconstruction:
      // lhs_ls is the left-hand side matrix for solving the least-squares
      // system using the normal equation approach, for each mesh element.
      // It is indexed as follows:
      // The first index is the element id;
      // the second index is the row id of the 3-by-3 matrix;
      // the third index is the column id of the 3-by-3 matrix.
      std::vector< std::array< std::array< tk::real, 3 >, 3 > >
        lhs_ls( U.nunk(), {{ {{0.0, 0.0, 0.0}},
                             {{0.0, 0.0, 0.0}},
                             {{0.0, 0.0, 0.0}} }} );
      // rhs_ls is the right-hand side vector for solving the least-squares
      // system using the normal equation approach, for each element.
      // It is indexed as follows:
      // The first index is the element id;
      // the second index is the scalar equation which is being reconstructed;
      // the third index is the row id of the rhs vector.
      // two rhs_ls vectors are needed for reconstructing conserved and
      // primitive quantites separately
      std::vector< std::vector< std::array< tk::real, 3 > > >
        rhsu_ls( U.nunk(), std::vector< std::array< tk::real, 3 > >
          ( m_ncomp,
            {{ 0.0, 0.0, 0.0 }} ) );
      std::vector< std::vector< std::array< tk::real, 3 > > >
        rhsp_ls( U.nunk(), std::vector< std::array< tk::real, 3 > >
          ( nprim(),
            {{ 0.0, 0.0, 0.0 }} ) );

      // reconstruct x,y,z-derivatives of unknowns. For multimat, conserved and
      // primitive quantities are reconstructed separately.
      // 0. get lhs matrix, which is only geometry dependent
      tk::lhsLeastSq_P0P1(fd, geoElem, geoFace, lhs_ls);

      // 1. internal face contributions
      tk::intLeastSq_P0P1( m_ncomp, m_offset, rdof, fd, geoElem, U, rhsu_ls );
      tk::intLeastSq_P0P1( nprim(), m_offset, rdof, fd, geoElem, P, rhsp_ls );

      // 2. boundary face contributions
      for (const auto& b : bctypes)
      {
        tk::bndLeastSqConservedVar_P0P1( m_system, m_ncomp, m_offset, rdof,
          b.first, fd, geoFace, geoElem, t, b.second, U, rhsu_ls, nprim() );
        tk::bndLeastSqPrimitiveVar_P0P1( m_system, nprim(), m_offset, rdof,
          b.first, fd, geoFace, geoElem, t, b.second, P, rhsp_ls, m_ncomp );
      }

      // 3. solve 3x3 least-squares system
      tk::solveLeastSq_P0P1( m_ncomp, m_offset, rdof, lhs_ls, rhsu_ls, U );
      tk::solveLeastSq_P0P1( nprim(), m_offset, rdof, lhs_ls, rhsp_ls, P );

      // 4. transform reconstructed derivatives to Dubiner dofs
      tk::transform_P0P1( m_ncomp, m_offset, rdof, fd.Esuel().size()/4, inpoel,
                          coord, U );
      tk::transform_P0P1( nprim(), m_offset, rdof, fd.Esuel().size()/4, inpoel,
                          coord, P );
    }

    //! Limit second-order solution, and primitive quantities separately
    //! \param[in] t Physical time
    //! \param[in] geoFace Face geometry array
    //! \param[in] geoElem Element geometry array
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] ndofel Vector of local number of degrees of freedome
    //! \param[in,out] U Solution vector at recent time step
    //! \param[in,out] P Vector of primitives at recent time step
    void limit( [[maybe_unused]] tk::real t,
                [[maybe_unused]] const tk::Fields& geoFace,
                [[maybe_unused]] const tk::Fields& geoElem,
                const inciter::FaceData& fd,
                const std::vector< std::size_t >& inpoel,
                const tk::UnsMesh::Coords& coord,
                const std::vector< std::size_t >& ndofel,
                tk::Fields& U,
                tk::Fields& P ) const
    {
      Assert( U.nunk() == P.nunk(), "Number of unknowns in solution "
              "vector and primitive vector at recent time step incorrect" );

      const auto limiter = g_inputdeck.get< tag::discr, tag::limiter >();
      const auto nmat =
        g_inputdeck.get< tag::param, tag::multimat, tag::nmat >()[m_system];

      if (limiter == ctr::LimiterType::SUPERBEEP1)
      {
        // limit solution vector
        Superbee_P1( fd.Esuel(), inpoel, ndofel, m_offset, coord, U, nmat );

        // limit vector of primitives
        Superbee_P1( fd.Esuel(), inpoel, ndofel, m_offset, coord, P );
      }
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
    //! \param[in] ndofel Vector of local number of degrees of freedome
    //! \param[in,out] R Right-hand side vector computed
    void rhs( tk::real t,
              const tk::Fields& geoFace,
              const tk::Fields& geoElem,
              const inciter::FaceData& fd,
              const std::vector< std::size_t >& inpoel,
              const tk::UnsMesh::Coords& coord,
              const tk::Fields& U,
              const tk::Fields& P,
              const std::vector< std::size_t >& ndofel,
              tk::Fields& R ) const
    {
      const auto ndof = g_inputdeck.get< tag::discr, tag::ndof >();
      const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();
      const auto nmat =
        g_inputdeck.get< tag::param, tag::multimat, tag::nmat >()[m_system];

      Assert( U.nunk() == P.nunk(), "Number of unknowns in solution "
              "vector and primitive vector at recent time step incorrect" );
      Assert( U.nunk() == R.nunk(), "Number of unknowns in solution "
              "vector and right-hand side at recent time step incorrect" );
      Assert( U.nprop() == rdof*m_ncomp, "Number of components in solution "
              "vector must equal "+ std::to_string(rdof*m_ncomp) );
      Assert( P.nprop() == rdof*nprim(), "Number of components in primitive "
              "vector must equal "+ std::to_string(rdof*nprim()) );
      Assert( R.nprop() == ndof*m_ncomp, "Number of components in right-hand "
              "side vector must equal "+ std::to_string(ndof*m_ncomp) );
      Assert( inpoel.size()/4 == U.nunk(), "Connectivity inpoel has incorrect "
              "size" );
      Assert( fd.Inpofa().size()/3 == fd.Esuf().size()/2,
              "Mismatch in inpofa size" );
      Assert( ndof == 1, "DGP1/2 not set up for multi-material" );

      // set rhs to zero
      R.fill(0.0);

      // allocate space for Riemann derivatives used in non-conservative terms
      std::vector< std::vector< tk::real > >
        riemannDeriv( 3*nmat+1, std::vector<tk::real>(U.nunk(),0.0) );

      // configure a no-op lambda for prescribed velocity
      auto velfn = [this]( ncomp_t, ncomp_t, tk::real, tk::real, tk::real ){
        return std::vector< std::array< tk::real, 3 > >( m_ncomp ); };

      // supported boundary condition types and associated state functions
      std::vector< std::pair< std::vector< bcconf_t >, tk::StateFn > > bctypes{{
        { m_bcdir, Dirichlet },
        { m_bcsym, Symmetry },
        { m_bcextrapolate, Extrapolate } }};

      // compute internal surface flux integrals
      tk::surfInt( m_system, nmat, m_offset, ndof, rdof, inpoel, coord,
                   fd, geoFace, AUSM::flux, velfn, U, P, ndofel, R,
                   riemannDeriv );

      // compute source term integrals
      tk::srcInt( m_system, m_ncomp, m_offset, t, ndof, inpoel, coord, geoElem,
                  Problem::src, ndofel, R );

      if(ndof > 1)
        // compute volume integrals
        tk::volInt( m_system, m_ncomp, m_offset, ndof, inpoel, coord, geoElem,
                    flux, velfn, U, ndofel, R );

      // compute boundary surface flux integrals
      for (const auto& b : bctypes)
        tk::bndSurfInt( m_system, nmat, m_offset, ndof, rdof, b.first,
                        fd, geoFace, inpoel, coord, t, AUSM::flux, velfn,
                        b.second, U, P, ndofel, R, riemannDeriv );

      Assert( riemannDeriv.size() == 3*nmat+1, "Size of Riemann derivative "
              "vector incorrect" );

      // get derivatives from riemannDeriv
      for (std::size_t k=0; k<riemannDeriv.size(); ++k)
      {
        Assert( riemannDeriv[k].size() == U.nunk(), "Riemann derivative vector "
                "for non-conservative terms has incorrect size" );
        for (std::size_t e=0; e<U.nunk(); ++e)
          riemannDeriv[k][e] /= geoElem(e, 0, 0);
      }

      // compute volume integrals of non-conservative terms
      tk::nonConservativeInt( m_system, nmat, m_offset, ndof, rdof, inpoel,
                              coord, geoElem, U, P, riemannDeriv, ndofel,
                              R );

      // compute finite pressure relaxation terms
      if (g_inputdeck.get< tag::param, tag::multimat, tag::prelax >()[m_system])
      {
        const auto ct = g_inputdeck.get< tag::param, tag::multimat,
                                         tag::prelax_timescale >()[m_system];
        tk::pressureRelaxationInt( m_system, m_ncomp, nmat, m_offset, ndof,
                                   rdof, geoElem, U, ndofel, ct, R );
      }
    }

    //! Compute the minimum time step size
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] geoFace Face geometry array
    //! \param[in] geoElem Element geometry array
    //! \param[in] U Solution vector at recent time step
//    //! \param[in] ndofel Vector of local number of degrees of freedom
    //! \return Minimum time step size
    tk::real dt( const std::array< std::vector< tk::real >, 3 >& coord,
                 const std::vector< std::size_t >& inpoel,
                 const inciter::FaceData& fd,
                 const tk::Fields& geoFace,
                 const tk::Fields& geoElem,
                 const std::vector< std::size_t >& /*ndofel*/,
                 const tk::Fields& U ) const
    {
      const auto ndof = g_inputdeck.get< tag::discr, tag::ndof >();
      const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();
      const auto nmat =
        g_inputdeck.get< tag::param, tag::multimat, tag::nmat >()[m_system];

      const auto& esuf = fd.Esuf();
      const auto& inpofa = fd.Inpofa();

      // Number of quadrature points for face integration
      auto ng = tk::NGfa(ndof);

      // arrays for quadrature points
      std::array< std::vector< tk::real >, 2 > coordgp;
      std::vector< tk::real > wgp;

      coordgp[0].resize( ng );
      coordgp[1].resize( ng );
      wgp.resize( ng );

      tk::real rhok, rho, u, v, w, p, a, vn, dSV_l, dSV_r;
      std::vector< tk::real > delt( U.nunk(), 0.0 );

      const auto& cx = coord[0];
      const auto& cy = coord[1];
      const auto& cz = coord[2];

      // get quadrature point weights and coordinates for triangle
      tk::GaussQuadratureTri( ng, coordgp, wgp );

      // compute internal surface maximum characteristic speed
      for (std::size_t f=0; f<esuf.size()/2; ++f)
      {
        std::size_t el = static_cast< std::size_t >(esuf[2*f]);
        auto er = esuf[2*f+1];

        // Extract the left element coordinates
        std::array< std::array< tk::real, 3>, 4 > coordel_l {{
          {{ cx[inpoel[4*el  ]], cy[inpoel[4*el  ]], cz[inpoel[4*el  ]] }},
          {{ cx[inpoel[4*el+1]], cy[inpoel[4*el+1]], cz[inpoel[4*el+1]] }},
          {{ cx[inpoel[4*el+2]], cy[inpoel[4*el+2]], cz[inpoel[4*el+2]] }},
          {{ cx[inpoel[4*el+3]], cy[inpoel[4*el+3]], cz[inpoel[4*el+3]] }} }};

        // Compute the determinant of Jacobian matrix
        auto detT_l =
          tk::Jacobian(coordel_l[0], coordel_l[1], coordel_l[2], coordel_l[3]);

        // Extract the face coordinates
        std::array< std::array< tk::real, 3>, 3 > coordfa {{
          {{ cx[ inpofa[3*f  ] ], cy[ inpofa[3*f  ] ], cz[ inpofa[3*f  ] ] }},
          {{ cx[ inpofa[3*f+1] ], cy[ inpofa[3*f+1] ], cz[ inpofa[3*f+1] ] }},
          {{ cx[ inpofa[3*f+2] ], cy[ inpofa[3*f+2] ], cz[ inpofa[3*f+2] ] }}
        }};

        dSV_l = 0.0;
        dSV_r = 0.0;

        // Gaussian quadrature
        for (std::size_t igp=0; igp<ng; ++igp)
        {
          // Compute the coordinates of quadrature point at physical domain
          auto gp = tk::eval_gp( igp, coordfa, coordgp );

          // Compute the basis function for the left element
          auto B_l = tk::eval_basis( ndof,
            tk::Jacobian(coordel_l[0], gp, coordel_l[2], coordel_l[3])/detT_l,
            tk::Jacobian(coordel_l[0], coordel_l[1], gp, coordel_l[3])/detT_l,
            tk::Jacobian(coordel_l[0], coordel_l[1], coordel_l[2], gp)/detT_l );

          auto wt = wgp[igp] * geoFace(f,0,0);

          std::array< std::vector< tk::real >, 2 > ugp;

          // left element
          for (ncomp_t c=0; c<m_ncomp; ++c)
          {
            auto mark = c*rdof;
            ugp[0].push_back( U(el, mark, m_offset) );
          }

          rho = 0.0;
          for (std::size_t k=0; k<nmat; ++k)
            rho += ugp[0][densityIdx(nmat, k)];

          u = ugp[0][momentumIdx(nmat, 0)]/rho;
          v = ugp[0][momentumIdx(nmat, 1)]/rho;
          w = ugp[0][momentumIdx(nmat, 2)]/rho;

          a = 0.0;
          for (std::size_t k=0; k<nmat; ++k)
          {
            rhok = ugp[0][densityIdx(nmat, k)]/ugp[0][volfracIdx(nmat, k)];
            p = eos_pressure< tag::multimat >( 0, rhok, u, v, w,
                                               ugp[0][energyIdx(nmat, k)]/ugp[0][volfracIdx(nmat, k)],
                                               k );
            a = std::max( a, eos_soundspeed< tag::multimat >( 0, rhok, p, k ) );
          }

          vn = u*geoFace(f,1,0) + v*geoFace(f,2,0) + w*geoFace(f,3,0);

          dSV_l = wt * (std::fabs(vn) + a);

          // right element
          if (er > -1) {

            // nodal coordinates of the right element
            std::size_t eR = static_cast< std::size_t >( er );

            // Extract the right element coordinates
            std::array< std::array< tk::real, 3>, 4 > coordel_r {{
              {{ cx[inpoel[4*eR  ]], cy[inpoel[4*eR  ]], cz[inpoel[4*eR  ]] }},
              {{ cx[inpoel[4*eR+1]], cy[inpoel[4*eR+1]], cz[inpoel[4*eR+1]] }},
              {{ cx[inpoel[4*eR+2]], cy[inpoel[4*eR+2]], cz[inpoel[4*eR+2]] }},
              {{ cx[inpoel[4*eR+3]], cy[inpoel[4*eR+3]], cz[inpoel[4*eR+3]] }}
            }};

            // Compute the determinant of Jacobian matrix
            auto detT_r =
              tk::Jacobian(coordel_r[0],coordel_r[1],coordel_r[2],coordel_r[3]);

            // Compute the coordinates of quadrature point at physical domain
            gp = tk::eval_gp( igp, coordfa, coordgp );

            // Compute the basis function for the right element
            auto B_r = tk::eval_basis( ndof,
             tk::Jacobian(coordel_r[0], gp, coordel_r[2], coordel_r[3])/detT_r,
             tk::Jacobian(coordel_r[0], coordel_r[1], gp, coordel_r[3])/detT_r,
             tk::Jacobian(coordel_r[0], coordel_r[1], coordel_r[2], gp)/detT_r);

            for (ncomp_t c=0; c<5; ++c)
            {
              auto mark = c*rdof;
              ugp[1].push_back( U(eR, mark, m_offset) );
            }

            rho = 0.0;
            for (std::size_t k=0; k<nmat; ++k)
              rho += ugp[1][densityIdx(nmat, k)];

            u = ugp[1][momentumIdx(nmat, 0)]/rho;
            v = ugp[1][momentumIdx(nmat, 1)]/rho;
            w = ugp[1][momentumIdx(nmat, 2)]/rho;

            a = 0.0;
            for (std::size_t k=0; k<nmat; ++k)
            {
              rhok = ugp[1][densityIdx(nmat, k)]/ugp[1][volfracIdx(nmat, k)];
              p = eos_pressure< tag::multimat >( 0, rhok, u, v, w,
                                                 ugp[1][energyIdx(nmat, k)]/ugp[1][volfracIdx(nmat, k)],
                                                 k );
              a = std::max( a, eos_soundspeed< tag::multimat >( 0, rhok, p, k ) );
            }

            vn = u*geoFace(f,1,0) + v*geoFace(f,2,0) + w*geoFace(f,3,0);

            dSV_r = wt * (std::fabs(vn) + a);
            delt[eR] += std::max( dSV_l, dSV_r );
          }

          delt[el] += std::max( dSV_l, dSV_r );
        }
      }

      tk::real mindt = std::numeric_limits< tk::real >::max();

      // compute allowable dt
      for (std::size_t e=0; e<U.nunk(); ++e)
      {
        mindt = std::min( mindt, geoElem(e,0,0)/delt[e] );
      }

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
      const auto rdof = g_inputdeck.get< tag::discr, tag::rdof >();
      const auto nmat =
        g_inputdeck.get< tag::param, tag::multimat, tag::nmat >()[0];

      std::array< std::array< tk::real, 4 >, 3 > v;
      v[0] = U.extract( momentumIdx(nmat, 0)*rdof, m_offset, N );
      v[1] = U.extract( momentumIdx(nmat, 1)*rdof, m_offset, N );
      v[2] = U.extract( momentumIdx(nmat, 2)*rdof, m_offset, N );

      std::vector< std::array< tk::real, 4 > > ar;
      ar.resize(nmat);
      for (std::size_t k=0; k<nmat; ++k)
        ar[k] = U.extract( densityIdx(nmat, k)*rdof, m_offset, N );

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

    //! \brief Query all side set IDs the user has configured for all components
    //!   in this PDE system
    //! \param[in,out] conf Set of unique side set IDs to add to
    void side( std::unordered_set< int >& conf ) const
    { Problem::side( conf ); }

    //! Return field names to be output to file
    //! \return Vector of strings labelling fields output in file
    std::vector< std::string > fieldNames() const
    { return Problem::fieldNames( m_ncomp ); }

    //! Return field output going to file
    //! \param[in] t Physical time
    //! \param[in] geoElem Element geometry array
    //! \param[in,out] U Solution vector at recent time step
    //! \return Vector of vectors to be output to file
    std::vector< std::vector< tk::real > >
    fieldOutput( tk::real t,
                 const tk::Fields& geoElem,
                 tk::Fields& U ) const
    {
      std::array< std::vector< tk::real >, 3 > coord;
      std::vector< tk::real > v;
      v        = geoElem.extract(0,0);
      coord[0] = geoElem.extract(1,0);
      coord[1] = geoElem.extract(2,0);
      coord[2] = geoElem.extract(3,0);

      return Problem::fieldOutput( m_system, m_ncomp, m_offset, t,
                                   0.0, v, coord, U );
    }

    //! Return nodal field output going to file
    std::vector< std::vector< tk::real > >
    avgElemToNode( const std::vector< std::size_t >& /*inpoel*/,
                   const tk::UnsMesh::Coords& /*coord*/,
                   const tk::Fields& /*geoElem*/,
                   const tk::Fields& /*U*/ ) const
    {
      std::vector< std::vector< tk::real > > out;
      return out;
    }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    std::vector< std::string > names() const
    { return Problem::names( m_ncomp ); }

    //! Return analytic solution (if defined by Problem) at xi, yi, zi, t
    //! \param[in] xi X-coordinate at which to evaluate the analytic solution
    //! \param[in] yi Y-coordinate at which to evaluate the analytic solution
    //! \param[in] zi Z-coordinate at which to evaluate the analytic solution
    //! \param[in] t Physical time at which to evaluate the analytic solution
    //! \return Vector of analytic solution at given location and time
    std::vector< tk::real >
    analyticSolution( tk::real xi, tk::real yi, tk::real zi, tk::real t ) const
    {
      auto s = Problem::solution( m_system, m_ncomp, xi, yi, zi, t );
      return std::vector< tk::real >( begin(s), end(s) );
    }

  private:
    //! Equation system index
    const ncomp_t m_system;
    //! Number of components in this PDE system
    const ncomp_t m_ncomp;
    //! Offset PDE system operates from
    const ncomp_t m_offset;
    ////! Riemann solver
    //RiemannSolver m_riemann;
    //! Dirichlet BC configuration
    const std::vector< bcconf_t > m_bcdir;
    //! Symmetric BC configuration
    const std::vector< bcconf_t > m_bcsym;
    //! Extrapolation BC configuration
    const std::vector< bcconf_t > m_bcextrapolate;

    //! Evaluate conservative part of physical flux function for this PDE system
    //! \param[in] system Equation system index
    //! \param[in] ncomp Number of scalar components in this PDE system
    //! \param[in] ugp Numerical solution at the Gauss point at which to
    //!   evaluate the flux
    //! \return Flux vectors for all components in this PDE system
    //! \note The function signature must follow tk::FluxFn
    static tk::FluxFn::result_type
    flux( ncomp_t system,
          [[maybe_unused]] ncomp_t ncomp,
          const std::vector< tk::real >& ugp,
          const std::vector< std::array< tk::real, 3 > >& )
    {
      Assert( ugp.size() == ncomp, "Size mismatch" );
      const auto nmat =
        g_inputdeck.get< tag::param, tag::multimat, tag::nmat >()[system];

      tk::real rho(0.0), p(0.0);
      for (std::size_t k=0; k<nmat; ++k)
        rho += ugp[densityIdx(nmat, k)];

      auto u = ugp[momentumIdx(nmat, 0)] / rho;
      auto v = ugp[momentumIdx(nmat, 1)] / rho;
      auto w = ugp[momentumIdx(nmat, 2)] / rho;

      std::vector< tk::real > pk( nmat, 0.0 );
      for (std::size_t k=0; k<nmat; ++k)
      {
        pk[k] = eos_pressure< tag::multimat >( system,
                                               ugp[densityIdx(nmat, k)]/ugp[volfracIdx(nmat, k)],
                                               u, v, w,
                                               ugp[energyIdx(nmat, k)]/ugp[volfracIdx(nmat, k)],
                                               k );
        p += ugp[volfracIdx(nmat, k)] * pk[k];
      }

      std::vector< std::array< tk::real, 3 > > fl( ugp.size() );

      // conservative part of momentum flux
      fl[momentumIdx(nmat, 0)][0] = ugp[momentumIdx(nmat, 0)] * u + p;
      fl[momentumIdx(nmat, 1)][0] = ugp[momentumIdx(nmat, 1)] * u;
      fl[momentumIdx(nmat, 2)][0] = ugp[momentumIdx(nmat, 2)] * u;

      fl[momentumIdx(nmat, 0)][1] = ugp[momentumIdx(nmat, 0)] * v;
      fl[momentumIdx(nmat, 1)][1] = ugp[momentumIdx(nmat, 1)] * v + p;
      fl[momentumIdx(nmat, 2)][1] = ugp[momentumIdx(nmat, 2)] * v;

      fl[momentumIdx(nmat, 0)][2] = ugp[momentumIdx(nmat, 0)] * w;
      fl[momentumIdx(nmat, 1)][2] = ugp[momentumIdx(nmat, 1)] * w;
      fl[momentumIdx(nmat, 2)][2] = ugp[momentumIdx(nmat, 2)] * w + p;

      for (std::size_t k=0; k<nmat; ++k)
      {
        // conservative part of volume-fraction flux
        fl[volfracIdx(nmat, k)][0] = 0.0;
        fl[volfracIdx(nmat, k)][1] = 0.0;
        fl[volfracIdx(nmat, k)][2] = 0.0;

        // conservative part of material continuity flux
        fl[densityIdx(nmat, k)][0] = u * ugp[densityIdx(nmat, k)];
        fl[densityIdx(nmat, k)][1] = v * ugp[densityIdx(nmat, k)];
        fl[densityIdx(nmat, k)][2] = w * ugp[densityIdx(nmat, k)];

        // conservative part of material total-energy flux
        auto hmat = ugp[energyIdx(nmat, k)] + ugp[volfracIdx(nmat, k)]*pk[k];
        fl[energyIdx(nmat, k)][0] = u * hmat;
        fl[energyIdx(nmat, k)][1] = v * hmat;
        fl[energyIdx(nmat, k)][2] = w * hmat;
      }

      // NEED TO RETURN m_ncomp flux vectors in fl, not 5

      return fl;
    }

    //! \brief Boundary state function providing the left and right state of a
    //!   face at Dirichlet boundaries
    //! \param[in] system Equation system index
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
    Dirichlet( ncomp_t system, ncomp_t ncomp, const std::vector< tk::real >& ul,
               tk::real x, tk::real y, tk::real z, tk::real t,
               const std::array< tk::real, 3 >& )
    {
      const auto nmat =
        g_inputdeck.get< tag::param, tag::multimat, tag::nmat >()[system];

      auto ur = Problem::solution( system, ncomp, x, y, z, t );
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
        tk::real rhomat = ur[densityIdx(nmat, k)] / ur[volfracIdx(nmat, k)];
        tk::real rhoemat = ur[energyIdx(nmat, k)] / ur[volfracIdx(nmat, k)];
        ur[ncomp+pressureIdx(nmat, k)] =
          eos_pressure< tag::multimat >( system, rhomat,
            ur[ncomp+velocityIdx(nmat, 0)], ur[ncomp+velocityIdx(nmat, 1)],
            ur[ncomp+velocityIdx(nmat, 2)], rhoemat, k );
      }

      Assert( ur.size() == ncomp+nmat+3, "Incorrect size for appended "
              "boundary state vector" );

      return {{ std::move(ul), std::move(ur) }};
    }

    //! \brief Boundary state function providing the left and right state of a
    //!   face at symmetry boundaries
    //! \param[in] system Equation system index
    //! \param[in] ncomp Number of scalar components in this PDE system
    //! \param[in] ul Left (domain-internal) state
    //! \param[in] fn Unit face normal
    //! \return Left and right states for all scalar components in this PDE
    //!   system
    //! \note The function signature must follow tk::StateFn. For multimat, the
    //!   left or right state is the vector of conserved quantities, followed by
    //!   the vector of primitive quantities appended to it.
    static tk::StateFn::result_type
    Symmetry( ncomp_t system, ncomp_t ncomp, const std::vector< tk::real >& ul,
              tk::real, tk::real, tk::real, tk::real,
              const std::array< tk::real, 3 >& fn )
    {
      const auto nmat =
        g_inputdeck.get< tag::param, tag::multimat, tag::nmat >()[system];

      Assert( ul.size() == ncomp+nmat+3, "Incorrect size for appended internal "
              "state vector" );

      tk::real rho(0.0);
      for (std::size_t k=0; k<nmat; ++k)
        rho += ul[densityIdx(nmat, k)];

      auto ur = ul;

      // Internal cell velocity components
      auto v1l = ul[ncomp+velocityIdx(nmat, 0)];
      auto v2l = ul[ncomp+velocityIdx(nmat, 1)];
      auto v3l = ul[ncomp+velocityIdx(nmat, 2)];
      // Normal component of velocity
      auto vnl = v1l*fn[0] + v2l*fn[1] + v3l*fn[2];
      // Ghost state velocity components
      auto v1r = v1l - 2.0*vnl*fn[0];
      auto v2r = v2l - 2.0*vnl*fn[1];
      auto v3r = v3l - 2.0*vnl*fn[2];
      // Boundary condition
      for (std::size_t k=0; k<nmat; ++k)
      {
        ur[volfracIdx(nmat, k)] = ul[volfracIdx(nmat, k)];
        ur[densityIdx(nmat, k)] = ul[densityIdx(nmat, k)];
        ur[energyIdx(nmat, k)] = ul[energyIdx(nmat, k)];
      }
      ur[momentumIdx(nmat, 0)] = rho * v1r;
      ur[momentumIdx(nmat, 1)] = rho * v2r;
      ur[momentumIdx(nmat, 2)] = rho * v3r;

      // Internal cell primitive quantities using the separately reconstructed
      // primitive quantities. This is used to get ghost state for primitive
      // quantities
      v1l = ul[ncomp+velocityIdx(nmat, 0)];
      v2l = ul[ncomp+velocityIdx(nmat, 1)];
      v3l = ul[ncomp+velocityIdx(nmat, 2)];
      // Normal component of velocity
      vnl = v1l*fn[0] + v2l*fn[1] + v3l*fn[2];
      // Ghost state velocity components
      v1r = v1l - 2.0*vnl*fn[0];
      v2r = v2l - 2.0*vnl*fn[1];
      v3r = v3l - 2.0*vnl*fn[2];

      // get primitives in boundary state
      // velocity
      ur[ncomp+velocityIdx(nmat, 0)] = v1r;
      ur[ncomp+velocityIdx(nmat, 1)] = v2r;
      ur[ncomp+velocityIdx(nmat, 2)] = v3r;
      // material pressures
      for (std::size_t k=0; k<nmat; ++k)
        ur[ncomp+pressureIdx(nmat, k)] = ul[ncomp+pressureIdx(nmat, k)];

      Assert( ur.size() == ncomp+nmat+3, "Incorrect size for appended boundary "
              "state vector" );

      return {{ std::move(ul), std::move(ur) }};
    }

    //! \brief Boundary state function providing the left and right state of a
    //!   face at extrapolation boundaries
    //! \param[in] ul Left (domain-internal) state
    //! \return Left and right states for all scalar components in this PDE
    //!   system
    //! \note The function signature must follow tk::StateFn. For multimat, the
    //!   left or right state is the vector of conserved quantities, followed by
    //!   the vector of primitive quantities appended to it.
    static tk::StateFn::result_type
    Extrapolate( ncomp_t, ncomp_t, const std::vector< tk::real >& ul,
                 tk::real, tk::real, tk::real, tk::real,
                 const std::array< tk::real, 3 >& )
    {
      return {{ ul, ul }};
    }
};

} // dg::

} // inciter::

#endif // MultiMatDG_h
