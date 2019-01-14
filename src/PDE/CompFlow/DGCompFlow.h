// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/DGCompFlow.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Compressible single-material flow using discontinuous Galerkin
     finite elements
  \details   This file implements calls to the physics operators governing
    compressible single-material flow using discontinuous Galerkin
    discretizations.
*/
// *****************************************************************************
#ifndef DGCompFlow_h
#define DGCompFlow_h

#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <map>

#include "Macro.h"
#include "Exception.h"
#include "Vector.h"
#include "ContainerUtil.h"
#include "UnsMesh.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "Integrate/Quadrature.h"
#include "Integrate/Initialize.h"
#include "Integrate/Mass.h"
#include "Integrate/Surface.h"
#include "Integrate/Boundary.h"
#include "Integrate/Volume.h"
#include "Integrate/Source.h"
#include "Integrate/Riemann/RiemannFactory.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

namespace dg {

//! \brief CompFlow used polymorphically with tk::DGPDE
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Physics - physics configuration, see PDE/CompFlow/Physics.h
//!   - Problem - problem configuration, see PDE/CompFlow/Problem.h
//! \note The default physics is Euler, set in inciter::deck::check_compflow()
template< class Physics, class Problem >
class CompFlow {

  private:
    using ncomp_t = kw::ncomp::info::expect::type;
    using bcconf_t = kw::sideset::info::expect::type;
    using eq = tag::compflow;

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
    explicit CompFlow( ncomp_t c ) :
      m_system( c ),
      m_ncomp( g_inputdeck.get< tag::component, eq >().at(c) ),
      m_offset( g_inputdeck.get< tag::component >().offset< eq >(c) ),
      m_riemann( tk::cref_find( RiemannSolvers(),
                   g_inputdeck.get< tag::discr, tag::flux >() ) ),
      m_bcdir( config< tag::bcdir >( c ) ),
      m_bcsym( config< tag::bcsym >( c ) ),
      m_bcextrapolate( config< tag::bcextrapolate >( c ) )
      //ErrChk( !m_bcdir.empty() || !m_bcsym.empty() || !m_bcextrapolate.empty(),
      //        "Boundary conditions not set in control file for DG CompFlow" );
    {}

    //! Initalize the compressible flow equations, prepare for time integration
    //! \param[in] L Block diagonal mass matrix
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in,out] unk Array of unknowns
    //! \param[in] t Physical time
    void initialize( const tk::Fields& L,
                     const std::vector< std::size_t >& inpoel,
                     const tk::UnsMesh::Coords& coord,
                     tk::Fields& unk,
                     tk::real t ) const
    {
      tk::initialize( m_system, m_ncomp, m_offset, L, inpoel, coord,
                      Problem::solution, unk, t );
    }

    //! Compute the left hand side block-diagonal mass matrix
    //! \param[in] geoElem Element geometry array
    //! \param[in,out] l Block diagonal mass matrix
    void lhs( const tk::Fields& geoElem, tk::Fields& l ) const {
      tk::mass( m_ncomp, m_offset, geoElem, l );
    }

    //! Compute right hand side
    //! \param[in] t Physical time
    //! \param[in] geoFace Face geometry array
    //! \param[in] geoElem Element geometry array
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] U Solution vector at recent time step
    //! \param[in] limFunc Limiter function for higher-order solution dofs
    //! \param[in,out] R Right-hand side vector computed
    void rhs( tk::real t,
              const tk::Fields& geoFace,
              const tk::Fields& geoElem,
              const inciter::FaceData& fd,
              const std::vector< std::size_t >& inpoel,
              const tk::UnsMesh::Coords& coord,
              const tk::Fields& U,
              const tk::Fields& limFunc,
              tk::Fields& R ) const
    {
      const auto ndof = g_inputdeck.get< tag::discr, tag::ndof >();

      Assert( U.nunk() == R.nunk(), "Number of unknowns in solution "
              "vector and right-hand side at recent time step incorrect" );
      Assert( U.nprop() == ndof*5 && R.nprop() == ndof*5,
              "Number of components in solution and right-hand side vector "
              "must equal "+ std::to_string(ndof*5) );
      Assert( inpoel.size()/4 == U.nunk(), "Connectivity inpoel has incorrect "
              "size" );
      Assert( fd.Inpofa().size()/3 == fd.Esuf().size()/2,
              "Mismatch in inpofa size" );

      // set rhs to zero
      R.fill(0.0);

      // configure Riemann flux function
      auto rieflxfn =
       [this]( const std::array< tk::real, 3 >& fn,
               const std::array< std::vector< tk::real >, 2 >& u,
               const std::vector< std::array< tk::real, 3 > >& v )
             { return m_riemann.flux( fn, u, v ); };
      // configure a no-op lambda for prescribed velocity
      auto velfn = [this]( ncomp_t, ncomp_t, tk::real, tk::real, tk::real ){
        return std::vector< std::array< tk::real, 3 > >( this->m_ncomp ); };

      // supported boundary condition types and associated state functions
      std::vector< std::pair< std::vector< bcconf_t >, tk::StateFn > > bctypes{{
        { m_bcdir, Dirichlet },
        { m_bcsym, Symmetry },
        { m_bcextrapolate, Extrapolate } }};

      // Number of integration points
      std::size_t NGfa;
      
      switch(ndof)
      {
        case 1:           //DG(P0)
          NGfa = 1;
          break;

        case 4:           //DG(P1)
          NGfa = 3;
          break;

        case 10:          //DG(P2)
          NGfa = 6;
          break;

        default:
          Throw( "dg::Compflow::rhs() not defined for NDOF=" +
               std::to_string(ndof) );
      }

      // arrays for quadrature points for face integration
      std::vector< std::vector< tk::real > > coordgpfa;
      std::vector< tk::real > wgpfa;
      coordgpfa.resize( NGfa, std::vector< tk::real >(2) );
      wgpfa.resize(NGfa);
      
      // get quadrature point weights and coordinates for triangle
      tk::GaussQuadratureTri( NGfa, coordgpfa, wgpfa );

      // compute internal surface flux integrals      
      for (std::size_t igp=0; igp<NGfa; ++igp)
      {
        tk::surfInt( m_system, m_ncomp, m_offset, coordgpfa[igp], wgpfa[igp],
                     inpoel, coord, fd, geoFace, rieflxfn, velfn, U, limFunc, R );
      }
      
      switch(ndof)
      {
        case 1:       // DG(P0)
          // compute source term intehrals
          tk::srcIntP0( m_system, m_ncomp, m_offset,
                        t, geoElem, Problem::src, R );
          // compute boundary surface flux integrals
          for (const auto& b : bctypes)
            tk::sidesetIntP0( m_system, m_ncomp, m_offset, b.first, fd,
                              geoFace, t, rieflxfn, velfn, b.second, U, R );
          break;

        case 4:       // DG(P1)
          // compute source term intehrals
          tk::srcIntP1( m_system, m_ncomp, m_offset,
                        t, inpoel, coord, geoElem, Problem::src, R );
          // compute volume integrals
          tk::volIntP1( m_system, m_ncomp, m_offset, inpoel, coord, geoElem, flux,
                        velfn, U, limFunc, R );
          // compute boundary surface flux integrals
          for (const auto& b : bctypes)
            tk::sidesetIntP1( m_system, m_ncomp, m_offset, b.first, fd, geoFace, inpoel, 
                              coord, t, rieflxfn, velfn, b.second, U, limFunc, R );
          break;

        case 10:      // DG(P2)
          // compute source term intehrals
          tk::srcIntP2( m_system, m_ncomp, m_offset, 
                        t, inpoel, coord, geoElem, Problem::src, R );
          // compute volume integrals
          tk::volIntP2( m_system, m_ncomp, m_offset, inpoel, coord, geoElem, flux,
                        velfn, U, R );
          // compute boundary surface flux integrals
          for (const auto& b : bctypes)
            tk::sidesetIntP2( m_system, m_ncomp, m_offset, b.first, fd, geoFace, inpoel,
                              coord, t, rieflxfn, velfn, b.second, U, R );
          break;

        default:
          Throw( "dg::Compflow::rhs() not defined for NDOF=" +
               std::to_string(ndof) );
      }
    }

    //! Compute the minimum time step size
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] geoFace Face geometry array
    //! \param[in] geoElem Element geometry array
    //! \param[in] limFunc Limiter function for higher-order solution dofs
    //! \param[in] U Solution vector at recent time step
    //! \return Minimum time step size
    tk::real dt( const std::array< std::vector< tk::real >, 3 >& coord,
                 const std::vector< std::size_t >& inpoel,
                 const inciter::FaceData& fd,
                 const tk::Fields& geoFace,
                 const tk::Fields& geoElem,
                 const tk::Fields& limFunc,
                 const tk::Fields& U ) const
    {
      const auto ndof = g_inputdeck.get< tag::discr, tag::ndof >();
      const tk::real g = g_inputdeck.get< tag::param, eq, tag::gamma >()[0];

      const auto& esuf = fd.Esuf();
      const auto& inpofa = fd.Inpofa();

      // Number of integration points
      constexpr std::size_t NG = 3;

      // arrays for quadrature points
      std::vector< std::vector< tk::real > > coordgp;
      std::vector< tk::real > wgp;
      coordgp.resize( NG, std::vector< tk::real >(2) );
      wgp.resize(NG);
      
      tk::real rho, u, v, w, rhoE, p, a, vn, dSV_l, dSV_r;
      std::vector< tk::real > delt( U.nunk(), 0.0 );

      const auto& cx = coord[0];
      const auto& cy = coord[1];
      const auto& cz = coord[2];

      // get quadrature point weights and coordinates for triangle
      tk::GaussQuadratureTri( NG, coordgp, wgp );

      // compute internal surface maximum characteristic speed
      for (std::size_t f=0; f<esuf.size()/2; ++f)
      {
        std::size_t el = static_cast< std::size_t >(esuf[2*f]);
        auto er = esuf[2*f+1];

        // nodal coordinates of the left element
        std::array< tk::real, 3 >
          p1_l{{ cx[ inpoel[4*el] ],
                 cy[ inpoel[4*el] ],
                 cz[ inpoel[4*el] ] }},
          p2_l{{ cx[ inpoel[4*el+1] ],
                 cy[ inpoel[4*el+1] ],
                 cz[ inpoel[4*el+1] ] }},
          p3_l{{ cx[ inpoel[4*el+2] ],
                 cy[ inpoel[4*el+2] ],
                 cz[ inpoel[4*el+2] ] }},
          p4_l{{ cx[ inpoel[4*el+3] ],
                 cy[ inpoel[4*el+3] ],
                 cz[ inpoel[4*el+3] ] }};

        auto detT_l = tk::Jacobian( p1_l, p2_l, p3_l, p4_l );


        auto x1 = cx[ inpofa[3*f]   ];
        auto y1 = cy[ inpofa[3*f]   ];
        auto z1 = cz[ inpofa[3*f]   ];

        auto x2 = cx[ inpofa[3*f+1] ];
        auto y2 = cy[ inpofa[3*f+1] ];
        auto z2 = cz[ inpofa[3*f+1] ];

        auto x3 = cx[ inpofa[3*f+2] ];
        auto y3 = cy[ inpofa[3*f+2] ];
        auto z3 = cz[ inpofa[3*f+2] ];

        dSV_l = 0.0;
        dSV_r = 0.0;

        // Gaussian quadrature
        for (std::size_t igp=0; igp<NG; ++igp)
        {
          // Barycentric coordinates for the triangular face
          auto shp1 = 1.0 - coordgp[igp][0] - coordgp[igp][1];
          auto shp2 = coordgp[igp][0];
          auto shp3 = coordgp[igp][1];

          // transformation of the quadrature point from the 2D reference/master
          // element to physical domain, to obtain its physical (x,y,z)
          // coordinates.
          auto xgp = x1*shp1 + x2*shp2 + x3*shp3;
          auto ygp = y1*shp1 + y2*shp2 + y3*shp3;
          auto zgp = z1*shp1 + z2*shp2 + z3*shp3;

          tk::real detT_gp(0);

          // transformation of the physical coordinates of the quadrature point
          // to reference space for the left element to be able to compute
          // basis functions on the left element.
          detT_gp = tk::Jacobian( p1_l, {{ xgp, ygp, zgp }}, p3_l, p4_l );
          auto xi_l = detT_gp / detT_l;
          detT_gp = tk::Jacobian( p1_l, p2_l, {{ xgp, ygp, zgp }}, p4_l );
          auto eta_l = detT_gp / detT_l;
          detT_gp = tk::Jacobian( p1_l, p2_l, p3_l, {{ xgp, ygp, zgp }} );
          auto zeta_l = detT_gp / detT_l;

          // basis functions at igp for the left element
          auto B2l = 2.0 * xi_l + eta_l + zeta_l - 1.0;
          auto B3l = 3.0 * eta_l + zeta_l - 1.0;
          auto B4l = 4.0 * zeta_l - 1.0;

          auto wt = wgp[igp] * geoFace(f,0,0);

          std::array< std::vector< tk::real >, 2 > ugp;

          // left element
          for (ncomp_t c=0; c<5; ++c)
          {
            auto mark = c*ndof;
            auto lmark = c*(ndof-1);
            ugp[0].push_back(  U(el, mark,   m_offset)
                             + limFunc(el, lmark+0, 0) * U(el, mark+1, m_offset) * B2l
                             + limFunc(el, lmark+1, 0) * U(el, mark+2, m_offset) * B3l
                             + limFunc(el, lmark+2, 0) * U(el, mark+3, m_offset) * B4l );
          }

          rho = ugp[0][0];
          u = ugp[0][1]/rho;
          v = ugp[0][2]/rho;
          w = ugp[0][3]/rho;
          rhoE = ugp[0][4];
          p = (g-1.0)*(rhoE - rho*(u*u + v*v + w*w)/2.0);

          a = std::sqrt(g * p / rho);

          vn = u*geoFace(f,1,0) + v*geoFace(f,2,0) + w*geoFace(f,3,0);

          dSV_l = wt * (std::fabs(vn) + a);

          // right element
          if (er > -1) {

            // nodal coordinates of the right element
            std::size_t eR = static_cast< std::size_t >( er );
            std::array< tk::real, 3 >
              p1_r{{ cx[ inpoel[4*eR] ],
                     cy[ inpoel[4*eR] ],
                     cz[ inpoel[4*eR] ] }},
              p2_r{{ cx[ inpoel[4*eR+1] ],
                     cy[ inpoel[4*eR+1] ],
                     cz[ inpoel[4*eR+1] ] }},
              p3_r{{ cx[ inpoel[4*eR+2] ],
                     cy[ inpoel[4*eR+2] ],
                     cz[ inpoel[4*eR+2] ] }},
              p4_r{{ cx[ inpoel[4*eR+3] ],
                     cy[ inpoel[4*eR+3] ],
                     cz[ inpoel[4*eR+3] ] }};

            auto detT_r = tk::Jacobian( p1_r, p2_r, p3_r, p4_r );

            // transformation of the physical coordinates of the quadrature
            // point to reference space for the right element
            detT_gp = tk::Jacobian( p1_r, {{ xgp, ygp, zgp }}, p3_r, p4_r );
            auto xi_r = detT_gp / detT_r;
            detT_gp = tk::Jacobian( p1_r, p2_r, {{ xgp, ygp, zgp }}, p4_r );
            auto eta_r = detT_gp / detT_r;
            detT_gp = tk::Jacobian( p1_r, p2_r, p3_r, {{ xgp, ygp, zgp }} );
            auto zeta_r = detT_gp / detT_r;

            // basis functions at igp for the right element
            auto B2r = 2.0 * xi_r + eta_r + zeta_r - 1.0;
            auto B3r = 3.0 * eta_r + zeta_r - 1.0;
            auto B4r = 4.0 * zeta_r - 1.0;

            for (ncomp_t c=0; c<5; ++c)
            {
              auto mark = c*ndof;
              auto lmark = c*(ndof-1);
              ugp[1].push_back(  U(eR, mark,   m_offset)
                               + limFunc(eR, lmark+0, 0) * U(eR, mark+1, m_offset) * B2r
                               + limFunc(eR, lmark+1, 0) * U(eR, mark+2, m_offset) * B3r
                               + limFunc(eR, lmark+2, 0) * U(eR, mark+3, m_offset) * B4r );
            }

            rho = ugp[1][0];
            u = ugp[1][1]/rho;
            v = ugp[1][2]/rho;
            w = ugp[1][3]/rho;
            rhoE = ugp[1][4];
            p = (g-1.0)*(rhoE - rho*(u*u + v*v + w*w)/2.0);

            a = std::sqrt(g * p / rho);

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

    //! Extract the velocity field at cell nodes
    //! \param[in] U Solution vector at recent time step
    //! \param[in] N Element node indices
    //! \return Array of the four values of the velocity field
    std::array< std::array< tk::real, 4 >, 3 >
    velocity( const tk::Fields& U,
              const std::array< std::vector< tk::real >, 3 >&,
              const std::array< std::size_t, 4 >& N ) const
    {
      std::array< std::array< tk::real, 4 >, 3 > v;
      v[0] = U.extract( 1, m_offset, N );
      v[1] = U.extract( 2, m_offset, N );
      v[2] = U.extract( 3, m_offset, N );
      auto r = U.extract( 0, m_offset, N );
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
    fieldOutput( const tk::Fields& /*L*/,
                 const std::vector< std::size_t >& /*inpoel*/,
                 const tk::UnsMesh::Coords& /*coord*/,
                 tk::real t,
                 const tk::Fields& geoElem,
                 tk::Fields& U ) const
    {
      std::array< std::vector< tk::real >, 3 > coord;
      std::vector< tk::real > v;
      v        = geoElem.extract(0,0);
      coord[0] = geoElem.extract(1,0);
      coord[1] = geoElem.extract(2,0);
      coord[2] = geoElem.extract(3,0);

      return Problem::fieldOutput( m_system, m_ncomp, m_offset, t, 0.0, v,
                                   coord, U );
    }

    //! Return nodal field output going to file
    std::vector< std::vector< tk::real > >
    avgElemToNode( const std::vector< std::size_t >& inpoel,
                   const tk::UnsMesh::Coords& coord,
                   const tk::Fields& /*geoElem*/,
                   const tk::Fields& limFunc,
                   const tk::Fields& U ) const
    {
      const auto ndof = g_inputdeck.get< tag::discr, tag::ndof >();
      const tk::real g = g_inputdeck.get< tag::param, eq, tag::gamma >()[0];

      const auto& cx = coord[0];
      const auto& cy = coord[1];
      const auto& cz = coord[2];

      std::vector< tk::real > count(cx.size(), 0);
      std::vector< std::vector< tk::real > >
        out( 6, std::vector< tk::real >( cx.size(), 0.0 ) );

      for (std::size_t e=0; e<inpoel.size()/4 ; ++e)
      {
        // nodal coordinates of the left element
        std::array< std::array< tk::real, 3 >, 4 >
          pi{{ {{ cx[ inpoel[4*e] ],
                 cy[ inpoel[4*e] ],
                 cz[ inpoel[4*e] ] }},
              {{ cx[ inpoel[4*e+1] ],
                 cy[ inpoel[4*e+1] ],
                 cz[ inpoel[4*e+1] ] }},
              {{ cx[ inpoel[4*e+2] ],
                 cy[ inpoel[4*e+2] ],
                 cz[ inpoel[4*e+2] ] }},
              {{ cx[ inpoel[4*e+3] ],
                 cy[ inpoel[4*e+3] ],
                 cz[ inpoel[4*e+3] ] }} }};
        auto detT = tk::Jacobian( pi[0], pi[1], pi[2], pi[3] );

        for (std::size_t i=0; i<4; ++i)
        {
          tk::real detT_gp;
          // transformation of the physical coordinates of the quadrature point
          // to reference space for the left element to be able to compute
          // basis functions on the left element.
          detT_gp = tk::Jacobian( pi[0], pi[i], pi[2], pi[3] );
          auto xi = detT_gp / detT;
          detT_gp = tk::Jacobian( pi[0], pi[1], pi[i], pi[3] );
          auto eta = detT_gp / detT;
          detT_gp = tk::Jacobian( pi[0], pi[1], pi[2], pi[i] );
          auto zeta = detT_gp / detT;

          auto B2 = 2.0 * xi + eta + zeta - 1.0;
          auto B3 = 3.0 * eta + zeta - 1.0;
          auto B4 = 4.0 * zeta - 1.0;

          std::vector< tk::real > ugp(5,0);

          for (ncomp_t c=0; c<5; ++c)
          {
            if (ndof == 1) {
              ugp[c] =  U(e, c, m_offset);
            } else {
              auto mark = c*ndof;
              auto lmark = c*(ndof-1);
              ugp[c] =  U(e, mark,   m_offset)
                      + limFunc(e, lmark+0, 0) * U(e, mark+1, m_offset) * B2
                      + limFunc(e, lmark+1, 0) * U(e, mark+2, m_offset) * B3
                      + limFunc(e, lmark+2, 0) * U(e, mark+3, m_offset) * B4;
            }
          }

          auto u = ugp[1] / ugp[0];
          auto v = ugp[2] / ugp[0];
          auto w = ugp[3] / ugp[0];
          auto p = (g - 1) * (ugp[4] - 0.5 * ugp[0] * (u*u + v*v + w*w) );

          out[0][ inpoel[4*e+i] ] += ugp[0];
          out[1][ inpoel[4*e+i] ] += u;
          out[2][ inpoel[4*e+i] ] += v;
          out[3][ inpoel[4*e+i] ] += w;
          out[4][ inpoel[4*e+i] ] += ugp[4]/ugp[0];
          out[5][ inpoel[4*e+i] ] += p;
          count[ inpoel[4*e+i] ] += 1.0;
        }
      }

      // average
      for (std::size_t i=0; i<cx.size(); ++i)
        for (std::size_t c=0; c<6; ++c)
          out[c][i] /= count[i];

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
    //! Riemann solver
    RiemannSolver m_riemann;
    //! Dirichlet BC configuration
    const std::vector< bcconf_t > m_bcdir;
    //! Symmetric BC configuration
    const std::vector< bcconf_t > m_bcsym;
    //! Extrapolation BC configuration
    const std::vector< bcconf_t > m_bcextrapolate;

    //! Evaluate physical flux function for this PDE system
    //! \param[in] system Equation system index
    //! \param[in] ncomp Number of scalar components in this PDE system
    //! \param[in] ugp Numerical solution at the Gauss point at which to
    //!   evaluate the flux
    //! \return Flux vectors for all components in this PDE system
    //! \note The function signature must follow tk::FluxFn
    static tk::FluxFn::result_type
    flux( ncomp_t system,
          ncomp_t ncomp,
          const std::vector< tk::real >& ugp,
          const std::vector< std::array< tk::real, 3 > >& )
    {
      Assert( ugp.size() == ncomp, "Size mismatch" );
      IGNORE(ncomp);

      const auto g = g_inputdeck.get< tag::param, eq, tag::gamma >()[ system ];

      auto u = ugp[1] / ugp[0];
      auto v = ugp[2] / ugp[0];
      auto w = ugp[3] / ugp[0];
      auto p = (g - 1) * (ugp[4] - 0.5 * ugp[0] * (u*u + v*v + w*w) );

      std::vector< std::array< tk::real, 3 > > fl( ugp.size() );

      fl[0][0] = ugp[1];
      fl[1][0] = ugp[1] * u + p;
      fl[2][0] = ugp[1] * v;
      fl[3][0] = ugp[1] * w;
      fl[4][0] = u * (ugp[4] + p);

      fl[0][1] = ugp[2];
      fl[1][1] = ugp[2] * u;
      fl[2][1] = ugp[2] * v + p;
      fl[3][1] = ugp[2] * w;
      fl[4][1] = v * (ugp[4] + p);

      fl[0][2] = ugp[3];
      fl[1][2] = ugp[3] * u;
      fl[2][2] = ugp[3] * v;
      fl[3][2] = ugp[3] * w + p;
      fl[4][2] = w * (ugp[4] + p);

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
    //! \note The function signature must follow tk::StateFn
    static tk::StateFn::result_type
    Dirichlet( ncomp_t system, ncomp_t ncomp, const std::vector< tk::real >& ul,
               tk::real x, tk::real y, tk::real z, tk::real t,
               const std::array< tk::real, 3 >& )
    {
      return {{ ul, Problem::solution( system, ncomp, x, y, z, t ) }};
    }

    //! \brief Boundary state function providing the left and right state of a
    //!   face at symmetry boundaries
    //! \param[in] ul Left (domain-internal) state
    //! \param[in] fn Unit face normal
    //! \return Left and right states for all scalar components in this PDE
    //!   system
    //! \note The function signature must follow tk::StateFn
    static tk::StateFn::result_type
    Symmetry( ncomp_t, ncomp_t, const std::vector< tk::real >& ul,
              tk::real, tk::real, tk::real, tk::real,
              const std::array< tk::real, 3 >& fn )
    {
      std::vector< tk::real > ur(5);
      // Internal cell velocity components
      auto v1l = ul[1]/ul[0];
      auto v2l = ul[2]/ul[0];
      auto v3l = ul[3]/ul[0];
      // Normal component of velocity
      auto vnl = v1l*fn[0] + v2l*fn[1] + v3l*fn[2];
      // Ghost state velocity components
      auto v1r = v1l - 2.0*vnl*fn[0];
      auto v2r = v2l - 2.0*vnl*fn[1];
      auto v3r = v3l - 2.0*vnl*fn[2];
      // Boundary condition
      ur[0] = ul[0];
      ur[1] = ur[0] * v1r;
      ur[2] = ur[0] * v2r;
      ur[3] = ur[0] * v3r;
      ur[4] = ul[4];
      return {{ std::move(ul), std::move(ur) }};
    }

    //! \brief Boundary state function providing the left and right state of a
    //!   face at extrapolation boundaries
    //! \param[in] ul Left (domain-internal) state
    //! \return Left and right states for all scalar components in this PDE
    //!   system
    //! \note The function signature must follow tk::StateFn
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

#endif // DGCompFlow_h
