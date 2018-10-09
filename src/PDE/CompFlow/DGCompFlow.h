// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/DGCompFlow.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Compressible single-material flow using discontinuous Galerkin
  \details   This file implements the physics operators governing compressible
    single-material flow using discontinuous Galerkin discretization.
*/
// *****************************************************************************
#ifndef DGCompFlow_h
#define DGCompFlow_h

#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <map>

#include <brigand/sequences/list.hpp>
#include <brigand/algorithms/for_each.hpp>

#include "Macro.h"
#include "Exception.h"
#include "Vector.h"
#include "ContainerUtil.h"
#include "RiemannSolver.h"
#include "Riemann/HLLC.h"
#include "Riemann/LaxFriedrichs.h"
#include "UnsMesh.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

namespace dg {

//! \brief CompFlow used polymorphically with tk::DGPDE
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Physics - physics configuration, see PDE/CompFlow/Physics.h
//!   - Problem - problem configuration, see PDE/CompFlow/Problems.h
//! \note The default physics is Euler, set in inciter::deck::check_compflow()
template< class Physics, class Problem >
class CompFlow {

  private:
    using ncomp_t = kw::ncomp::info::expect::type;
    using bcconf_t = kw::sideset::info::expect::type;

    //! Factory for Riemann solvers
    //! \details This factory is used to store the constructors as a
    //!   std::function of specific Riemann solvers that can be invoked at a
    //!   later time compared to the point where the map is populated. The key
    //!   is an enum, uniquely idenfitying a specific Riemann solver. The value
    //!   is std::function storing a constructor to be invoked. The type of
    //!   object stored in std::function is a generic (base) class constructor,
    //!   which provides a polymorphyic interface (overridable functions) that
    //!   specific (child) Riemann solvers override.
    using RiemannFactory =
      std::map< ctr::FluxType, std::function< RiemannSolver() > >;

    //! Register a Riemann solver into the Riemann solver factory
    struct registerRiemannSolver {
      //! Factory to which to register the Riemann solver
      RiemannFactory& factory;
      //! Constructor
      //! \param[in] f Factory
      registerRiemannSolver( RiemannFactory& f ) : factory( f ) {}
      //! \brief Function call operator templated on the type that implements
      //!   a specific Riemann solver
      template< typename U > void operator()( brigand::type_<U> ) {
         // Function object holding the (default) constructor to be called later
         // without bound arguments, since all specific Riemann solvers'
         // constructors are compiler-generated (default) constructors, and thus
         // taking no arguments.
         std::function< U() > c = boost::value_factory< U >();
         // Associate constructor function object to flux type in factory
         factory.emplace( U::type(),
           boost::bind(boost::value_factory< RiemannSolver >(), std::move(c)) );
      }
    };

    // Register all supported Riemann solvers into a factory
    //! \return Riemann solver factory
    RiemannFactory RiemannSolvers() {
      namespace mpl = boost::mpl;
      using RiemannSolverList = brigand::list< HLLC, LaxFriedrichs >;
      RiemannFactory r;
      brigand::for_each< RiemannSolverList >( registerRiemannSolver( r ) );
      return r;
    }

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
      const auto& v = g_inputdeck.get< tag::param, tag::compflow, bctag >();
      if (v.size() > c) bc = v[c];
      return bc;
    }

  public:
    //! Constructor
    explicit CompFlow( ncomp_t c ) :
      m_ncomp( g_inputdeck.get< tag::component >().get< tag::transport >().at(c) ),
      m_offset( g_inputdeck.get< tag::component >().offset< tag::compflow >(c) ),
      m_riemann( tk::cref_find( RiemannSolvers(),
                   g_inputdeck.get< tag::discr, tag::flux >() ) ),
      m_bcdir( config< tag::bcdir >( c ) ),
      m_bcsym( config< tag::bcsym >( c ) ),
      m_bcextrapolate( config< tag::bcextrapolate >( c ) ),
      m_ndof( 4 )
      //ErrChk( !m_bcdir.empty() || !m_bcsym.empty() || !m_bcextrapolate.empty(),
      //        "Boundary conditions not set in control file for DG CompFlow" );
    {
      //Problem::errchk( m_c, m_ncomp );
    }

    //! Initalize the compressible flow equations for DG
    //! \param[in] L Element mass matrix
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
      const auto ndof = g_inputdeck.get< tag::discr, tag::ndof >();
      if (ndof == 1)
        initializep0( L, inpoel, coord, unk, t );
      else if (ndof == 4)
        initializep1( L, inpoel, coord, unk, t );
      else Throw( "DGCompFlow::initialize() not defined" );
    }

    //! Compute the left hand side block-diagonal mass matrix
    //! \param[in] geoElem Element geometry array
    //! \param[in,out] l Block diagonal mass matrix
    void lhs( const tk::Fields& geoElem, tk::Fields& l ) const
    {
      Assert( geoElem.nunk() == l.nunk(), "Size mismatch" );
      std::size_t nelem = geoElem.nunk();

      for (std::size_t e=0; e<nelem; ++e)
      {
        for (ncomp_t c=0; c<5; ++c)
          l(e, c, m_offset) = geoElem(e,0,0);
      }
    }

    //! Compute the left hand side P1 block-diagonal mass matrix
    //! \param[in] geoElem Element geometry array
    //! \param[in,out] l Block diagonal mass matrix
    void lhsp1( const tk::Fields& geoElem, tk::Fields& l ) const
    {
      Assert( geoElem.nunk() == l.nunk(), "Size mismatch" );
      std::size_t nelem = geoElem.nunk();

      for (std::size_t e=0; e<nelem; ++e)
      {
        for (ncomp_t c=0; c<5; ++c)
        {
          auto mark = c*m_ndof;
          l(e, mark,   m_offset) = geoElem(e,0,0);
          l(e, mark+1, m_offset) = geoElem(e,0,0) / 10.0;
          l(e, mark+2, m_offset) = geoElem(e,0,0) * 3.0/10.0;
          l(e, mark+3, m_offset) = geoElem(e,0,0) * 3.0/5.0;
        }
      }
    }

    //! Compute right hand side
    //! \param[in] t Physical time
    //! \param[in] geoElem Element geometry array
    //! \param[in] geoFace Face geometry array
    //! \param[in] fd Face connectivity data object
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] U Solution vector at recent time step
    //! \param[in,out] R Right-hand side vector computed
    void rhs( tk::real t,
              const tk::Fields& geoFace,
              const tk::Fields& geoElem,
              const inciter::FaceData& fd,
              const std::vector< std::size_t >& inpoel,
              const tk::UnsMesh::Coords& coord,
              const tk::Fields& U,
              tk::Fields& R ) const
    {
      Assert( U.nunk() == R.nunk(), "Number of unknowns in solution "
              "vector and right-hand side at recent time step incorrect" );
      Assert( U.nunk() == geoElem.nunk(), "Number of unknowns in solution "
              "vector and element-geometry at recent time step incorrect" );
      Assert( U.nprop() == 5 && R.nprop() == 5,
              "Number of components in solution and right-hand side vector " 
              "must equal "+ std::to_string(5) );

      IGNORE(coord);
      IGNORE(inpoel);

      // set rhs to zero
      R.fill(0.0);

      const auto& esuf = fd.Esuf();

      // compute internal surface flux integrals
      for (auto f=fd.Nbfac(); f<esuf.size()/2; ++f)
      {
        std::size_t el = static_cast< std::size_t >(esuf[2*f]);
        std::size_t er = static_cast< std::size_t >(esuf[2*f+1]);
        auto farea = geoFace(f,0,0);

        auto flux =
          m_riemann.flux( f, geoFace, {{U.extract(el), U.extract(er)}} );

        for (ncomp_t c=0; c<5; ++c) {
          R(el, c, m_offset) -= farea * flux[c];
          R(er, c, m_offset) += farea * flux[c];
        }
      }

      // compute boundary surface flux integrals
      bndIntegral< Dir >( m_bcdir, fd, geoFace, t, U, R );
      bndIntegral< Sym >( m_bcsym, fd, geoFace, t, U, R );
      bndIntegral< Extrapolate >( m_bcextrapolate, fd, geoFace, t, U, R );

      // Add source term to right hand side
      for (std::size_t e=0; e<geoElem.nunk(); ++e) {
        auto vole = geoElem(e,0,0);
        auto xc = geoElem(e,1,0);
        auto yc = geoElem(e,2,0);
        auto zc = geoElem(e,3,0);
        auto s = Problem::src(0, xc, yc, zc, t);
        for (ncomp_t c=0; c<5; ++c)
          R(e, c, m_offset) += vole * s[c];
      }
    }

    //! Compute P1 right hand side
    //! \param[in] t Physical time
    //! \param[in] geoElem Element geometry array
    //! \param[in] geoFace Face geometry array
    //! \param[in] fd Face connectivity data object
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] U Solution vector at recent time step
    //! \param[in,out] R Right-hand side vector computed
    void rhsp1( tk::real t,
                const tk::Fields& geoFace,
                const tk::Fields& geoElem,
                const inciter::FaceData& fd,
                const std::vector< std::size_t >& inpoel,
                const tk::UnsMesh::Coords& coord,
                const tk::Fields& U,
                tk::Fields& R ) const
    {
      Assert( U.nunk() == R.nunk(), "Number of unknowns in solution "
              "vector and right-hand side at recent time step incorrect" );
      Assert( U.nunk() == geoElem.nunk(), "Number of unknowns in solution "
              "vector and element-geometry at recent time step incorrect" );
      Assert( U.nprop() == 5 && R.nprop() == 5,
              "Number of components in solution and right-hand side vector "
              "must equal "+ std::to_string(5) );

      const auto& bface = fd.Bface();
      const auto& esuf = fd.Esuf();
      const auto& inpofa = fd.Inpofa();

      // set rhs to zero
      R.fill(0.0);

      // arrays for quadrature points
      std::array< std::vector< tk::real >, 3 > coordgp;
      std::vector< tk::real > wgp;

      coordgp[0].resize(3,0);
      coordgp[1].resize(3,0);
      coordgp[2].resize(3,0);

      wgp.resize(3,0);

      // get quadrature point weights and coordinates for triangle
      GaussQuadrature( 2, coordgp, wgp );

      // compute internal surface flux integrals
      for (auto f=fd.Nbfac(); f<esuf.size()/2; ++f)
      {
        std::size_t el = static_cast< std::size_t >(esuf[2*f]);
        std::size_t er = static_cast< std::size_t >(esuf[2*f+1]);

        // nodal coordinates of the left element
        std::array< tk::real, 3 >
          p1_l{{ coord[0][ inpoel[4*el] ],
                 coord[1][ inpoel[4*el] ],
                 coord[2][ inpoel[4*el] ] }},
          p2_l{{ coord[0][ inpoel[4*el+1] ],
                 coord[1][ inpoel[4*el+1] ],
                 coord[2][ inpoel[4*el+1] ] }},
          p3_l{{ coord[0][ inpoel[4*el+2] ],
                 coord[1][ inpoel[4*el+2] ],
                 coord[2][ inpoel[4*el+2] ] }},
          p4_l{{ coord[0][ inpoel[4*el+3] ],
                 coord[1][ inpoel[4*el+3] ],
                 coord[2][ inpoel[4*el+3] ] }};

        // nodal coordinates of the right element
        std::array< tk::real, 3 >
          p1_r{{ coord[0][ inpoel[4*er] ],
                 coord[1][ inpoel[4*er] ],
                 coord[2][ inpoel[4*er] ] }},
          p2_r{{ coord[0][ inpoel[4*er+1] ],
                 coord[1][ inpoel[4*er+1] ],
                 coord[2][ inpoel[4*er+1] ] }},
          p3_r{{ coord[0][ inpoel[4*er+2] ],
                 coord[1][ inpoel[4*er+2] ],
                 coord[2][ inpoel[4*er+2] ] }},
          p4_r{{ coord[0][ inpoel[4*er+3] ],
                 coord[1][ inpoel[4*er+3] ],
                 coord[2][ inpoel[4*er+3] ] }};

        auto detT_l = getJacobian( p1_l, p2_l, p3_l, p4_l );
        auto detT_r = getJacobian( p1_r, p2_r, p3_r, p4_r );

        auto x1 = coord[0][ inpofa[3*f]   ];
        auto y1 = coord[1][ inpofa[3*f]   ];
        auto z1 = coord[2][ inpofa[3*f]   ];

        auto x2 = coord[0][ inpofa[3*f+1] ];
        auto y2 = coord[1][ inpofa[3*f+1] ];
        auto z2 = coord[2][ inpofa[3*f+1] ];

        auto x3 = coord[0][ inpofa[3*f+2] ];
        auto y3 = coord[1][ inpofa[3*f+2] ];
        auto z3 = coord[2][ inpofa[3*f+2] ];

        // Gaussian quadrature
        for (std::size_t igp=0; igp<3; ++igp)
        {
          // Barycentric coordinates for the triangular face
          auto shp1 = 1.0 - coordgp[0][igp] - coordgp[1][igp];
          auto shp2 = coordgp[0][igp];
          auto shp3 = coordgp[1][igp];

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
          detT_gp = getJacobian( p1_l, {{ xgp, ygp, zgp }}, p3_l, p4_l );
          auto xi_l = detT_gp / detT_l;
          detT_gp = getJacobian( p1_l, p2_l, {{ xgp, ygp, zgp }}, p4_l );
          auto eta_l = detT_gp / detT_l;
          detT_gp = getJacobian( p1_l, p2_l, p3_l, {{ xgp, ygp, zgp }} );
          auto zeta_l = detT_gp / detT_l;

          // basis functions at igp for the left element
          auto B2l = 2.0 * xi_l + eta_l + zeta_l - 1.0;
          auto B3l = 3.0 * eta_l + zeta_l - 1.0;
          auto B4l = 4.0 * zeta_l - 1.0;

          // transformation of the physical coordinates of the quadrature point
          // to reference space for the right element
          detT_gp = getJacobian( p1_r, {{ xgp, ygp, zgp }}, p3_r, p4_r );
          auto xi_r = detT_gp / detT_r;
          detT_gp = getJacobian( p1_r, p2_r, {{ xgp, ygp, zgp }}, p4_r );
          auto eta_r = detT_gp / detT_r;
          detT_gp = getJacobian( p1_r, p2_r, p3_r, {{ xgp, ygp, zgp }} );
          auto zeta_r = detT_gp / detT_r;
          
          // basis functions at igp for the right element
          auto B2r = 2.0 * xi_r + eta_r + zeta_r - 1.0;
          auto B3r = 3.0 * eta_r + zeta_r - 1.0;
          auto B4r = 4.0 * zeta_r - 1.0;

          auto wt = wgp[igp] * geoFace(f,0,0);

          std::array< std::vector< tk::real >, 2 > ugp;

          for (ncomp_t c=0; c<m_ncomp; ++c)
          {
            auto mark = c*m_ndof;
            ugp[0].push_back(  U(el, mark,   m_offset)
                             + U(el, mark+1, m_offset) * B2l
                             + U(el, mark+2, m_offset) * B3l
                             + U(el, mark+3, m_offset) * B4l );
            ugp[1].push_back(  U(er, mark,   m_offset)
                             + U(er, mark+1, m_offset) * B2r
                             + U(er, mark+2, m_offset) * B3r
                             + U(er, mark+3, m_offset) * B4r );
          }

          auto flux = m_riemann.flux( f, geoFace, {{ugp[0], ugp[1]}} );

          for (ncomp_t c=0; c<m_ncomp; ++c)
          {
            auto mark = c*m_ndof;

            R(el, mark,   m_offset) -= wt * flux[c];
            R(el, mark+1, m_offset) -= wt * flux[c] * B2l;
            R(el, mark+2, m_offset) -= wt * flux[c] * B3l;
            R(el, mark+3, m_offset) -= wt * flux[c] * B4l;

            R(er, mark,   m_offset) += wt * flux[c];
            R(er, mark+1, m_offset) += wt * flux[c] * B2r;
            R(er, mark+2, m_offset) += wt * flux[c] * B3r;
            R(er, mark+3, m_offset) += wt * flux[c] * B4r;
          }
        }
      }

      // compute boundary surface flux integrals
      bndIntegralp1< Dir >( m_bcdir, bface, esuf, geoFace, inpoel, inpofa, coord, t, U, R );
      bndIntegralp1< Sym >( m_bcsym, bface, esuf, geoFace, inpoel, inpofa, coord, t, U, R );
      bndIntegralp1< Extrapolate >( m_bcextrapolate, bface, esuf, geoFace, inpoel, inpofa, 
				    			    coord, t, U, R );

      // resize quadrature point arrays
      coordgp[0].resize(5,0);
      coordgp[1].resize(5,0);
      coordgp[2].resize(5,0);

      wgp.resize(5,0);

      // get quadrature point weights and coordinates for tetrahedron
      GaussQuadrature( 3, coordgp, wgp );

      std::array< std::array< tk::real, 3 >, 3 > jacInv;

      // compute volume integrals
      for (std::size_t e=0; e<U.nunk(); ++e)
      {
        auto x1 = coord[0][ inpoel[4*e]   ];
        auto y1 = coord[1][ inpoel[4*e]   ];
        auto z1 = coord[2][ inpoel[4*e]   ];

        auto x2 = coord[0][ inpoel[4*e+1] ];
        auto y2 = coord[1][ inpoel[4*e+1] ];
        auto z2 = coord[2][ inpoel[4*e+1] ];

        auto x3 = coord[0][ inpoel[4*e+2] ];
        auto y3 = coord[1][ inpoel[4*e+2] ];
        auto z3 = coord[2][ inpoel[4*e+2] ];

        auto x4 = coord[0][ inpoel[4*e+3] ];
        auto y4 = coord[1][ inpoel[4*e+3] ];
        auto z4 = coord[2][ inpoel[4*e+3] ];

        jacInv = getJacInverse( {{x1, y1, z1}},
                                {{x2, y2, z2}},
                                {{x3, y3, z3}},
                                {{x4, y4, z4}} );

        auto db2dxi1 = 2.0;
        auto db2dxi2 = 1.0;
        auto db2dxi3 = 1.0;

        auto db3dxi1 = 0.0;
        auto db3dxi2 = 3.0;
        auto db3dxi3 = 1.0;

        auto db4dxi1 = 0.0;
        auto db4dxi2 = 0.0;
        auto db4dxi3 = 4.0;

        auto db2dx =  db2dxi1 * jacInv[0][0]
                    + db2dxi2 * jacInv[1][0]
                    + db2dxi3 * jacInv[2][0];

        auto db2dy =  db2dxi1 * jacInv[0][1]
                    + db2dxi2 * jacInv[1][1]
                    + db2dxi3 * jacInv[2][1];

        auto db2dz =  db2dxi1 * jacInv[0][2]
                    + db2dxi2 * jacInv[1][2]
                    + db2dxi3 * jacInv[2][2];

        auto db3dx =  db3dxi1 * jacInv[0][0]
                    + db3dxi2 * jacInv[1][0]
                    + db3dxi3 * jacInv[2][0];

        auto db3dy =  db3dxi1 * jacInv[0][1]
                    + db3dxi2 * jacInv[1][1]
                    + db3dxi3 * jacInv[2][1];

        auto db3dz =  db3dxi1 * jacInv[0][2]
                    + db3dxi2 * jacInv[1][2]
                    + db3dxi3 * jacInv[2][2];

        auto db4dx =  db4dxi1 * jacInv[0][0]
                    + db4dxi2 * jacInv[1][0]
                    + db4dxi3 * jacInv[2][0];        
	
	    auto db4dy =  db4dxi1 * jacInv[0][1]
                    + db4dxi2 * jacInv[1][1]
                    + db4dxi3 * jacInv[2][1];

        auto db4dz =  db4dxi1 * jacInv[0][2]
                    + db4dxi2 * jacInv[1][2]
                    + db4dxi3 * jacInv[2][2];

        // Gaussian quadrature
        for (std::size_t igp=0; igp<5; ++igp)
        {
          tk::real u, v, w, p;
          std::vector< tk::real > ugp;
          std::array< std::vector< tk::real >, 3 > flux;

          tk::real g = g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[e];

          ugp.resize(m_ncomp,0);

          flux[0].resize(m_ncomp,0);
          flux[1].resize(m_ncomp,0);
          flux[2].resize(m_ncomp,0);

          auto B2 = 2.0 * coordgp[0][igp] + coordgp[1][igp] + coordgp[2][igp] - 1.0;
          auto B3 = 3.0 * coordgp[1][igp] + coordgp[2][igp] - 1.0;
          auto B4 = 4.0 * coordgp[2][igp] - 1.0;

          auto wt = wgp[igp] * geoElem(e, 0, 0);
	  	  
	  	  for (ncomp_t c=0; c<m_ncomp; ++c)
	  	  {
			auto mark = c*m_ndof;

      	    ugp[c] =  U(e, mark,   m_offset)
      	            + U(e, mark+1, m_offset) * B2
      	            + U(e, mark+2, m_offset) * B3
      	            + U(e, mark+3, m_offset) * B4;
	  	  }

	  	  u = ugp[1] / ugp[0];
	  	  v = ugp[2] / ugp[0];
	  	  w = ugp[3] / ugp[0];
	  	  p = (g - 1) * (ugp[4] - 0.5 * ugp[0] * (u*u + v*v + w*w) );

	  	  flux[0][0] = ugp[1];
	  	  flux[0][1] = ugp[1] * u + p;
	  	  flux[0][2] = ugp[1] * v;
	  	  flux[0][3] = ugp[1] * w;
	  	  flux[0][4] = u * (ugp[4] + p);

	  	  flux[1][0] = ugp[2];
      	  flux[1][1] = ugp[2] * u;
      	  flux[1][2] = ugp[2] * v + p;
      	  flux[1][3] = ugp[2] * w;
      	  flux[1][4] = v * (ugp[4] + p);

	  	  flux[2][0] = ugp[3];
      	  flux[2][1] = ugp[3] * u;
      	  flux[2][2] = ugp[3] * v;
      	  flux[2][3] = ugp[3] * w + p;
      	  flux[2][4] = w * (ugp[4] + p);

      	  for (ncomp_t c=0; c<m_ncomp; ++c)
      	  {
	  	    auto mark = c*m_ndof;
	  	    
	  	    R(e, mark+1, m_offset) += wt * (flux[0][c]*db2dx + flux[1][c]*db2dy + flux[2][c]*db2dz);
	  	    R(e, mark+2, m_offset) += wt * (flux[0][c]*db3dx + flux[1][c]*db3dy + flux[2][c]*db3dz);
	  	    R(e, mark+3, m_offset) += wt * (flux[0][c]*db4dx + flux[1][c]*db4dy + flux[2][c]*db4dz);
	  	  }
		}
      }
    }

    //! Compute the minimum time step size
//     //! \param[in] U Solution vector at recent time step
//     //! \param[in] coord Mesh node coordinates
//     //! \param[in] inpoel Mesh element connectivity
    //! \return Minimum time step size
    tk::real dt( const std::array< std::vector< tk::real >, 3 >& /*coord*/,
                 const std::vector< std::size_t >& /*inpoel*/,
                 const tk::Fields& /*U*/ ) const
    {
      tk::real mindt = std::numeric_limits< tk::real >::max();
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
    { return Problem::fieldNames(); }

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

      return Problem::fieldOutput( 0, m_offset, t, 0.0, v, coord, U );
    }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    std::vector< std::string > names() const
    { return Problem::names(); }

    //! Return analytic solution (if defined by Problem) at xi, yi, zi, t
    //! \param[in] xi X-coordinate
    //! \param[in] yi Y-coordinate
    //! \param[in] zi Z-coordinate
    //! \param[in] t Physical time
    //! \return Vector of analytic solution at given location and time
    std::vector< tk::real >
    analyticSolution( tk::real xi, tk::real yi, tk::real zi, tk::real t ) const
    {
      auto s = Problem::solution( 0, xi, yi, zi, t );
      return std::vector< tk::real >( begin(s), end(s) );
    }

  private:
    //!< Number of components in this PDE
    const ncomp_t m_ncomp;
    //! Offset PDE operates from
    const ncomp_t m_offset;
    //! Riemann solver
    RiemannSolver m_riemann;
    //! Dirichlet BC configuration
    const std::vector< bcconf_t > m_bcdir;
    //! Symmetric BC configuration
    const std::vector< bcconf_t > m_bcsym;
    //! Extrapolation BC configuration
    const std::vector< bcconf_t > m_bcextrapolate;
    const uint8_t m_ndof;

    //! Initalize the compressible flow equations using problem policy
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in,out] unk Array of unknowns
    //! \param[in] t Physical time
    void initializep0( const tk::Fields&,
                       const std::vector< std::size_t >& inpoel,
                       const tk::UnsMesh::Coords& coord,
                       tk::Fields& unk,
                       tk::real t ) const
    {
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      for (std::size_t e=0; e<unk.nunk(); ++e) {
        // node ids
        const auto A = inpoel[e*4+0];
        const auto B = inpoel[e*4+1];
        const auto C = inpoel[e*4+2];
        const auto D = inpoel[e*4+3];
        // compute centroid
        auto xcc = (x[A]+x[B]+x[C]+x[D])/4.0;
        auto ycc = (y[A]+y[B]+y[C]+y[D])/4.0;
        auto zcc = (z[A]+z[B]+z[C]+z[D])/4.0;
        // evaluate solution at centroid
        const auto s = Problem::solution( 0, xcc, ycc, zcc, t );
        // initialize unknown vector with solution at centroids
        unk(e, 0, m_offset) = s[0];
        unk(e, 1, m_offset) = s[1];
        unk(e, 2, m_offset) = s[2];
        unk(e, 3, m_offset) = s[3];
        unk(e, 4, m_offset) = s[4];
      }
    }

    //! Initalize the compressible flow equations, prepare for time integration
//    //! \param[in] L Element mass matrix
//    //! \param[in] inpoel Element-node connectivity
//    //! \param[in] coord Array of nodal coordinates
//    //! \param[in,out] unk Array of unknowns
//    //! \param[in] t Physical time
    void initializep1( const tk::Fields& /*L*/,
                       const std::vector< std::size_t >& /*inpoel*/,
                       const tk::UnsMesh::Coords& /*coord*/,
                       tk::Fields& /*unk*/,
                       tk::real /*t*/ ) const
    {
      //Assert( L.nunk() == unk.nunk(), "Size mismatch" );
      //std::size_t nelem = unk.nunk();
    }

    //! \brief State policy class providing the left and right state of a face
    //!   at Dirichlet boundaries
    struct Dir {
      static std::array< std::vector< tk::real >, 2 >
      LR( const std::vector< tk::real >& U,
          tk::real xc, tk::real yc, tk::real zc,
          std::array< tk::real, 3 > /*fn*/,
          tk::real t ) {
        auto ul = U;
        auto ur = ul;
        const auto urbc = Problem::solution(0, xc, yc, zc, t);
        for (ncomp_t c=0; c<5; ++c)
          ur[c] = urbc[c];
        return {{ std::move(ul), std::move(ur) }};
      }
    };

    //! \brief State policy class providing the left and right state of a face
    //!   at symmetric boundaries
    struct Sym {
      static std::array< std::vector< tk::real >, 2 >
      LR( const std::vector< tk::real >& U,
          tk::real /*xc*/, tk::real /*yc*/, tk::real /*zc*/,
          std::array< tk::real, 3 > fn,
          tk::real /*t*/ ) {
        auto ul = U;
        auto ur = ul;
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
    };

    //! \brief State policy class providing the left and right state of a face
    //!   at extrapolation boundaries
    struct Extrapolate {
      static std::array< std::vector< tk::real >, 2 >
      LR( const std::vector< tk::real >& U,
          tk::real /*xc*/, tk::real /*yc*/, tk::real /*zc*/,
          std::array< tk::real, 3 > /*fn*/,
          tk::real /*t*/ ) {
        return {{ U, U }};
      }
    };

    //! Compute boundary surface integral for a number of faces
    //! \param[in] faces Face IDs at which to compute surface integral
    //! \param[in] esuf Elements surrounding face, see tk::genEsuf()
    //! \param[in] geoFace Face geometry array
    //! \param[in] t Physical time
    //! \param[in] U Solution vector at recent time step
    //! \param[in,out] R Right-hand side vector computed
    //! \tparam State Policy class providing the left and right state at
    //!   boundaries by its member function State::LR()
    template< class State >
    void surfInt( const std::vector< std::size_t >& faces,
                  const std::vector< int >& esuf,
                  const tk::Fields& geoFace,
                  tk::real t,
                  const tk::Fields& U,
                  tk::Fields& R ) const
    {
      for (const auto& f : faces) {
        std::size_t el = static_cast< std::size_t >(esuf[2*f]);
        Assert( esuf[2*f+1] == -1, "outside boundary element not -1" );
        auto farea = geoFace(f,0,0);
        auto xc = geoFace(f,4,0);
        auto yc = geoFace(f,5,0);
        auto zc = geoFace(f,6,0);
        std::array< tk::real, 3 > fn {{ geoFace(f,1,0),
                                        geoFace(f,2,0),
                                        geoFace(f,3,0) }};

        //--- fluxes
        auto flux = m_riemann.flux( f, geoFace, State::LR(U.extract( el ),xc,yc,zc,fn,t) );

        for (ncomp_t c=0; c<5; ++c)
          R(el, c, m_offset) -= farea * flux[c];
      }
    }

    //! Compute boundary surface flux integrals for a given boundary type
    //! \tparam BCType Specifies the type of boundary condition to apply
    //! \param bcconfig BC configuration vector for multiple side sets
    //! \param[in] fd Face connectivity data object
    //! \param[in] geoFace Face geometry array
    //! \param[in] t Physical time
    //! \param[in] U Solution vector at recent time step
    //! \param[in,out] R Right-hand side vector computed
    template< class BCType >
    void
    bndIntegral( const std::vector< bcconf_t >& bcconfig,
                 const inciter::FaceData& fd,
                 const tk::Fields& geoFace,
                 tk::real t,
                 const tk::Fields& U,
                 tk::Fields& R ) const
    {
      const auto& bface = fd.Bface();
      const auto& esuf = fd.Esuf();

      for (const auto& s : bcconfig) {       // for all bc sidesets
        auto bc = bface.find( std::stoi(s) );// faces for side set
        if (bc != end(bface))
          surfInt< BCType >( bc->second, esuf, geoFace, t, U, R );
      }
    }

    //! Compute boundary surface integral for a number of faces
    //! \param[in] faces Face IDs at which to compute surface integral
    //! \param[in] esuf Elements surrounding face, see tk::genEsuf()
    //! \param[in] geoFace Face geometry array
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] inpofa Face-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] t Physical time
    //! \param[in] U Solution vector at recent time step
    //! \param[in,out] R Right-hand side vector computed
    //! \tparam State Policy class providing the left and right state at
    //!   boundaries by its member function State::LR()
    template< class State >
    void surfIntp1( const std::vector< std::size_t >& faces,
                    const std::vector< int >& esuf,
                    const tk::Fields& geoFace,
                    const std::vector< std::size_t >& inpoel,
                    const std::vector< std::size_t >& inpofa,
                    const tk::UnsMesh::Coords& coord,
                    tk::real t,
                    const tk::Fields& U,
                    tk::Fields& R ) const
      {
      // arrays for quadrature points
      std::array< std::vector< tk::real >, 3 > coordgp;
      std::vector< tk::real > wgp;

      coordgp[0].resize(3,0);
      coordgp[1].resize(3,0);
      coordgp[2].resize(3,0);

      wgp.resize(3,0);

      // get quadrature point weights and coordinates for triangle
      GaussQuadrature( 2, coordgp, wgp );

      for (const auto& f : faces) {
        std::size_t el = static_cast< std::size_t >(esuf[2*f]);
        Assert( esuf[2*f+1] == -1, "outside boundary element not -1" );
		//auto xc = geoFace(f,4,0);
        //auto yc = geoFace(f,5,0);
        //auto zc = geoFace(f,6,0); 
        std::array< tk::real, 3 > fn {{ geoFace(f,1,0),
                                        geoFace(f,2,0),
                                        geoFace(f,3,0) }};

        // nodal coordinates of the left element
        std::array< tk::real, 3 >
          p1_l{{ coord[0][ inpoel[4*el] ],
                 coord[1][ inpoel[4*el] ],
                 coord[2][ inpoel[4*el] ] }},
          p2_l{{ coord[0][ inpoel[4*el+1] ],
                 coord[1][ inpoel[4*el+1] ],
                 coord[2][ inpoel[4*el+1] ] }},
          p3_l{{ coord[0][ inpoel[4*el+2] ],
                 coord[1][ inpoel[4*el+2] ],
                 coord[2][ inpoel[4*el+2] ] }},
          p4_l{{ coord[0][ inpoel[4*el+3] ],
                 coord[1][ inpoel[4*el+3] ],
                 coord[2][ inpoel[4*el+3] ] }};

        auto detT_l = getJacobian( p1_l, p2_l, p3_l, p4_l );

        auto x1 = coord[0][ inpofa[3*f]   ];
        auto y1 = coord[1][ inpofa[3*f]   ];
        auto z1 = coord[2][ inpofa[3*f]   ];

        auto x2 = coord[0][ inpofa[3*f+1] ];
        auto y2 = coord[1][ inpofa[3*f+1] ];
        auto z2 = coord[2][ inpofa[3*f+1] ];

        auto x3 = coord[0][ inpofa[3*f+2] ];
        auto y3 = coord[1][ inpofa[3*f+2] ];
        auto z3 = coord[2][ inpofa[3*f+2] ];

        // Gaussian quadrature
        for (std::size_t igp=0; igp<3; ++igp)
        {
          // Barycentric coordinates for the triangular face
          auto shp1 = 1.0 - coordgp[0][igp] - coordgp[1][igp];
          auto shp2 = coordgp[0][igp];
          auto shp3 = coordgp[1][igp];

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
          detT_gp = getJacobian( p1_l, {{ xgp, ygp, zgp }}, p3_l, p4_l );
          auto xi_l = detT_gp / detT_l;
          detT_gp = getJacobian( p1_l, p2_l, {{ xgp, ygp, zgp }}, p4_l );
          auto eta_l = detT_gp / detT_l;
          detT_gp = getJacobian( p1_l, p2_l, p3_l, {{ xgp, ygp, zgp }} );
          auto zeta_l = detT_gp / detT_l;

          // basis functions at igp for the left element
          auto B2l = 2.0 * xi_l + eta_l + zeta_l - 1.0;
          auto B3l = 3.0 * eta_l + zeta_l - 1.0;
          auto B4l = 4.0 * zeta_l - 1.0;

          auto wt = wgp[igp] * geoFace(f,0,0);

	  	  std::vector< tk::real > ugp;

          for (ncomp_t c=0; c<m_ncomp; ++c)
          {
            auto mark = c*m_ndof;
            ugp.push_back (  U(el, mark,   m_offset)
                           + U(el, mark+1, m_offset) * B2l
                           + U(el, mark+2, m_offset) * B3l
                           + U(el, mark+3, m_offset) * B4l );
          }

          //--- fluxes
          auto flux = m_riemann.flux( f, geoFace, State::LR(ugp,xgp,ygp,zgp,fn,t) );

          for (ncomp_t c=0; c<m_ncomp; ++c)
          {
            auto mark = c*m_ndof;

            R(el, mark,   m_offset) -= wt * flux[c];
            R(el, mark+1, m_offset) -= wt * flux[c] * B2l;
            R(el, mark+2, m_offset) -= wt * flux[c] * B3l;
            R(el, mark+3, m_offset) -= wt * flux[c] * B4l;
          }
        }
      }
    }

    //! Compute boundary surface flux integrals for a given boundary type
    //! \tparam BCType Specifies the type of boundary condition to apply
    //! \param bcconfig BC configuration vector for multiple side sets
    //! \param[in] bface Boundary faces side-set information
    //! \param[in] esuf Elements surrounding faces
    //! \param[in] geoFace Face geometry array
    //! \param[in] inpoel Element-node connectivity
    //! \param[in] inpofa Face-node connectivity
    //! \param[in] coord Array of nodal coordinates
    //! \param[in] t Physical time
    //! \param[in] U Solution vector at recent time step
    //! \param[in,out] R Right-hand side vector computed
    template< class BCType >
    void
    bndIntegralp1( const std::vector< bcconf_t >& bcconfig,
                   const std::map< int, std::vector< std::size_t > >& bface,
                   const std::vector< int >& esuf,
                   const tk::Fields& geoFace,
                   const std::vector< std::size_t >& inpoel,
                   const std::vector< std::size_t >& inpofa,
                   const tk::UnsMesh::Coords& coord,
                   tk::real t,
                   const tk::Fields& U,
                   tk::Fields& R ) const
    {
      for (const auto& s : bcconfig) {       // for all bc sidesets
        auto bc = bface.find( std::stoi(s) );// faces for side set
        if (bc != end(bface))
          surfIntp1< BCType >( bc->second, esuf, geoFace, inpoel, inpofa,
                               coord, t, U, R );
      }
    }

    //! Gaussian quadrature points locations and weights
    //! \param[in] ndimn Dimension of integration domain
    //! \param[in,out] coordgp Coordinates of quadrature points
    //! \param[in,out] wgp Weights of quadrature points
    void
    GaussQuadrature( uint8_t ndimn,
                     std::array< std::vector< tk::real >, 3 >& coordgp,
                     std::vector< tk::real >& wgp ) const
    {
      if (ndimn == 3)
        {
          coordgp[0][0] = 0.25;
          coordgp[1][0] = 0.25;
          coordgp[2][0] = 0.25;
          wgp[0]        = -12.0/15.0;

          coordgp[0][1] = 1.0/6.0;
          coordgp[1][1] = 1.0/6.0;
          coordgp[2][1] = 1.0/6.0;
          wgp[1]        = 9.0/20.0;

          coordgp[0][2] = 0.5;
          coordgp[1][2] = 1.0/6.0;
          coordgp[2][2] = 1.0/6.0;
          wgp[2]        = 9.0/20.0;

          coordgp[0][3] = 1.0/6.0;
          coordgp[1][3] = 0.5;
          coordgp[2][3] = 1.0/6.0;
          wgp[3]        = 9.0/20.0;

          coordgp[0][4] = 1.0/6.0;
          coordgp[1][4] = 1.0/6.0;
          coordgp[2][4] = 0.5;
          wgp[4]        = 9.0/20.0;
        }
      else if (ndimn == 2)
        {
          coordgp[0][0] = 2.0/3.0;
          coordgp[1][0] = 1.0/6.0;
          wgp[0]        = 1.0/3.0;

          coordgp[0][1] = 1.0/6.0;
          coordgp[1][1] = 2.0/3.0;
          wgp[1]        = 1.0/3.0;

          coordgp[0][2] = 1.0/6.0;
          coordgp[1][2] = 1.0/6.0;
          wgp[2]        = 1.0/3.0;
        }
      else
        {
          std::cout << "Incorrect dimensionality input to GaussQuadrature\n";
        }
    }

    //! Determinant of Jacobian of transformation
    //! \param[in] p1 (x,y,z) coordinates of 1st local node in the tetrahedron
    //! \param[in] p2 (x,y,z) coordinates of 2nd local node in the tetrahedron
    //! \param[in] p3 (x,y,z) coordinates of 3rd local node in the tetrahedron
    //! \param[in] p4 (x,y,z) coordinates of 4th local node in the tetrahedron
    //! \return Determinant of the Jacobian of transformation of physical
    //!   tetrahedron to reference (xi, eta, zeta) space
    tk::real
    getJacobian( const std::array< tk::real, 3 >& p1,
                 const std::array< tk::real, 3 >& p2,
                 const std::array< tk::real, 3 >& p3,
                 const std::array< tk::real, 3 >& p4 ) const
    {
      tk::real detT;

      detT = (p2[0]-p1[0])
              * ((p3[1]-p1[1])*(p4[2]-p1[2]) - (p4[1]-p1[1])*(p3[2]-p1[2]))
            -(p3[0]-p1[0])
              * ((p2[1]-p1[1])*(p4[2]-p1[2]) - (p4[1]-p1[1])*(p2[2]-p1[2]))
            +(p4[0]-p1[0])
              * ((p2[1]-p1[1])*(p3[2]-p1[2]) - (p3[1]-p1[1])*(p2[2]-p1[2]));

      return detT;
    }

    //! Inverse of Jacobian of transformation
    //! \param[in] p1 (x,y,z) coordinates of 1st local node in the tetrahedron
    //! \param[in] p2 (x,y,z) coordinates of 2nd local node in the tetrahedron
    //! \param[in] p3 (x,y,z) coordinates of 3rd local node in the tetrahedron
    //! \param[in] p4 (x,y,z) coordinates of 4th local node in the tetrahedron
    //! \return Inverse of the Jacobian of transformation of physical
    //!   tetrahedron to reference (xi, eta, zeta) space
    std::array< std::array< tk::real, 3 >, 3 >
    getJacInverse( const std::array< tk::real, 3 >& p1,
                   const std::array< tk::real, 3 >& p2,
                   const std::array< tk::real, 3 >& p3,
                   const std::array< tk::real, 3 >& p4 ) const
    {
      std::array< std::array< tk::real, 3 >, 3 > jacInv;

      auto detJ = getJacobian( p1, p2, p3, p4 );

      jacInv[0][0] =  (  (p3[1]-p1[1])*(p4[2]-p1[2])
                       - (p4[1]-p1[1])*(p3[2]-p1[2])) / detJ;
      jacInv[1][0] = -(  (p2[1]-p1[1])*(p4[2]-p1[2])
                       - (p4[1]-p1[1])*(p2[2]-p1[2])) / detJ;
      jacInv[2][0] =  (  (p2[1]-p1[1])*(p3[2]-p1[2])
                       - (p3[1]-p1[1])*(p2[2]-p1[2])) / detJ;

      jacInv[0][1] = -(  (p3[0]-p1[0])*(p4[2]-p1[2])
                       - (p4[0]-p1[0])*(p3[2]-p1[2])) / detJ;
      jacInv[1][1] =  (  (p2[0]-p1[0])*(p4[2]-p1[2])
                       - (p4[0]-p1[0])*(p2[2]-p1[2])) / detJ;
      jacInv[2][1] = -(  (p2[0]-p1[0])*(p3[2]-p1[2])
                       - (p3[0]-p1[0])*(p2[2]-p1[2])) / detJ;

      jacInv[0][2] =  (  (p3[0]-p1[0])*(p4[1]-p1[1])
                       - (p4[0]-p1[0])*(p3[1]-p1[1])) / detJ;
      jacInv[1][2] = -(  (p2[0]-p1[0])*(p4[1]-p1[1])
                       - (p4[0]-p1[0])*(p2[1]-p1[1])) / detJ;
      jacInv[2][2] =  (  (p2[0]-p1[0])*(p3[1]-p1[1])
                       - (p3[0]-p1[0])*(p2[1]-p1[1])) / detJ;

      return jacInv;
    }
};

} // dg::

} // inciter::

#endif // DGCompFlow_h
