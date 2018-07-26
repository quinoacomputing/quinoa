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
      m_offset( g_inputdeck.get< tag::component >().offset< tag::compflow >(c) ),
      m_riemann( tk::cref_find( RiemannSolvers(),
                   g_inputdeck.get< tag::discr, tag::flux >() ) ),
      m_bcdir( config< tag::bcdir >( c ) ),
      m_bcsym( config< tag::bcsym >( c ) ),
      m_bcextrapolate( config< tag::bcextrapolate >( c ) ),
      m_ndof( 4 )
      //ErrChk( !m_bcdir.empty() || !m_bcsym.empty() || !m_bcextrapolate.empty(),
      //        "Boundary conditions not set in control file for DG CompFlow" );
    {}

    //! Initalize the compressible flow equations, prepare for time integration
    //! \param[in] geoElem Element geometry array
    //! \param[in,out] unk Array of unknowns
    //! \param[in] t Physical time
    void initialize( const tk::Fields& geoElem,
                     tk::Fields& unk,
                     tk::real t ) const
    {
      Assert( geoElem.nunk() == unk.nunk(), "Size mismatch" );
      std::size_t nelem = unk.nunk();

      for (std::size_t e=0; e<nelem; ++e)
      {
        auto xcc = geoElem(e,1,0);
        auto ycc = geoElem(e,2,0);
        auto zcc = geoElem(e,3,0);

        const auto s = Problem::solution( 0, xcc, ycc, zcc, t );
        unk(e, 0, m_offset) = s[0];
        unk(e, 1, m_offset) = s[1];
        unk(e, 2, m_offset) = s[2];
        unk(e, 3, m_offset) = s[3];
        unk(e, 4, m_offset) = s[4];
      }
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
    //! \param[in] U Solution vector at recent time step
    //! \param[in,out] R Right-hand side vector computed
    void rhs( tk::real t,
              const tk::Fields& geoFace,
              const tk::Fields& geoElem,
              const inciter::FaceData& fd,
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
    void rhsp1( tk::real /*t*/,
                const tk::Fields& /*geoFace*/,
                const tk::Fields& /*geoElem*/,
                const inciter::FaceData& /*fd*/,
                const std::vector< std::size_t >& /*inpoel*/,
                const tk::UnsMesh::Coords& /*coord*/,
                const tk::Fields& /*U*/,
                tk::Fields& /*R*/ ) const
    {
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
    //! \param[in] V Total mesh volume
    //! \param[in] geoElem Element geometry array
    //! \param[in,out] U Solution vector at recent time step
    //! \return Vector of vectors to be output to file
    std::vector< std::vector< tk::real > >
    fieldOutput( tk::real t,
                 tk::real V,
                 const tk::Fields& geoElem,
                 tk::Fields& U ) const
    {
      std::array< std::vector< tk::real >, 3 > coord;
      std::vector< tk::real > v;
      v        = geoElem.extract(0,0);
      coord[0] = geoElem.extract(1,0);
      coord[1] = geoElem.extract(2,0);
      coord[2] = geoElem.extract(3,0);

      return Problem::fieldOutput( 0, m_offset, t, V, v, coord, U );
    }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    std::vector< std::string > names() const
    { return Problem::names(); }

  private:
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

    //! \brief State policy class providing the left and right state of a face
    //!   at Dirichlet boundaries
    struct Dir {
      static std::array< std::vector< tk::real >, 2 >
      LR( const tk::Fields& U, std::size_t e,
          tk::real xc, tk::real yc, tk::real zc,
          std::array< tk::real, 3 > /*fn*/,
          tk::real t ) {
        auto ul = U.extract( e );
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
      LR( const tk::Fields& U, std::size_t e,
          tk::real /*xc*/, tk::real /*yc*/, tk::real /*zc*/,
          std::array< tk::real, 3 > fn,
          tk::real /*t*/ ) {
        auto ul = U.extract( e );
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
      LR( const tk::Fields& U, std::size_t e,
          tk::real /*xc*/, tk::real /*yc*/, tk::real /*zc*/,
          std::array< tk::real, 3 > /*fn*/,
          tk::real /*t*/ ) {
        return {{ U.extract( e ), U.extract( e ) }};
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
        auto flux = m_riemann.flux( f, geoFace, State::LR(U,el,xc,yc,zc,fn,t) );

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
};

} // dg::

} // inciter::

#endif // DGCompFlow_h
