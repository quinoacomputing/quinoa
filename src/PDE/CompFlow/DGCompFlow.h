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
#include <unordered_map>

#include "Macro.h"
#include "Exception.h"
#include "Vector.h"

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

    //! Extract BC configuration ignoring if BC not specified
    //! \param[in] c Equation system index (among multiple systems configured)
    //! \return Vector of BC config of type bcconf_t used to apply BCs for all
    //!   scalar components this Transport eq system is configred for
    //! \note A more preferable way of catching errors such as this function
    //!   hides is during parsing, so that we don't even get here if BCs are not
    //!   correctly specified. For now we simply ignore if BCs are not
    //!   specified by allowing empty BC vectors from the user input.
    template< typename bctag >
    std::vector< bcconf_t >
    config( ncomp_t c ) {
      std::vector< bcconf_t > bc;
      const auto& v = g_inputdeck.get< tag::param, tag::compflow, bctag >();
      ErrChk( !v.empty(), "Boundary conditions not set in control file for DG " 
              "CompFlow" );
      if (v.size() > c) bc = v[c];
      return bc;
    }

  public:
    //! \brief Constructor
    explicit CompFlow( ncomp_t c ) :
      m_offset( 0 ),
      m_bcdir( config< tag::bcdir >( c ) )
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

      const auto& bface = fd.Bface();
      const auto& esuf = fd.Esuf();

      // set rhs to zero
      R.fill(0.0);

      // compute internal surface flux integrals
      for (auto f=fd.Nbfac(); f<esuf.size()/2; ++f)
      {
        std::size_t el = static_cast< std::size_t >(esuf[2*f]);
        std::size_t er = static_cast< std::size_t >(esuf[2*f+1]);
        auto farea = geoFace(f,0,0);

        // Fluxes
        auto flux = numericalFluxFunc( f, geoFace, {{U.extract(el), U.extract(er)}} );

        for (ncomp_t c=0; c<5; ++c) {
          R(el, c, m_offset) -= farea * flux[c];
          R(er, c, m_offset) += farea * flux[c];
        }
      }

      // compute boundary surface flux integrals
      for (const auto& s : m_bcdir) {    // for all dirbc sidesets
        auto bc = bface.find( std::stoi(s) );  // faces for dir bc side set
        if (bc != end(bface))
          surfInt< Dir >( bc->second, esuf, geoFace, U, R, t );
      }

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
    const ncomp_t m_offset;             //!< Offset PDE operates from
    //! Dirichlet BC configuration
    const std::vector< bcconf_t > m_bcdir;

    //! \brief State policy class providing the left and right state of a face
    //!   at symmetric boundaries
    struct Sym {
      static std::array< std::vector< tk::real >, 2 >
      LR( const tk::Fields& U, std::size_t e,
          tk::real /*xc*/, tk::real /*yc*/, tk::real /*zc*/,
          tk::real /*t*/ ) {
        return {{ U.extract( e ), U.extract( e ) }};
      }
    };

    //! \brief State policy class providing the left and right state of a face
    //!   at Dirichlet boundaries
    struct Dir {
      static std::array< std::vector< tk::real >, 2 >
      LR( const tk::Fields& U, std::size_t e,
          tk::real xc, tk::real yc, tk::real zc,
          tk::real t ) {
        auto ul = U.extract( e );
        auto ur = ul;
        const auto urbc = Problem::solution(0, xc, yc, zc, t);
        for (ncomp_t c=0; c<5; ++c)
          ur[c] = urbc[c];
        return {{ std::move(ul), std::move(ur) }};
      }
    };

    //! Compute boundary surface integral
    //! \param[in] faces Face IDs at which to compute surface integral
    //! \param[in] esuf Elements surrounding face, see tk::genEsuf()
    //! \param[in] geoFace Face geometry array
    //! \param[in] U Solution vector at recent time step
    //! \param[in,out] R Right-hand side vector computed
    //! \param[in] t Physical time
    //! \tparam State Policy class providing the left and right state at
    //!   boundaries by its member function State::LR()
    template< class State >
    void surfInt( const std::vector< std::size_t >& faces,
                  const std::vector< int >& esuf,
                  const tk::Fields& geoFace,
                  const tk::Fields& U,
                  tk::Fields& R,
                  tk::real t ) const
    {
      for (const auto& f : faces) {
        std::size_t el = static_cast< std::size_t >(esuf[2*f]);
        Assert( esuf[2*f+1] == -1, "outside boundary element not -1" );
        auto farea = geoFace(f,0,0);
        auto xc = geoFace(f,4,0);
        auto yc = geoFace(f,5,0);
        auto zc = geoFace(f,6,0);

        //--- fluxes
        auto flux = numericalFluxFunc( f, geoFace, State::LR(U,el,xc,yc,zc,t) );

        for (ncomp_t c=0; c<5; ++c)
          R(el, c, m_offset) -= farea * flux[c];
      }
    }

    //! HLLC approximate Riemann solver
    //! \param[in] f Face ID
    //! \param[in] geoFace Face geometry array
    //! \param[in] u Left and right unknown/state vector
    //! \return Riemann solution using central difference method
    std::vector< tk::real >
    numericalFluxFunc( std::size_t f,
                     const tk::Fields& geoFace,
                     const std::array< std::vector< tk::real >, 2 >& u )
                   const
    {
      std::vector< tk::real > flux( u[0].size(), 0 );

      std::array< tk::real, 3 > fn {{ geoFace(f,1,0),
                                      geoFace(f,2,0),
                                      geoFace(f,3,0) }};

      // ratio of specific heats
      auto g = g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[0];

      // Primitive variables
      auto rhol = u[0][0];
      auto rhor = u[1][0];

      auto pl = (g-1.0)*(u[0][4] - (u[0][1]*u[0][1] +
                                    u[0][2]*u[0][2] +
                                    u[0][3]*u[0][3]) / (2.0*rhol));

      auto pr = (g-1.0)*(u[1][4] - (u[1][1]*u[1][1] +
                                    u[1][2]*u[1][2] +
                                    u[1][3]*u[1][3]) / (2.0*rhor));

      auto al = sqrt(g * pl / rhol);
      auto ar = sqrt(g * pr / rhor);

      // Face-normal velocities
      auto ul = u[0][1]/rhol;
      auto vl = u[0][2]/rhol;
      auto wl = u[0][3]/rhol;

      tk::real vnl = ul*fn[0] + vl*fn[1] + wl*fn[2];

      auto ur = u[1][1]/rhor;
      auto vr = u[1][2]/rhor;
      auto wr = u[1][3]/rhor;

      tk::real vnr = ur*fn[0] + vr*fn[1] + wr*fn[2];

      // Roe-averaged variables
      auto rlr = sqrt(rhor/rhol);
      auto rlr1 = 1.0 + rlr;

      auto vnroe = (vnr*rlr + vnl)/rlr1 ;
      auto aroe = (ar*rlr + al)/rlr1 ;

      // Signal velocities
      auto Sl = fmin(vnl-al, vnroe-aroe);
      auto Sr = fmax(vnr+ar, vnroe+aroe);
      auto Sm = ( rhor*vnr*(Sr-vnr) - rhol*vnl*(Sl-vnl) + pl-pr )
               /( rhor*(Sr-vnr) - rhol*(Sl-vnl) );

      // Middle-zone (star) variables
      auto pStar = rhol*(vnl-Sl)*(vnl-Sm) + pl;
      auto uStar = u;

      uStar[0][0] = (Sl-vnl) * rhol/ (Sl-Sm);
      uStar[0][1] = ((Sl-vnl) * u[0][1] + (pStar-pl)*fn[0]) / (Sl-Sm);
      uStar[0][2] = ((Sl-vnl) * u[0][2] + (pStar-pl)*fn[1]) / (Sl-Sm);
      uStar[0][3] = ((Sl-vnl) * u[0][3] + (pStar-pl)*fn[2]) / (Sl-Sm);
      uStar[0][4] = ((Sl-vnl) * u[0][4] - pl*vnl + pStar*Sm) / (Sl-Sm);

      uStar[1][0] = (Sr-vnr) * rhor/ (Sr-Sm);
      uStar[1][1] = ((Sr-vnr) * u[1][1] + (pStar-pr)*fn[0]) / (Sr-Sm);
      uStar[1][2] = ((Sr-vnr) * u[1][2] + (pStar-pr)*fn[1]) / (Sr-Sm);
      uStar[1][3] = ((Sr-vnr) * u[1][3] + (pStar-pr)*fn[2]) / (Sr-Sm);
      uStar[1][4] = ((Sr-vnr) * u[1][4] - pr*vnr + pStar*Sm) / (Sr-Sm);

      // Numerical fluxes
      if (Sl > 0.0) {
        flux[0] = u[0][0] * vnl;
        flux[1] = u[0][1] * vnl + pl*fn[0];
        flux[2] = u[0][2] * vnl + pl*fn[1];
        flux[3] = u[0][3] * vnl + pl*fn[2];
        flux[4] = ( u[0][4] + pl ) * vnl;
      }
      else if (Sl <= 0.0 && Sm > 0.0) {
        flux[0] = uStar[0][0] * Sm;
        flux[1] = uStar[0][1] * Sm + pStar*fn[0];
        flux[2] = uStar[0][2] * Sm + pStar*fn[1];
        flux[3] = uStar[0][3] * Sm + pStar*fn[2];
        flux[4] = ( uStar[0][4] + pStar ) * Sm;
      }
      else if (Sm <= 0.0 && Sr >= 0.0) {
        flux[0] = uStar[1][0] * Sm;
        flux[1] = uStar[1][1] * Sm + pStar*fn[0];
        flux[2] = uStar[1][2] * Sm + pStar*fn[1];
        flux[3] = uStar[1][3] * Sm + pStar*fn[2];
        flux[4] = ( uStar[1][4] + pStar ) * Sm;
      }
      else {
        flux[0] = u[1][0] * vnr;
        flux[1] = u[1][1] * vnr + pr*fn[0];
        flux[2] = u[1][2] * vnr + pr*fn[1];
        flux[3] = u[1][3] * vnr + pr*fn[2];
        flux[4] = ( u[1][4] + pr ) * vnr;
      }

      return flux;


      //// Lax-Friedrichs flux function

      //std::vector< tk::real > flux( u[0].size(), 0 );
      //                        fluxl( u[0].size(), 0 ),
      //                        fluxr( u[0].size(), 0 );

      //// ratio of specific heats
      //auto g = g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[0];

      //// Primitive variables
      //auto rhol = u[0][0];
      //auto rhor = u[1][0];

      //auto pl = (g-1.0)*(u[0][4] - (u[0][1]*u[0][1] +
      //                              u[0][2]*u[0][2] +
      //                              u[0][3]*u[0][3]) / (2.0*rhol));

      //auto pr = (g-1.0)*(u[1][4] - (u[1][1]*u[1][1] +
      //                              u[1][2]*u[1][2] +
      //                              u[1][3]*u[1][3]) / (2.0*rhor));

      //auto al = sqrt(g * pl / rhol);
      //auto ar = sqrt(g * pr / rhor);

      //// Face-normal velocities
      //auto ul = u[0][1]/rhol;
      //auto vl = u[0][2]/rhol;
      //auto wl = u[0][3]/rhol;

      //tk::real vnl = ul*fn[0] + vl*fn[1] + wl*fn[2];

      //auto ur = u[1][1]/rhor;
      //auto vr = u[1][2]/rhor;
      //auto wr = u[1][3]/rhor;

      //tk::real vnr = ur*fn[0] + vr*fn[1] + wr*fn[2];

      //// Flux functions
      //fluxl[0] = u[0][0] * vnl;
      //fluxl[1] = u[0][1] * vnl + pl*fn[0];
      //fluxl[2] = u[0][2] * vnl + pl*fn[1];
      //fluxl[3] = u[0][3] * vnl + pl*fn[2];
      //fluxl[4] = ( u[0][4] + pl ) * vnl;

      //fluxr[0] = u[1][0] * vnr;
      //fluxr[1] = u[1][1] * vnr + pr*fn[0];
      //fluxr[2] = u[1][2] * vnr + pr*fn[1];
      //fluxr[3] = u[1][3] * vnr + pr*fn[2];
      //fluxr[4] = ( u[1][4] + pr ) * vnr;

      //auto lambda = fmax(al,ar) + fmax(fabs(vnl),fabs(vnr));
    
      //// Numerical flux function
      //for(ncomp_t c=0; c<5; ++c)
      //{
      //  flux[c] = 0.5 * ( fluxl[c] + fluxr[c]
      //                   - lambda * (u[1][c] - u[0][c]) );
      //}
    
      //return flux;
    }
};

} // dg::

} // inciter::

#endif // DGCompFlow_h
