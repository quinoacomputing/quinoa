// *****************************************************************************
/*!
  \file      src/PDE/TransportProblem.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     Problem configurations for transport equations
  \details   This file defines policy classes for a transport equation, defined
    in PDE/Transport.h.

    General requirements on transport equation problem policy classes:

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::ProblemType type() noexcept {
          return ctr::ProblemType::SHEAR_DIFF;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.
*/
// *****************************************************************************
#ifndef TransportProblem_h
#define TransportProblem_h

#include <vector>
#include <array>

#include <cmath>

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "Inciter/Options/Problem.h"

namespace inciter {

/*! Transport PDE problem: diffusion of a shear layer
    \details This class implements the analytical solutions for the test
    problem, adopted from Okubo Akira Karweit Michael J. , (1969),
    [DIFFUSION FROM A CONTINUOUS SOURCE IN A UNIFORM SHEAR
    FLOW](http://onlinelibrary.wiley.com/doi/10.4319/lo.1969.14.4.0514/abstract),
    Limnology and Oceanography, 14, doi: 10.4319/lo.1969.14.4.0514. In essence,
    this is a test problem for the advection-diffusion equation in 3D where the
    analytical solution is known in a closed form as the solution evolves in
    time. The initial solution is a Gaussian that is advected and diffused in
    time with an imposed constant-in-time velocity field that features
    advection and shear. Also, the diffusion coefficients can be different in
    the three coordinate directions. Note that t0 as well as all three
    components of the diffusion must be larger than zero at t=t0 to have a
    well-defined initial condition.
   
    In a nutshell, the equation solved is
    \f[
      \frac{\partial S}{\partial t} + \left(V_0 + \Omega_y y + \Omega_z z
      \right) \frac{\partial S}{\partial x} =
      A_x \frac{\partial^2S}{\partial x^2} +
      A_y \frac{\partial^2S}{\partial y^2} +
      A_z \frac{\partial^2S}{\partial z^2}
    \f]
    whose solution is given by
    \f[
      S(t,x,y,z,) = \frac{1}{8\pi^{3/2}(A_xA_yA_z)^{1/2}t^{3/2}
                             (1+\phi_3^2t^2)^{1/2}}
                    \exp\left[ -\frac{x-V_0t-(\Omega_yy+\Omega_zz)^2/2}
                                     {4A_xt(1+\phi_3^2t^2}
                               -\frac{y^2}{4A_yt}
                               -\frac{z^2}{4A_zt} \right]
    \f]
    where \f$ \phi_3^2 = (\Omega_y^2A_y/A_x + \Omega_z^2A_z/A_x)/12\f$.
    See also the paper.
*/
class TransportProblemShearDiff {

  private:
    //! Evaluate analytical solution at (x,y,z,t) for all components
    //! \param[in] e Equation system index, i.e., which transport equation
    //!   system we operate on among the systems of PDEs
    //! \param[in] x X coordinate where to evaluate the solution
    //! \param[in] y Y coordinate where to evaluate the solution
    //! \param[in] z Z coordinate where to evaluate the solution
    //! \param[in] t Time where to evaluate the solution
    //! \return Values of all components evaluated at (x,y,z,t)
    static std::vector< tk::real >
    solution( tk::ctr::ncomp_type e,
              tk::ctr::ncomp_type ncomp,
              tk::real x, tk::real y, tk::real z, tk::real t )
    {
      using tag::param; using tag::transport;
      const auto& u0 = g_inputdeck.get< param, transport, tag::u0 >()[e];
      const auto& d = g_inputdeck.get< param, transport, tag::diffusivity >()[e];
      const auto& l = g_inputdeck.get< param, transport, tag::lambda >()[e];
      std::vector< tk::real > r( ncomp );
      for (ncomp_t c=0; c<ncomp; ++c) {
        const auto li = 2*c;
        const auto di = 3*c;
        const auto phi3s = (l[li+0]*l[li+0]*d[di+1]/d[di+0] +
                            l[li+1]*l[li+1]*d[di+2]/d[di+0]) / 12.0;
        r[c] =
            1.0 / ( 8.0 * std::pow(M_PI,3.0/2.0) *
                    std::sqrt(d[di+0]*d[di+1]*d[di+2]) *
                    std::pow(t,3.0/2.0) * std::sqrt(1.0+phi3s*t*t) ) *
            exp( -std::pow( x - u0[c]*t -
                            0.5*(l[li+0]*y + l[li+1]*z)*t, 2.0 ) /
                  ( 4.0 * d[di+0] * t * (1.0 + phi3s*t*t) )
                 -y*y / ( 4.0 * d[di+1] * t )
                 -z*z / ( 4.0 * d[di+2] * t ) );
      }
      return r;
    }

    //! \brief Evaluate the increment from t to t+dt of the analytical solution
    //!   at (x,y,z) for all components
    //! \param[in] e Equation system index, i.e., which transporte equation
    //!   system we operate on among the systems of PDEs
    //! \param[in] x X coordinate where to evaluate the solution
    //! \param[in] y Y coordinate where to evaluate the solution
    //! \param[in] z Z coordinate where to evaluate the solution
    //! \param[in] t Time where to evaluate the solution increment starting from
    //! \param[in] dt Time increment at which evaluate the solution increment to
    //! \return Increment in values of all components evaluated at (x,y,z,t+dt)
    static std::vector< tk::real >
    solinc( tk::ctr::ncomp_type e,
            tk::ctr::ncomp_type ncomp,
            tk::real x, tk::real y, tk::real z, tk::real t, tk::real dt )
    {
      auto st1 = solution( e, ncomp, x, y, z, t );
      auto st2 = solution( e, ncomp, x, y, z, t+dt );
      std::transform( begin(st1), end(st1), begin(st2), begin(st2),
                      []( tk::real s, tk::real& d ){ return d -= s; } );
      return st2;
    }

  public:
    //! Do error checking on PDE parameters
    //! \param[in] e Equation system index, i.e., which transport equation
    //!   system we operate on among the systems of PDEs
    //! \param[in] ncomp Number of components in this transport equation
    static void errchk( tk::ctr::ncomp_type e, tk::ctr::ncomp_type ncomp )
    {
      using tag::param; using tag::transport;
      const auto& u0 = g_inputdeck.get< param, transport, tag::u0 >()[e];
      ErrChk( ncomp == u0.size(),
        "Wrong number of advection-diffusion PDE parameters 'u0'" );
      const auto& lambda = g_inputdeck.get< param, transport, tag::lambda >()[e];
      ErrChk( 2*ncomp == lambda.size(),
        "Wrong number of advection-diffusion PDE parameters 'lambda'" );
      const auto& d = g_inputdeck.get< param, transport, tag::diffusivity >()[e];
      ErrChk( 3*ncomp == d.size(),
        "Wrong number of advection-diffusion PDE parameters 'diffusivity'" );
    }

    //! Set initial conditions for dispersion in simple shear flow
    //! \param[in] coord Mesh node coordinates
    //! \param[in,out] unk Array of unknowns
    //! \param[in] e Equation system index, i.e., which transport equation
    //!   system we operate on among the systems of PDEs
    //! \param[in] ncomp Number of components in this transport equation
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    //! \param[in] t Physical time
    static void init( const std::array< std::vector< tk::real >, 3 >& coord,
                      tk::Fields& unk,
                      tk::ctr::ncomp_type e,
                      tk::ctr::ncomp_type ncomp,
                      tk::ctr::ncomp_type offset,
                      tk::real t )
    {
      Assert( coord[0].size() == unk.nunk(), "Size mismatch" );
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];
      for (ncomp_t i=0; i<x.size(); ++i) {
        const auto s = solution( e, ncomp, x[i], y[i], z[i], t );
        for (ncomp_t c=0; c<ncomp; ++c)
          unk(i,c,offset) = s[c];
      }
    }

    //! \brief Query all side set IDs the user has configured for all components
    //!   in this PDE system
    //! \param[in,out] conf Set of unique side set IDs to add to
    static void side( std::unordered_set< int >& conf ) {
      using tag::param; using tag::transport; using tag::bcdir;
      for (const auto& s : g_inputdeck.get< param, transport, bcdir >())
        for (const auto& i : s)
          conf.insert( std::stoi(i) );
    }

    //! Assign prescribed shear velocity to nodes of tetrahedron element
    //! \param[in] y y coordinate at which to assign velocity
    //! \param[in] z Z coordinate at which to assign velocity
    //! \param[in] e Equation system index, i.e., which transport equation
    //!   system we operate on among the systems of PDEs
    //! \param[in] ncomp Number of components in this transport equation
    //! \return Velocity assigned to all vertices of a tetrehedron, size:
    //!   ncomp * ndim = [ncomp][3]
    static std::vector< std::array< tk::real, 3 > >
    prescribedVelocity( tk::real,
                        tk::real y,
                        tk::real z,
                        tk::ctr::ncomp_type e,
                        tk::ctr::ncomp_type ncomp )
    {
      using tag::param; using tag::transport;
      const auto& u0 = g_inputdeck.get< param, transport, tag::u0 >()[e];
      const auto& l = g_inputdeck.get< param, transport, tag::lambda >()[e];
      std::vector< std::array< tk::real, 3 > > vel( ncomp );
      for (ncomp_t c=0; c<ncomp; ++c)
        vel[c] = {{ u0[c] + l[2*c+0]*y + l[2*c+1]*z, 0.0, 0.0 }};
      return vel;
    }

    //! Return the velocity field at cell nodes
    //! \return Array of the four values of the three velocity coordinates
    static std::array< std::array< tk::real, 4 >, 3 >
    velocity( const tk::Fields&,
              const std::array< std::vector< tk::real >, 3 >&,
              const std::array< std::size_t, 4 >& )
    { return {{ {{0.0, 0.0, 0.0, 0.0}},
                {{0.0, 0.0, 0.0, 0.0}},
                {{0.0, 0.0, 0.0, 0.0}} }}; }

    //! \brief Query Dirichlet boundary condition value on a given side set for
    //!    all components in this PDE system
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] ncomp Number of components in this transport equation
    //! \param[in] t Physical time
    //! \param[in] deltat Time step size
    //! \param[in] side Pair of side set ID and node IDs on the side set
    //! \param[in] coord Mesh node coordinates
    //! \return Vector of pairs of bool and boundary condition value associated
    //!   to mesh node IDs at which Dirichlet boundary conditions are set. Note
    //!   that instead of the actual boundary condition value, we return the
    //!   increment between t+dt and t, since that is what the solution requires
    //!   as we solve for the soution increments and not the solution itself.
    static std::unordered_map< std::size_t,
                               std::vector< std::pair< bool, tk::real > > >
    dirbc( tk::ctr::ncomp_type e,
           tk::ctr::ncomp_type ncomp,
           tk::real t,
           tk::real deltat,
           const std::pair< const int, std::vector< std::size_t > >& side,
           const std::array< std::vector< tk::real >, 3 >& coord )
    {
      using tag::param; using tag::transport; using tag::bcdir;
      using NodeBC = std::vector< std::pair< bool, tk::real > >;
      std::unordered_map< std::size_t, NodeBC > bc;
      const auto& ubc = g_inputdeck.get< param, transport, bcdir >();
      if (!ubc.empty()) {
        Assert( ubc.size() > e, "Indexing out of Dirichlet BC eq-vector" );
        const auto& x = coord[0];
        const auto& y = coord[1];
        const auto& z = coord[2];
        for (const auto& b : ubc[e])
          if (std::stoi(b) == side.first)
            for (auto n : side.second) {
              Assert( x.size() > n, "Indexing out of coordinate array" );
              auto s = solinc( e, ncomp, x[n], y[n], z[n], t, deltat );
              auto& nbc = bc[n] = NodeBC( ncomp );
              for (ncomp_t c=0; c<ncomp; ++c)
                nbc[c] = { true, s[c] };
            }
      }
      return bc;
    }

    //! Problem type enum accessor
    //! \return Problem type as enum
    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::SHEAR_DIFF; }
};

//! Transport PDE problem: rotation of Zalesak's slotted cylinder
//! \see Zalesak, S. (1978) Fully Multidimensional Flux-Corrected Transport
//!   Algorithms for Fluids. Journal of Computational Physics 31, 335-362.
//! \see Leveque R.J. (1996) High-Resolution Conservative Algorithms for
//!   Advection in Incompressible Flow. SIAM Journal on Numerical Analyis 33,
//!   627-665.
class TransportProblemSlotCyl {

  private:
    //! Evaluate analytical solution at (x,y,t) for all components
    //! \param[in] x X coordinate where to evaluate the solution
    //! \param[in] y Y coordinate where to evaluate the solution
    //! \param[in] t Time where to evaluate the solution
    //! \return Values of all components evaluated at (x,y,t)
    static std::vector< tk::real >
    solution( tk::ctr::ncomp_type ncomp,
              tk::real x, tk::real y, tk::real t )
    {
      using std::sin; using std::cos;
      std::vector< tk::real > s( ncomp, 0.0 );
      for (ncomp_t c=0; c<ncomp; ++c) {
        auto T = t + 2.0*M_PI/ncomp * c;
        const tk::real R0 = 0.15;
        // center of the cone
        tk::real x0 = 0.5;
        tk::real y0 = 0.25;
        tk::real r = std::sqrt((x0-0.5)*(x0-0.5) + (y0-0.5)*(y0-0.5));
        tk::real kx = 0.5 + r*sin( T );
        tk::real ky = 0.5 - r*cos( T );
        // center of the hump
        x0 = 0.25;
        y0 = 0.5;
        r = std::sqrt((x0-0.5)*(x0-0.5) + (y0-0.5)*(y0-0.5));
        tk::real hx = 0.5 + r*sin( T-M_PI/2.0 ),
                 hy = 0.5 - r*cos( T-M_PI/2.0 );
        // center of the slotted cylinder
        x0 = 0.5;
        y0 = 0.75;
        r = std::sqrt((x0-0.5)*(x0-0.5) + (y0-0.5)*(y0-0.5));
        tk::real cx = 0.5 + r*sin( T+M_PI ),
                 cy = 0.5 - r*cos( T+M_PI );
        // end points of the cylinder slot
        tk::real i1x = 0.525, i1y = cy - r*cos( std::asin(0.025/r) ),
                 i2x = 0.525, i2y = 0.8,
                 i3x = 0.475, i3y = 0.8;
        // rotate end points of cylinder slot
        tk::real ri1x = 0.5 + cos(T)*(i1x-0.5) - sin(T)*(i1y-0.5),
                 ri1y = 0.5 + sin(T)*(i1x-0.5) + cos(T)*(i1y-0.5),
                 ri2x = 0.5 + cos(T)*(i2x-0.5) - sin(T)*(i2y-0.5),
                 ri2y = 0.5 + sin(T)*(i2x-0.5) + cos(T)*(i2y-0.5),
                 ri3x = 0.5 + cos(T)*(i3x-0.5) - sin(T)*(i3y-0.5),
                 ri3y = 0.5 + sin(T)*(i3x-0.5) + cos(T)*(i3y-0.5);
        // direction of slot sides
        tk::real v1x = ri2x-ri1x, v1y = ri2y-ri1y,
                 v2x = ri3x-ri2x, v2y = ri3y-ri2y;
        // lengths of direction of slot sides vectors
        tk::real v1 = std::sqrt(v1x*v1x + v1y*v1y),
                 v2 = std::sqrt(v2x*v2x + v2y*v2y);
        // cone
        r = std::sqrt((x-kx)*(x-kx) + (y-ky)*(y-ky)) / R0;
        if (r<1.0) s[c] = 0.6*(1.0-r);
        // hump
        r = std::sqrt((x-hx)*(x-hx) + (y-hy)*(y-hy)) / R0;
        if (r<1.0) s[c] = 0.2*(1.0+cos(M_PI*std::min(r,1.0)));
        // cylinder
        r = std::sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy)) / R0;
        const std::array< tk::real, 2 > r1{{ v1x, v1y }},
                                        r2{{ x-ri1x, y-ri1y }};
        const auto d1 = (r1[0]*r2[1] - r2[0]*r1[1]) / v1;
        const std::array< tk::real, 2 > r3{{ v2x, v2y }},
                                        r4{{ x-ri2x, y-ri2y }};
        const auto d2 = (r3[0]*r4[1] - r4[0]*r3[1]) / v2;
        if (r<1.0 && (d1>0.05 || d1<0.0 || d2<0.0)) s[c] = 0.6;
      }
      return s;
    }

    //! \brief Evaluate the increment from t to t+dt of the analytical solution
    //!   at (x,y,z) for all components
    //! \param[in] x X coordinate where to evaluate the solution
    //! \param[in] y Y coordinate where to evaluate the solution
    //! \param[in] t Time where to evaluate the solution increment starting from
    //! \param[in] dt Time increment at which evaluate the solution increment to
    //! \return Increment in values of all components evaluated at (x,y,t+dt)
    static std::vector< tk::real >
    solinc( tk::ctr::ncomp_type ncomp,
            tk::real x, tk::real y, tk::real t, tk::real dt )
    {
      auto st1 = solution( ncomp, x, y, t );
      auto st2 = solution( ncomp, x, y, t+dt );
      std::transform( begin(st1), end(st1), begin(st2), begin(st2),
                      []( tk::real s, tk::real& d ){ return d -= s; } );
      return st2;
    }

  public:
    //! Do error checking on PDE parameters
    static void errchk( tk::ctr::ncomp_type, tk::ctr::ncomp_type ) {}

    //! Set initial conditions for Zalesak's slotted cylinder case
    //! \param[in] coord Mesh node coordinates
    //! \param[in,out] unk Array of unknowns
    //! \param[in] ncomp Number of components in this transport equation
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    //! \param[in] t Physical time
    static void init( const std::array< std::vector< tk::real >, 3 >& coord,
                      tk::Fields& unk,
                      tk::ctr::ncomp_type,
                      tk::ctr::ncomp_type ncomp,
                      tk::ctr::ncomp_type offset,
                      tk::real t )
    {
      Assert( coord[0].size() == unk.nunk(), "Size mismatch" );
      const auto& x = coord[0];
      const auto& y = coord[1];
      for (ncomp_t i=0; i<x.size(); ++i) {
        const auto s = solution( ncomp, x[i], y[i], t );
        for (ncomp_t c=0; c<ncomp; ++c)
          unk(i,c,offset) = s[c];
      }
    }

    //! \brief Query all side set IDs the user has configured for all components
    //!   in this PDE system
    //! \param[in,out] conf Set of unique side set IDs to add to
    static void side( std::unordered_set< int >& conf ) {
      using tag::param; using tag::transport; using tag::bcdir;
      for (const auto& s : g_inputdeck.get< param, transport, bcdir >())
        for (const auto& i : s)
          conf.insert( std::stoi(i) );
    }

    //! \brief Query Dirichlet boundary condition value on a given side set for
    //!    all components in this PDE system
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] ncomp Number of components in this transport equation
    //! \param[in] t Physical time
    //! \param[in] deltat Time step size
    //! \param[in] side Pair of side set ID and node IDs on the side set
    //! \param[in] coord Mesh node coordinates
    //! \return Vector of pairs of bool and boundary condition value associated
    //!   to mesh node IDs at which Dirichlet boundary conditions are set. Note
    //!   that instead of the actual boundary condition value, we return the
    //!   increment between t+dt and t, since that is what the solution requires
    //!   as we solve for the soution increments and not the solution itself.
    static std::unordered_map< std::size_t,
                               std::vector< std::pair< bool, tk::real > > >
    dirbc( tk::ctr::ncomp_type e,
           tk::ctr::ncomp_type ncomp,
           tk::real t,
           tk::real deltat,
           const std::pair< const int, std::vector< std::size_t > >& side,
           const std::array< std::vector< tk::real >, 3 >& coord )
    {
      using tag::param; using tag::transport; using tag::bcdir;
      using NodeBC = std::vector< std::pair< bool, tk::real > >;
      std::unordered_map< std::size_t, NodeBC > bc;
      const auto& ubc = g_inputdeck.get< param, transport, bcdir >();
      if (!ubc.empty()) {
        Assert( ubc.size() > e, "Indexing out of Dirichlet BC eq-vector" );
        const auto& x = coord[0];
        const auto& y = coord[1];
        for (const auto& b : ubc[e])
          if (std::stoi(b) == side.first)
            for (auto n : side.second) {
              Assert( x.size() > n, "Indexing out of coordinate array" );
              auto s = solinc( ncomp, x[n], y[n], t, deltat );
              auto& nbc = bc[n] = NodeBC( ncomp );
              for (ncomp_t c=0; c<ncomp; ++c)
                nbc[c] = { true, s[c] };
            }
      }
      return bc;
    }

    //! Assign prescribed shear velocity to nodes of tetrahedron element
    //! \param[in] x X coordinate at which to assign velocity
    //! \param[in] y y coordinate at which to assign velocity
    //! \param[in] ncomp Number of components in this transport equation
    //! \return Velocity assigned to all vertices of a tetrehedron, size:
    //!   ncomp * ndim = [ncomp][3]
    static std::vector< std::array< tk::real, 3 > >
    prescribedVelocity( tk::real x,
                        tk::real y,
                        tk::real,
                        tk::ctr::ncomp_type,
                        tk::ctr::ncomp_type ncomp )
    {
      std::vector< std::array< tk::real, 3 > > vel( ncomp );
      for (ncomp_t c=0; c<ncomp; ++c)
        vel[c] = {{ 0.5-y, x-0.5, 0.0 }};
      return vel;
    }

    //! Return the velocity field at cell nodes
    //! \return Array of the four values of the three velocity coordinates
    static std::array< std::array< tk::real, 4 >, 3 >
    velocity( const tk::Fields&,
              const std::array< std::vector< tk::real >, 3 >& coord,
              const std::array< std::size_t, 4 >& N )
    {
      const auto& x = coord[0];
      const auto& y = coord[1];
      return {{ {{ 0.5-y[N[0]], 0.5-y[N[1]], 0.5-y[N[2]], 0.5-y[N[3]] }},
                {{ x[N[0]]-0.5, x[N[1]]-0.5, x[N[2]]-0.5, x[N[3]]-0.5 }},
                {{ 0.0,         0.0,         0.0,         0.0 }} }};
    }

    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::SLOT_CYL; }
};

//! List of all transport equation problem policies
using TransportProblems = boost::mpl::vector< TransportProblemShearDiff
                                            , TransportProblemSlotCyl
                                            >;

} // inciter::

#endif // TransportProblem_h
