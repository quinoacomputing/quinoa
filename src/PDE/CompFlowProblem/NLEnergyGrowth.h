// *****************************************************************************
/*!
  \file      src/PDE/CompFlowProblem/NLEnergyGrowth.h
  \author    F. Gonzalez
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Problem configurations for the compressible flow equations
  \details   This file defines a policy classe for the compressible flow
    equations, defined in PDE/CompFlow.h. See PDE/CompFlow.h for general
    requirements on flow equations problem policy classes.
*/
// *****************************************************************************
#ifndef CompFlowProblemNLEnergyGrowth_h
#define CompFlowProblemNLEnergyGrowth_h

#include <string>
#include <unordered_set>

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Inciter/Options/Problem.h"

namespace inciter {

//! CompFlow system of PDEs problem: nonlinear energy growth
//! \see Waltz, et. al, "Manufactured solutions for the three-dimensional Euler
//!   equations with relevance to Inertial Confinement Fusion", Journal of
//!   Computational Physics 267 (2014) 196-209.
class CompFlowProblemNLEnergyGrowth {

  private:
    //! Compute internal energy parameter
    //! \param[in] bx Parameter betax
    //! \param[in] by Parameter betay
    //! \param[in] bz Parameter betaz
    //! \param[in] x X coordinate to evaluate at
    //! \param[in] y Y coordinate to evaluate at
    //! \param[in] z Z coordinate to evaluate at
    //! \return Internal energy parameter
    static tk::real hx( tk::real bx, tk::real by, tk::real bz,
                        tk::real x, tk::real y, tk::real z )
    { return std::cos(bx*M_PI*x) * std::cos(by*M_PI*y) * std::cos(bz*M_PI*z); }

    //! Compute a power of the negative cube of the internal energy
    //! \param[in] ce Internal energy parameter
    //! \param[in] kappa Internal energy parameter
    //! \param[in] t Physical time
    //! \param[in] h Internal energy parameter
    //! \param[in] p Power
    //! \return Negative cube of the internal energy raised to power p
    static tk::real ec( tk::real ce, tk::real kappa, tk::real t, tk::real h,
                        tk::real p )
    { return std::pow( -3.0*(ce + kappa*h*h*t), p ); }

  public:

    //! Set initial conditions
    //! \param[in] coord Mesh node coordinates
    //! \param[in,out] unk Array of unknowns
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    //! \param[in] t Physical time
    static void init( const std::array< std::vector< tk::real >, 3 >& coord,
                      const std::vector< std::size_t >&,
                      const std::unordered_map< std::size_t,
                              std::vector< std::pair< bool, tk::real > > >&,
                      tk::Fields& unk,
                      tk::ctr::ncomp_type e,
                      tk::ctr::ncomp_type offset,
                      tk::real t )
    {
      // manufactured solution parameters
      const auto& ce =
        g_inputdeck.get< tag::param, tag::compflow, tag::ce >()[e];
      const auto& r0 =
        g_inputdeck.get< tag::param, tag::compflow, tag::r0 >()[e];
      const auto& alpha =
        g_inputdeck.get< tag::param, tag::compflow, tag::alpha >()[e];
      const auto& k =
        g_inputdeck.get< tag::param, tag::compflow, tag::kappa >()[e];
      const auto& bx =
        g_inputdeck.get< tag::param, tag::compflow, tag::betax >()[e];
      const auto& by =
        g_inputdeck.get< tag::param, tag::compflow, tag::betay >()[e];
      const auto& bz =
        g_inputdeck.get< tag::param, tag::compflow, tag::betaz >()[e];
      // temporal component of density field
      tk::real ft = std::exp( -alpha*t );
      // set initial and boundary conditions
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];
      for (ncomp_t i=0; i<x.size(); ++i) {
        auto& r  = unk(i,0,offset); // rho
        auto& ru = unk(i,1,offset); // rho * u
        auto& rv = unk(i,2,offset); // rho * v
        auto& rw = unk(i,3,offset); // rho * w
        auto& re = unk(i,4,offset); // rho * e, e:total=kinetic+internal energy
        // spatial component of density field
        tk::real gx = 1.0 - (x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
        // internal energy parameter
        const auto h = hx( bx, by, bz, x[i], y[i], z[i] );
        // solution at t
        r = r0 + ft*gx;
        ru = 0.0;
        rv = 0.0;
        rw = 0.0;
        re = r * ec(ce,k,t,h,-1.0/3.0);
      }
    }

    //! Add source term to rhs for nonlinear energy growth manufactured solution
    //! \param[in] t Physical time
    //! \param[in] coord Mesh node coordinates
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] dt Size of time step
    //! \param[in] N Element node indices
    //! \param[in] mass Element mass matrix, nnode*nnode [4][4]
    //! \param[in] r Pointers to right hand side at component and offset
    //! \param[in,out] R Right-hand side vector contributing to
    static void
    sourceRhs( tk::real t,
               const std::array< std::vector< tk::real >, 3 >& coord,
               tk::ctr::ncomp_type e,
               tk::real dt,
               const std::array< std::size_t, 4 >& N,
               const std::array< std::array< tk::real, 4 >, 4 >& mass,
               const std::array< const tk::real*, 5 >& r,
               tk::Fields& R )
    {
      // manufactured solution parameters
      const auto& a =
        g_inputdeck.get< tag::param, tag::compflow, tag::alpha >()[e];
      const auto& bx =
        g_inputdeck.get< tag::param, tag::compflow, tag::betax >()[e];
      const auto& by =
        g_inputdeck.get< tag::param, tag::compflow, tag::betay >()[e];
      const auto& bz =
        g_inputdeck.get< tag::param, tag::compflow, tag::betaz >()[e];
      const auto& ce =
        g_inputdeck.get< tag::param, tag::compflow, tag::ce >()[e];
      const auto& kappa =
        g_inputdeck.get< tag::param, tag::compflow, tag::kappa >()[e];
      const auto& r0 =
        g_inputdeck.get< tag::param, tag::compflow, tag::r0 >()[e];
      // ratio of specific heats
      tk::real g =
        g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[e];

      // mesh node coordinates
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      // spatial component of density field
      std::array< tk::real, 4 >
        gx{{ 1.0 - (x[N[0]]*x[N[0]] + y[N[0]]*y[N[0]] + z[N[0]]*z[N[0]]),
             1.0 - (x[N[1]]*x[N[1]] + y[N[1]]*y[N[1]] + z[N[1]]*z[N[1]]),
             1.0 - (x[N[2]]*x[N[2]] + y[N[2]]*y[N[2]] + z[N[2]]*z[N[2]]),
             1.0 - (x[N[3]]*x[N[3]] + y[N[3]]*y[N[3]] + z[N[3]]*z[N[3]]) }};

      // derivative of spatial component of density field
      std::array< std::array< tk::real, 4 >, 3 >
        dg{{  {{ -2.0*x[N[0]], -2.0*x[N[1]], -2.0*x[N[2]], -2.0*x[N[3]] }},
              {{ -2.0*y[N[0]], -2.0*y[N[1]], -2.0*y[N[2]], -2.0*y[N[3]] }},
              {{ -2.0*z[N[0]], -2.0*z[N[1]], -2.0*z[N[2]], -2.0*z[N[3]] }} }};

      // spatial component of energy field
      // NOTE: look into making beta_{x,y,z} into a vector
      std::array< tk::real, 4 > h{{
        hx( bx, by, bz, x[N[0]], y[N[0]], z[N[0]] ),
        hx( bx, by, bz, x[N[1]], y[N[1]], z[N[1]] ),
        hx( bx, by, bz, x[N[2]], y[N[2]], z[N[2]] ),
        hx( bx, by, bz, x[N[3]], y[N[3]], z[N[3]] ) }};

      // derivative of spatial component of energy field
      std::array< std::array< tk::real, 4 >, 3 >
        dh{{ {{ -bx*M_PI*std::sin(bx*M_PI*x[N[0]])*std::cos(by*M_PI*y[N[0]])*
                  std::cos(bz*M_PI*z[N[0]]),
                -bx*M_PI*std::sin(bx*M_PI*x[N[1]])*std::cos(by*M_PI*y[N[1]])*
                  std::cos(bz*M_PI*z[N[1]]),
                -bx*M_PI*std::sin(bx*M_PI*x[N[2]])*std::cos(by*M_PI*y[N[2]])*
                  std::cos(bz*M_PI*z[N[2]]),
                -bx*M_PI*std::sin(bx*M_PI*x[N[3]])*std::cos(by*M_PI*y[N[3]])*
                  std::cos(bz*M_PI*z[N[3]]) }},
             {{ -by*M_PI*std::cos(bx*M_PI*x[N[0]])*std::sin(by*M_PI*y[N[0]])*
                  std::cos(bz*M_PI*z[N[0]]),
                -by*M_PI*std::cos(bx*M_PI*x[N[1]])*std::sin(by*M_PI*y[N[1]])*
                  std::cos(bz*M_PI*z[N[1]]),
                -by*M_PI*std::cos(bx*M_PI*x[N[2]])*std::sin(by*M_PI*y[N[2]])*
                  std::cos(bz*M_PI*z[N[2]]),
                -by*M_PI*std::cos(bx*M_PI*x[N[3]])*std::sin(by*M_PI*y[N[3]])*
                  std::cos(bz*M_PI*z[N[3]]) }},
             {{ -bz*M_PI*std::cos(bx*M_PI*x[N[0]])*std::cos(by*M_PI*y[N[0]])*
                  std::sin(bz*M_PI*z[N[0]]),
                -bz*M_PI*std::cos(bx*M_PI*x[N[1]])*std::cos(by*M_PI*y[N[1]])*
                  std::sin(bz*M_PI*z[N[1]]),
                -bz*M_PI*std::cos(bx*M_PI*x[N[2]])*std::cos(by*M_PI*y[N[2]])*
                  std::sin(bz*M_PI*z[N[2]]),
                -bz*M_PI*std::cos(bx*M_PI*x[N[3]])*std::cos(by*M_PI*y[N[3]])*
                  std::sin(bz*M_PI*z[N[3]]) }} }};

      // density source
      std::array< tk::real, 4 >
        Sr{{ -a*std::exp(-a*t)*gx[0],
             -a*std::exp(-a*t)*gx[1],
             -a*std::exp(-a*t)*gx[2],
             -a*std::exp(-a*t)*gx[3] }};

      std::array< tk::real, 4 > ec4{{
        ec(ce,kappa,t,h[0],-4.0/3.0),
        ec(ce,kappa,t,h[1],-4.0/3.0),
        ec(ce,kappa,t,h[2],-4.0/3.0),
        ec(ce,kappa,t,h[3],-4.0/3.0) }};

      std::array< tk::real, 4 > ec1{{
        ec(ce,kappa,t,h[0],-1.0/3.0),
        ec(ce,kappa,t,h[1],-1.0/3.0),
        ec(ce,kappa,t,h[2],-1.0/3.0),
        ec(ce,kappa,t,h[3],-1.0/3.0) }};

      // energy source
      std::array< tk::real, 4 > Se{{
        (r0 + std::exp(-a*t)*gx[0])*kappa*h[0]*h[0]*ec4[0]
            - a*std::exp(-a*t)*gx[0]*ec1[0],
        (r0 + std::exp(-a*t)*gx[1])*kappa*h[1]*h[1]*ec4[1]
            - a*std::exp(-a*t)*gx[1]*ec1[1],
        (r0 + std::exp(-a*t)*gx[2])*kappa*h[2]*h[2]*ec4[2]
            - a*std::exp(-a*t)*gx[2]*ec1[2],
        (r0 + std::exp(-a*t)*gx[3])*kappa*h[3]*h[3]*ec4[3]
            - a*std::exp(-a*t)*gx[3]*ec1[3] }};

      // momentum source
      std::array< std::array< tk::real, 4 >, 3 > Sm{{
        {{ 2.0*kappa*h[0]*t*(r0+std::exp(-a*t)*gx[0])*(g-1)*ec4[0]*dh[0][0] +
             ec1[0]*(g-1)*std::exp(-a*t)*dg[0][0],
           2.0*kappa*h[1]*t*(r0+std::exp(-a*t)*gx[1])*(g-1)*ec4[1]*dh[0][1] +
             ec1[1]*(g-1)*std::exp(-a*t)*dg[0][1],
           2.0*kappa*h[2]*t*(r0+std::exp(-a*t)*gx[2])*(g-1)*ec4[2]*dh[0][2] +
             ec1[2]*(g-1)*std::exp(-a*t)*dg[0][2],
           2.0*kappa*h[3]*t*(r0+std::exp(-a*t)*gx[3])*(g-1)*ec4[3]*dh[0][3] +
             ec1[3]*(g-1)*std::exp(-a*t)*dg[0][3] }},
        {{ 2.0*kappa*h[0]*t*(r0+std::exp(-a*t)*gx[0])*(g-1)*ec4[0]*dh[1][0] +
             ec1[0]*(g-1)*std::exp(-a*t)*dg[1][0],
           2.0*kappa*h[1]*t*(r0+std::exp(-a*t)*gx[1])*(g-1)*ec4[1]*dh[1][1] +
             ec1[1]*(g-1)*std::exp(-a*t)*dg[1][1],
           2.0*kappa*h[2]*t*(r0+std::exp(-a*t)*gx[2])*(g-1)*ec4[2]*dh[1][2] +
             ec1[2]*(g-1)*std::exp(-a*t)*dg[1][2],
           2.0*kappa*h[3]*t*(r0+std::exp(-a*t)*gx[3])*(g-1)*ec4[3]*dh[1][3] +
             ec1[3]*(g-1)*std::exp(-a*t)*dg[1][3] }},
        {{ 2.0*kappa*h[0]*t*(r0+std::exp(-a*t)*gx[0])*(g-1)*ec4[0]*dh[2][0] +
             ec1[0]*(g-1)*std::exp(-a*t)*dg[2][0],
           2.0*kappa*h[1]*t*(r0+std::exp(-a*t)*gx[1])*(g-1)*ec4[1]*dh[2][1] +
             ec1[1]*(g-1)*std::exp(-a*t)*dg[2][1],
           2.0*kappa*h[2]*t*(r0+std::exp(-a*t)*gx[2])*(g-1)*ec4[2]*dh[2][2] +
             ec1[2]*(g-1)*std::exp(-a*t)*dg[2][2],
           2.0*kappa*h[3]*t*(r0+std::exp(-a*t)*gx[3])*(g-1)*ec4[3]*dh[2][3] +
             ec1[3]*(g-1)*std::exp(-a*t)*dg[2][3] }} }};

      // add momentum and energy source at element nodes
      for (std::size_t j=0; j<4; ++j)
        for (std::size_t k=0; k<4; ++k) {
          // source contribution to mass rhs
          R.var(r[0],N[j]) += dt * mass[j][k] * Sr[k];
          // source contribution to momentum rhs
          for (std::size_t l=0; l<3; ++l)
            R.var(r[l+1],N[j]) += dt * mass[j][k] * Sm[l][k];
          // source contribution to enerhy rhs
          R.var(r[4],N[j]) += dt * mass[j][k] * Se[k];
        }
    }

    //! \brief Query all side set IDs the user has configured for all components
    //!   in this PDE system
    //! \param[in,out] conf Set of unique side set IDs to add to
    static void side( std::unordered_set< int >& conf ) {
      using tag::param; using tag::compflow; using tag::bcdir;
      for (const auto& s : g_inputdeck.get< param, compflow, bcdir >())
        for (const auto& i : s)
          conf.insert( std::stoi(i) );
    }

    //! \brief Query Dirichlet boundary condition value on a given side set for
    //!    all components in this PDE system
    //! \param[in] sideset Side set ID
    //! \return Vector of pairs of bool and BC value for all components
    static std::vector< std::pair< bool, tk::real > > dirbc( int sideset ) {
      using tag::param; using tag::compflow; using tag::bcdir;
      std::vector< std::pair< bool, tk::real > > bc( 5, { false, 0.0 } );
      for (const auto& s : g_inputdeck.get< param, compflow, bcdir >())
        for (const auto& i : s)
          if (std::stoi(i) == sideset)
            for (auto& b : bc)
               b = { true, 0.0 };
      return bc;
    }

    //! Return field names to be output to file
    //! \return Vector of strings labelling fields output in file
    static std::vector< std::string > fieldNames() {
      std::vector< std::string > n;
      n.push_back( "density numerical" );
      n.push_back( "density analytical" );
      n.push_back( "x-velocity numerical" );
      n.push_back( "x-velocity analytical" );
      n.push_back( "y-velocity numerical" );
      n.push_back( "y-velocity analytical" );
      n.push_back( "z-velocity numerical" );
      n.push_back( "z-velocity analytical" );
      n.push_back( "specific total energy numerical" );
      n.push_back( "specific total energy analytical" );
      return n;
    }

    //! Return field output going to file
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    //! \param[in] t Physical time
    //! \param[in] coord Mesh node coordinates
    //! \param[in] U Solution vector at recent time step stage
    //! \return Vector of vectors to be output to file
    static std::vector< std::vector< tk::real > >
    fieldOutput( tk::ctr::ncomp_type e,
                 tk::ctr::ncomp_type offset,
                 tk::real t,
                 tk::real,
                 const std::vector< tk::real >&,
                 const std::array< std::vector< tk::real >, 3 >& coord,
                 const tk::Fields& U )
    {
      // manufactured solution parameters
      const auto& a =
        g_inputdeck.get< tag::param, tag::compflow, tag::alpha >()[e];
      const auto& bx =
        g_inputdeck.get< tag::param, tag::compflow, tag::betax >()[e];
      const auto& by =
        g_inputdeck.get< tag::param, tag::compflow, tag::betay >()[e];
      const auto& bz =
        g_inputdeck.get< tag::param, tag::compflow, tag::betaz >()[e];
      const auto& ce =
        g_inputdeck.get< tag::param, tag::compflow, tag::ce >()[e];
      const auto& kappa =
        g_inputdeck.get< tag::param, tag::compflow, tag::kappa >()[e];
      const auto& r0 =
        g_inputdeck.get< tag::param, tag::compflow, tag::r0 >()[e];

      std::vector< std::vector< tk::real > > out;
      const auto r  = U.extract( 0, offset );
      const auto ru = U.extract( 1, offset );
      const auto rv = U.extract( 2, offset );
      const auto rw = U.extract( 3, offset );
      const auto re = U.extract( 4, offset );

      // mesh node coordinates
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      std::vector< tk::real > rho = r;
      out.push_back( rho );
      tk::real ft = std::exp(-a*t);
      for (std::size_t i=0; i<rho.size(); ++i) {
        tk::real gx = 1.0 - (x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
        rho[i] = r0 + ft*gx;
      }
      out.push_back( rho );

      std::vector< tk::real > u = ru;
      std::transform( r.begin(), r.end(), u.begin(), u.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( u );
      for (std::size_t i=0; i<u.size(); ++i) u[i] = 0.0;
      out.push_back( u );

      std::vector< tk::real > v = rv;
      std::transform( r.begin(), r.end(), v.begin(), v.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( v );
      for (std::size_t i=0; i<v.size(); ++i) v[i] = 0.0;
      out.push_back( v );

      std::vector< tk::real > w = rw;
      std::transform( r.begin(), r.end(), w.begin(), w.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( w );
      for (std::size_t i=0; i<w.size(); ++i) w[i] = 0.0;
      out.push_back( w );

      std::vector< tk::real > E = re;
      std::transform( r.begin(), r.end(), E.begin(), E.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( E );
      for (std::size_t i=0; i<E.size(); ++i)
        E[i] = ec(ce, kappa, t, hx(bx,by,bz,x[i],y[i],z[i]), -1.0/3.0);
      out.push_back( E );

      return out;
   }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    static std::vector< std::string > names()
    { return { "r", "ru", "rv", "rw", "re" }; }

    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::NL_ENERGY_GROWTH; }
};

} // inciter::

#endif // CompFlowProblemNLEnergyGrowth_h
