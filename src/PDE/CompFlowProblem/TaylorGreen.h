// *****************************************************************************
/*!
  \file      src/PDE/CompFlowProblem/TaylorGreen.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Problem configurations for the compressible flow equations
  \details   This file defines policy classes for the compressible flow
    equations, defined in PDE/CompFlow.h. See PDE/CompFlow.h for general
    requirements on flow equations problem policy classes.
*/
// *****************************************************************************
#ifndef CompFlowProblemTaylorGreen_h
#define CompFlowProblemTaylorGreen_h

#include <string>
#include <unordered_set>

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Inciter/Options/Problem.h"

namespace inciter {

//! CompFlow system of PDEs problem: Taylor-Green
//! \see G.I. Taylor, A.E. Green, "Mechanism of the Production of Small Eddies
//!   from Large Ones", Proc. R. Soc. Lond. A 1937 158 499-521; DOI:
//!   10.1098/rspa.1937.0036. Published 3 February 1937
//! \see Waltz, et. al, "Verification of a three-dimensional unstructured finite
//!   element method using analytic and manufactured solutions", Computers and
//!   Fluids, 2013, Vol.81, pp.57-67.
class CompFlowProblemTaylorGreen {
  public:

    //! Set initial conditions
    //! \param[in] coord Mesh node coordinates
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in,out] unk Array of unknowns
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    static void init( const std::array< std::vector< tk::real >, 3 >& coord,
                      const std::vector< std::size_t >&,
                      tk::Fields& unk,
                      tk::ctr::ncomp_type e,
                      tk::ctr::ncomp_type offset,
                      tk::real )
    {
      Assert( coord[0].size() == unk.nunk(), "Size mismatch" );
      // dynamic = kinematic viscosity, since rho assumed 1.0
      // ratio of specific heats
      tk::real g =
        g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[e];
      // set initial and boundary conditions
      const auto& x = coord[0];
      const auto& y = coord[1];
      for (ncomp_t i=0; i<x.size(); ++i) {
        auto& r  = unk(i,0,offset); // rho
        auto& ru = unk(i,1,offset); // rho * u
        auto& rv = unk(i,2,offset); // rho * v
        auto& rw = unk(i,3,offset); // rho * w
        auto& re = unk(i,4,offset); // rho * e
        r = 1.0;
        ru = std::sin(M_PI*x[i]) * std::cos(M_PI*y[i]);
        rv = -std::cos(M_PI*x[i]) * std::sin(M_PI*y[i]);
        rw = 0.0;
        tk::real p = 10.0 +
          r/4.0*( std::cos(2.0*M_PI*x[i]) + std::cos(2.0*M_PI*y[i]) );
        re = p/(g-1.0)/r + 0.5*(ru*ru + rv*rv + rw*rw)/r/r;
      }
    }

    //! Add source term to rhs
    //! \param[in] coord Mesh node coordinates
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] dt Size of time step
    //! \param[in] N Element node indices
    //! \param[in] mass Element mass matrix, nnode*nnode [4][4]
    //! \param[in] r Pointers to right hand side at component and offset
    //! \param[in,out] R Right-hand side vector contributing to
    static void
    sourceRhs( tk::real,
               const std::array< std::vector< tk::real >, 3 >& coord,
               tk::ctr::ncomp_type,
               tk::real dt,
               const std::array< std::size_t, 4 >& N,
               const std::array< std::array< tk::real, 4 >, 4 >& mass,
               const std::array< const tk::real*, 5 >& r,
               tk::Fields& R )
    {
      // mesh node coordinates
      const auto& x = coord[0];
      const auto& y = coord[1];

      // compute energy source
      std::array< tk::real, 4 > Se{{
        3.0*M_PI/8.0*(std::cos(3.0*M_PI*x[N[0]])*std::cos(M_PI*y[N[0]]) -
                      std::cos(M_PI*x[N[0]])*std::cos(3.0*M_PI*y[N[0]])),
        3.0*M_PI/8.0*(std::cos(3.0*M_PI*x[N[1]])*std::cos(M_PI*y[N[1]]) -
                      std::cos(M_PI*x[N[1]])*std::cos(3.0*M_PI*y[N[1]])),
        3.0*M_PI/8.0*(std::cos(3.0*M_PI*x[N[2]])*std::cos(M_PI*y[N[2]]) -
                      std::cos(M_PI*x[N[2]])*std::cos(3.0*M_PI*y[N[2]])),
        3.0*M_PI/8.0*(std::cos(3.0*M_PI*x[N[3]])*std::cos(M_PI*y[N[3]]) -
                      std::cos(M_PI*x[N[3]])*std::cos(3.0*M_PI*y[N[3]])) }};

      std::array< tk::real, 4 > p;
      for (std::size_t i=0; i<4; ++i)
         p[i] = 10.0 +
           1.0/4.0*( std::cos(2.0*M_PI*x[N[i]]) + std::cos(2.0*M_PI*y[N[i]]) );

      // add source term at element nodes
      for (std::size_t alpha=0; alpha<4; ++alpha)
        for (std::size_t beta=0; beta<4; ++beta) {
          // source contribution to energy rhs
          R.var(r[4],N[alpha]) += dt * mass[alpha][beta] * Se[beta];
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
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] side Pair of side set ID and node IDs on the side set
    //! \return Vector of pairs of bool and boundary condition value associated
    //!   to mesh node IDs at which Dirichlet boundary conditions are set. Note
    //!   that instead of the actual boundary condition value, we return the
    //!   increment between t+dt and t, since that is what the solution requires
    //!   as we solve for the soution increments and not the solution itself.
    static std::unordered_map< std::size_t,
                               std::vector< std::pair< bool, tk::real > > >
    dirbc( tk::ctr::ncomp_type e,
           tk::real,
           tk::real,
           const std::pair< const int, std::vector< std::size_t > >& side,
           const std::array< std::vector< tk::real >, 3 >& )
    {
      using tag::param; using tag::compflow; using tag::bcdir;
      using NodeBC = std::vector< std::pair< bool, tk::real > >;
      std::unordered_map< std::size_t, NodeBC > bc;
      const auto& ubc = g_inputdeck.get< param, compflow, bcdir >();
      Assert( ubc.size() > e, "Indexing out of Dirichlet BC eq-vector" );
      for (const auto& b : ubc[e])
        if (std::stoi(b) == side.first)
        for (auto n : side.second)
          bc[n] = {{ {true,0.0}, {true,0.0}, {true,0.0}, {true,0.0},
                     {true,0.0} }};
      return bc;
    }

    //! Return field names to be output to file
    //! \return Vector of strings labelling fields output in file
    static std::vector< std::string > fieldNames() {
      std::vector< std::string > n;
      n.push_back( "density_numerical" );
      n.push_back( "density_analytical" );
      n.push_back( "x-velocity_numerical" );
      n.push_back( "x-velocity_analytical" );
      n.push_back( "err(u)" );
      n.push_back( "y-velocity_numerical" );
      n.push_back( "y-velocity_analytical" );
      n.push_back( "err(v)" );
      n.push_back( "z-velocity_numerical" );
      n.push_back( "z-velocity_analytical" );
      n.push_back( "specific_total_energy_numerical" );
      n.push_back( "specific_total_energy_analytical" );
      n.push_back( "err(E)" );
      n.push_back( "pressure_numerical" );
      n.push_back( "pressure_analytical" );
      return n;
    }

    //! Return field output going to file
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    //! \param[in] V Total mesh volume
    //! \param[in] vol Nodal mesh volumes
    //! \param[in] coord Mesh node coordinates
    //! \param[in] U Solution vector at recent time step stage
    //! \return Vector of vectors to be output to file
    static std::vector< std::vector< tk::real > >
    fieldOutput( tk::ctr::ncomp_type e,
                 tk::ctr::ncomp_type offset,
                 tk::real,
                 tk::real V,
                 const std::vector< tk::real >& vol,
                 const std::array< std::vector< tk::real >, 3 >& coord,
                 const tk::Fields& U )
    {
      // ratio of specific heats
      tk::real g =
        g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[e];

      std::vector< std::vector< tk::real > > out;
      const auto r  = U.extract( 0, offset );
      const auto ru = U.extract( 1, offset );
      const auto rv = U.extract( 2, offset );
      const auto rw = U.extract( 3, offset );
      const auto re = U.extract( 4, offset );

      // mesh node coordinates
      const auto& x = coord[0];
      const auto& y = coord[1];

      out.push_back( r );
      out.push_back( std::vector< tk::real >( r.size(), 1.0 ) );

      std::vector< tk::real > u = ru;
      std::transform( r.begin(), r.end(), u.begin(), u.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( u );
      std::vector< tk::real > ua = ru;
      for (std::size_t i=0; i<ua.size(); ++i)
        ua[i] = std::sin(M_PI*x[i]) * std::cos(M_PI*y[i]);
      out.push_back( ua );

      // error in x-velocity
      auto err = u;
      for (std::size_t i=0; i<u.size(); ++i)
         err[i] = std::pow( ua[i] - u[i], 2.0 ) * vol[i] / V;
       out.push_back( err );

      std::vector< tk::real > v = rv;
      std::vector< tk::real > va = rv;
      std::transform( r.begin(), r.end(), v.begin(), v.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( v );
      for (std::size_t i=0; i<va.size(); ++i)
        va[i] = -std::cos(M_PI*x[i]) * std::sin(M_PI*y[i]);
      out.push_back( va );

      // error in v-velocity
      for (std::size_t i=0; i<v.size(); ++i)
        err[i] = std::pow( va[i] - v[i], 2.0 ) * vol[i] / V;
      out.push_back( err );

      std::vector< tk::real > w = rw;
      std::vector< tk::real > wa = rw;
      std::transform( r.begin(), r.end(), w.begin(), w.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( w );
      for (std::size_t i=0; i<wa.size(); ++i)
        wa[i] = 0.0;
      out.push_back( wa );

      std::vector< tk::real > E = re;
      std::vector< tk::real > Ea = re;
      std::vector< tk::real > Pa( r.size(), 0.0 );
      std::transform( r.begin(), r.end(), E.begin(), E.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( E );
      for (std::size_t i=0; i<Ea.size(); ++i) {
        Pa[i] = 10.0 +
          r[i]/4.0*(std::cos(2.0*M_PI*x[i]) + std::cos(2.0*M_PI*y[i]));
        Ea[i] = Pa[i]/(g-1.0)/r[i] +
                0.5*(ua[i]*ua[i] + va[i]*va[i] + wa[i]*wa[i])/r[i];
      }
      out.push_back( Ea );

      // error in total specific energy
      for (std::size_t i=0; i<v.size(); ++i)
        err[i] = std::pow( Ea[i] - E[i], 2.0 ) * vol[i] / V;
      out.push_back( err );

      std::vector< tk::real > P( r.size(), 0.0 );
      for (std::size_t i=0; i<P.size(); ++i)
        P[i] = (g-1.0)*(E[i] - (u[i]*u[i] + v[i]*v[i] + w[i]*w[i])/2.0/r[i]);
      out.push_back( P );
      out.push_back( Pa );

      return out;
   }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    static std::vector< std::string > names()
    { return { "r", "ru", "rv", "rw", "re" }; }

    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::TAYLOR_GREEN; }
};

} // inciter::

#endif // CompFlowProblemTaylorGreen_h
