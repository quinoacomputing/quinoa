// *****************************************************************************
/*!
  \file      src/PDE/CompFlowProblem.h
  \author    J. Bakosi
  \date      Thu 25 Aug 2016 11:58:18 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Problem configurations for the compressible flow equations
  \details   This file defines policy classes for the compressible flow
    equations, defined in PDE/CompFlow.h.

    General requirements on flow equations problem policy classes:

    - Must define the static function _type()_, returning the enum value of the
      policy option. Example:
      \code{.cpp}
        static ctr::ProblemType type() noexcept {
          return ctr::ProblemType::USER_DEFINED;
        }
      \endcode
      which returns the enum value of the option from the underlying option
      class, collecting all possible options for coefficients policies.
*/
// *****************************************************************************
#ifndef CompFlowProblem_h
#define CompFlowProblem_h

#include <cmath>

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Inciter/Options/Problem.h"

namespace inciter {

//! CompFlow system of PDEs problem: user defined
class CompFlowProblemUserDefined {
  public:

    //! Set initial conditions
    //! \param[in] coord Mesh node coordinates
    //! \param[in,out] unk Array of unknowns
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    static void init( const ctr::InputDeck&,
                      const std::array< std::vector< tk::real >, 3 >& coord,
                      tk::Fields& unk,
                      tk::ctr::ncomp_type e,
                      tk::ctr::ncomp_type offset )
    {
      IGNORE(e);
      const auto& x = coord[0];
      for (ncomp_t i=0; i<x.size(); ++i) {
         unk( i, 0, offset ) = 1.0;     // density
         unk( i, 1, offset ) = 0.0;     // density * velocity
         unk( i, 2, offset ) = 0.0;
         unk( i, 3, offset ) = -1.0;
         unk( i, 4, offset ) = 293.0;     // density * specific total energy
      }
    }

    //! Add source term to rhs
    //! \details No-op for user-defined problem
    static void
    sourceRhs( const std::array< std::vector< tk::real >, 3 >&,
               tk::ctr::ncomp_type,
               tk::real,
               tk::real,
               tk::real,
               const std::array< std::size_t, 4 >&,
               const std::array< std::array< tk::real, 3 >, 4 >&,
               const std::array< std::array< tk::real, 4 >, 4 >&,
               const std::array< const tk::real*, 5 >&,
               std::array< std::array< tk::real, 4 >, 5 >&,
               tk::Fields&,
               tk::Fields& ) {}

    //! Return field names to be output to file
    //! \return Vector of strings labelling fields output in file
    static std::vector< std::string > names() {
      std::vector< std::string > n;
      n.push_back( "density" );
      n.push_back( "x-velocity" );
      n.push_back( "y-velocity" );
      n.push_back( "z-velocity" );
      n.push_back( "specific total energy" );
      n.push_back( "pressure" );
      n.push_back( "temperature" );
      return n;
    }

    //! Return field output going to file
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    //! \param[in] U Solution vector at recent time step stage
    //! \return Vector of vectors to be output to file
    static std::vector< std::vector< tk::real > >
    output( tk::ctr::ncomp_type e,
            tk::ctr::ncomp_type offset,
            tk::real,
            const std::array< std::vector< tk::real >, 3 >&,
            const tk::Fields& U )
    {
      IGNORE(e);
      std::vector< std::vector< tk::real > > out;
      const auto r = U.extract( 0, offset );
      const auto ru = U.extract( 1, offset );
      const auto rv = U.extract( 2, offset );
      const auto rw = U.extract( 3, offset );
      const auto re = U.extract( 4, offset );
      out.push_back( r );
      std::vector< tk::real > u = ru;
      std::transform( r.begin(), r.end(), u.begin(), u.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( u );
      std::vector< tk::real > v = rv;
      std::transform( r.begin(), r.end(), v.begin(), v.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( v );
      std::vector< tk::real > w = rw;
      std::transform( r.begin(), r.end(), w.begin(), w.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( w );
      std::vector< tk::real > E = re;
      std::transform( r.begin(), r.end(), E.begin(), E.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( E );
      std::vector< tk::real > p = r;
      tk::real g = g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[0];
      for (std::size_t i=0; i<p.size(); ++i)
        p[i] = (g-1.0)*r[i]*(E[i] - (u[i]*u[i] + v[i]*v[i] + w[i]*w[i])/2.0);
      out.push_back( p );
      std::vector< tk::real > T = r;
      tk::real cv = g_inputdeck.get< tag::param, tag::compflow, tag::cv >()[0];
      for (std::size_t i=0; i<T.size(); ++i)
        T[i] = cv*(E[i] - (u[i]*u[i] + v[i]*v[i] + w[i]*w[i])/2.0);
      out.push_back( p );
      return out;
   }

   static ctr::ProblemType type() noexcept
   { return ctr::ProblemType::USER_DEFINED; }
};

//! CompFlow system of PDEs problem: vortical flow
//! \see Waltz, et. al, "Manufactured solutions for the three-dimensional Euler
//!   equations with relevance to Inertial Confinement Fusion", Journal of
//!   Computational Physics 267 (2014) 196-209.
class CompFlowProblemVorticalFlow {
  public:

    //! Set initial conditions
    //! \param[in] coord Mesh node coordinates
    //! \param[in,out] unk Array of unknowns
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    static void init( const ctr::InputDeck&,
                      const std::array< std::vector< tk::real >, 3 >& coord,
                      tk::Fields& unk,
                      tk::ctr::ncomp_type e,
                      tk::ctr::ncomp_type offset )
    {
      IGNORE(e);

      // manufactured solution parameters
      tk::real g =
        g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[e];
      const auto& p0 =
        g_inputdeck.get< tag::param, tag::compflow, tag::p0 >()[e];

      const auto& x = coord[0];
      for (ncomp_t i=0; i<x.size(); ++i) {
         unk( i, 0, offset ) = 1.0;
         unk( i, 1, offset ) = 0.0;
         unk( i, 2, offset ) = 0.0;
         unk( i, 3, offset ) = 0.0;
         unk( i, 4, offset ) = p0/(g-1);
      }
    }

    //! Add source term to rhs for vortical flow manufactured solution
    //! \param[in] coord Mesh node coordinates
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] mult Multiplier differentiating the different stages in
    //!    multi-stage time stepping
    //! \param[in] dt Size of time step
    //! \param[in] J Element Jacobi determinant
    //! \param[in] N Element node indices
    //! \param[in] grad Shape function derivatives, nnode*ndim [4][3]
    //! \param[in] mass Element mass matrix, nnode*nnode [4][4]
    //! \param[in] r Pointers to right hand side at component and offset
    //! \param[in,out] s Solution at element nodes at recent time step stage
    //! \param[in,out] R Right-hand side vector contributing to
    static void
    sourceRhs( const std::array< std::vector< tk::real >, 3 >& coord,
               tk::ctr::ncomp_type e,
               tk::real mult,
               tk::real dt,
               tk::real J,
               const std::array< std::size_t, 4 >& N,
               const std::array< std::array< tk::real, 3 >, 4 >& grad,
               const std::array< std::array< tk::real, 4 >, 4 >& mass,
               const std::array< const tk::real*, 5 >& r,
               std::array< std::array< tk::real, 4 >, 5 >& s,
               tk::Fields& R,
               tk::Fields& U )
    {
      // manufactured solution parameters
      const auto& a =
        g_inputdeck.get< tag::param, tag::compflow, tag::alpha >()[e];
      const auto& b =
        g_inputdeck.get< tag::param, tag::compflow, tag::beta >()[e];
      const auto& p0 =
        g_inputdeck.get< tag::param, tag::compflow, tag::p0 >()[e];
      // ratio of specific heats
      tk::real g =
        g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[e];

      // mesh node coordinates
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      // set boundary conditions (directly on the solution)
      const auto eps = 1.0e-8;//std::numeric_limits< tk::real >::epsilon();
      for (std::size_t j=0; j<4; ++j)
        if ( std::fabs(x[N[j]]+0.5) < eps || std::fabs(x[N[j]]-0.5) < eps ||
             std::fabs(y[N[j]]+0.5) < eps || std::fabs(y[N[j]]-0.5) < eps ||
             std::fabs(z[N[j]]+0.5) < eps || std::fabs(z[N[j]]-0.5) < eps )
        {
          s[1][j] = U(N[j],1,0) = a*x[N[j]] - b*y[N[j]];
          s[2][j] = U(N[j],2,0) = b*x[N[j]] + a*y[N[j]];
          s[3][j] = U(N[j],3,0) = -2.0*a*z[N[j]];
          s[4][j] = U(N[j],4,0) = (s[1][j]*s[1][j] +
                                   s[2][j]*s[2][j] +
                                   s[3][j]*s[3][j])/2.0 +
                                  (p0 - 2.0*a*a*z[N[j]]*z[N[j]])/(g-1.0);
        }

      // add density source at element nodes
      tk::real c = mult * dt * J/24.0;
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<4; ++j)
          for (std::size_t k=0; k<4; ++k)
            R.var(r[0],N[j]) += c * grad[k][i] * s[i+1][k];

      // momentum source
      std::array< std::array< tk::real, 4 >, 2 > Sm{{
        {{ s[0][0]*((a*a-b*b)*x[N[0]] - 2.0*a*b*y[N[0]]),
           s[0][1]*((a*a-b*b)*x[N[1]] - 2.0*a*b*y[N[1]]),
           s[0][2]*((a*a-b*b)*x[N[2]] - 2.0*a*b*y[N[2]]),
           s[0][3]*((a*a-b*b)*x[N[3]] - 2.0*a*b*y[N[3]]) }},
        {{ s[0][0]*((a*a-b*b)*y[N[0]] + 2.0*a*b*x[N[0]]),
           s[0][1]*((a*a-b*b)*y[N[1]] + 2.0*a*b*x[N[1]]),
           s[0][2]*((a*a-b*b)*y[N[2]] + 2.0*a*b*x[N[2]]),
           s[0][3]*((a*a-b*b)*y[N[3]] + 2.0*a*b*x[N[3]]) }} }};

      // energy source
      std::array< tk::real, 4 > Se{{
        Sm[0][0]*s[1][0]/s[0][0] + Sm[1][0]*s[2][0]/s[0][0] +
          8.0*s[0][0]*a*a*a*z[N[0]]*z[N[0]]/(g-1.0),
        Sm[0][1]*s[1][1]/s[0][1] + Sm[1][1]*s[2][1]/s[0][1] +
          8.0*s[0][1]*a*a*a*z[N[1]]*z[N[1]]/(g-1.0),
        Sm[0][2]*s[1][2]/s[0][2] + Sm[1][2]*s[2][2]/s[0][2] +
          8.0*s[0][2]*a*a*a*z[N[2]]*z[N[2]]/(g-1.0),
        Sm[0][3]*s[1][3]/s[0][3] + Sm[1][3]*s[2][3]/s[0][3] +
          8.0*s[0][3]*a*a*a*z[N[3]]*z[N[3]]/(g-1.0) }};

      c = mult * dt;
      for (std::size_t j=0; j<4; ++j)
        for (std::size_t k=0; k<4; ++k) {
          // source contribution to momentum rhs
          for (std::size_t l=0; l<2; ++l)
            R.var(r[l+1],N[j]) += c * mass[j][k] * Sm[l][k];
          // source contribution to enerhy rhs
          R.var(r[4],N[j]) += c * mass[j][k] * Se[k];
        }
    }

    //! Return field names to be output to file
    //! \return Vector of strings labelling fields output in file
    static std::vector< std::string > names() {
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
    //! \param[in] coord Mesh node coordinates
    //! \param[in] U Solution vector at recent time step stage
    //! \return Vector of vectors to be output to file
    static std::vector< std::vector< tk::real > >
    output( tk::ctr::ncomp_type e,
            tk::ctr::ncomp_type offset,
            tk::real,
            const std::array< std::vector< tk::real >, 3 >& coord,
            const tk::Fields& U )
    {
      // manufactured solution parameters
      const auto& a =
        g_inputdeck.get< tag::param, tag::compflow, tag::alpha >()[e];
      const auto& b =
        g_inputdeck.get< tag::param, tag::compflow, tag::beta >()[e];
      const auto& p0 =
        g_inputdeck.get< tag::param, tag::compflow, tag::p0 >()[e];
      // ratio of specific heats
      tk::real g =
        g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[e];

      std::vector< std::vector< tk::real > > out;
      const auto r = U.extract( 0, offset );
      const auto ru = U.extract( 1, offset );
      const auto rv = U.extract( 2, offset );
      const auto rw = U.extract( 3, offset );
      const auto re = U.extract( 4, offset );

      // mesh node coordinates
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      out.push_back( r );
      out.push_back( std::vector< tk::real >( r.size(), 1.0 ) );

      std::vector< tk::real > u = ru;
      std::transform( r.begin(), r.end(), u.begin(), u.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( u );
      for (std::size_t i=0; i<u.size(); ++i) u[i] = a*x[i] - b*y[i];
      out.push_back( u );

      std::vector< tk::real > v = rv;
      std::transform( r.begin(), r.end(), v.begin(), v.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( v );
      for (std::size_t i=0; i<v.size(); ++i) v[i] = b*x[i] + a*y[i];
      out.push_back( v );

      std::vector< tk::real > w = rw;
      std::transform( r.begin(), r.end(), w.begin(), w.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( w );
      for (std::size_t i=0; i<w.size(); ++i) w[i] = -2.0*a*z[i];
      out.push_back( w );

      std::vector< tk::real > E = re;
      std::transform( r.begin(), r.end(), E.begin(), E.begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      out.push_back( E );
      for (std::size_t i=0; i<E.size(); ++i)
         E[i] = 0.5*(u[i]*u[i] + v[i]*v[i] + w[i]*w[i]) +
                (p0 - 2.0*r[i]*a*a*z[i]*z[i])/r[i]/(g-1.0);
      out.push_back( E );

      return out;
   }

    static ctr::ProblemType type() noexcept
    { return ctr::ProblemType::VORTICAL_FLOW; }
};

//! List of all CompFlow problem policies
using CompFlowProblems = boost::mpl::vector< CompFlowProblemUserDefined
                                           , CompFlowProblemVorticalFlow >;

} // inciter::

#endif // CompFlowProblem_h
