// *****************************************************************************
/*!
  \file      src/PDE/CompFlowProblem.h
  \author    J. Bakosi
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
#include <cstdlib>
#include <string>
#include <unordered_set>

#include <boost/mpl/vector.hpp>

#include "Types.h"
#include "Inciter/Options/Problem.h"

namespace inciter {

//! CompFlow system of PDEs problem: user defined
class CompFlowProblemUserDefined {
  public:

    //! Set initial conditions
    //! \param[in] coord Mesh node coordinates
    //! \param[in] gid Global node IDs of owned elements
    //! \param[in] bc Vector of pairs of bool and boundary condition value
    //!   associated to mesh node IDs at which to set Dirichlet boundary
    //!   conditions 
    //! \param[in,out] unk Array of unknowns
    //! \param[in] e Equation system index, i.e., which compressible
    //!   flow equation system we operate on among the systems of PDEs
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    static void init( const std::array< std::vector< tk::real >, 3 >& coord,
                      const std::vector< std::size_t >& gid,
                      const std::unordered_map< std::size_t,
                              std::vector< std::pair< bool, tk::real > > >& bc,
                      tk::Fields& unk,
                      tk::ctr::ncomp_type e,
                      tk::ctr::ncomp_type offset )
    {
      IGNORE(e);
      const auto& x = coord[0];
      for (ncomp_t i=0; i<x.size(); ++i) {
        // domain points
        unk(i,0,offset) = 1.0;        // density
        unk(i,1,offset) = 0.0;        // density * velocity
        unk(i,2,offset) = 0.0;
        unk(i,3,offset) = 1.0;
        unk(i,4,offset) = 293.0;      // density * specific total energy
        // boundary conditions
        const auto b = bc.find( gid[i] );
        if (b != end(bc)) {
          const auto& v = b->second;
          Assert( v.size() == 5, "Incorrect BC vector size" );
          for (std::size_t c=0; c<v.size(); ++c)
            if (v[c].first) unk(i,c,offset) = v[c].second;
        }
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
               const std::array< std::array< tk::real, 4 >, 4 >&,
               const std::array< std::array< tk::real, 3 >, 4 >&,
               const std::array< const tk::real*, 5 >&,
               std::array< std::array< tk::real, 4 >, 5 >&,
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

    //! \brief Query all side set IDs the user has configured for all components
    //!   in this PDE system
    //! \param[in,out] conf Set of unique side set IDs to add to
    static void side( std::unordered_set< int >& conf ) {
      using tag::param; using tag::compflow; using tag::bcdir;
      for (const auto& s : g_inputdeck.get< param, compflow, bcdir >())
        conf.insert( std::stoi(s[0]) );
    }

    //! \brief Query Dirichlet boundary condition value on a given side set for
    //!    all components in this PDE system
    //! \param[in] sideset Side set ID
    //! \return Vector of pairs of bool and BC value for all components
    static std::vector< std::pair< bool, tk::real > >
    dirbc( int sideset ) {
      using tag::param; using tag::compflow; using tag::bcdir;
      std::vector< std::pair< bool, tk::real > > b( 5, { false, 0.0 } );
      for (const auto& s : g_inputdeck.get< param, compflow, bcdir >()) {
        Assert( s.size() == 3, "Side set vector size incorrect" );
        if (std::stoi(s[0]) == sideset)
          b[ static_cast<std::size_t>(std::stol(s[1])-1) ] =
            { true, std::atof(s[2].c_str()) };
      }
      return b;
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
    static void init( const std::array< std::vector< tk::real >, 3 >& coord,
                      const std::vector< std::size_t >&,
                      const std::unordered_map< std::size_t,
                              std::vector< std::pair< bool, tk::real > > >&,
                      tk::Fields& unk,
                      tk::ctr::ncomp_type e,
                      tk::ctr::ncomp_type offset )
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
      // set initial and boundary conditions
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];
      for (ncomp_t i=0; i<x.size(); ++i) {
        auto& r  = unk(i,0,offset); // rho
        auto& ru = unk(i,1,offset); // rho * u
        auto& rv = unk(i,2,offset); // rho * v
        auto& rw = unk(i,3,offset); // rho * w
        auto& re = unk(i,4,offset); // rho * e
        r = 1.0;
        ru = a*x[i] - b*y[i];
        rv = b*x[i] + a*y[i];
        rw = -2.0*a*z[i];
        re = (ru*ru + rv*rv + rw*rw)/2.0 + (p0-2.0*a*a*z[i]*z[i])/(g-1.0);
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
    //! \param[in] mass Element mass matrix, nnode*nnode [4][4]
    //! \param[in] grad Shape function derivatives, nnode*ndim [4][3]
    //! \param[in] r Pointers to right hand side at component and offset
    //! \param[in,out] u Solution at element nodes at recent time step stage
    //! \param[in,out] R Right-hand side vector contributing to
    static void
    sourceRhs( const std::array< std::vector< tk::real >, 3 >& coord,
               tk::ctr::ncomp_type e,
               tk::real mult,
               tk::real dt,
               tk::real J,
               const std::array< std::size_t, 4 >& N,
               const std::array< std::array< tk::real, 4 >, 4 >& mass,
               const std::array< std::array< tk::real, 3 >, 4 >& grad,
               const std::array< const tk::real*, 5 >& r,
               std::array< std::array< tk::real, 4 >, 5 >& u,
               tk::Fields& R )
    {
      // manufactured solution parameters
      const auto& a =
        g_inputdeck.get< tag::param, tag::compflow, tag::alpha >()[e];
      const auto& b =
        g_inputdeck.get< tag::param, tag::compflow, tag::beta >()[e];
      // ratio of specific heats
      tk::real g =
        g_inputdeck.get< tag::param, tag::compflow, tag::gamma >()[e];

      // mesh node coordinates
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      std::array< tk::real, 4 > ru{{ a*x[N[0]] - b*y[N[0]],
                                     a*x[N[1]] - b*y[N[1]],
                                     a*x[N[2]] - b*y[N[2]],
                                     a*x[N[3]] - b*y[N[3]] }};
      std::array< tk::real, 4 > rv{{ b*x[N[0]] + a*y[N[0]],
                                     b*x[N[1]] + a*y[N[1]],
                                     b*x[N[2]] + a*y[N[2]],
                                     b*x[N[3]] + a*y[N[3]] }};

      // compute momentum source
      std::array< std::array< tk::real, 4 >, 3 >
        Sm{{ {{ a*ru[0] - b*rv[0],
                a*ru[1] - b*rv[1],
                a*ru[2] - b*rv[2],
                a*ru[3] - b*rv[3] }},
             {{ b*ru[0] + a*rv[0],
                b*ru[1] + a*rv[1],
                b*ru[2] + a*rv[2],
                b*ru[3] + a*rv[3] }} }};

      // compute energy source
      std::array< tk::real, 4 > Se{{
        Sm[0][0]*ru[0] + Sm[1][0]*rv[0] + 8.0*a*a*a*z[N[0]]*z[N[0]]/(g-1.0),
        Sm[0][1]*ru[1] + Sm[1][1]*rv[1] + 8.0*a*a*a*z[N[1]]*z[N[1]]/(g-1.0),
        Sm[0][2]*ru[2] + Sm[1][2]*rv[2] + 8.0*a*a*a*z[N[2]]*z[N[2]]/(g-1.0),
        Sm[0][3]*ru[3] + Sm[1][3]*rv[3] + 8.0*a*a*a*z[N[3]]*z[N[3]]/(g-1.0) }};

      // add momentum and energy source at element nodes
      tk::real c = mult * dt;
      for (std::size_t j=0; j<4; ++j)
        for (std::size_t k=0; k<4; ++k) {
          // source contribution to mass rhs
          for (std::size_t i=0; i<3; ++i)
            R.var(r[0],N[j]) += c * J/24.0 * grad[k][i] * u[i+1][k];
          // source contribution to momentum rhs
          for (std::size_t l=0; l<2; ++l)
            R.var(r[l+1],N[j]) += c * mass[j][k] * Sm[l][k];
          // source contribution to enerhy rhs
          R.var(r[4],N[j]) += c * mass[j][k] * Se[k];
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
      const auto r  = U.extract( 0, offset );
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
                (p0 - 2.0*a*a*z[i]*z[i])/(g-1.0);
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
