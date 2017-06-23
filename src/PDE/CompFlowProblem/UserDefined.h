// *****************************************************************************
/*!
  \file      src/PDE/CompFlowProblem/UserDefined.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Problem configurations for the compressible flow equations
  \details   This file defines a policy classe for the compressible flow
    equations, defined in PDE/CompFlow.h. See PDE/CompFlow.h for general
    requirements on flow equations problem policy classes.
*/
// *****************************************************************************
#ifndef CompFlowProblemUserDefined_h
#define CompFlowProblemUserDefined_h

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
    //! \param[in,out] unk Array of unknowns
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    static void init( const std::array< std::vector< tk::real >, 3 >& coord,
                      const std::vector< std::size_t >&,
                      tk::Fields& unk,
                      tk::ctr::ncomp_type,
                      tk::ctr::ncomp_type offset,
                      tk::real )
    {
      Assert( coord[0].size() == unk.nunk(), "Size mismatch" );
      const auto& x = coord[0];
      for (ncomp_t i=0; i<x.size(); ++i) {
        // domain points
        unk(i,0,offset) = 1.0;        // density
        unk(i,1,offset) = 0.0;        // density * velocity
        unk(i,2,offset) = 0.0;
        unk(i,3,offset) = 1.0;
        unk(i,4,offset) = 293.0;      // density * specific total energy
      }
    }

    //! Add source term to rhs
    //! \details No-op for user-defined problem
    static void
    sourceRhs( tk::real,
               const std::array< std::vector< tk::real >, 3 >&,
               tk::ctr::ncomp_type,
               tk::real,
               const std::array< std::size_t, 4 >&,
               const std::array< std::array< tk::real, 4 >, 4 >&,
               const std::array< const tk::real*, 5 >&,
               tk::Fields& ) {}

    //! Return field names to be output to file
    //! \return Vector of strings labelling fields output in file
    static std::vector< std::string > fieldNames() {
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

//     //! \brief Query Dirichlet boundary condition value on a given side set for
//     //!    all components in this PDE system
//     //! \param[in] sideset Side set ID
//     //! \return Vector of pairs of bool and BC value for all components
//     static std::vector< std::pair< bool, tk::real > >
//     dirbc( int sideset ) {
//       using tag::param; using tag::compflow; using tag::bcdir;
//       std::vector< std::pair< bool, tk::real > > b( 5, { false, 0.0 } );
//       for (const auto& s : g_inputdeck.get< param, compflow, bcdir >()) {
//         Assert( s.size() == 3, "Side set vector size incorrect" );
//         if (std::stoi(s[0]) == sideset)
//           b[ static_cast<std::size_t>(std::stol(s[1])-1) ] =
//             { true, std::atof(s[2].c_str()) };
//       }
//       return b;
//     }

    //! \brief Query Dirichlet boundary condition value on a given side set for
    //!    all components in this PDE system
    //! \return Vector of pairs of bool and boundary condition value associated
    //!   to mesh node IDs at which Dirichlet boundary conditions are set. Note
    //!   that instead of the actual boundary condition value, we return the
    //!   increment between t+dt and t, since that is what the solution requires
    //!   as we solve for the soution increments and not the solution itself.
    static std::unordered_map< std::size_t,
                               std::vector< std::pair< bool, tk::real > > >
    dirbc( tk::ctr::ncomp_type,
           tk::real,
           tk::real,
           const std::pair< const int, std::vector< std::size_t > >&,
           const std::array< std::vector< tk::real >, 3 >& )
    {
      using tag::param; using tag::compflow; using tag::bcdir;
      using NodeBC = std::vector< std::pair< bool, tk::real > >;
      std::unordered_map< std::size_t, NodeBC > bc;
      // TODO: include functionality from above dirbc() commented
      return bc;
    }

    //! Return field output going to file
    //! \param[in] offset System offset specifying the position of the system of
    //!   PDEs among other systems
    //! \param[in] U Solution vector at recent time step stage
    //! \return Vector of vectors to be output to file
    static std::vector< std::vector< tk::real > >
    fieldOutput( tk::ctr::ncomp_type,
                 tk::ctr::ncomp_type offset,
                 tk::real,
                 tk::real,
                 const std::vector< tk::real >&,
                 const std::array< std::vector< tk::real >, 3 >&,
                 const tk::Fields& U )
    {
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
      out.push_back( T );
      return out;
    }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    static std::vector< std::string > names()
    { return { "r", "ru", "rv", "rw", "re" }; }

   static ctr::ProblemType type() noexcept
   { return ctr::ProblemType::USER_DEFINED; }
};
} // inciter::

#endif // CompFlowProblemUserDefined_h
