// *****************************************************************************
/*!
  \file      src/PDE/Poisson.h
  \author    J. Bakosi
  \date      Wed 17 Aug 2016 08:45:25 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Poisson equation
  \details   This file implements the Poisson equation.
*/
// *****************************************************************************
#ifndef Poisson_h
#define Poisson_h

#include <cmath>

#include "PoissonPhysics.h"
#include "PoissonProblem.h"
#include "Vector.h"
#include "DerivedData.h"
#include "Exception.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! \brief Poisso equation used polymorphically with tk::PDE
//! \details The template argument(s) specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Physics - physics configuration, see PDE/PoissonPhysics.h
//!   - Problem - problem configuration, see PDE/PoissonProblem.h
template< class Physics, class Problem >
class Poisson {

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! Constructor
    //! \param[in] c Equation system index (among multiple systems configured)
    //! \author J. Bakosi
    explicit Poisson( ncomp_t c ) :
      m_c( c ),
      m_ncomp( g_inputdeck.get< tag::component >().get< tag::poisson >().at(c) ),
      m_offset( g_inputdeck.get< tag::component >().offset< tag::poisson >(c) )
    {
      Problem::template errchk< tag::poisson >( g_inputdeck, m_c, m_ncomp );
    }

    //! Initalize the Poisson equation(s) using problem policy
    //! \param[in,out] unk Array of unknowns
    //! \param[in] coord Mesh node coordinates
    //! \param[in] t Physical time
    //! \author J. Bakosi
    void initialize( const std::array< std::vector< tk::real >, 3 >& coord,
                     tk::MeshNodes& unk,
                     tk::real t ) const
    {
      Problem::template
        init< tag::poisson >
            ( g_inputdeck, coord, unk, m_c, m_ncomp, m_offset, t );
    }

    //! Compute the left hand side sparse matrix
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] psup Linked lists storing IDs of points surrounding points
    //! \param[in,out] lhsd Diagonal of the sparse matrix storing nonzeros
    //! \param[in,out] lhso Off-diagonal of the sparse matrix storing nonzeros
    //! \details Sparse matrix storing the nonzero matrix values at rows and
    //!   columns given by psup. The format is similar to compressed row
    //!   storage, but the diagonal and off-diagonal data are stored in separate
    //!   vectors. For the off-diagonal data the local row and column indices,
    //!   at which values are nonzero, are stored by psup (psup1 and psup2,
    //!   where psup2 holds the indices at which psup1 holds the point ids
    //!   surrounding points, see also tk::genPsup()). Note that the number of
    //!   mesh points (our chunk) npoin = psup.second.size()-1.
    //! \author J. Bakosi
    void lhs( const std::array< std::vector< tk::real >, 3 >& coord,
              const std::vector< std::size_t >& inpoel,
              const std::pair< std::vector< std::size_t >,
                               std::vector< std::size_t > >& psup,
              tk::MeshNodes& lhsd,
              tk::MeshNodes& lhso ) const
    {
      Assert( psup.second.size()-1 == coord[0].size(),
              "Number of mesh points and number of global IDs unequal" );
      Assert( lhsd.nunk() == psup.second.size()-1, "Number of unknowns in "
              "diagonal sparse matrix storage incorrect" );
      Assert( lhso.nunk() == psup.first.size(), "Number of unknowns in "
              "off-diagonal sparse matrix storage incorrect" );

      // Lambda to compute the sparse matrix vector index for row and column
      // indices. Used only for off-diagonal entries.
      auto spidx = [ &psup ]( std::size_t r, std::size_t c ) -> std::size_t {
        Assert( r != c, "Only for computing the off-diagonal indices" );
        for (auto i=psup.second[r]+1; i<=psup.second[r+1]; ++i)
          if (c == psup.first[i]) return i;
        Throw( "Cannot find row, column: " + std::to_string(r) + ',' +
               std::to_string(c) + " in sparse matrix" );
      };

      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      // Zero matrix for all components
      for (ncomp_t c=0; c<m_ncomp; ++c) {
        lhsd.fill( c, m_offset, 0.0 );
        lhso.fill( c, m_offset, 0.0 );
      }

      for (std::size_t e=0; e<inpoel.size()/4; ++e) {
        const auto A = inpoel[e*4+0];
        const auto B = inpoel[e*4+1];
        const auto C = inpoel[e*4+2];
        const auto D = inpoel[e*4+3];
        std::array< tk::real, 3 > ba{{ x[B]-x[A], y[B]-y[A], z[B]-z[A] }},
                                  ca{{ x[C]-x[A], y[C]-y[A], z[C]-z[A] }},
                                  da{{ x[D]-x[A], y[D]-y[A], z[D]-z[A] }};
        const auto J = tk::triple( ba, ca, da );

        // shape function derivatives
        std::array< std::array< tk::real, 3 >, 4 > grad;  // nnode*ndim [4][3]
        grad[1] = tk::crossdiv( ca, da, J );
        grad[2] = tk::crossdiv( da, ba, J );
        grad[3] = tk::crossdiv( ba, ca, J );
        for (std::size_t i=0; i<3; ++i)
          grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];

        const auto iab = spidx(A,B);
        const auto iac = spidx(A,C);
        const auto iad = spidx(A,D);
        const auto iba = spidx(B,A);
        const auto ibc = spidx(B,C);
        const auto ibd = spidx(B,D);
        const auto ica = spidx(C,A);
        const auto icb = spidx(C,B);
        const auto icd = spidx(C,D);
        const auto ida = spidx(D,A);
        const auto idb = spidx(D,B);
        const auto idc = spidx(D,C);

        for (ncomp_t c=0; c<m_ncomp; ++c) {
          const auto r = lhsd.cptr( c, m_offset );
          const auto s = lhso.cptr( c, m_offset );

          for (std::size_t i=0; i<4; ++i)
            for (std::size_t j=0; j<4; ++j)
              for (std::size_t k=0; k<3; ++k) {
                const auto val = J/6.0 * grad[i][k] * grad[j][k];
                if (i==j) {
                  lhsd.var( r, A ) -= val;
                  lhsd.var( r, B ) -= val;
                  lhsd.var( r, C ) -= val;
                  lhsd.var( r, D ) -= val;
                } else {
                  lhso.var( s, iab ) -= val;
                  lhso.var( s, iac ) -= val;
                  lhso.var( s, iad ) -= val;
                  lhso.var( s, iba ) -= val;
                  lhso.var( s, ibc ) -= val;
                  lhso.var( s, ibd ) -= val;
                  lhso.var( s, ica ) -= val;
                  lhso.var( s, icb ) -= val;
                  lhso.var( s, icd ) -= val;
                  lhso.var( s, ida ) -= val;
                  lhso.var( s, idb ) -= val;
                  lhso.var( s, idc ) -= val;
                }
              }
        }
      }
    }

    //! Compute right hand side
    //! \param[in] mult Multiplier differentiating the different stages in
    //!    multi-stage time stepping
    //! \param[in] dt Size of time step
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] U Solution vector at recent time step stage
    //! \param[in] Un Solution vector at previous time step
    //! \param[in,out] R Right-hand side vector computed
    //! \author J. Bakosi
    void rhs( tk::real mult,
              tk::real dt,
              const std::array< std::vector< tk::real >, 3 >& coord,
              const std::vector< std::size_t >& inpoel,
              const tk::MeshNodes& U,
              const tk::MeshNodes& Un,
              tk::MeshNodes& R ) const
    {
      IGNORE(mult);
      IGNORE(dt);
      IGNORE(coord);
      IGNORE(inpoel);
      IGNORE(U);
      IGNORE(Un);

      Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
              "vector at recent time step incorrect" );
      Assert( Un.nunk() == coord[0].size(), "Number of unknowns in solution "
              "vector at previous time step incorrect" );
      Assert( R.nunk() == coord[0].size(), "Number of unknowns in right-hand "
              "side vector incorrect" );

      // zero right hand side for all components
      for (ncomp_t c=0; c<m_ncomp; ++c) R.fill( c, m_offset, 0.0 );
    }

    //! Extract the velocity field at cell nodes
    //! \param[in] U Solution vector at recent time step stage
    //! \param[in] A Index of 1st cell node to query
    //! \param[in] B Index of 2nd cell node to query
    //! \param[in] C Index of 3rd cell node to query
    //! \param[in] D Index of 4th cell node to query
    //! \return Array of the four values of the three velocity coordinates
    std::vector< std::array< tk::real, 4 > >
    velocity( const tk::MeshNodes& U,
              ncomp_t A, ncomp_t B, ncomp_t C, ncomp_t D ) const
    {
      IGNORE(U); IGNORE(A); IGNORE(B); IGNORE(C); IGNORE(D);
      std::vector< std::array< tk::real, 4 > > v;
      return v;
    }

    //! \brief Query if a Dirichlet boundary condition has set by the user on
    //!   a given side set for any component in the PDE system
    //! \param[in] sideset Side set ID
    //! \return True if the user has set a Dirichlet boundary condition on the
    //!   side set given for any component in the PDE system.
    bool anydirbc( int sideset ) const {
      const auto& bc =
        g_inputdeck.get< tag::param, tag::poisson, tag::bc_dirichlet >();
      for (const auto& s : bc)
        if (static_cast<int>(std::round(s[0])) == sideset)
          return true;
      return false;
    }

    //! \brief Query Dirichlet boundary condition value set by the user on a
    //!   given side set for all components in this PDE system
    //! \param[in] sideset Side set ID
    //! \return Vector of pairs of bool and BC value for all components
    std::vector< std::pair< bool, tk::real > > dirbc( int sideset ) const {
      const auto& bc =
        g_inputdeck.get< tag::param, tag::poisson, tag::bc_dirichlet >();
      std::vector< std::pair< bool, tk::real > > b( m_ncomp, { false, 0.0 } );
      for (const auto& s : bc) {
        Assert( s.size() == 3, "Side set vector size incorrect" );
        if (static_cast<int>(std::round(s[0])) == sideset)
          b[ static_cast<std::size_t>(std::round(s[1]))-1 ] = { true, s[2] };
      }
      return b;
    }

    //! Return field names to be output to file
    //! \return Vector of strings labelling fields output in file
    //! \details This functions should be written in conjunction with output(),
    //!   which provides the vector of fields to be output
    std::vector< std::string > names() const {
      std::vector< std::string > n;
      const auto& depvar =
        g_inputdeck.get< tag::param, tag::poisson, tag::depvar >().at(m_c);
      // will output numerical solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        n.push_back( depvar + std::to_string(c) + "_numerical" );
//       // will output analytic solution for all components
//       for (ncomp_t c=0; c<m_ncomp; ++c)
//         n.push_back( depvar + std::to_string(c) + "_analytic" );
      return n;
    }

    //! Return field output going to file
    //! \param[in] t Physical time
    //! \param[in] coord Mesh node coordinates
    //! \param[in,out] U Solution vector at recent time step stage
    //! \return Vector of vectors to be output to file
    //! \details This functions should be written in conjunction with names(),
    //!   which provides the vector of field names
    //! \note U is overwritten
    std::vector< std::vector< tk::real > >
    output( tk::real t,
            const std::array< std::vector< tk::real >, 3 >& coord,
            tk::MeshNodes& U ) const
    {
      IGNORE(t);
      IGNORE(coord);
      std::vector< std::vector< tk::real > > out;
      // will output numerical solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        out.push_back( U.extract( c, m_offset ) );
//       // evaluate analytic solution at time t
//       initialize( coord, U, t );
//       // will output analytic solution for all components
//       for (ncomp_t c=0; c<m_ncomp; ++c)
//         out.push_back( U.extract( c, m_offset ) );
      return out;
    }

  private:
    const ncomp_t m_c;                  //!< Equation system index
    const ncomp_t m_ncomp;              //!< Number of components in this PDE
    const ncomp_t m_offset;             //!< Offset this PDE operates from
};

} // inciter::

#endif // Poisson_h
