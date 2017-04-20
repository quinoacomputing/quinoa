// *****************************************************************************
/*!
  \file      src/PDE/Transport.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Transport equation
  \details   This file implements the time integration of a transport equation
     of set of scalars.
*/
// *****************************************************************************
#ifndef Transport_h
#define Transport_h

#include <vector>
#include <array>
#include <limits>
#include <cmath>
#include <unordered_set>
#include <unordered_map>

#include "TransportPhysics.h"
#include "TransportProblem.h"
#include "Vector.h"
#include "DerivedData.h"
#include "Exception.h"
#include "Macro.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! \brief Transport equation used polymorphically with tk::PDE
//! \details The template argument(s) specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Physics - physics configuration, see PDE/TransportPhysics.h
//!   - Problem - problem configuration, see PDE/TransportProblem.h
//! \note The default physics is Advection, set in
//!    inciter::deck::check_transport()
template< class Physics, class Problem >
class Transport {

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! Constructor
    //! \param[in] c Equation system index (among multiple systems configured)
    //! \author J. Bakosi
    explicit Transport( ncomp_t c ) :
      m_c( c ),
      m_ncomp(
        g_inputdeck.get< tag::component >().get< tag::transport >().at(c) ),
      m_offset(
        g_inputdeck.get< tag::component >().offset< tag::transport >(c) )
    {
      Problem::template errchk< tag::transport >( m_c, m_ncomp );
    }

    //! Initalize the transport equations using problem policy
    //! \param[in] coord Mesh node coordinates
    //! \param[in,out] unk Array of unknowns
    //! \param[in] t Physical time
    //! \author J. Bakosi
    void initialize( const std::array< std::vector< tk::real >, 3 >& coord,
                     tk::Fields& unk,
                     tk::real t,
                     const std::vector< std::size_t >&,
                     const std::unordered_map< std::size_t,
                       std::vector< std::pair< bool, tk::real > > >& )
    const {
      Problem::template
       init< tag::transport >( coord, unk, m_c, m_ncomp, m_offset, t );
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
              tk::Fields& lhsd,
              tk::Fields& lhso ) const
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
        const auto J = tk::triple( ba, ca, da ) / 120.0;
        Assert( J > 0, "Element Jacobian non-positive" );

        for (ncomp_t c=0; c<m_ncomp; ++c) {
          const auto r = lhsd.cptr( c, m_offset );
          lhsd.var( r, A ) += 2.0 * J;
          lhsd.var( r, B ) += 2.0 * J;
          lhsd.var( r, C ) += 2.0 * J;
          lhsd.var( r, D ) += 2.0 * J;

          const auto s = lhso.cptr( c, m_offset );
          lhso.var( s, spidx(A,B) ) += J;
          lhso.var( s, spidx(A,C) ) += J;
          lhso.var( s, spidx(A,D) ) += J;

          lhso.var( s, spidx(B,A) ) += J;
          lhso.var( s, spidx(B,C) ) += J;
          lhso.var( s, spidx(B,D) ) += J;

          lhso.var( s, spidx(C,A) ) += J;
          lhso.var( s, spidx(C,B) ) += J;
          lhso.var( s, spidx(C,D) ) += J;

          lhso.var( s, spidx(D,A) ) += J;
          lhso.var( s, spidx(D,B) ) += J;
          lhso.var( s, spidx(D,C) ) += J;
        }
      }
    }

    //! Compute right hand side
    //! \param[in] mult Multiplier differentiating the different stages in
    //!    multi-stage time stepping
    //! \param[in] deltat Size of time step
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] U Solution vector at recent time step stage
    //! \param[in,out] R Right-hand side vector computed
    //! \author J. Bakosi
    void rhs( tk::real mult,
              tk::real deltat,
              const std::array< std::vector< tk::real >, 3 >& coord,
              const std::vector< std::size_t >& inpoel,
              const tk::Fields& U,
              tk::Fields& R ) const
    {
      Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
              "vector at recent time step incorrect" );
      Assert( R.nunk() == coord[0].size() && R.nprop() == m_ncomp,
              "Number of unknowns and/or number of components in right-hand "
              "side vector incorrect" );

      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      // zero right hand side for all components
      for (ncomp_t c=0; c<m_ncomp; ++c) R.fill( c, m_offset, 0.0 );

      for (std::size_t e=0; e<inpoel.size()/4; ++e) {
        const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                               inpoel[e*4+2], inpoel[e*4+3] }};
        // compute element Jacobi determinant
        const std::array< tk::real, 3 >
          ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
          ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
          da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
        const auto J = tk::triple( ba, ca, da );
        Assert( J > 0, "Element Jacobian non-positive" );

        // construct tetrahedron element-level matrices

        // consistent mass
        std::array< std::array< tk::real, 4 >, 4 > mass;  // nnode*nnode [4][4]
        mass[0][0] = mass[1][1] = mass[2][2] = mass[3][3] = J/60.0;  // diagonal
        mass[0][1] = mass[0][2] = mass[0][3] =                   // off-diagonal
        mass[1][0] = mass[1][2] = mass[1][3] =
        mass[2][0] = mass[2][1] = mass[2][3] =
        mass[3][0] = mass[3][1] = mass[3][2] = J/120.0;

        // shape function derivatives
        std::array< std::array< tk::real, 3 >, 4 > grad;  // nnode*ndim [4][3]
        grad[1] = tk::crossdiv( ca, da, J );
        grad[2] = tk::crossdiv( da, ba, J );
        grad[3] = tk::crossdiv( ba, ca, J );
        for (std::size_t i=0; i<3; ++i)
          grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];

        // access solution at element nodes at recent time step stage
        std::vector< std::array< tk::real, 4 > > u( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c) u[c] = U.extract( c, m_offset, N );
        // access pointer to right hand side at component and offset
        std::vector< const tk::real* > r( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c) r[c] = R.cptr( c, m_offset );

        // get velocity for problem
        const auto vel =
          Problem::template
            prescribedVelocity< tag::transport >( N, coord, m_c, m_ncomp );

        // add advection contribution to right hand side
        tk::real a = mult * deltat;
        for (ncomp_t c=0; c<m_ncomp; ++c)
          for (std::size_t i=0; i<4; ++i)
            for (std::size_t j=0; j<4; ++j)
              for (std::size_t k=0; k<3; ++k)
                for (std::size_t l=0; l<4; ++l)
                  R.var(r[c],N[j]) -= a * mass[j][i] * vel[c][k][i]
                                        * grad[l][k] * u[c][l];

        // add diffusion contribution to right hand side
        Physics::diffusionRhs( m_c, m_ncomp, mult, deltat, J, N, grad, u, r, R );
      }
    }

    //! Compute the minimum time step size
    //! \param[in] U Solution vector at recent time step stage
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \return Minimum time step size
    tk::real dt( const std::array< std::vector< tk::real >, 3 >& coord,
                 const std::vector< std::size_t >& inpoel,
                 const tk::Fields& U ) const
    {
      Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
              "vector at recent time step incorrect" );
      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];
      // compute the minimum dt across all elements we own
      tk::real mindt = std::numeric_limits< tk::real >::max();
      for (std::size_t e=0; e<inpoel.size()/4; ++e) {
        const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                               inpoel[e*4+2], inpoel[e*4+3] }};
        // compute cubic root of element volume as the characteristic length
        const std::array< tk::real, 3 >
          ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
          ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
          da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
        const auto L = std::cbrt( tk::triple( ba, ca, da ) / 6.0 );
        // access solution at element nodes at recent time step stage
        std::vector< std::array< tk::real, 4 > > u( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c) u[c] = U.extract( c, m_offset, N );
        // get velocity for problem
        const auto vel =
          Problem::template
            prescribedVelocity< tag::transport >( N, coord, m_c, m_ncomp );
        // compute the maximum length of the characteristic velocity (advection
        // velocity) across the four element nodes
        tk::real maxvel = 0.0;
        for (ncomp_t c=0; c<m_ncomp; ++c)
          for (std::size_t i=0; i<4; ++i) {
            auto v = std::sqrt( vel[c][0][i]*vel[c][0][i] + 
                                vel[c][1][i]*vel[c][1][i] +
                                vel[c][2][i]*vel[c][2][i] );
            if (v > maxvel) maxvel = v;
          }
        // compute element dt for the advection
        auto advection_dt = L / maxvel;
        // compute element dt based on diffusion
        auto diffusion_dt = Physics::diffusion_dt( m_c, m_ncomp, L, u );
        // compute minimum element dt
        auto elemdt = std::min( advection_dt, diffusion_dt );
        // find minimum dt across all elements
        if (elemdt < mindt) mindt = elemdt;
      }
      return mindt;
    }

    //! Extract the transport velocity field at cell nodes
    //! \param[in] U Solution vector at recent time step stage
    //! \param[in] coord Mesh node coordinates
    //! \param[in] N Element node indices    
    //! \return Array of the four values of the transport velocity
    std::array< std::array< tk::real, 4 >, 3 >
    velocity( const tk::Fields& U,
              const std::array< std::vector< tk::real >, 3 >& coord,
              const std::array< std::size_t, 4 >& N ) const
    { return Problem::velocity( U, coord, N ); }

    //! \brief Query all side set IDs the user has configured for all components
    //!   in this PDE system
    static void side( std::unordered_set< int >& ) {}

    //! \brief Query Dirichlet boundary condition value set by the user on a
    //!   given side set for all components in this PDE system
    //! \param[in] sideset Side set ID
    //! \return Vector of pairs of bool and BC value for all components
    std::vector< std::pair< bool, tk::real > > dirbc( int sideset ) const {
      std::vector< std::pair< bool, tk::real > > b( m_ncomp, { false, 0.0 } );
      IGNORE(sideset);
      return b;
    }

    //! Return field names to be output to file
    //! \return Vector of strings labelling fields output in file
    //! \details This functions should be written in conjunction with output(),
    //!   which provides the vector of fields to be output
    std::vector< std::string > names() const {
      std::vector< std::string > n;
      const auto& depvar =
        g_inputdeck.get< tag::param, tag::transport, tag::depvar >().at(m_c);
      // will output numerical solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        n.push_back( depvar + std::to_string(c) + "_numerical" );
      // will output analytic solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        n.push_back( depvar + std::to_string(c) + "_analytic" );
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
            tk::Fields& U ) const
    {
      std::vector< std::vector< tk::real > > out;
      // will output numerical solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        out.push_back( U.extract( c, m_offset ) );
      // evaluate analytic solution at time t
      initialize( coord, U, t, {{}}, {{}} );
      // will output analytic solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        out.push_back( U.extract( c, m_offset ) );
      return out;
    }

  private:
    const ncomp_t m_c;                  //!< Equation system index
    const ncomp_t m_ncomp;              //!< Number of components in this PDE
    const ncomp_t m_offset;             //!< Offset this PDE operates from
};

} // inciter::

#endif // Transport_h
