// *****************************************************************************
/*!
  \file      src/PDE/Transport/DGTransport.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Scalar transport using disccontinous Galerkin discretization
  \details   This file implements the physics operators governing transported
     scalars using disccontinuous Galerkin discretization.
*/
// *****************************************************************************
#ifndef DGTransport_h
#define DGTransport_h

#include <vector>
#include <array>
#include <limits>
#include <cmath>
#include <unordered_set>
#include <unordered_map>

#include "Macro.h"
#include "Exception.h"
#include "Vector.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

namespace dg {

//! \brief Transport equation used polymorphically with tk::DGPDE
//! \details The template argument(s) specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Physics - physics configuration, see PDE/Transport/Physics.h
//!   - Problem - problem configuration, see PDE/Transport/Problem.h
//! \note The default physics is DGAdvection, set in
//!    inciter::deck::check_transport()
template< class Physics, class Problem >
class Transport {

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! Constructor
    //! \param[in] c Equation system index (among multiple systems configured)
    explicit Transport( ncomp_t c ) :
      m_c( c ),
      m_ncomp(
        g_inputdeck.get< tag::component >().get< tag::transport >().at(c) ),
      m_offset(
        g_inputdeck.get< tag::component >().offset< tag::transport >(c) )
    {
      Problem::errchk( m_c, m_ncomp );
    }

    //! Initalize the transport equations using problem policy
//     //! \param[in] coord Mesh node coordinates
//     //! \param[in,out] unk Array of unknowns
//     //! \param[in] t Physical time
    void initialize( const std::array< std::vector< tk::real >, 3 >& /*coord*/,
                     tk::Fields& /*unk*/,
                     tk::real /*t*/,
                     const std::vector< std::size_t >& ) const
    {
      // Call Problem::solution() in a loop over all elements assigning the
      // initial conditions for cell centers. See cg::Transport::initialize()
      // for an example for a loop over all nodes using the node coordinates.
      // Instead of coord, we probably want to pass in a const-ref to
      // DG::m_geoElem and work with the cell centroid coordinates.
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
    //! \param[in] deltat Size of time step
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \param[in] U Solution vector at recent time step
    //! \param[in,out] R Right-hand side vector computed
    void rhs( tk::real,
              tk::real deltat,
              const std::array< std::vector< tk::real >, 3 >& coord,
              const std::vector< std::size_t >& inpoel,
              const tk::Fields& U,
              tk::Fields& Ue,
              tk::Fields& R ) const
    {
      using tag::transport;
      Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
              "vector at recent time step incorrect" );
      Assert( R.nunk() == coord[0].size() && R.nprop() == m_ncomp,
              "Number of unknowns and/or number of components in right-hand "
              "side vector incorrect" );
      Assert( U.nprop() == m_ncomp, "Number of components in solution vector "
              "must equal " + std::to_string(m_ncomp) );
      Assert( R.nprop() == m_ncomp, "Number of components in rhs must equal " +
              std::to_string(m_ncomp) );

      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      // 1st stage: update element values from node values (gather-add)
      for (std::size_t e=0; e<inpoel.size()/4; ++e) {

        // access node IDs
        const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                               inpoel[e*4+2], inpoel[e*4+3] }};
        // compute element Jacobi determinant
        const std::array< tk::real, 3 >
          ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
          ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
          da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
        const auto J = tk::triple( ba, ca, da );        // J = 6V
        Assert( J > 0, "Element Jacobian non-positive" );

        // shape function derivatives, nnode*ndim [4][3]
        std::array< std::array< tk::real, 3 >, 4 > grad;
        grad[1] = tk::crossdiv( ca, da, J );
        grad[2] = tk::crossdiv( da, ba, J );
        grad[3] = tk::crossdiv( ba, ca, J );
        for (std::size_t i=0; i<3; ++i)
          grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];

        // access solution at element nodes
        std::vector< std::array< tk::real, 4 > > u( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c) u[c] = U.extract( c, m_offset, N );
        // access solution at elements
        std::vector< const tk::real* > ue( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c) ue[c] = Ue.cptr( c, m_offset );

        // sum nodal averages to element
        for (ncomp_t c=0; c<m_ncomp; ++c) {
          Ue.var(ue[c],e) = 0.0;
          for (std::size_t a=0; a<4; ++a)
            Ue.var(ue[c],e) += u[c][a]/4.0;
        }

        // get prescribed velocity
        const std::array< std::vector<std::array<tk::real,3>>, 4 > vel{{
         Problem::prescribedVelocity(x[N[0]], y[N[0]], z[N[0]], m_c, m_ncomp),
         Problem::prescribedVelocity(x[N[0]], y[N[1]], z[N[1]], m_c, m_ncomp),
         Problem::prescribedVelocity(x[N[0]], y[N[2]], z[N[2]], m_c, m_ncomp),
         Problem::prescribedVelocity(x[N[0]], y[N[3]], z[N[3]], m_c, m_ncomp)}};

        // sum flux (advection) contributions to element
        tk::real d = deltat/2.0;
        for (std::size_t c=0; c<m_ncomp; ++c)
          for (std::size_t j=0; j<3; ++j)
            for (std::size_t a=0; a<4; ++a)
              Ue.var(ue[c],e) -= d * grad[a][j] * vel[a][c][j]*u[c][a];

      }


      // zero right hand side for all components
      for (ncomp_t c=0; c<m_ncomp; ++c) R.fill( c, m_offset, 0.0 );

      // 2nd stage: form rhs from element values (scatter-add)
      for (std::size_t e=0; e<inpoel.size()/4; ++e) {

        // access node IDs
        const std::array< std::size_t, 4 > N{{ inpoel[e*4+0], inpoel[e*4+1],
                                               inpoel[e*4+2], inpoel[e*4+3] }};
        // compute element Jacobi determinant
        const std::array< tk::real, 3 >
          ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
          ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
          da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
        const auto J = tk::triple( ba, ca, da );        // J = 6V
        Assert( J > 0, "Element Jacobian non-positive" );

        // shape function derivatives, nnode*ndim [4][3]
        std::array< std::array< tk::real, 3 >, 4 > grad;
        grad[1] = tk::crossdiv( ca, da, J );
        grad[2] = tk::crossdiv( da, ba, J );
        grad[3] = tk::crossdiv( ba, ca, J );
        for (std::size_t i=0; i<3; ++i)
          grad[0][i] = -grad[1][i]-grad[2][i]-grad[3][i];

        // access solution at elements
        std::vector< tk::real > ue( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c) ue[c] = Ue( e, c, m_offset );
        // access pointer to right hand side at component and offset
        std::vector< const tk::real* > r( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c) r[c] = R.cptr( c, m_offset );
        // access solution at nodes of element
        std::vector< std::array< tk::real, 4 > > u( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c) u[c] = U.extract( c, m_offset, N );

        // get prescribed velocity
        auto xc = (x[N[0]] + x[N[1]] + x[N[2]] + x[N[3]]) / 4.0;
        auto yc = (y[N[0]] + y[N[1]] + y[N[2]] + y[N[3]]) / 4.0;
        auto zc = (z[N[0]] + z[N[1]] + z[N[2]] + z[N[3]]) / 4.0;
        const auto vel =
          Problem::prescribedVelocity( xc, yc, zc, m_c, m_ncomp );

        // scatter-add flux contributions to rhs at nodes
        tk::real d = deltat * J/6.0;
        for (std::size_t c=0; c<m_ncomp; ++c)
          for (std::size_t j=0; j<3; ++j)
            for (std::size_t a=0; a<4; ++a)
              R.var(r[c],N[a]) += d * grad[a][j] * vel[c][j]*ue[c];

        // add (optional) diffusion contribution to right hand side
        Physics::diffusionRhs( m_c, m_ncomp, deltat, J, grad, N, u, r, R );

      }
    }

    //! Compute the minimum time step size
    //! \param[in] U Solution vector at recent time step
    //! \param[in] coord Mesh node coordinates
    //! \param[in] inpoel Mesh element connectivity
    //! \return Minimum time step size
    tk::real dt( const std::array< std::vector< tk::real >, 3 >& coord,
                 const std::vector< std::size_t >& inpoel,
                 const tk::Fields& U ) const
    {
      using tag::transport;
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
        // access solution at element nodes at recent time step
        std::vector< std::array< tk::real, 4 > > u( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c) u[c] = U.extract( c, m_offset, N );
        // get velocity for problem
        const std::array< std::vector<std::array<tk::real,3>>, 4 > vel{{
         Problem::prescribedVelocity(x[N[0]], y[N[0]], z[N[0]], m_c, m_ncomp),
         Problem::prescribedVelocity(x[N[0]], y[N[1]], z[N[1]], m_c, m_ncomp),
         Problem::prescribedVelocity(x[N[0]], y[N[2]], z[N[2]], m_c, m_ncomp),
         Problem::prescribedVelocity(x[N[0]], y[N[3]], z[N[3]], m_c, m_ncomp)}};
        // compute the maximum length of the characteristic velocity (advection
        // velocity) across the four element nodes
        tk::real maxvel = 0.0;
        for (ncomp_t c=0; c<m_ncomp; ++c)
          for (std::size_t i=0; i<4; ++i) {
            auto v = std::sqrt( vel[i][c][0]*vel[i][c][0] +
                                vel[i][c][1]*vel[i][c][1] +
                                vel[i][c][2]*vel[i][c][2] );
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

    //! \brief Query all side set IDs the user has configured for all components
    //!   in this PDE system
    //! \param[in,out] conf Set of unique side set IDs to add to
    void side( std::unordered_set< int >& conf ) const
    { Problem::side( conf ); }

    //! \brief Query Dirichlet boundary condition value on a given side set for
    //!    all components in this PDE system
//     //! \param[in] t Physical time
//     //! \param[in] deltat Time step size
//     //! \param[in] sides Pair of side set ID and face IDs on the side set
//     //! \param[in] coord Mesh face coordinates
    //! \return Vector of pairs of bool and boundary condition value associated
    //!   to mesh face IDs at which Dirichlet boundary conditions are set. Note
    //!   that instead of the actual boundary condition value, we return the
    //!   increment between t+dt and t, since that is what the solution requires
    //!   as we solve for the soution increments and not the solution itself.
    std::unordered_map< std::size_t, std::vector< std::pair<bool,tk::real> > >
    dirbc( tk::real /*t*/,
           tk::real /*deltat*/,
           const std::pair< const int, std::vector< std::size_t > >& /*sides*/,
           const std::array< std::vector< tk::real >, 3 >& /*coord*/ ) const
    {
      // Call Problem::solinc() within a search for all face centers of a side
      // set given in sides (key=setid, value=faceids). See
      // cg::Transport::dirbc() for an example search for all nodes of a side
      // set given in sides (key=setid, value=nodeids). Instead of coord, we
      // probably want to pass in a const-ref to DG::m_geoFace and work with the
      // face centroid coordinates.
      using FaceBC = std::vector< std::pair< bool, tk::real > >;
      std::unordered_map< std::size_t, FaceBC > bc;
      return bc;
    }

    //! Return field names to be output to file
    //! \return Vector of strings labelling fields output in file
    //! \details This functions should be written in conjunction with
    //!   fieldOutput(), which provides the vector of fields to be output
    std::vector< std::string > fieldNames() const {
      std::vector< std::string > n;
      const auto& depvar =
        g_inputdeck.get< tag::param, tag::transport, tag::depvar >().at(m_c);
      // will output numerical solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        n.push_back( depvar + std::to_string(c) + "_numerical" );
      // will output analytic solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        n.push_back( depvar + std::to_string(c) + "_analytic" );
      // will output error for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        n.push_back( depvar + std::to_string(c) + "_error" );
      return n;
    }

    //! Return field output going to file
    //! \param[in] t Physical time
    //! \param[in] V Total mesh volume
    //! \param[in] coord Mesh node coordinates
    //! \param[in] v Nodal volumes
    //! \param[in,out] U Solution vector at recent time step
    //! \return Vector of vectors to be output to file
    //! \details This functions should be written in conjunction with names(),
    //!   which provides the vector of field names
    //! \note U is overwritten
    std::vector< std::vector< tk::real > >
    fieldOutput( tk::real t,
                 tk::real V,
                 const std::array< std::vector< tk::real >, 3 >& coord,
                 const std::vector< tk::real >& v,
                 tk::Fields& U ) const
    {
      std::vector< std::vector< tk::real > > out;
      // will output numerical solution for all components
      auto E = U;
      for (ncomp_t c=0; c<m_ncomp; ++c)
        out.push_back( U.extract( c, m_offset ) );
      // evaluate analytic solution at time t
      initialize( coord, U, t, {{}} );
      // will output analytic solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        out.push_back( U.extract( c, m_offset ) );
      // will output error for all components
      for (ncomp_t c=0; c<m_ncomp; ++c) {
        auto u = U.extract( c, m_offset );
        auto e = E.extract( c, m_offset );
        Assert( u.size() == e.size(), "Size mismatch" );
        Assert( u.size() == v.size(), "Size mismatch" );
        for (std::size_t i=0; i<u.size(); ++i)
          e[i] = std::pow( e[i] - u[i], 2.0 ) * v[i] / V;
        out.push_back( e );
      }
      return out;
    }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    std::vector< std::string > names() const {
      std::vector< std::string > n;
      const auto& depvar =
        g_inputdeck.get< tag::param, tag::transport, tag::depvar >().at(m_c);
      // construct the name of the numerical solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        n.push_back( depvar + std::to_string(c) );
      return n;
    }

  private:
    const ncomp_t m_c;                  //!< Equation system index
    const ncomp_t m_ncomp;              //!< Number of components in this PDE
    const ncomp_t m_offset;             //!< Offset this PDE operates from
};

} // dg::
} // inciter::

#endif // DGTransport_h
