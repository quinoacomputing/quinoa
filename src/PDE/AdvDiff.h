//******************************************************************************
/*!
  \file      src/PDE/AdvDiff.h
  \author    J. Bakosi
  \date      Mon 02 May 2016 04:52:35 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Advection-diffusion equation of a transported scalar
  \details   This file implements the time integration of the
    advection-diffusion equation of a single scalar.
*/
//******************************************************************************
#ifndef AdvDiff_h
#define AdvDiff_h

#include "AdvDiffProblem.h"
#include "Vector.h"
#include "DerivedData.h"
#include "Exception.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! \brief Advection-diffusion equation used polymorphically with tk::PDE
//! \details The template argument(s) specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Problem - problem configuration, see PDE/AdvDiffProblem.h
template< class Problem >
class AdvDiff {

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! Constructor
    //! \param[in] c Equation system index (among multiple systems configured)
    //! \author J. Bakosi
    explicit AdvDiff( ncomp_t c ) :
      m_c( c ),
      m_ncomp( g_inputdeck.get< tag::component >().get< tag::advdiff >().at(c) ),
      m_offset( g_inputdeck.get< tag::component >().offset< tag::advdiff >(c) )
    {
      Problem::template errchk< tag::advdiff >( g_inputdeck, m_c, m_ncomp );
    }

    //! Initalize the advection-diffusion equations using problem policy
    //! \param[in,out] unk Array of unknowns
    //! \param[in] coord Mesh node coordinates
    //! \param[in] t Physical time
    //! \author J. Bakosi
    void initialize( const std::array< std::vector< tk::real >, 3 >& coord,
                     tk::MeshNodes& unk,
                     tk::real t )
    {
      Problem::template
        init< tag::advdiff >
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
              tk::MeshNodes& lhso )
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
              tk::MeshNodes& R )
    {
      Assert( U.nunk() == coord[0].size(), "Number of unknowns in solution "
              "vector at recent time step incorrect" );
      Assert( Un.nunk() == coord[0].size(), "Number of unknowns in solution "
              "vector at previous time step incorrect" );
      Assert( R.nunk() == coord[0].size(), "Number of unknowns in right-hand "
              "side vector incorrect" );

      const auto& x = coord[0];
      const auto& y = coord[1];
      const auto& z = coord[2];

      // get reference to diffusivities for all components
      const auto& diff =
        g_inputdeck.get< tag::param, tag::advdiff, tag::diffusivity >().at(m_c);

      // zero right hand side for all components
      for (ncomp_t c=0; c<m_ncomp; ++c) R.fill( c, m_offset, 0.0 );

      for (std::size_t e=0; e<inpoel.size()/4; ++e) {
        const auto A = inpoel[e*4+0];
        const auto B = inpoel[e*4+1];
        const auto C = inpoel[e*4+2];
        const auto D = inpoel[e*4+3];

        // compute element Jacobi determinant
        const std::array< tk::real, 3 > ba{{ x[B]-x[A], y[B]-y[A], z[B]-z[A] }},
                                        ca{{ x[C]-x[A], y[C]-y[A], z[C]-z[A] }},
                                        da{{ x[D]-x[A], y[D]-y[A], z[D]-z[A] }};
        const auto J = tk::triple( ba, ca, da );

        // construct tetrahedron element-level matrices

        // consistent mass
        std::array< std::array< tk::real, 4 >, 4 > mass;  // nnode*nnode [4][4]
        // diagonal
        mass[0][0] = mass[1][1] = mass[2][2] = mass[3][3] = J/60.0;
        // off-diagonal
        mass[0][1] = mass[0][2] = mass[0][3] =
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

        // access solution at element nodes at time n
        std::vector< std::array< tk::real, 4 > > u( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c)
          u[c] = Un.extract( c, m_offset, A, B, C, D );
        // access solution at element nodes at recent time step stage
        std::vector< std::array< tk::real, 4 > > s( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c)
          s[c] = U.extract( c, m_offset, A, B, C, D );
        // access pointer to right hand side at component and offset
        std::vector< const tk::real* > r( m_ncomp );
        for (ncomp_t c=0; c<m_ncomp; ++c) r[c] = R.cptr( c, m_offset );

        // get velocity for problem
        const auto vel =
          Problem::template
            velocity< tag::advdiff >
                    ( g_inputdeck, A, B, C, D, coord, m_c, m_ncomp );

        // add mass contribution to right hand side
        for (ncomp_t c=0; c<m_ncomp; ++c)
          for (std::size_t j=0; j<4; ++j) {
            R.var(r[c],A) += mass[0][j] * u[c][j];
            R.var(r[c],B) += mass[1][j] * u[c][j];
            R.var(r[c],C) += mass[2][j] * u[c][j];
            R.var(r[c],D) += mass[3][j] * u[c][j];
          }

        // add advection contribution to right hand side
        for (ncomp_t c=0; c<m_ncomp; ++c)
          for (std::size_t j=0; j<4; ++j)
            for (std::size_t k=0; k<3; ++k)
              for (std::size_t l=0; l<4; ++l) {
                R.var(r[c],A) -= mult * dt * mass[0][j] * vel[c][k][j]
                                      * grad[l][k] * s[c][l];
                R.var(r[c],B) -= mult * dt * mass[1][j] * vel[c][k][j]
                                      * grad[l][k] * s[c][l];
                R.var(r[c],C) -= mult * dt * mass[2][j] * vel[c][k][j]
                                      * grad[l][k] * s[c][l];
                R.var(r[c],D) -= mult * dt * mass[3][j] * vel[c][k][j]
                                      * grad[l][k] * s[c][l];
              }

        // add diffusion contribution to right hand side
        for (ncomp_t c=0; c<m_ncomp; ++c)
          for (std::size_t j=0; j<4; ++j)
            for (std::size_t k=0; k<3; ++k) {
              R.var(r[c],A) -= mult * dt * diff[c] * J * grad[0][k]
                                    * grad[j][k] * s[c][j];
              R.var(r[c],B) -= mult * dt * diff[c] * J * grad[1][k]
                                    * grad[j][k] * s[c][j];
              R.var(r[c],C) -= mult * dt * diff[c] * J * grad[2][k]
                                    * grad[j][k] * s[c][j];
              R.var(r[c],D) -= mult * dt * diff[c] * J * grad[3][k]
                                    * grad[j][k] * s[c][j];
            }
      }
    }

    //! Advance unknowns according to the Euler equations
    //! \param[in,out] unk Array of unknowns
    //! \param[in] dt Time step size
    //! \param[in] t Physical time
    //! \author J. Bakosi
    //void advance( tk::MeshNodes& unk, tk::real dt, tk::real t ) {}
    void advance( tk::MeshNodes&, tk::real, tk::real ) {}

    //! Return field names to be output to file
    //! \return Vector of strings labelling fields output in file
    //! \details This functions should be written in conjunction with output(),
    //!   which provides the vector of fields to be output
    std::vector< std::string > names() {
      std::vector< std::string > n;
      const auto& depvar =
        g_inputdeck.get< tag::param, tag::advdiff, tag::depvar >().at(m_c);
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
            tk::MeshNodes& U )
    {
      std::vector< std::vector< tk::real > > out;
      // will output numerical solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        out.push_back( U.extract( c, m_offset ) );
      // evaluate analytic solution at time t
      initialize( coord, U, t );
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

#endif // AdvDiff_h
