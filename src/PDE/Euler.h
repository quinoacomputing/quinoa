// *****************************************************************************
/*!
  \file      src/PDE/Euler.h
  \author    J. Bakosi
  \date      Wed 17 Aug 2016 07:46:41 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Euler equations describing compressible flow
  \details   This file implements the time integration of the Euler equations
    governing compressible fluid flow.
*/
// *****************************************************************************
#ifndef Euler_h
#define Euler_h

#include <cmath>

#include "Macro.h"
#include "EulerProblem.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! \brief Euler equations used polymorphically with tk::PDE
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Problem - problem configuration, see PDE/EulerProblem.h
template< class Problem >
class Euler {

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! \brief Constructor
    //! \author J. Bakosi
    explicit Euler( ncomp_t c ) :
      m_c( c ),
      m_ncomp( g_inputdeck.get< tag::component >().get< tag::euler >().at(c) ),
      m_offset( g_inputdeck.get< tag::component >().offset< tag::euler >(c) )
    {}

    //! Initalize the Euler equations, prepare for time integration
    //! \param[in,out] unk Array of unknowns
    //! \author J. Bakosi
    void initialize( const std::array< std::vector< tk::real >, 3 >& coord,
                     tk::MeshNodes& unk,
                     tk::real t ) const
    {
      IGNORE(coord);
      IGNORE(unk);
      IGNORE(t);
      //! Set initial conditions using problem configuration policy
      //Problem::template init< tag::euler >( g_inputdeck, unk, m_offset );
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
      IGNORE(coord);
      IGNORE(inpoel);
      IGNORE(psup);
      IGNORE(lhsd);
      IGNORE(lhso);
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
      IGNORE(R);
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
      std::vector< std::array< tk::real, 4 > > v( 3 );
      v[0] = U.extract( 1, m_offset, A, B, C, D );
      v[1] = U.extract( 2, m_offset, A, B, C, D );
      v[2] = U.extract( 3, m_offset, A, B, C, D );
      auto r = U.extract( 0, m_offset, A, B, C, D );
      std::transform( r.begin(), r.end(), v[0].begin(), v[0].begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      std::transform( r.begin(), r.end(), v[1].begin(), v[1].begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      std::transform( r.begin(), r.end(), v[2].begin(), v[2].begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      return v;
    }

    //! \brief Query if a Dirichlet boundary condition has set by the user on
    //!   any side set for any component in the PDE system
    //! \param[in] sideset Side set ID
    //! \return True if the user has set a Dirichlet boundary condition on any
    //!   of the side sets for any component in the PDE system.
    bool anydirbc( int sideset ) const {
      const auto& bc =
        g_inputdeck.get< tag::param, tag::euler, tag::bc_dirichlet >();
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
        g_inputdeck.get< tag::param, tag::euler, tag::bc_dirichlet >();
      std::vector< std::pair< bool, tk::real > > b( m_ncomp, { false, 0.0 } );
      IGNORE(sideset);
      IGNORE(bc);
      return b;
    }

    //! Return field names to be output to file
    //! \return Vector of strings labelling fields output in file
    std::vector< std::string > names() const {
      std::vector< std::string > n( m_ncomp );
      // ...
      return n;
    }

    //! Return field output going to file
    //! \param[in] t Physical time
    //! \param[in] coord Mesh node coordinates
    //! \param[in] U Solution vector at recent time step stage
    //! \return Vector of vectors to be output to file
    std::vector< std::vector< tk::real > >
    output( tk::real t,
            const std::array< std::vector< tk::real >, 3 >& coord,
            const tk::MeshNodes& U ) const
    {
      IGNORE(t);
      IGNORE(coord);
      IGNORE(U);
      std::vector< std::vector< tk::real > > out;
      // ...
      return out;
   }

  private:
    const ncomp_t m_c;                  //!< Equation system index
    const ncomp_t m_ncomp;              //!< Number of components
    const ncomp_t m_offset;             //!< Offset PDE operates from
};

} // inciter::

#endif // Euler_h
