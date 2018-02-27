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
                     tk::real /*t*/ ) const
    {
      // Call Problem::solution() in a loop over all elements assigning the
      // initial conditions for cell centers. See cg::Transport::initialize()
      // for an example for a loop over all nodes using the node coordinates.
      // Instead of coord, we probably want to pass in a const-ref to
      // DG::m_geoElem and work with the cell centroid coordinates.
      std::cout << "type: " <<
        std::to_string( static_cast< uint8_t >( Problem::type() )) << '\n';
    }

    //! Compute the left hand side sparse matrix
//     //! \param[in] coord Mesh node coordinates
//     //! \param[in] inpoel Mesh element connectivity
//     //! \param[in] psup Linked lists storing IDs of points surrounding points
//     //! \param[in,out] lhsd Diagonal of the sparse matrix storing nonzeros
//     //! \param[in,out] lhso Off-diagonal of the sparse matrix storing nonzeros
    //! \details Sparse matrix storing the nonzero matrix values at rows and
    //!   columns given by psup. The format is similar to compressed row
    //!   storage, but the diagonal and off-diagonal data are stored in separate
    //!   vectors. For the off-diagonal data the local row and column indices,
    //!   at which values are nonzero, are stored by psup (psup1 and psup2,
    //!   where psup2 holds the indices at which psup1 holds the point ids
    //!   surrounding points, see also tk::genPsup()). Note that the number of
    //!   mesh points (our chunk) npoin = psup.second.size()-1.
    void lhs( const std::array< std::vector< tk::real >, 3 >& /*coord*/,
              const std::vector< std::size_t >& /*inpoel*/,
              const std::pair< std::vector< std::size_t >,
                               std::vector< std::size_t > >& /*psup*/,
              tk::Fields& /*lhsd*/,
              tk::Fields& /*lhso*/ ) const
    {
    }

    //! Compute right hand side
//     //! \param[in] deltat Size of time step
//     //! \param[in] coord Mesh node coordinates
//     //! \param[in] inpoel Mesh element connectivity
//     //! \param[in] U Solution vector at recent time step
//     //! \param[in,out] R Right-hand side vector computed
    void rhs( tk::real,
              tk::real /*deltat*/,
              const std::array< std::vector< tk::real >, 3 >& /*coord*/,
              const std::vector< std::size_t >& /*inpoel*/,
              const tk::Fields& /*U*/,
              tk::Fields& /*Ue*/,
              tk::Fields& /*R*/ ) const
    {
    }

    //! Compute the minimum time step size
//     //! \param[in] U Solution vector at recent time step
//     //! \param[in] coord Mesh node coordinates
//     //! \param[in] inpoel Mesh element connectivity
    //! \return Minimum time step size
    tk::real dt( const std::array< std::vector< tk::real >, 3 >& /*coord*/,
                 const std::vector< std::size_t >& /*inpoel*/,
                 const tk::Fields& /*U*/ ) const
    {
      tk::real mindt = std::numeric_limits< tk::real >::max();
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
      initialize( coord, U, t );
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
