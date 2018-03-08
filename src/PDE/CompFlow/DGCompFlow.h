// *****************************************************************************
/*!
  \file      src/PDE/CompFlow/DGCompFlow.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Compressible single-material flow using discontinuous Galerkin
  \details   This file implements the physics operators governing compressible
    single-material flow using discontinuous Galerkin discretization.
*/
// *****************************************************************************
#ifndef DGCompFlow_h
#define DGCompFlow_h

#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>

#include "Macro.h"
#include "Exception.h"
#include "Vector.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

namespace dg {

//! \brief CompFlow used polymorphically with tk::DGPDE
//! \details The template arguments specify policies and are used to configure
//!   the behavior of the class. The policies are:
//!   - Physics - physics configuration, see PDE/CompFlow/Physics.h
//!   - Problem - problem configuration, see PDE/CompFlow/Problems.h
//! \note The default physics is Euler, set in inciter::deck::check_compflow()
template< class Physics, class Problem >
class CompFlow {

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! \brief Constructor
    explicit CompFlow( ncomp_t ) : m_offset( 0 ) {}

    //! Initalize the compressible flow equations, prepare for time integration
//     //! \param[in] geoElem Element geometry array
//     //! \param[in,out] unk Array of unknowns
//     //! \param[in] t Physical time
    void initialize( const tk::Fields& /*geoElem*/,
                     tk::Fields& /*unk*/,
                     tk::real /*t*/ ) const
    {
      // Loop over all cells and call Problem::solution() to query the initial
      // conditions or the analytical solutions (whichever Problem::solution()
      // evaluates for the given problem). See cg::CompFlow::initialize() for a
      // node-based example.
    }

    //! Compute the left hand side block-diagonal mass matrix
//     //! \param[in] geoElem Element geometry array
//     //! \param[in] inpoel Mesh element connectivity
//     //! \param[in] psup Linked lists storing IDs of points surrounding points
//     //! \param[in,out] lhs Left-hand side block-diagonal mass matrix
    void lhs( const tk::Fields& /*geoElem*/,
              tk::Fields& /*lhs*/ ) const
    {
    }

    //! Compute right hand side
//     //! \param[in] t Physical time
//     //! \param[in] geoElem Element geometry array
//     //! \param[in] fd Face connectivity data object
//     //! \param[in] U Solution vector at recent time step
//     //! \param[in,out] R Right-hand side vector computed
    void rhs( tk::real /*t*/,
              const tk::Fields& /*geoElem*/,
              const inciter::FaceData& /*fd*/,
              const tk::Fields& /*U*/,
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

    //! Extract the velocity field at cell nodes
    //! \param[in] U Solution vector at recent time step
    //! \param[in] N Element node indices
    //! \return Array of the four values of the velocity field
    std::array< std::array< tk::real, 4 >, 3 >
    velocity( const tk::Fields& U,
              const std::array< std::vector< tk::real >, 3 >&,
              const std::array< std::size_t, 4 >& N ) const
    {
      std::array< std::array< tk::real, 4 >, 3 > v;
      v[0] = U.extract( 1, m_offset, N );
      v[1] = U.extract( 2, m_offset, N );
      v[2] = U.extract( 3, m_offset, N );
      auto r = U.extract( 0, m_offset, N );
      std::transform( r.begin(), r.end(), v[0].begin(), v[0].begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      std::transform( r.begin(), r.end(), v[1].begin(), v[1].begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      std::transform( r.begin(), r.end(), v[2].begin(), v[2].begin(),
                      []( tk::real s, tk::real& d ){ return d /= s; } );
      return v;
    }

    //! \brief Query all side set IDs the user has configured for all components
    //!   in this PDE system
    //! \param[in,out] conf Set of unique side set IDs to add to
    void side( std::unordered_set< int >& conf ) const
    { Problem::side( conf ); }

    //! Return field names to be output to file
    //! \return Vector of strings labelling fields output in file
    std::vector< std::string > fieldNames() const
    { return Problem::fieldNames(); }

    //! Return field output going to file
    //! \param[in] t Physical time
    //! \param[in] V Total mesh volume
    //! \param[in] geoElem Element geometry array
    //! \param[in,out] U Solution vector at recent time step
    //! \return Vector of vectors to be output to file
    std::vector< std::vector< tk::real > >
    fieldOutput( tk::real t,
                 tk::real V,
                 const tk::Fields& geoElem,
                 tk::Fields& U ) const
    {
      std::array< std::vector< tk::real >, 3 > coord;
      std::vector< tk::real > v;
      v        = geoElem.extract(0,0);
      coord[0] = geoElem.extract(1,0);
      coord[1] = geoElem.extract(2,0);
      coord[2] = geoElem.extract(3,0);

      return Problem::fieldOutput( 0, m_offset, t, V, v, coord, U );
    }

    //! Return names of integral variables to be output to diagnostics file
    //! \return Vector of strings labelling integral variables output
    std::vector< std::string > names() const
    { return Problem::names(); }

  private:
    const ncomp_t m_offset;             //!< Offset PDE operates from
};

} // dg::

} // inciter::

#endif // DGCompFlow_h
