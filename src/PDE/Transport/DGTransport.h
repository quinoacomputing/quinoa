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
    //! \param[in] geoElem Element geometry array
    //! \param[in,out] unk Array of unknowns
    //! \param[in] t Physical time
    void initialize( const tk::Fields& geoElem,
                     tk::Fields& unk,
                     tk::real t ) const
    {
      std::size_t nelem = unk.nunk();

      for (std::size_t e=0; e<nelem; ++e)
      {
        auto xcc = geoElem(e,1,0);
        auto ycc = geoElem(e,2,0);
        auto zcc = geoElem(e,3,0);

        const auto s = Problem::solution( m_c, m_ncomp, xcc, ycc, zcc, t );
        for (ncomp_t c=0; c<m_ncomp; ++c)
          unk(e, c, m_offset)   = s[c];
      }
    }

    //! Compute the left hand side mass matrix
    //! \param[in] geoElem Element geometry array
    //! \param[in,out] lhs Block diagonal mass matrix matrix
    void lhs( const tk::Fields& geoElem,
              tk::Fields& lhs ) const
    {
      std::size_t nelem = geoElem.nunk();

      for (std::size_t e=0; e<nelem; ++e)
      {
        for (ncomp_t c=0; c<m_ncomp; ++c)
          lhs(e, c, m_offset) = geoElem(e,0,0);
      }
    }

    //! Compute right hand side
    //! \param[in] geoFace Face geometry array
    //! \param[in] fd Face connectivity and boundary conditions object
    //! \param[in] U Solution vector at recent time step
    //! \param[in,out] R Right-hand side vector computed
    void rhs( tk::real,
              const tk::Fields& geoFace,
              const inciter::FaceData& fd,
              const tk::Fields& U,
              tk::Fields& R ) const
    {
      Assert( U.nunk() == R.nunk(), "Number of unknowns in solution "
              "vector and right-hand side at recent time step incorrect" );
      Assert( U.nprop() == m_ncomp && R.nprop() == m_ncomp,
              "Number of components in solution and right-hand side vector " 
              "must equal "+ std::to_string(m_ncomp) );

      auto& esuf = fd.Esuf();
      auto& bface = fd.Bface();

      // set rhs to zero
      R.fill(0.0);

      // compute internal surface flux integrals
      for (auto f=fd.Nbfac(); f<fd.Ntfac(); ++f)
      {
        std::size_t el = static_cast< std::size_t >(esuf[2*f]);
        std::size_t er = static_cast< std::size_t >(esuf[2*f+1]);

        auto farea = geoFace(f,0,0);

        std::array< tk::real, 3 > fn {{ geoFace(f,1,0),
                                        geoFace(f,2,0),
                                        geoFace(f,3,0) }};

        auto xc = geoFace(f,4,0);
        auto yc = geoFace(f,5,0);
        auto zc = geoFace(f,6,0);

        auto ul = U.extract(el);
        auto ur = U.extract(er);

        //--- upwind fluxes
        auto flux = upwindFlux(xc, yc, zc, ul, ur, fn);

        for (ncomp_t c=0; c<m_ncomp; ++c)
        {
          R(el, c, m_offset) -= farea * flux[c];
          R(er, c, m_offset) += farea * flux[c];
        }
      }

      // compute boundary surface flux integrals

      // symmetry boundary condition
      auto bc = bface.find(1);

      if (bc != bface.end())
      {
        for (const auto& f : bc->second)
        {
          std::size_t el = static_cast< std::size_t >(esuf[2*f]);

          Assert( esuf[2*f+1] == -1,
                  "outside boundary element not -1" );

          auto farea = geoFace(f,0,0);

          std::array< tk::real, 3 > fn {{ geoFace(f,1,0),
                                          geoFace(f,2,0),
                                          geoFace(f,3,0) }};

          auto xc = geoFace(f,4,0);
          auto yc = geoFace(f,5,0);
          auto zc = geoFace(f,6,0);

          auto ul = U.extract(el);
          auto ur = U.extract(el);

          //--- upwind fluxes
          auto flux = upwindFlux(xc, yc, zc, ul, ur, fn);

          for (ncomp_t c=0; c<m_ncomp; ++c)
          {
            R(el, c, m_offset) -= farea * flux[c];
          }
        }
      }

      // inlet boundary condition
      bc = bface.find(2);

      if (bc != bface.end())
      {
        for (const auto& f : bc->second)
        {
          std::size_t el = static_cast< std::size_t >(esuf[2*f]);

          Assert( esuf[2*f+1] == -1,
                  "outside boundary element not -1" );

          auto farea = geoFace(f,0,0);

          std::array< tk::real, 3 > fn {{ geoFace(f,1,0),
                                          geoFace(f,2,0),
                                          geoFace(f,3,0) }};

          auto xc = geoFace(f,4,0);
          auto yc = geoFace(f,5,0);
          auto zc = geoFace(f,6,0);

          auto ul = U.extract(el);
          std::vector< tk::real > ur(ul.size(),0);

          //--- upwind fluxes
          auto flux = upwindFlux(xc, yc, zc, ul, ur, fn);

          for (ncomp_t c=0; c<m_ncomp; ++c)
          {
            R(el, c, m_offset) -= farea * flux[c];
          }
        }
      }

      // outlet boundary condition
      bc = bface.find(3);

      if (bc != bface.end())
      {
        for (const auto& f : bc->second)
        {
          std::size_t el = static_cast< std::size_t >(esuf[2*f]);

          Assert( esuf[2*f+1] == -1,
                  "outside boundary element not -1" );

          auto farea = geoFace(f,0,0);

          std::array< tk::real, 3 > fn {{ geoFace(f,1,0),
                                          geoFace(f,2,0),
                                          geoFace(f,3,0) }};

          auto xc = geoFace(f,4,0);
          auto yc = geoFace(f,5,0);
          auto zc = geoFace(f,6,0);

          auto ul = U.extract(el);
          auto ur = U.extract(el);

          //--- upwind fluxes
          auto flux = upwindFlux(xc, yc, zc, ul, ur, fn);

          for (ncomp_t c=0; c<m_ncomp; ++c)
          {
            R(el, c, m_offset) -= farea * flux[c];
          }
        }
      }
    }

    std::vector< tk::real >
    upwindFlux( tk::real xc,
                tk::real yc, 
                tk::real zc,
                std::vector< tk::real > ul,
                std::vector< tk::real > ur,
                std::array< tk::real, 3 > fn ) const
    // *****************************************************************************
    // Riemann solver using upwind method
    //! \param[in] xc X coordinate at which to assign advection velocity
    //! \param[in] yc Y coordinate at which to assign advection velocity
    //! \param[in] zc Z coordinate at which to assign advection velocity
    //! \param[in] ul Left unknown/state vector
    //! \param[in] ur Right unknown/state vector
    //! \param[in] fn Face unit normal vector
    //! \return Riemann solution using upwind method
    // *****************************************************************************
    {
        std::vector< tk::real > flux(ul.size(),0);

        const auto vel = Problem::prescribedVelocity( xc, yc, zc, m_c, m_ncomp );
    
        for(ncomp_t c=0; c<m_ncomp; ++c)
        {
          auto ax = vel[c][0];
          auto ay = vel[c][1];
          auto az = vel[c][2];

          // wave speed
          tk::real swave = ax*fn[0] + ay*fn[1] + az*fn[2];
    
          // upwinding
          tk::real splus  = 0.5 * (swave + fabs(swave));
          tk::real sminus = 0.5 * (swave - fabs(swave));
    
          flux[c] = splus * ul[c] + sminus * ur[c];
        }
    
        return flux;
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
    //! \param[in] geoElem Element geometry array
    //! \param[in,out] U Solution vector at recent time step
    //! \return Vector of vectors to be output to file
    //! \details This functions should be written in conjunction with names(),
    //!   which provides the vector of field names
    //! \note U is overwritten
    std::vector< std::vector< tk::real > >
    fieldOutput( tk::real t,
                 tk::real /*V*/,
                 const tk::Fields& geoElem,
                 tk::Fields& U ) const
    {
      std::vector< std::vector< tk::real > > out;
      // will output numerical solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        out.push_back( U.extract( c, m_offset ) );
      // evaluate analytic solution at time t
      auto E = U;
      initialize( geoElem, E, t );
      // will output analytic solution for all components
      for (ncomp_t c=0; c<m_ncomp; ++c)
        out.push_back( E.extract( c, m_offset ) );
      // will output error for all components
      for (ncomp_t c=0; c<m_ncomp; ++c) {
        auto u = U.extract( c, m_offset );
        auto e = E.extract( c, m_offset );
        for (std::size_t i=0; i<u.size(); ++i)
          e[i] = std::pow( e[i] - u[i], 2.0 ) * geoElem(i,0,0);
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
