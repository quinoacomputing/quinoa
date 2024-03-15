// *****************************************************************************
/*!
  \file      src/Inciter/FluxCorrector.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     FluxCorrector performs limiting for transport equations
  \details   FluxCorrector performs limiting for transport equations. Each
    FluxCorrector object performs the limiting procedure, according to a
    flux-corrected transport algorithm, on a chunk of the full load (part of the
    mesh).
*/
// *****************************************************************************
#ifndef FluxCorrector_h
#define FluxCorrector_h

#include <utility>
#include <unordered_map>

#include "NoWarning/pup.hpp"

#include "Keywords.hpp"
#include "Fields.hpp"
#include "Inciter/InputDeck/New2InputDeck.hpp"

namespace inciter {

extern ctr::New2InputDeck g_newinputdeck;

//! FluxCorrector is used to perform flux-corrected transport
//! \see Löhner, R., Morgan, K., Peraire, J. and Vahdati, M. (1987), Finite
//!   element flux-corrected transport (FEM–FCT) for the Euler and Navier–Stokes
//!   equations. Int. J. Numer. Meth. Fluids, 7: 1093–1109.
//!   doi:10.1002/fld.1650071007
class FluxCorrector {

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

  public:
    //! Constructor
    //! \param[in] is Size of the mesh element connectivity vector (inpoel size)
    explicit FluxCorrector( std::size_t is = 0 ) :
      m_aec( is, g_newinputdeck.get< newtag::ncomp >() ),
      m_sys( g_newinputdeck.get< newtag::sysfctvar >() ),
      m_vel( findvel() ) {}

    //! Find components of a velocity for equation systems
    //! \return List of 3 component indices to treat as a velocity
    //! \warning Currently, this is only a punt for single-material flow: we
    //!   simply take the components 1,2,3 as the velocity for each system of
    //!   type Eq
    std::vector< std::array< ncomp_t, 3 > >
    findvel() {
      std::vector< std::array< ncomp_t, 3 > > vel;
      vel.push_back( {1, 2, 3} );
      return vel;
    }

    //! Resize state (e.g., after mesh refinement)
    void resize( std::size_t is ) { m_aec.resize( is ); }

    //! Compute antidiffusive element contributions (AEC)
    void aec(
      const std::array< std::vector< tk::real >, 3 >& coord,
      const std::vector< std::size_t >& inpoel,
      const std::vector< tk::real >& vol,
      const std::unordered_map< std::size_t,
              std::vector< std::pair< bool, tk::real > > >& bcdir,
      const std::unordered_map< int,
        std::unordered_set< std::size_t > >& symbcnodemap,
      const std::unordered_map< int,
        std::unordered_map< std::size_t, std::array< tk::real, 4 > > >& bnorm,
      const tk::Fields& Un,
      tk::Fields& P );

    //! Verify the assembled antidiffusive element contributions
    bool verify( std::size_t nchare,
                 const std::vector< std::size_t >& inpoel,
                 const tk::Fields& dUh,
                 const tk::Fields& dUl ) const;

    //! Compute mass diffusion contribution to the rhs of the low order system
    tk::Fields diff( const std::array< std::vector< tk::real >, 3 >& coord,
                     const std::vector< std::size_t >& inpoel,
                     const tk::Fields& Un ) const;

    //! \brief Compute the maximum and minimum unknowns of all elements
    //!   surrounding nodes
    void alw( const std::vector< std::size_t >& inpoel,
              const tk::Fields& Un,
              const tk::Fields& Ul,
              tk::Fields& Q ) const;

    //! Compute limited antiffusive element contributions and apply to mesh nodes
    void lim( const std::vector< std::size_t >& inpoel,
              const std::unordered_map< std::size_t,
                      std::vector< std::pair< bool, tk::real > > >& bcdir,
              const tk::Fields& P,
              const tk::Fields& Ul,
              tk::Fields& Q,
              tk::Fields& A ) const;

    // Collect mesh output fields from FCT
    std::tuple< std::vector< std::string >,
                std::vector< std::vector< tk::real > > >
    fields( const std::vector< std::size_t >& inpoel ) const;

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er& p ) {
      p | m_aec;
      p | m_sys;
      p | m_vel;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i FluxCorrector object reference
    friend void operator|( PUP::er& p, FluxCorrector& i ) { i.pup(p); }
    //@}

  private:
   //! Antidiffusive element contributions for all scalar components
   tk::Fields m_aec;
   //! Component indices to treat as a system
   std::vector< ncomp_t > m_sys;
   //! Component indices to treat as a velocity vector for multiple systems
   std::vector< std::array< ncomp_t, 3 > > m_vel;
};

} // inciter::

#endif // FluxCorrector_h
