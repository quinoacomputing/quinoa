// *****************************************************************************
/*!
  \file      src/Inciter/FluxCorrector.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
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

#include "NoWarning/pup.h"

#include "Keywords.h"
#include "Fields.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! FluxCorrector is used to perform flux-corrected transport
//! \see Löhner, R., Morgan, K., Peraire, J. and Vahdati, M. (1987), Finite
//!   element flux-corrected transport (FEM–FCT) for the Euler and Navier–Stokes
//!   equations. Int. J. Numer. Meth. Fluids, 7: 1093–1109.
//!   doi:10.1002/fld.1650071007
class FluxCorrector {

  public:
    //! Constructor
    //! \param[in] is Size of the mesh element connectivity vector (inpoel size)
    explicit FluxCorrector( std::size_t is = 0 ) :
      m_aec( is, g_inputdeck.get< tag::component >().nprop() ) {}

    //! Compute antidiffusive element contributions (AEC)
    void aec( const std::array< std::vector< tk::real >, 3 >& coord,
              const std::vector< std::size_t >& inpoel,
              const std::vector< tk::real >& vol,
              const std::unordered_map< std::size_t,
                      std::vector< std::pair< bool, tk::real > > >& bc,
              const std::vector< std::size_t >& gid,
              const tk::Fields& dUh,
              const tk::Fields& Un,
              tk::Fields& P );

    //! Verify the assembled antidiffusive element contributions
    bool verify( std::size_t nchare,
                 const std::vector< std::size_t >& inpoel,
                 const tk::Fields& dUh,
                 const tk::Fields& dUl ) const;

    //! Compute lumped mass matrix lhs for low order system
    tk::Fields lump( const std::array< std::vector< tk::real >, 3 >& coord,
                     const std::vector< std::size_t >& inpoel ) const;

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
              const tk::Fields& P,
              const tk::Fields& Ul,
              tk::Fields& Q,
              tk::Fields& A ) const;

    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er& p ) {
      p | m_aec;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i FluxCorrector object reference
    //! \author J. Bakosi
    friend void operator|( PUP::er& p, FluxCorrector& i ) { i.pup(p); }
    //@}

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

   //! Antidiffusive element contributions for all scalar components
   tk::Fields m_aec;
};

} // inciter::

#endif // FluxCorrector_h
