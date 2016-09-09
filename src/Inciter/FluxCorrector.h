// *****************************************************************************
/*!
  \file      src/Inciter/FluxCorrector.h
  \author    J. Bakosi
  \date      Thu 08 Sep 2016 07:27:39 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     FluxCorrector performs limiting for transport equations
  \details   FluxCorrector performs limiting for transport equations. There is a
    FluxCorrector Charm++ array element bound to each Carrier array element..
    Each FluxCorrector object performs the limiting procedure, according to a
    flux-corrected transport algorithm, on a chunk of the full load (part of the
    mesh).
*/
// *****************************************************************************
#ifndef FluxCorrector_h
#define FluxCorrector_h

#include <utility>

#include "NoWarning/pup.h"

#include "Keywords.h"
#include "MeshNodes.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! FluxCorrector Charm++ chare array used to perform flux-corrected transport
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
    tk::MeshNodes
    aec( const std::array< std::vector< tk::real >, 3 >& coord,
         const std::vector< std::size_t >& inpoel,
         const tk::MeshNodes& Un,
         const tk::MeshNodes& Uh );

    //! Compute lumped mass matrix lhs for low order system
    tk::MeshNodes
    lump( const std::array< std::vector< tk::real >, 3 >& coord,
          const std::vector< std::size_t >& inpoel ) const;

    //! Compute mass diffusion contribution to the rhs of the low order system
    tk::MeshNodes
    diff( const std::array< std::vector< tk::real >, 3 >& coord,
          const std::vector< std::size_t >& inpoel,
          const tk::MeshNodes& Un ) const;

    //! Compute the allowed solution increments and decrements at mesh nodes
    tk::MeshNodes
    allowed( const std::vector< std::size_t >& inpoel,
             const tk::MeshNodes& Un,
             const tk::MeshNodes& Ul ) const;

    //! Perform limiting
    void limit( const std::vector< std::size_t >& inpoel,
                const tk::MeshNodes& P,
                tk::MeshNodes& Q,
                tk::MeshNodes& U ) const;

    ///@{
    //! \brief Pack/Unpack serialize member function
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
   tk::MeshNodes m_aec;

};

} // inciter::

#endif // FluxCorrector_h
