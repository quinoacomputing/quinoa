// *****************************************************************************
/*!
  \file      src/Inciter/FluxCorrector.h
  \author    J. Bakosi
  \date      Fri 02 Sep 2016 09:56:07 AM MDT
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

namespace inciter {

//! FluxCorrector Charm++ chare array used to perform flux-corrected transport
//! \see Löhner, R., Morgan, K., Peraire, J. and Vahdati, M. (1987), Finite
//!   element flux-corrected transport (FEM–FCT) for the Euler and Navier–Stokes
//!   equations. Int. J. Numer. Meth. Fluids, 7: 1093–1109.
//!   doi:10.1002/fld.1650071007
class FluxCorrector {

  public:
    //! Constructor
    explicit FluxCorrector(
      std::pair< std::vector< std::size_t >,
                 std::vector< std::size_t > >&& esup = {} );

    //! Const-ref accessor to elements surrounding points
    const std::pair< std::vector< std::size_t >,
                     std::vector< std::size_t > >& esup() const
    { return m_esup; }

    //! Compute antidiffusive element contributions
    tk::MeshNodes
    aec( const std::array< std::vector< tk::real >, 3 >& coord,
         const std::vector< std::size_t >& inpoel,
         const tk::MeshNodes& u,
         const tk::MeshNodes& uh );

    //! Compute mass diffusion to be added to the low order solution rhs
    std::pair< tk::MeshNodes, tk::MeshNodes >
    low( const std::array< std::vector< tk::real >, 3 >& coord,
         const std::vector< std::size_t >& inpoel,
         const tk::MeshNodes& Un );

    ///@{
    //! \brief Pack/Unpack serialize member function
    void pup( PUP::er& p ) {
      p | m_esup;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i FluxCorrector object reference
    //! \author J. Bakosi
    friend void operator|( PUP::er& p, FluxCorrector& i ) { i.pup(p); }
    //@}

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

    //! Elements surrounding points of the mesh chunk we operate on
    std::pair< std::vector< std::size_t >,
               std::vector< std::size_t > > m_esup;
};

} // inciter::

#endif // FluxCorrector_h
