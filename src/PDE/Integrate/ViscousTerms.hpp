// *****************************************************************************
/*!
  \file      src/PDE/Integrate/ViscousTerms.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions for computing viscous surface terms of a system
     of PDEs in DG methods
  \details   This file contains functionality for computing viscous surface
     terms of a system of PDEs used in discontinuous Galerkin methods for
     various orders of numerical representation.
*/
// *****************************************************************************
#ifndef ViscousTerms_h
#define ViscousTerms_h

#include <array>
#include <vector>

#include "Fields.hpp"
#include "EoS/EOS.hpp"
#include "PUPUtil.hpp"

namespace tk {

//! Viscous terms for multispecies PDEs using the modified-gradient approach
class MultiSpeciesViscousTermsP0P1 {
  public:
    //! Constructor
    MultiSpeciesViscousTermsP0P1( std::size_t nspec, std::size_t rdof );

    //! Return local polynomial order
    std::size_t localDof( std::size_t ) const { return m_rdof; }

    //! Reconstruct conserved variables and primitives at a face point
    std::vector< tk::real >
    stateAt(
      const std::vector< inciter::EOS >& mat_blk,
      const Fields& U,
      const Fields& P,
      std::size_t e,
      std::size_t ndof,
      const std::vector< tk::real >& B ) const;

    //! Compute gradients of quantities for an interior element
    void
    gradientIntElem(
      const Fields& U,
      const Fields& P,
      std::size_t elem,
      const std::array< std::vector< tk::real >, 3 >& dBdx,
      std::array< std::array< tk::real, 3 >, 4>& grad ) const;

    //! Compute the multispecies viscous flux at an interior face
    void
    interiorFlux(
      const std::vector< inciter::EOS >& mat_blk,
      std::size_t ncomp,
      const std::array< std::vector< tk::real >, 2 >& state,
      const std::array< std::vector< tk::real >, 2 >& cellAvgState,
      const std::array< tk::real, 3 >& fn,
      const std::array< std::array< tk::real, 3 >, 2 >& centroid,
      const std::array< std::array< std::array< tk::real, 3 >, 4>, 2 >& grad,
      std::vector< tk::real >& fl ) const;

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er& p ) {
      p | m_nspec;
      p | m_rdof;
    }

    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] v Viscous terms object reference
    friend void operator|( PUP::er& p, MultiSpeciesViscousTermsP0P1& v )
    { v.pup(p); }
    //@}

  private:
    std::size_t m_nspec;
    std::size_t m_rdof;
};


//! Viscous terms for multispecies PDEs using high-order DG using the DDG approach
class MultiSpeciesViscousTermsDGP1 {
  public:
    //! Constructor
    MultiSpeciesViscousTermsDGP1( std::size_t nspec, std::size_t rdof );

    //! Return local polynomial order
    std::size_t localDof( std::size_t ) const { return m_rdof; }

    //! DDG Functions will go here

    //! Reconstruct conserved variables and primitives at a face point
    std::vector< tk::real >
    stateAt(
      const std::vector< inciter::EOS >& mat_blk,
      const Fields& U,
      const Fields& P,
      std::size_t e,
      std::size_t ndof,
      const std::vector< tk::real >& B ) const;

    //! Compute gradients of quantities for an interior element
    void
    gradientIntElem(
      const Fields& U,
      const Fields& P,
      std::size_t elem,
      const std::array< std::vector< tk::real >, 3 >& dBdx,
      std::array< std::array< tk::real, 3 >, 4>& grad ) const;

    //! Compute the multispecies viscous flux at an interior face
    void
    interiorFlux(
      const std::vector< inciter::EOS >& mat_blk,
      std::size_t ncomp,
      const std::array< std::vector< tk::real >, 2 >& state,
      const std::array< tk::real, 3 >& fn,
      const std::array< std::array< std::array< tk::real, 3 >, 4>, 2 >& grad,
      std::vector< tk::real >& fl ) const;

    //! Compute the viscous volume flux on a tet element
    void
    volumeFlux(
    const std::vector< inciter::EOS >& mat_blk,
    std::size_t ncomp,
    const std::array< std::vector< tk::real >, 2 >& state,
    const std::vector< std::array< tk::real, 3 > >& visc_fl) const;

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er& p ) {
      p | m_nspec;
      p | m_rdof;
    }

    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] v Viscous terms object reference
    friend void operator|( PUP::er& p, MultiSpeciesViscousTermsDGP1& v )
    { v.pup(p); }
    //@}

  private:
    std::size_t m_nspec;
    std::size_t m_rdof;
};

} // tk::

#endif // ViscousTerms_h
