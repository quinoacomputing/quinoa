// *****************************************************************************
/*!
  \file      src/Inciter/ElemDiagnostics.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     ElemDiagnostics class for collecting diagnostics
  \details   ElemDiagnostics class for collecting diagnostics, e.g., residuals,
    and various norms of errors while solving partial differential equations.
*/
// *****************************************************************************
#ifndef ElemDiagnostics_h
#define ElemDiagnostics_h

#include "Discretization.hpp"
#include "PUPUtil.hpp"
#include "Diagnostics.hpp"

namespace inciter {

//! ElemDiagnostics class used to compute diagnostics while integrating PDEs
class ElemDiagnostics {

  public:
    //! Configure Charm++ custom reduction types initiated from this class
    static void registerReducers();

    //! Compute diagnostics, e.g., residuals, norms of errors, etc.
    bool compute( Discretization& d,
                  const std::size_t nchGhost,
                  const tk::Fields& geoElem,
                  const std::vector< std::size_t >& ndofel,
                  const tk::Fields& u,
                  const tk::Fields& un ) const;

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    void pup( PUP::er& ) {}
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] d Diagnostics object reference
    friend void operator|( PUP::er& p, ElemDiagnostics& d ) { d.pup(p); }
    //@}

  private:
    //! Compute diagnostics for DG
    void compute_diag( const Discretization& d,
                       const std::size_t ndofel,
                       const std::size_t nchGhost,
                       const tk::Fields& geoElem,
                       const std::vector< std::size_t >& pIndex,
                       const tk::Fields& u,
                       const tk::Fields& un,
                       std::vector< std::vector< tk::real > >& diag ) const;
};

} // inciter::

#endif // ElemDiagnostics_h
