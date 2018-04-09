// *****************************************************************************
/*!
  \file      src/Inciter/ElemDiagnostics.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     ElemDiagnostics class for collecting diagnostics
  \details   ElemDiagnostics class for collecting diagnostics, e.g., residuals,
    and various norms of errors while solving partial differential equations.
*/
// *****************************************************************************
#ifndef ElemDiagnostics_h
#define ElemDiagnostics_h

#include <unordered_set>

#include "Discretization.h"
#include "PUPUtil.h"
#include "Diagnostics.h"

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
                  const tk::Fields& u );

    ///@{
    //! \brief Pack/Unpack serialize member function
    void pup( PUP::er& ) {}
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] d Diagnostics object reference
    friend void operator|( PUP::er& p, ElemDiagnostics& d ) { d.pup(p); }
    //@}
};

} // inciter::

#endif // ElemDiagnostics_h
