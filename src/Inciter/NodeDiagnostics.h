// *****************************************************************************
/*!
  \file      src/Inciter/NodeDiagnostics.h
  \copyright 2012-2015, J. Bakosi, 2016-2019, Triad National Security, LLC.
  \brief     NodeDiagnostics class for collecting diagnostics
  \details   NodeDiagnostics class for collecting diagnostics, e.g., residuals,
    and various norms of errors while solving partial differential equations.
*/
// *****************************************************************************
#ifndef NodeDiagnostics_h
#define NodeDiagnostics_h

#include <unordered_set>

#include "Discretization.h"
#include "PUPUtil.h"
#include "Diagnostics.h"

namespace inciter {

//! NodeDiagnostics class used to compute diagnostics while integrating PDEs
class NodeDiagnostics {

  public:
    //! Configure Charm++ custom reduction types initiated from this class
    static void registerReducers();

    //! Compute diagnostics, e.g., residuals, norms of errors, etc.
    bool compute( Discretization& d, const tk::Fields& u ) const;

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    void pup( PUP::er & ) {}
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] d Diagnostics object reference
    friend void operator|( PUP::er& p, NodeDiagnostics& d ) { d.pup(p); }
    //@}
};

} // inciter::

#endif // NodeDiagnostics_h
