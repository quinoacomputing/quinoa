// *****************************************************************************
/*!
  \file      src/Inciter/NodeDiagnostics.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
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
    //! Constructor
    explicit NodeDiagnostics( const Discretization& d );

    //! Configure Charm++ custom reduction types initiated from this class
    static void registerReducers();

    //! Compute diagnostics, e.g., residuals, norms of errors, etc.
    bool compute( Discretization& d, const tk::Fields& u );

    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      p | m_slave;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] d Diagnostics object reference
    friend void operator|( PUP::er& p, NodeDiagnostics& d ) { d.pup(p); }
    //@}

  private:
    //! Slave mesh node local IDs
    //! \details Local IDs of those mesh nodes to which we contribute to but do
    //!   not own. Ownership here is defined by having a lower chare ID than any
    //!   other chare that also contributes to the node.
    std::unordered_set< std::size_t > m_slave;
};

} // inciter::

#endif // NodeDiagnostics_h
