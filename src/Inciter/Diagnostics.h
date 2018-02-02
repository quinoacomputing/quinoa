// *****************************************************************************
/*!
  \file      src/Inciter/Diagnostics.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Diagnostics class for collecting diagnostics
  \details   Diagnostics class for collecting diagnostics, e.g., residuals, and
    various norms of errors while solving partial differential equations.
*/
// *****************************************************************************
#ifndef Diagnostics_h
#define Diagnostics_h

#include <unordered_set>

#include "Discretization.h"
#include "PUPUtil.h"

namespace inciter {

//! Number of entries in diagnostics vector (of vectors), see compute()
const std::size_t NUMDIAG = 6;

//! Diagnostics labels
enum Diag { L2SOL=0,    //!< L2 norm of numerical solution
            L2ERR,      //!< L2 norm of numerical-analytic solution
            LINFERR,    //!< L_inf norm of numerical-analytic solution
            ITER,       //!< Iteration count
            TIME,       //!< Physical time
            DT };       //!< Time step size

//! Diagnostics class used to compute diagnostics while integrating PDEs
class Diagnostics {

  public:
    //! Constructor
    explicit Diagnostics( const Discretization& d );

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
    friend void operator|( PUP::er& p, Diagnostics& d ) { d.pup(p); }
    //@}

  private:
    //! Slave mesh node local IDs
    //! \details Local IDs of those mesh nodes to which we contribute to but do
    //!   not own. Ownership here is defined by having a lower chare ID than any
    //!   other chare that also contributes to the node.
    std::unordered_set< std::size_t > m_slave;
};

} // inciter::

#endif // Diagnostics_h
