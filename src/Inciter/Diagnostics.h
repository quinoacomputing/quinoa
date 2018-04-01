// *****************************************************************************
/*!
  \file      src/Inciter/Diagnostics.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Common data for collecting diagnostics
  \details   Common data for collecting (node-, elem-, etc) diagnostics, e.g.,
    residuals, and various norms of errors while solving partial differential
    equations.
*/
// *****************************************************************************
#ifndef Diagnostics_h
#define Diagnostics_h

namespace inciter {

//! Number of entries in diagnostics vector (of vectors)
const std::size_t NUMDIAG = 6;

//! Diagnostics labels
enum Diag { L2SOL=0,    //!< L2 norm of numerical solution
            L2ERR,      //!< L2 norm of numerical-analytic solution
            LINFERR,    //!< L_inf norm of numerical-analytic solution
            ITER,       //!< Iteration count
            TIME,       //!< Physical time
            DT };       //!< Time step size

} // inciter::

#endif // Diagnostics_h
