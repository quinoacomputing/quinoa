// *****************************************************************************
/*!
  \file      src/Inciter/Diagnostics.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
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
const std::size_t NUMDIAG = 9;

//! Diagnostics labels
enum Diag { L2SOL=0,    //!< L2 norm of numerical solution
            L2ERR,      //!< L2 norm of numerical-analytic solution
            L2RES,      //!< L2 norm of the residual
            LINFERR,    //!< L_inf norm of numerical-analytic solution
            TOTALSOL,   //!< Sum of conserved solution over entire domain
            RESFORCE,   //!< Resultant force vector on mesh boundaries
            ITER,       //!< Iteration count
            TIME,       //!< Physical time
            DT };       //!< Time step size

} // inciter::

#endif // Diagnostics_h
