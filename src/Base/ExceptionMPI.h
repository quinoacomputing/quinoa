// *****************************************************************************
/*!
  \file      src/Base/ExceptionMPI.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Exception macros interoperating with MPI
  \details   Exception macros interoperating with MPI.
*/
// *****************************************************************************
#ifndef ExceptionMPI_h
#define ExceptionMPI_h

#include "NoWarning/mpi.h"

#include "Exception.h"

namespace tk {

//! \brief Assert macro that only throws an exception if expr fails.
//! \details If NDEBUG is defined (e.g. cmake's RELEASE or OPTIMIZED mode), do
//!    nothing, expr is not evaluated. If NDEBUG is not defined, evaluate expr.
//!    If expr is true, do nothing. If expr is false, throw Exception with
//!    arguments passed in. The behavior is similar to libc's assert macro, but
//!    throwing an Exception instead will also generate a nice call-trace and
//!    will attempt to free memory. This macro should be used to detect
//!    programmer errors.
//! \author J. Bakosi
#ifdef NDEBUG
#  define AssertMPI(expr, ...) (static_cast<void>(0))
#else  // NDEBUG
#  define AssertMPI(expr, ...) \
   ((expr) ? static_cast<void>(0) : Throw(__VA_ARGS__))
#endif // NDEBUG

//! \brief ErrChkMPI macro that only throws an exception if expr fails.
//! \details The behavior of this macro is the same whether NDEBUG is defined or
//!    not: expr is always evaluated. If expr is true, do nothing. If expr is
//!    false, throw Exception with arguments passed in. This macro should be
//!    used to detect user or runtime errors. The main difference compared to
//!    the vanilla ErrChk macro, defined in Base/Exception.h, is that this macro
//!    always does an MPI_Allreduce on the outcome of expr. Since MPI_Allreduce
//!    combines values of err from all MPI ranks and distributes the result back
//!    to all ranks, the result of expr (err) appear on all ranks, thus if any
//!    of the ranks had expr false, all will throw. Thus this macro should be
//!    used to detect user or runtime errors from code sections that are
//!    executed on multiple MPI ranks. Within asynchrounous Charm++ chares, the
//!    simpler ErrChk macro suffices.
//! \author J. Bakosi
#ifdef QUINOA_CONFIG_MPI_ENABLED
#define ErrChkMPI(expr, ...) \
{ \
  int err = (expr) ? 0 : 1; \
  int globalerr; \
  MPI_Allreduce( &err, &globalerr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD ); \
  if ( globalerr ) Throw( __VA_ARGS__); \
}
#else
#define ErrChkMPI(expr, ...)\
{\
  int QUINOA_ErrChkMPI_err = (expr) ? 0 : 1; \
  if ( QUINOA_ErrChkMPI_err ) Throw( __VA_ARGS__); \
}
#endif

} // tk::

#endif // ExceptionMPI_h
