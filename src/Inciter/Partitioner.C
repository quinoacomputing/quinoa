//******************************************************************************
/*!
  \file      src/Inciter/Partitioner.C
  \author    J. Bakosi
  \date      Tue 19 Jan 2016 03:44:23 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Charm++ chare partitioner group used to perform mesh partitioning
  \details   Charm++ chare partitioner group used to parform mesh partitioning.
*/
//******************************************************************************

#include "Partitioner.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "conductor.decl.h"
#include "performer.decl.h"
#include "linsysmerger.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace inciter {

//! \brief Charm++ global mesh node IDs merger reducer
//! \details This variable is defined here in the .C file and declared as extern
//!   in Partitioner.h. If instead one defines it in the header (as static),
//!   a new version of the variable is created any time the header file is
//!   included, yielding no compilation nor linking errors. However, that leads
//!   to runtime errors, since Partitioner::registerNodesMerger(), a Charm++
//!   "initnode" entry method, *may* fill one while contribute() may use the
//!   other (unregistered) one. Result: undefined behavior, segfault, and
//!   formatting the internet ...
CkReduction::reducerType NodesMerger;

}

// Some compilers (e.g., GNU and Intel) do not find some of the SDAG code
// generated for Charm++ entry methods defined entirely inside .ci files, such
// as Partitioner<>::wait4owned(), thus we must explicitly spell out all
// possible instantiations of Partitioner here, similar to that in the .ci file
// to instantiate registration and delivery of code for the individual
// specializations. See also
// https://isocpp.org/wiki/faq/templates#separate-template-class-defn-from-decl.
template class inciter::Partitioner<
                 inciter::CProxy_Conductor,
                 inciter::CProxy_Performer,
                 tk::CProxy_LinSysMerger< inciter::CProxy_Conductor,
                                          inciter::CProxy_Performer > >;

#include "partitioner.def.h"
