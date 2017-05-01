// *****************************************************************************
/*!
  \file      src/LinSys/LinSysMerger.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Linear system merger
  \details   Linear system merger.
*/
// *****************************************************************************

#include "LinSysMerger.h"

#include "NoWarning/carrier.decl.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wweak-template-vtables"
#endif

namespace tk {

//! \brief Charm++ reducers used by LinSysMerger
//! \details These variables are defined here in the .C file and declared as
//!   extern in LinSysMerger.h. If instead one defines them in the header (as
//!   static), a new version of any of these variables is created any time the
//!   header file is included, yielding no compilation nor linking errors.
//!   However, that leads to runtime errors, since
//!   LinSysMerger::registerBCMerger(), a Charm++ "initnode" entry method, *may*
//!   fill one while contribute() may use the other (unregistered) one. Result:
//!   undefined behavior, segfault, and formatting the internet ...
CkReduction::reducerType BCVectorMerger;
CkReduction::reducerType BCMapMerger;
CkReduction::reducerType BCValMerger;

}

// Some compilers (e.g., GNU and Intel) do not find some of the SDAG code
// generated for Charm++ entry methods defined entirely inside .ci files, such
// as LinSysnMerger<>::wait4sol(), thus we must explicitly spell out all
// possible instantiations of LinSysMerger here, similar to that in the .ci file
// to instantiate registration and delivery of code for the individual
// specializations. See also
// https://isocpp.org/wiki/faq/templates#separate-template-class-defn-from-decl.
template class tk::LinSysMerger< inciter::CProxy_Transporter,
                                 inciter::CProxy_Carrier,
                                 inciter::AuxSolverLumpMassDiff >;

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wunused-parameter"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
  #pragma clang diagnostic ignored "-Wreorder"
  #pragma clang diagnostic ignored "-Wundef"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wswitch-default"
#endif

#include "linsysmerger.def.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif
