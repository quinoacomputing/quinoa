//******************************************************************************
/*!
  \file      src/LinSys/LinSysMerger.C
  \author    J. Bakosi
  \date      Tue 03 May 2016 11:56:55 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Linear system merger
  \details   Linear system merger.
*/
//******************************************************************************

#include "LinSysMerger.h"

#include "NoWarning/performer.decl.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wweak-template-vtables"
#endif

namespace tk {

// Some compilers (e.g., GNU and Intel) do not find some of the SDAG code
// generated for Charm++ entry methods defined entirely inside .ci files, such
// as LinSysnMerger<>::wait4init(), thus we must explicitly spell out all
// possible instantiations of LinSysMerger here, similar to that in the .ci file
// to instantiate registration and delivery of code for the individual
// specializations. See also
// https://isocpp.org/wiki/faq/templates#separate-template-class-defn-from-decl.
template class LinSysMerger< inciter::CProxy_Conductor,
                             inciter::CProxy_Performer >;

} // tk::

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
#endif

#include "linsysmerger.def.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif
