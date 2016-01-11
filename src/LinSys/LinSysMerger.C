//******************************************************************************
/*!
  \file      src/LinSys/LinSysMerger.C
  \author    J. Bakosi
  \date      Fri 08 Jan 2016 06:10:42 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Linear system merger
  \details   Linear system merger.
*/
//******************************************************************************

#include "LinSysMerger.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "performer.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

// Some compilers (e.g., GNU and Intel) do not find some of the SDAG code
// generated for Charm++ entry methods defined entirely inside .ci files, such
// as LinSysnMerger<>::wait4init(), thus we must explicitly spell out all
// possible instantiations of LinSysMerger here, similar to that in the .ci file
// to instantiate registration and delivery of code for the individual
// specializations. See also
// https://isocpp.org/wiki/faq/templates#separate-template-class-defn-from-decl.
template class tk::LinSysMerger< inciter::CProxy_Conductor,
                                 inciter::CProxy_Performer >;

#include "linsysmerger.def.h"
