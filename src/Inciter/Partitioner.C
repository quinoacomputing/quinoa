//******************************************************************************
/*!
  \file      src/Inciter/Partitioner.C
  \author    J. Bakosi
  \date      Mon 21 Dec 2015 09:20:32 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Charm++ chare partitioner group used to perform mesh partitioning
  \details   Charm++ chare partitioner group used to parform mesh partitioning.
*/
//******************************************************************************

#include "Partitioner.h"

// Some compilers (e.g., GNU and Intel) do not find some of the SDAG code
// generated for Charm++ entry methods defined entirely inside .ci files, such
// as Partitioner<>::wait4owned(), thus we must explicitly spell out all
// possible instantiations of Partitioner here, similar to that in the .ci file
// to instantiate registration and delivery of code for the individual
// specializations. See also
// https://isocpp.org/wiki/faq/templates#separate-template-class-defn-from-decl.
template class inciter::Partitioner< inciter::CProxy_Conductor >;

#include "partitioner.def.h"
