// *****************************************************************************
/*!
  \file      src/Inciter/Partitioner.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Charm++ chare partitioner group used to perform mesh partitioning
  \details   Charm++ chare partitioner group used to parform mesh partitioning.
*/
// *****************************************************************************

#include "Partitioner.h"
#include "AuxSolver.h"

#include "NoWarning/carrier.decl.h"
#include "NoWarning/particlewriter.decl.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wweak-template-vtables"
#endif

// Some compilers (e.g., GNU and Intel) do not find some of the SDAG code
// generated for Charm++ entry methods defined entirely inside .ci files, such
// as Partitioner<>::wait4owned(), thus we must explicitly spell out all
// possible instantiations of Partitioner here, similar to that in the .ci file
// to instantiate registration and delivery of code for the individual
// specializations. See also
// https://isocpp.org/wiki/faq/templates#separate-template-class-defn-from-decl.
template class inciter::Partitioner<
                 inciter::CProxy_Transporter,
                 inciter::CProxy_Carrier,
                 tk::CProxy_LinSysMerger< inciter::CProxy_Transporter,
                                          inciter::CProxy_Carrier,
                                          inciter::AuxSolverLumpMassDiff >,
                 tk::CProxy_ParticleWriter< inciter::CProxy_Transporter > >;

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#include "NoWarning/partitioner.def.h"
