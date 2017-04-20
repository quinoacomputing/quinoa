// *****************************************************************************
/*!
  \file      src/NoWarning/Zoltan2_PartitioningProblem.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include Zoltan2_PartitioningProblem.hpp with turning off specific
             compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_Zoltan2_PartitioningProblem_h
#define nowarning_Zoltan2_PartitioningProblem_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wunused-parameter"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wsign-compare"
  #pragma clang diagnostic ignored "-Wdocumentation"
  #pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wunused-exception-parameter"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
  #pragma clang diagnostic ignored "-Wswitch-enum"
  #pragma clang diagnostic ignored "-Wdouble-promotion"
  #pragma clang diagnostic ignored "-Wfloat-equal"
  #pragma clang diagnostic ignored "-Wconversion"
  #pragma clang diagnostic ignored "-Wfloat-conversion"
  #pragma clang diagnostic ignored "-Wdeprecated"
  #pragma clang diagnostic ignored "-Wundefined-func-template"
  #pragma clang diagnostic ignored "-Wcast-align"
  #pragma clang diagnostic ignored "-Wconditional-uninitialized"
  #pragma clang diagnostic ignored "-Wunused-local-typedef"
  #pragma clang diagnostic ignored "-Wundef"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-parameter"
  #pragma GCC diagnostic ignored "-Wcast-qual"
  #pragma GCC diagnostic ignored "-Wfloat-equal"
  #pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
  #pragma GCC diagnostic ignored "-Wswitch-default"
  #pragma GCC diagnostic ignored "-Wunused-local-typedefs"
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
  #pragma GCC diagnostic ignored "-Wshadow"
#endif

#include <Zoltan2_PartitioningProblem.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_Zoltan2_PartitioningProblem_h
