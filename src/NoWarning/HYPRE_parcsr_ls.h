//******************************************************************************
/*!
  \file      src/NoWarning/HYPRE_parcsr_ls.h
  \author    J. Bakosi
  \date      Thu 05 May 2016 01:42:23 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include HYPRE_parcsr_ls.h with turning off specific compiler
             warnings
*/
//******************************************************************************
#ifndef nowarning_HYPRE_parcsr_ls_h
#define nowarning_HYPRE_parcsr_ls_h

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Weffc++"
#endif

#include <HYPRE_parcsr_ls.h>

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_HYPRE_parcsr_ls_h
