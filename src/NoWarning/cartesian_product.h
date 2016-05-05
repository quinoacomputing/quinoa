//******************************************************************************
/*!
  \file      src/NoWarning/cartesian_product.h
  \author    J. Bakosi
  \date      Mon 02 May 2016 08:01:26 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include cartesian_product.h with turning off specific compiler
             warnings
*/
//******************************************************************************
#ifndef nowarning_cartesian_product_h
#define nowarning_cartesian_product_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wold-style-cast"
#endif

#include <boost/mpl/cartesian_product.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_cartesian_product_h
