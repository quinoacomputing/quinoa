//******************************************************************************
/*!
  \file      src/UnitTest/TestArray.h
  \author    J. Bakosi
  \date      Thu 12 Mar 2015 10:14:36 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Simple test Charm++ array for testing arrays
  \details   Simple test Charm++ array for testing arrays.
*/
//******************************************************************************
#ifndef TestArray_h
#define TestArray_h

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <testarray.decl.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace tut {

//! Charm++ chare array definition for testing arrays
struct TestArray : CBase_TestArray {
  explicit TestArray() {}
  explicit TestArray( CkMigrateMessage* ) {}
};

} // tut::

#endif // TestArray_h
