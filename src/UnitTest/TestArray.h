// *****************************************************************************
/*!
  \file      src/UnitTest/TestArray.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Simple test Charm++ array for testing arrays
  \details   Simple test Charm++ array for testing arrays.
*/
// *****************************************************************************
#ifndef TestArray_h
#define TestArray_h

#include "NoWarning/testarray.decl.h"

namespace tut {

//! Charm++ chare array definition for testing arrays
class TestArray : public CBase_TestArray {};

} // tut::

#endif // TestArray_h
