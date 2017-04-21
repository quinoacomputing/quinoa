// *****************************************************************************
/*!
  \file      src/UnitTest/TestArray.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Simple test Charm++ array for testing arrays
  \details   Simple test Charm++ array for testing arrays.
*/
// *****************************************************************************
#ifndef TestArray_h
#define TestArray_h

namespace tut {

//! Charm++ chare array definition for testing arrays
class TestArray : public CBase_TestArray {
  public:
    explicit TestArray() {}
    explicit TestArray( CkMigrateMessage* ) {}
};

} // tut::

#endif // TestArray_h
