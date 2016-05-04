//******************************************************************************
/*!
  \file      src/UnitTest/TestArray.h
  \author    J. Bakosi
  \date      Tue 03 May 2016 11:06:48 AM MDT
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     Simple test Charm++ array for testing arrays
  \details   Simple test Charm++ array for testing arrays.
*/
//******************************************************************************
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
