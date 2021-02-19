// *****************************************************************************
/*!
  \file      src/UnitTest/TUTConfig.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Template Unit Test unit test configuration for all tests
  \details   Template Unit Test unit test configuration for all tests.
*/
// *****************************************************************************
#ifndef TUTConfig_h
#define TUTConfig_h

namespace tut {

//! \brief Maximum number of tests in every test group to attempt to run
//! \details If any of the unit test groups have more tests than this number,
//!   this should be increased. All test groups use this value to override the
//!   default template argument for tut::test_group<>.
const int MAX_TESTS_IN_GROUP = 80;

} // tut::

#endif // TUTConfig_h
