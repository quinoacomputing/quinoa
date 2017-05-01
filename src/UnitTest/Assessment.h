// *****************************************************************************
/*!
  \file      src/UnitTest/Assessment.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit test suite assessment
  \details   Unit test suite assessment.
*/
// *****************************************************************************
#ifndef Assessment_h
#define Assessment_h

#include <vector>
#include <cstddef>
#include <iosfwd>

namespace tk { class Print; }

namespace unittest {

//! Evaluate a single unit test
void evaluate( std::vector< std::string > status,
               std::size_t& ncomplete,
               std::size_t& nwarn,
               std::size_t& nskip,
               std::size_t& nexcp,
               std::size_t& nfail );

//! Echo final assessment after the full unit test suite has finished
bool assess( const tk::Print& print,
             std::string suite,
             std::size_t nfail,
             std::size_t nwarn,
             std::size_t nskip,
             std::size_t nexcp,
             std::size_t ncomplete );

} // unittest::

#endif // Assessment_h
