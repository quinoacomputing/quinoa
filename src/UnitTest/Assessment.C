// *****************************************************************************
/*!
  \file      src/UnitTest/Assessment.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit test suite assessment
  \details   Unit test suite assessment.
*/
// *****************************************************************************

#include <string>

#include "NoWarning/format.h"

#include "Print.h"
#include "Assessment.h"

namespace unittest {

void
evaluate( std::vector< std::string > status,
          std::size_t& ncomplete,
          std::size_t& nwarn,
          std::size_t& nskip,
          std::size_t& nexcp,
          std::size_t& nfail )
// *****************************************************************************
//  Evaluate a single unit test
//! \param[in] status Vector of strings containing the test results. See
//!   unittest::TUTTest constructor for the expected structure of status.
//! \param[inout] ncomplete Number of completed tests
//! \param[inout] nwarn Number of tests with a warning
//! \param[inout] nskip Number of skipped tests
//! \param[inout] nexcp Number of tests with an exception
//! \param[inout] nfail Number of failed tests
//! \author J. Bakosi
// *****************************************************************************
{
  if (status[2] != "8") {             // only care about non-dummy tests
    ++ncomplete;                      // count number of tests completed
    if (status[2] == "3")             // count number of tests with a warning
      ++nwarn;
    else if (status[2] == "7")        // count number of skipped tests
      ++nskip;
    else if (status[2] == "2")        // count number of tests throwing
      ++nexcp;
    else if (status[2] != "0")        // count number of failed tests
      ++nfail;
  }
}

bool
assess( const tk::Print& print,
        std::string suite,
        std::size_t nfail,
        std::size_t nwarn,
        std::size_t nskip,
        std::size_t nexcp,
        std::size_t ncomplete )
// *****************************************************************************
// Echo final assessment after the full unit test suite has finished
//! \param[in] print Pretty printer
//! \param[in] suite Name of the test suite
//! \param[in] nfail Number of failed tests
//! \param[in] nwarn Number of tests with a warning
//! \param[in] nskip Number of skipped tests
//! \param[in] nexcp Number of tests with an exception
//! \param[in] ncomplete Number of completed tests
//! \return True of all tests passed, false if there was at least a failure or
//!   an exception
//! \author J. Bakosi
// *****************************************************************************
{
  if (!nfail && !nwarn && !nskip && !nexcp) {
    print.note< tk::QUIET >
      ( "All " + std::to_string(ncomplete) + " " + suite + " tests passed" );
  } else {
    std::string skip, warn, fail, excp;
    if (nwarn) warn = "finished with a warning: " + std::to_string(nwarn);
    if (nskip) skip = std::string(nwarn ? ", " : "") +
                      "skipped: " + std::to_string(nskip);
    if (nexcp) excp = std::string(nskip || nwarn ? ", " : "") +
                      "threw exception: " + std::to_string(nexcp);
    if (nfail) fail = std::string(nexcp || nskip || nwarn ?
                      ", " : "") + "failed: " + std::to_string(nfail);
    print.note< tk::QUIET >
              ( "Of " + std::to_string(ncomplete) + " " + suite +
                " tests total: " + warn + skip + excp + fail );
  }

  return nfail || nexcp ? false : true;
}

} // unittest::
