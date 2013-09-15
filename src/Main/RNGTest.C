//******************************************************************************
/*!
  \file      src/Main/RNGTest.C
  \author    J. Bakosi
  \date      Sat 14 Sep 2013 08:08:42 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa random number generator test suite
  \details   Quinoa random number generator test suite
*/
//******************************************************************************

#include <iostream>

#include <QuinoaConfig.h>
#include <Exception.h>
#include <RNGTestDriver.h>

namespace rngtest {

static void echoName()
//******************************************************************************
//  Echo Name
//! \author  J. Bakosi
//******************************************************************************
{
  std::cout << "==========================================\n";
  std::cout << "Quinoa: Random number generator test suite\n";
  std::cout << "==========================================" << std::endl;
}

static void echoBuildInfo()
//******************************************************************************
//  Echo build environment
//! \details Echo information read from [build]/Base/QuinoaConfig.h filled by
//!          CMake based on src/MainQuinoaConfig.h.in
//! \author  J. Bakosi
//******************************************************************************
{
  std::cout << "\nBuild environment:"
               "\n------------------\n";
  std::cout << " * Executable                  : " << QUINOA_RNGTEST_EXECUTABLE
            << "\n";
  std::cout << " * Version                     : " << QUINOA_VERSION << "\n";
  std::cout << " * Release                     : " << QUINOA_RELEASE << "\n";
  std::cout << " * Revision                    : " << QUINOA_GIT_COMMIT << "\n";
  std::cout << " * CMake build type            : " << QUINOA_BUILD_TYPE << "\n";
  std::cout << " * MPI C++ compiler            : " << QUINOA_MPI_COMPILER<<"\n";
  std::cout << " * MPI underlying C++ compiler : " << QUINOA_COMPILER << "\n";
  std::cout << " * Build date                  : " << QUINOA_BUILD_DATE << "\n";
#ifdef NDEBUG
  std::cout << " * Built without asserts" << "\n";
#else  // NDEBUG
  std::cout << " * Built with asserts" << "\n";
#endif // NDEBUG
  std::cout << std::endl;
}

} // namespace rngtest

using namespace quinoa;
using namespace rngtest;

int main(int argc, char* argv[])
//******************************************************************************
//  Quinoa random number generator test suite
//! \author  J. Bakosi
//******************************************************************************
{
  ErrCode error = ErrCode::SUCCESS;
  try {

    // Echo program name
    echoName();
    // Echo build environment
    echoBuildInfo();

    // Create driver
    RNGTestDriver driver(argc, argv);

    // Run RNG tests
    driver.execute();

  } // Catch and handle Quina::Exceptions
    catch (Exception& qe) {
      error = qe.handleException();
    }
    // Catch std::exceptions and transform them into Quinoa::Exceptions without
    // file:line:func information
    catch (std::exception& se) {
      Exception qe(ExceptType::RUNTIME, se.what());
      error = qe.handleException();
    }
    // Catch uncaught exceptions and still do cleanup
    catch (...) {
      Exception qe(ExceptType::UNCAUGHT, "Non-standard exception");
      error = qe.handleException();
    }

  // Return error code
  return static_cast<int>(error);
}
