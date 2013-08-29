//******************************************************************************
/*!
  \file      src/Random/MKLTest.h
  \author    J. Bakosi
  \date      Thu Aug 29 17:15:52 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL random number generator tests
  \details   MKL random number generator tests
*/
//******************************************************************************
#ifndef MKLTest_h
#define MKLTest_h

#include <vector>

#include <mkl_vsl.h>

#include <RNGTestControl.h>
#include <Option.h>
#include <RNGOptions.h>

namespace rngtest {

using quinoa::control::Option;

//! MKL random number generator tests
class MKLTest {

  public:
    //! Constructor
    explicit MKLTest(RNGTestControl* const control);

    //! Destructor
    ~MKLTest();

  private:
    //! Don't permit copy constructor
    MKLTest(const MKLTest&) = delete;
    //! Don't permit copy assigment
    MKLTest& operator=(const MKLTest&) = delete;
    //! Don't permit move constructor
    MKLTest(MKLTest&&) = delete;
    //! Don't permit move assigment
    MKLTest& operator=(MKLTest&&) = delete;

    //! Special error handler for MKL calls
    void MKLErrChk(int vslerr) const;

    const Option<select::RNG> m_rng;               //!< Available RNGs
    const std::vector<select::RNGType> m_testrng;  //!< RNGs to test

    std::vector<VSLStreamStatePtr> m_stream;       //!< RNG streams
};

} // namespace rngtest

#endif // MKLTest_h
