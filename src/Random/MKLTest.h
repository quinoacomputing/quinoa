//******************************************************************************
/*!
  \file      src/Random/MKLTest.h
  \author    J. Bakosi
  \date      Mon 30 Sep 2013 10:18:19 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL random number generator tests
  \details   MKL random number generator tests
*/
//******************************************************************************
#ifndef MKLTest_h
#define MKLTest_h

#include <vector>

#include <mkl_vsl.h>

#include <RNGTest/InputDeck/InputDeck.h>
#include <Option.h>
#include <Quinoa/Options/RNG.h>

namespace rngtest {

//! MKL random number generator tests
class MKLTest {

  public:
    //! Constructor
    explicit MKLTest(const InputDeck& control);

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

    const quinoa::ctr::Option<quinoa::sel::RNG> m_rng;    //!< Available RNGs
    const std::vector<quinoa::sel::RNGType> m_testrng;    //!< RNGs to test

    std::vector<VSLStreamStatePtr> m_stream;              //!< RNG streams
};

} // namespace rngtest

#endif // MKLTest_h
