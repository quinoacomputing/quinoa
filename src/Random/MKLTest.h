//******************************************************************************
/*!
  \file      src/Random/MKLTest.h
  \author    J. Bakosi
  \date      Sat 24 Aug 2013 07:29:46 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL random number generator tests
  \details   MKL random number generator tests
*/
//******************************************************************************
#ifndef MKLTest_h
#define MKLTest_h

#include <vector>

#include <mkl_vsl.h>

#include <Option.h>
#include <RNGOptions.h>

namespace Quinoa {

class Control;

//! MKL random number generator tests
class MKLTest {

  public:
    //! Constructor
    explicit MKLTest(Control* const control);

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

    const control::Option<select::RNG> m_rng;      //!< Available RNGs
    const std::vector<select::RNGType> m_testrng;  //!< RNGs to test

    std::vector<VSLStreamStatePtr> m_stream;       //!< RNG streams
};

} // namespace Quinoa

#endif // MKLTest_h
