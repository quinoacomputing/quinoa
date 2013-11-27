//******************************************************************************
/*!
  \file      src/RNG/TestU01.h
  \author    J. Bakosi
  \date      Wed 27 Nov 2013 12:43:10 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TestU01 statistical tests
  \details   TestU01 statistical tests
*/
//******************************************************************************
#ifndef TestU01_h
#define TestU01_h

#include <StatTest.h>

namespace rngtest {

//! TestU01 : StatTest
template< typename Result,                      //!< Results type
          Result* (*Creator)(void),             //!< Results creator function
          void (*Deleter)(Result *),            //!< Results deleter function
          double (*Run)(unif01_Gen*, Result*),  //!< Test runner function
          typename Info >                       //!< Test info
class TestU01 : public StatTest {

  public:
    //! Constructor
    explicit TestU01( const unif01_Gen* const gen ) :
      m_gen( gen ), m_res( ResultPtr( Creator() ) ) {};

    //! Destructor
    ~TestU01() noexcept override = default;

    //! Run
    double run() override {
      // Pretty awful that TestU01 does not guarantee the constness of gen
      return Run( const_cast<unif01_Gen*>(m_gen), m_res.get() );
    }

    //! Test name accessor
    const char* name() const override { return Info::name(); }

  private:
    //! Don't permit copy constructor
    TestU01(const TestU01&) = delete;
    //! Don't permit copy assigment
    TestU01& operator=(const TestU01&) = delete;
    //! Don't permit move constructor
    TestU01(TestU01&&) = delete;
    //! Don't permit move assigment
    TestU01& operator=(TestU01&&) = delete;

    const unif01_Gen* const m_gen;          //!< Raw ptr to TestU01 generator

    //! TestU01 results type with a custom deleter by TestU01
    using ResultPtr = TestU01Ptr< Result, Deleter >;
    //! TestU01 results struct (wrapped to std::unique_ptr)
    ResultPtr m_res;
};

} // rngtest::

#endif // TestU01_h
