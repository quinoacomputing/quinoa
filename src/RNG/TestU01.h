//******************************************************************************
/*!
  \file      src/RNG/TestU01.h
  \author    J. Bakosi
  \date      Mon 25 Nov 2013 10:48:38 PM MST
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
template< typename ResPtr, ResPtr* (*Creator)(void), void (*Deleter)(ResPtr *),
          double (*Run)(unif01_Gen*,
                        typename TestU01Ptr<ResPtr, Deleter>::element_type*) >
class TestU01 : public StatTest {

  public:
    //! Constructor
    explicit TestU01(const unif01_Gen* const gen) :
      m_gen( gen ), m_res( ResultPtr( Creator() ) ) {};

    //! Destructor
    ~TestU01() noexcept override = default;

    //! Run
    double run() override {
      return Run( const_cast<unif01_Gen*>(m_gen), m_res.get() );
    }

    //! Test name accessor
    const char* name() const override { return m_res->name; }

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
    using ResultPtr = TestU01Ptr< ResPtr, Deleter >;
    //! TestU01 results struct (wrapped to std::unique_ptr)
    ResultPtr m_res;
};

} // rngtest::

#endif // TestU01_h
