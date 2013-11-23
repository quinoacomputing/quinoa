//******************************************************************************
/*!
  \file      src/RNG/SmallCrush.h
  \author    J. Bakosi
  \date      Fri 22 Nov 2013 06:31:07 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SmallCrush battery
  \details   SmallCrush battery
*/
//******************************************************************************
#ifndef SmallCrush_h
#define SmallCrush_h

extern "C" {
  #include <unif01.h>
}

#include <TestU01.h>
#include <StatTest.h>

namespace rngtest {

//! SmallCrush : TestU01
class SmallCrush : public TestU01 {

  public:
    //! Constructor
    explicit SmallCrush(const Base& base);

    //! Destructor
    ~SmallCrush() noexcept override = default;

    //! Run battery of RNG tests
    void run() override;

  private:
    //! Don't permit copy constructor
    SmallCrush(const SmallCrush&) = delete;
    //! Don't permit copy assigment
    SmallCrush& operator=(const SmallCrush&) = delete;
    //! Don't permit move constructor
    SmallCrush(SmallCrush&&) = delete;
    //! Don't permit move assigment
    SmallCrush& operator=(SmallCrush&&) = delete;

    std::string m_rngname;

    //!< TestU01 external generator custom deleter provided by TestU01
    struct Gen01Deleter {
      void operator()( unif01_Gen* gen ) {
        unif01_DeleteExternGen01( gen );
      }
    };
    //!< TestU01 external generator type with a custom deleter by TestU01
    using Gen01Pointer = std::unique_ptr< unif01_Gen, Gen01Deleter >;

    //!< TestU01 external generator
    Gen01Pointer m_gen;

    std::vector< std::unique_ptr< StatTest > > m_tests;
};

} // rngtest::

#endif // SmallCrush_h
