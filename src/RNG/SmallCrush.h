//******************************************************************************
/*!
  \file      src/RNG/SmallCrush.h
  \author    J. Bakosi
  \date      Thu 21 Nov 2013 02:48:49 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SmallCrush battery
  \details   SmallCrush battery
*/
//******************************************************************************
#ifndef SmallCrush_h
#define SmallCrush_h

#include <TestU01.h>

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

    std::unique_ptr< tk::RNG > m_rng;           //!< Random number generator
};

} // rngtest::

#endif // SmallCrush_h
