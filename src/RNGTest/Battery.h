//******************************************************************************
/*!
  \file      src/RNGTest/Battery.h
  \author    J. Bakosi
  \date      Wed 09 Apr 2014 04:08:01 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Battery base
  \details   Battery base
*/
//******************************************************************************
#ifndef Battery_h
#define Battery_h

#include <Base.h>
#include <StatTest.h>

namespace rngtest {

//! Battery
class Battery {

  public:
    //! Constructor
    explicit Battery(const Base& base) : m_base(base) {}

    //! Destructor
    virtual ~Battery() noexcept = default;

    //! Run battery of RNG tests
    virtual void run() = 0;

    //! Print list of registered statistical tests
    virtual void print() const = 0;

    //! Container type for statistical tests
    using TestContainer = std::vector< std::unique_ptr< StatTest > >;
    using Tsize = TestContainer::size_type;

    //! Return number of statistical tests in battery
    virtual Tsize ntest() const = 0;

    //! Return number of statistics produced by battery
    virtual StatTest::Psize nstat() const = 0;

  protected:
    const Base& m_base;                   //!< Essentials

  private:
    //! Don't permit copy constructor
    Battery(const Battery&) = delete;
    //! Don't permit copy assigment
    Battery& operator=(const Battery&) = delete;
    //! Don't permit move constructor
    Battery(Battery&&) = delete;
    //! Don't permit move assigment
    Battery& operator=(Battery&&) = delete;
};

} // rngtest::

#endif // Battery_h
