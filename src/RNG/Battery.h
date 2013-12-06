//******************************************************************************
/*!
  \file      src/RNG/Battery.h
  \author    J. Bakosi
  \date      Thu 05 Dec 2013 11:43:28 PM MST
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
