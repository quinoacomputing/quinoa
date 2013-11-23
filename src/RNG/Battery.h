//******************************************************************************
/*!
  \file      src/RNG/Battery.h
  \author    J. Bakosi
  \date      Sat 23 Nov 2013 03:28:28 PM MST
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

  protected:
    const Base& m_base;                   //!< Essentials

    //! Statistical test type
    struct Test {
      std::unique_ptr< StatTest > ptr;
      double pval;
      // Constructor that moves
      template< class TestType >
      Test(std::unique_ptr< TestType >&& Ptr, double&& Pval) :
          ptr(std::move(Ptr)), pval(std::move(Pval)) {}
    };

    //! Container type for statistical tests and resulting p-values
    using TestContainer = std::vector< Test >;

    //! Add statistical test to battery
    template< class TestType, class GenPtrType >
    void add( TestContainer& tests, const GenPtrType& gen ) {
      std::unique_ptr< TestType > ptr( new TestType( gen.get() ) );
      tests.push_back( std::move( Test( std::move(ptr), -1.0 ) ) );
    }

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
