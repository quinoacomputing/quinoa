//******************************************************************************
/*!
  \file      src/RNG/BirthdaySpacings.h
  \author    J. Bakosi
  \date      Sat 23 Nov 2013 04:43:26 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Statistical tests suggested by George Marsaglia
  \details   Statistical tests suggested by George Marsaglia
*/
//******************************************************************************
#ifndef BirthdaySpacings_h
#define BirthdaySpacings_h

#include <TestU01Wrap.h>
#include <Marsaglia.h>

namespace rngtest {

//! BirthdaySpacings : Marsaglia
class BirthdaySpacings : public Marsaglia {

  public:
    //! Constructor
    explicit BirthdaySpacings(const unif01_Gen* const gen);

    //! Destructor
    ~BirthdaySpacings() noexcept override = default;

    //! Run
    double run() override;

    //! Test name accessor
    const char* name() const override { return m_res->name; }

  private:
    //! Don't permit copy constructor
    BirthdaySpacings(const BirthdaySpacings&) = delete;
    //! Don't permit copy assigment
    BirthdaySpacings& operator=(const BirthdaySpacings&) = delete;
    //! Don't permit move constructor
    BirthdaySpacings(BirthdaySpacings&&) = delete;
    //! Don't permit move assigment
    BirthdaySpacings& operator=(BirthdaySpacings&&) = delete;

    const unif01_Gen* const m_gen;          //!< Raw ptr to TestU01 generator

    //! TestU01 Poisson results type with a custom deleter by TestU01
    using PoissonResPtr = TestU01Ptr< sres_Poisson, sres_DeletePoisson >;
    //! TestU01 Poisson results struct (wrapped to std::unique_ptr)
    PoissonResPtr m_res;
};

} // rngtest::

#endif // BirthdaySpacings_h
