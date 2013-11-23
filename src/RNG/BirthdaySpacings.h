//******************************************************************************
/*!
  \file      src/RNG/BirthdaySpacings.h
  \author    J. Bakosi
  \date      Fri 22 Nov 2013 05:46:37 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Statistical tests suggested by George Marsaglia
  \details   Statistical tests suggested by George Marsaglia
*/
//******************************************************************************
#ifndef BirthdaySpacings_h
#define BirthdaySpacings_h

#include <Marsaglia.h>

namespace rngtest {

//! BirthdaySpacings : Marsaglia
class BirthdaySpacings : public Marsaglia {

  public:
    //! Constructor
    explicit BirthdaySpacings(const unif01_Gen* const gen) : m_gen(gen) {}

    //! Destructor
    ~BirthdaySpacings() noexcept override = default;

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
};

} // rngtest::

#endif // BirthdaySpacings_h
