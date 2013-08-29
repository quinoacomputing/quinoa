//******************************************************************************
/*!
  \file      src/Random/Random.h
  \author    J. Bakosi
  \date      Thu Aug 29 15:24:55 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator base
  \details   Random number generator base
*/
//******************************************************************************
#ifndef Random_h
#define Random_h

#include <QuinoaTypes.h>

namespace quinoa {

//! Random number generator base
class Random {

  protected:
    //! Constructor: Default, compiler generated
    explicit Random() = default;

    //! Destructor: Default, compiler generated
    virtual ~Random() noexcept = default;

  private:
    //! Don't permit copy constructor
    Random(const Random&) = delete;
    //! Don't permit copy assigment
    Random& operator=(const Random&) = delete;
    //! Don't permit move constructor
    Random(Random&&) = delete;
    //! Don't permit move assigment
    Random& operator=(Random&&) = delete;
};

} // namespace quinoa

#endif // Random_h
