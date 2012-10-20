//******************************************************************************
/*!
  \file      src/Random/Random.h
  \author    J. Bakosi
  \date      Fri 19 Oct 2012 04:31:52 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator base
  \details   Random number generator base
*/
//******************************************************************************
#ifndef Random_h
#define Random_h

#include <QuinoaTypes.h>

namespace Quinoa {

//! Random number generator base
class Random {

  public:
    //! Constructor: Default, compiler generated
    Random() = default;

    //! Destructor: Default, compiler generated
    ~Random() = default;

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

} // namespace Quinoa

#endif // Random_h
