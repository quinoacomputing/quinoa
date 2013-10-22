//******************************************************************************
/*!
  \file      src/RNG/RNG.h
  \author    J. Bakosi
  \date      Tue Oct 22 15:41:06 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator base
  \details   Random number generator base
*/
//******************************************************************************
#ifndef RNG_h
#define RNG_h

#include <Types.h>

namespace tk {

//! Random number generator base
class RNG {

  protected:
    //! Constructor: Default, compiler generated
    explicit RNG() = default;

    //! Destructor: Default, compiler generated
    virtual ~RNG() noexcept = default;

  private:
    //! Don't permit copy constructor
    RNG(const RNG&) = delete;
    //! Don't permit copy assigment
    RNG& operator=(const RNG&) = delete;
    //! Don't permit move constructor
    RNG(RNG&&) = delete;
    //! Don't permit move assigment
    RNG& operator=(RNG&&) = delete;
};

} // namespace tk

#endif // RNG_h
