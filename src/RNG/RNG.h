//******************************************************************************
/*!
  \file      src/RNG/RNG.h
  \author    J. Bakosi
  \date      Thu 21 Nov 2013 02:58:20 PM MST
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

  public:
    //! Constructor: Default, compiler generated
    explicit RNG() = default;

    //! Destructor: Default, compiler generated
    virtual ~RNG() noexcept = default;

    //! Uniform RNG interface
    virtual void uniform(int tid, int num, double* r) const = 0;

    //! Gaussian RNG interface
    virtual void gaussian(int tid, int num, double* r) const = 0;

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
