//******************************************************************************
/*!
  \file      src/SDE/SkewNormal.h
  \author    J. Bakosi
  \date      Wed 01 Jan 2014 02:02:00 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Skew-normal SDE
  \details   Skew-normal SDE
*/
//******************************************************************************
#ifndef SkewNormal_h
#define SkewNormal_h

#include <SDE.h>

namespace quinoa {

//! SkewNormal : SDE
class SkewNormal : public SDE {

  public:
    //! Constructor
    explicit SkewNormal() = default;

    //! Destructor
    ~SkewNormal() noexcept override = default;

  private:
    //! Don't permit copy constructor
    SkewNormal(const SkewNormal&) = delete;
    //! Don't permit copy assigment
    SkewNormal& operator=(const SkewNormal&) = delete;
    //! Don't permit move constructor
    SkewNormal(SkewNormal&&) = delete;
    //! Don't permit move assigment
    SkewNormal& operator=(SkewNormal&&) = delete;
};

} // quinoa::

#endif // SkewNormal_h
