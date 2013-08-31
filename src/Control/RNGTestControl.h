//******************************************************************************
/*!
  \file      src/Control/RNGTestControl.h
  \author    J. Bakosi
  \date      Fri 30 Aug 2013 06:15:44 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator test suite control
  \details   Random number generator test suite control
*/
//******************************************************************************
#ifndef RNGTestControl_h
#define RNGTestControl_h

#include <Control.h>
#include <RNGTestControlTypes.h>

namespace rngtest {

//! RNGTestControl : Control< specialized to RNGTest's control >
class RNGTestControl : public quinoa::Control< control::Bundle,
                                               control::BoolBundle,
                                               control::BundlePosition > {

  public:
    //! Constructor
    explicit RNGTestControl() noexcept
      : quinoa::Control< control::Bundle,
                         control::BoolBundle,
                         control::BundlePosition >(control::defaults) {}

    //! Destructor
    virtual ~RNGTestControl() noexcept = default;

  private:
    //! Don't permit copy constructor
    RNGTestControl(const RNGTestControl&) = delete;
    //! Don't permit copy assigment
    RNGTestControl& operator=(const RNGTestControl&) = delete;
    //! Don't permit move constructor
    RNGTestControl(RNGTestControl&&) = delete;
    //! Don't permit move assigment
    RNGTestControl& operator=(RNGTestControl&&) = delete;
};

} // namespace rngtest

#endif // RNGTestControl_h
