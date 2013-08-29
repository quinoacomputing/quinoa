//******************************************************************************
/*!
  \file      src/Control/RNGTestControl.h
  \author    J. Bakosi
  \date      Wed 28 Aug 2013 09:15:37 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator test suite control
  \details   Random number generator test suite control
*/
//******************************************************************************
#ifndef RNGTestControl_h
#define RNGTestControl_h

#include <Control.h>
#include <RNGTestControlTypes.h>

namespace Quinoa {

//! RNGTestControl : Control< specialized to RNGTest's control >
class RNGTestControl : public Control< control::Bundle,
                                       control::BoolBundle,
                                       control::BundlePosition > {

  public:
    //! Constructor
    explicit RNGTestControl() noexcept : Control(control::defaults) {}

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

} // namespace Quinoa

#endif // RNGTestControl_h
