//******************************************************************************
/*!
  \file      src/Control/RNGTestControl.h
  \author    J. Bakosi
  \date      Tue 03 Sep 2013 10:49:45 PM MDT
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

//! RNGTestControl : Control<specialized to RNGTest>, see RNGTestControlTypes.h
class RNGTestControl :
  public quinoa::Control< // tag              type
                          control::title,     std::string,
                          control::suite,     select::RNGTestType,
                          control::generator, std::vector<select::RNGType> > {

  public:
    //! Constructor: set defaults
    explicit RNGTestControl() = default;

// //! Default bundle for RNGTest's control
// const Bundle defaults(
//   "",                                  //!< Title
//   select::RNGTestType::NO_RNGTEST,     //!< RNG test suite
//   std::vector<select::RNGType>()       //!< Random number generators
// );

    //! Destructor
    ~RNGTestControl() noexcept override = default;

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
