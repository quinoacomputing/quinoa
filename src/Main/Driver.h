//******************************************************************************
/*!
  \file      src/Main/Driver.h
  \author    J. Bakosi
  \date      Mon 26 May 2014 04:50:01 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver base
  \details   Driver base
*/
//******************************************************************************
#ifndef Driver_h
#define Driver_h

#include <Exception.h>

namespace tk {

//! Driver base
class Driver {

  protected:
    //! Constructor
    explicit Driver() = default;

    //! Destructor
    virtual ~Driver() = default;

    //! Execute
    virtual void execute() const {
      Throw( ExceptType::WARNING,
             "Driver::execute() called; override not required and undefined" );
    }

  private:
    //! Don't permit copy constructor
    Driver(const Driver&) = delete;
    //! Don't permit assigment constructor
    Driver& operator=(const Driver&) = delete;
    //! Don't permit move constructor
    Driver(Driver&&) = delete;
    //! Don't permit move assignment
    Driver& operator=(Driver&&) = delete;
};

} // namespace tk

#endif // Driver_h
