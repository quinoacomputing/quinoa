//******************************************************************************
/*!
  \file      src/SDE/SDE.h
  \author    J. Bakosi
  \date      Wed 01 Jan 2014 01:05:24 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SDE base
  \details   SDE base
*/
//******************************************************************************
#ifndef SDE_h
#define SDE_h

#include <Types.h>
#include <Exception.h>

namespace quinoa {

//! SDE base
class SDE {

  protected:
    //! Constructor: protected, designed to be base-only
    explicit SDE() = default;

    //! Destructor: protected, designed to be freed via children-only
    virtual ~SDE() noexcept = default;

  private:
    //! Don't permit copy constructor
    SDE(const SDE&) = delete;
    //! Don't permit copy assigment
    SDE& operator=(const SDE&) = delete;
    //! Don't permit move constructor
    SDE(SDE&&) = delete;
    //! Don't permit move assigment
    SDE& operator=(SDE&&) = delete;
};

} // quinoa::

#endif // SDE_h
