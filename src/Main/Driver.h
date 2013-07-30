//******************************************************************************
/*!
  \file      src/Main/Driver.h
  \author    J. Bakosi
  \date      Mon 29 Jul 2013 09:30:56 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver base class declaration
  \details   Driver base class declaration
*/
//******************************************************************************
#ifndef Driver_h
#define Driver_h

namespace Quinoa {

class Control;

//! Driver base class
class Driver {

  public:
    //! Constructor
    Driver(int argc, char** argv);

    //! Destructor
    virtual ~Driver() noexcept;

    //! Solve
    virtual void execute() const = 0;

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    virtual void finalize() noexcept;

    //! Const control object accessor
    Control* control() const noexcept { return m_control; }

  private:
    //! Don't permit copy constructor
    Driver(const Driver&) = delete;
    //! Don't permit assigment constructor
    Driver& operator=(const Driver&) = delete;
    //! Don't permit move constructor
    Driver(Driver&&) = delete;
    //! Don't permit move assignment
    Driver& operator=(Driver&&) = delete;

    Control* m_control;               //!< Control object
};

} // namespace Quinoa

#endif // Driver_h
