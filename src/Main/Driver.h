//******************************************************************************
/*!
  \file      src/Main/Driver.h
  \author    J. Bakosi
  \date      Thu Aug 29 15:17:09 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver base
  \details   Driver base
*/
//******************************************************************************
#ifndef Driver_h
#define Driver_h

namespace quinoa {

class Timer;

//! Driver base class
class Driver {

  public:
    //! Constructor
    explicit Driver();

    //! Destructor
    virtual ~Driver() noexcept;

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    virtual void finalize() noexcept;

    //! Execute
    virtual void execute() const = 0;

  protected:
    Timer* m_timer;                   //!< Timer object

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

} // namespace quinoa

#endif // Driver_h
