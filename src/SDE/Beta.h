//******************************************************************************
/*!
  \file      src/SDE/Beta.h
  \author    J. Bakosi
  \date      Wed 01 Jan 2014 01:06:58 PM MST
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Beta SDE
  \details   Beta SDE
*/
//******************************************************************************
#ifndef Beta_h
#define Beta_h

#include <SDE.h>

namespace quinoa {

//! Beta : SDE
class Beta : public SDE {

  public:
    //! Constructor
    explicit Beta() = default;

    //! Destructor
    ~Beta() noexcept override = default;

  private:
    //! Don't permit copy constructor
    Beta(const Beta&) = delete;
    //! Don't permit copy assigment
    Beta& operator=(const Beta&) = delete;
    //! Don't permit move constructor
    Beta(Beta&&) = delete;
    //! Don't permit move assigment
    Beta& operator=(Beta&&) = delete;
};

} // quinoa::

#endif // Beta_h
