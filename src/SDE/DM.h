//******************************************************************************
/*!
  \file      src/SDE/DM.h
  \author    J. Bakosi
  \date      Wed 01 Jan 2014 01:33:10 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet mix model
  \details   Dirichlet mix model
*/
//******************************************************************************
#ifndef DM_h
#define DM_h

#include <Mix.h>

namespace quinoa {

//! DM : Mix
class DM : public Mix {

  public:
    //! Constructor
    explicit DM() {}

    //! Destructor
    ~DM() noexcept override = default;

  private:
    //! Don't permit copy constructor
    DM(const DM&) = delete;
    //! Don't permit copy assigment
    DM& operator=(const DM&) = delete;
    //! Don't permit move constructor
    DM(DM&&) = delete;
    //! Don't permit move assigment
    DM& operator=(DM&&) = delete;
};

} // quinoa::

#endif // DM_h
