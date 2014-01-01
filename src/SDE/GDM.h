//******************************************************************************
/*!
  \file      src/SDE/GDM.h
  \author    J. Bakosi
  \date      Wed 01 Jan 2014 01:32:59 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Generalized Dirichlet mix model
  \details   Generalized Dirichlet mix model
*/
//******************************************************************************
#ifndef GDM_h
#define GDM_h

#include <Mix.h>

namespace quinoa {

//! GDM : Mix
class GDM : public Mix {

  public:
    //! Constructor
    explicit GDM() {}

    //! Destructor
    ~GDM() noexcept override = default;

  private:
    //! Don't permit copy constructor
    GDM(const GDM&) = delete;
    //! Don't permit copy assigment
    GDM& operator=(const GDM&) = delete;
    //! Don't permit move constructor
    GDM(GDM&&) = delete;
    //! Don't permit move assigment
    GDM& operator=(GDM&&) = delete;
};

} // quinoa::

#endif // GDM_h
