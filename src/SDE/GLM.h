//******************************************************************************
/*!
  \file      src/SDE/GLM.h
  \author    J. Bakosi
  \date      Wed 01 Jan 2014 01:56:20 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Generalized Langevin hydrodynamics model
  \details   Generalized Langevin hydrodynamics model
*/
//******************************************************************************
#ifndef GLM_h
#define GLM_h

#include <Hydro.h>

namespace quinoa {

//! GLM : Hydro
class GLM : public Hydro {

  public:
    //! Constructor
    explicit GLM() = default;

    //! Destructor
    ~GLM() noexcept override = default;

  private:
    //! Don't permit copy constructor
    GLM(const GLM&) = delete;
    //! Don't permit copy assigment
    GLM& operator=(const GLM&) = delete;
    //! Don't permit move constructor
    GLM(GLM&&) = delete;
    //! Don't permit move assigment
    GLM& operator=(GLM&&) = delete;
};

} // quinoa::

#endif // GLM_h
