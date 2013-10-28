//******************************************************************************
/*!
  \file      src/Model/Hydro/GLM/GLM.h
  \author    J. Bakosi
  \date      Mon Oct 28 08:55:37 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Generalized Langevin hydrodynamics model
  \details   Generalized Langevin hydrodynamics model
*/
//******************************************************************************
#ifndef GLM_h
#define GLM_h

#include <Hydro/Hydro.h>

namespace quinoa {

//! GLM : Hydro
class GLM : public Hydro {

  public:
    //! Constructor
    explicit GLM() {}
//     explicit GLM(const Base& base, tk::real* const particles) :
//       Hydro(base, particles) {
//     }

    //! Destructor
    ~GLM() noexcept override = default;

//     //! Initialize particles
//     void init() override;
// 
//     //! Advance particles
//     void advance(int p, int tid, tk::real dt) override;

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
