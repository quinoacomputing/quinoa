//******************************************************************************
/*!
  \file      src/Model/Mix/Dirichlet/Dirichlet.h
  \author    J. Bakosi
  \date      Thu 08 Nov 2012 06:20:49 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet mix model
  \details   Dirichlet mix model
*/
//******************************************************************************
#ifndef Dirichlet_h
#define Dirichlet_h

#include <MixModel.h>
#include <DirichletSetup.h>

namespace Quinoa {

//! Dirichlet : MixModel
class Dirichlet : public MixModel {

  public:
    //! Constructor
    Dirichlet() = default;

    //! Destructor
    ~Dirichlet() = default;

    //! Set initial conditions
    void setIC() {}

  private:
    //! Don't permit copy constructor
    Dirichlet(const Dirichlet&) = delete;
    //! Don't permit copy assigment
    Dirichlet& operator=(const Dirichlet&) = delete;
    //! Don't permit move constructor
    Dirichlet(Dirichlet&&) = delete;
    //! Don't permit move assigment
    Dirichlet& operator=(Dirichlet&&) = delete;

    DirichletSetup m_setup;
};

} // namespace Quinoa

#endif // Dirichlet_h
