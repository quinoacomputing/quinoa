//******************************************************************************
/*!
  \file      src/Models/Mix/Dirichlet/Dirichlet.h
  \author    J. Bakosi
  \date      Wed 07 Nov 2012 08:33:35 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet mix model
  \details   Dirichlet mix model
*/
//******************************************************************************
#ifndef Dirichlet_h
#define Dirichlet_h

#include <MixModel.h>

namespace Quinoa {

//! Dirichlet : MixModel
class Dirichlet : MixModel {

  public:
    //! Constructor
    Dirichlet() = default;

    //! Destructor
    ~Dirichlet() = default;

  private:
    //! Don't permit copy constructor
    Dirichlet(const Dirichlet&) = delete;
    //! Don't permit copy assigment
    Dirichlet& operator=(const Dirichlet&) = delete;
    //! Don't permit move constructor
    Dirichlet(Dirichlet&&) = delete;
    //! Don't permit move assigment
    Dirichlet& operator=(Dirichlet&&) = delete;
};

} // namespace Quinoa

#endif // Dirichlet_h
