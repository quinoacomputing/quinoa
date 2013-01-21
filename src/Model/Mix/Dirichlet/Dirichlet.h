//******************************************************************************
/*!
  \file      src/Model/Mix/Dirichlet/Dirichlet.h
  \author    J. Bakosi
  \date      Mon 21 Jan 2013 11:28:21 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet mix model
  \details   Dirichlet mix model
*/
//******************************************************************************
#ifndef Dirichlet_h
#define Dirichlet_h

#include <QuinoaTypes.h>
#include <Mix.h>

namespace Quinoa {

//! Dirichlet : Mix
class Dirichlet : public Mix {

  public:
    //! Constructor
    Dirichlet(const int& nscalar);

    //! Destructor
    virtual ~Dirichlet() {}

    //! Echo information on Dirichlet model
    virtual void echo();

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
