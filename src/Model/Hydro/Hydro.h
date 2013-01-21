//******************************************************************************
/*!
  \file      src/Model/Hydro/Hydro.h
  \author    J. Bakosi
  \date      Mon 21 Jan 2013 11:58:36 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Hydro base
  \details   Hydro base
*/
//******************************************************************************
#ifndef Hydro_h
#define Hydro_h

#include <string>

#include <QuinoaTypes.h>
#include <Model.h>

namespace Quinoa {

using namespace std;

//! Hydro model base
class Hydro : public Model {

  public:
    //! Constructor
    Hydro(const string& name);

    //! Destructor
    virtual ~Hydro() {}

    //! Interface for echo information on mix model
    virtual void echo() = 0;

  protected:
    const string m_name;           //!< Name of hydro model

  private:
    //! Don't permit copy constructor
    Hydro(const Hydro&) = delete;
    //! Don't permit copy assigment
    Hydro& operator=(const Hydro&) = delete;
    //! Don't permit move constructor
    Hydro(Hydro&&) = delete;
    //! Don't permit move assigment
    Hydro& operator=(Hydro&&) = delete;
};

} // namespace Quinoa

#endif // Hydro_h
