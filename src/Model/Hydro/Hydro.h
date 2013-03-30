//******************************************************************************
/*!
  \file      src/Model/Hydro/Hydro.h
  \author    J. Bakosi
  \date      Sat 30 Mar 2013 01:18:22 PM MDT
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

  //! Number of particle properties for hydrodynamics: position + velocity
  const int NPROP = 6;

  public:
    //! Constructor
    Hydro(Memory* const memory,
          Paradigm* const paradigm,
          Control* const control,
          const string& name);

    //! Destructor
    virtual ~Hydro() {}

    //! Interface for initializing particles
    virtual void init() = 0;

    //! Interface for advancing particles in hydro model
    virtual void advance(const real dt) = 0;

    //! Interface for echo information on hydro model
    virtual void echo() = 0;

    //! Accessor to number of particle properties
    virtual int nprop() const { return m_nprop; }

  protected:
    const int m_nprop;              //!< Number of hydrodynamics properties
    const int m_npar;               //!< Number of particles

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
