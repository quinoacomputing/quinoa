//******************************************************************************
/*!
  \file      src/Model/Hydro/Hydro.h
  \author    J. Bakosi
  \date      Fri May 10 17:22:14 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Hydro base
  \details   Hydro base
*/
//******************************************************************************
#ifndef Hydro_h
#define Hydro_h

#include <QuinoaTypes.h>
#include <Model.h>
#include <Control.h>

namespace Quinoa {

using namespace std;

//! Hydro model base
class Hydro : public Model {

  public:
    //! Constructor
    explicit Hydro(Memory* const memory,
                   Paradigm* const paradigm,
                   Control* const control) : 
      Model(memory, paradigm, control, control->get<control::NPAR>()),
      m_nvelocity(control->get<control::NVELOCITY>()) {
      ErrChk(m_nvelocity > 0, FATAL,
             "Wrong number of particle velocity components"); }

    //! Destructor
    virtual ~Hydro() noexcept = default;

    //! Interface for initializing particles
    virtual void init() = 0;

    //! Interface for advancing particles in hydro model
    virtual void advance(const real dt) = 0;

    //! Interface for echo information on hydro model
    virtual void echo() const = 0;

  protected:
    const int m_nvelocity;          //!< Number of hydrodynamics properties

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
