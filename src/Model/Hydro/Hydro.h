//******************************************************************************
/*!
  \file      src/Model/Hydro/Hydro.h
  \author    J. Bakosi
  \date      Fri May 10 17:51:30 2013
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

//! Hydro model base for CRTP
//! See: http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
template< typename HydroType >
class Hydro : public Model {

  public:
    //! Constructor
    explicit Hydro(Memory* const memory,
                   Paradigm* const paradigm,
                   Control* const control,
                   const int nvelocity,
                   real* const velocities) : 
      Model(memory, paradigm, control, control->get<control::NPAR>()),
      m_nvelocity(nvelocity),
      m_velocities(velocities) {
      ErrChk(m_nvelocity > 0, FATAL,
             "Wrong number of particle velocity components");
      Assert(m_velocities != nullptr, FATAL, "Velocity pointer null?");
    }

    //! Destructor
    virtual ~Hydro() noexcept = default;

    //! Initialize particles
    void init() { static_cast<HydroType*>(this)->init(); }

    //! Advance particles in hydrodynamics model
    void advance(const real& dt) { static_cast<HydroType*>(this)->advance(dt); }

    //! Echo information on hydrodynamics model
    void echo() { static_cast<HydroType*>(this)->echo(); }

  protected:
    const int m_nvelocity;          //!< Number of hydrodynamics properties
    real* const m_velocities;       //!< Raw pointer to particle velocities

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
