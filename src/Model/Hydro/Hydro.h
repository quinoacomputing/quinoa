//******************************************************************************
/*!
  \file      src/Model/Hydro/Hydro.h
  \author    J. Bakosi
  \date      Fri Jul 26 15:21:34 2013
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
#include <MKLRandom.h>
#include <MKLRndStream.h>

namespace Quinoa {

//! Hydro model base for CRTP
//! See: http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
template< typename HydroType >
class Hydro : public Model {

  public:
    //! Constructor
    explicit Hydro(Memory* const memory,
                   Paradigm* const paradigm,
                   Control* const control,
                   real* const particles) :
        Model(memory,
              paradigm,
              control,
              particles,
              control->get<control::NPAR>(),
              control->nprop()),
        m_offset(control->velocityOffset()),
        m_nvelocity(control->get<control::NVELOCITY>()) {
        ErrChk(m_nvelocity > 0, ExceptType::FATAL,
               "Wrong number of velocities");
      }

    //! Destructor
    virtual ~Hydro() noexcept = default;

    //! CRTP interface: Initialize particles
    void init() { static_cast<HydroType*>(this)->init(); }

    //! CRTP interface: Advance particles in hydro model
    void advance(const real& dt) { static_cast<HydroType*>(this)->advance(dt); }

  protected:
    const int m_offset;             //!< Velocity-offset relative to base
    const int m_nvelocity;          //!< Number of particle velocities

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
