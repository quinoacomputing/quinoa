//******************************************************************************
/*!
  \file      src/Model/Hydro/Hydro.h
  \author    J. Bakosi
  \date      Sun 15 Sep 2013 05:27:13 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Hydro base
  \details   Hydro base
*/
//******************************************************************************
#ifndef Hydro_h
#define Hydro_h

#include <QuinoaTypes.h>
#include <Model.h>
#include <QuinoaControl.h>
#include <MKLRandom.h>
#include <MKLRndStream.h>

namespace quinoa {

//! Hydro model base for CRTP
//! See: http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
template< typename HydroType >
class Hydro : public Model {

  public:
    //! Constructor
    explicit Hydro(const Base& base, real* const particles) :
        Model(base,
              particles,
              base.control.get<control::component, control::npar>(),
              base.control.nprop()),
        m_offset(base.control.velocityOffset()),
        m_nvelocity(base.control.get<control::component, control::nvelocity>()) {
        ErrChk(m_nvelocity > 0, ExceptType::FATAL,
               "Wrong number of velocities");
      }

    //! Destructor
    ~Hydro() noexcept override = default;

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

} // namespace quinoa

#endif // Hydro_h
