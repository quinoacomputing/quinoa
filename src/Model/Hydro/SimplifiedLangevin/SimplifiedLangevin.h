//******************************************************************************
/*!
  \file      src/Model/Hydro/SimplifiedLangevin/SimplifiedLangevin.h
  \author    J. Bakosi
  \date      Mon 27 May 2013 01:07:55 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Simplified Langevin hydrodynamics model
  \details   Simplified Langevin hydrodynamics model
*/
//******************************************************************************
#ifndef SimplifiedLangevin_h
#define SimplifiedLangevin_h

#include <Memory.h>
#include <Hydro.h>

namespace Quinoa {

class Memory;
class Paradigm;
class JPDF;

//! SimplifiedLangevin : Hydro<SimplifiedLangevin> child for CRTP
//! See: http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
class SimplifiedLangevin : public Hydro<SimplifiedLangevin> {

  public:
    //! Constructor
    explicit SimplifiedLangevin(Memory* const memory,
                                Paradigm* const paradigm,
                                Control* const control,
                                real* const particles) :
      Hydro<SimplifiedLangevin>(memory, paradigm, control, particles),
      m_C0(control->get<control::C0>()) {
      //ErrChk on m_C0
    }

    //! Destructor
    virtual ~SimplifiedLangevin() noexcept = default;

    //! Return hydro model identification
    select::Hydro::Enum id() noexcept { return select::Hydro::Enum::SLM; }

    //! Initialize particles
    void init();

    //! Advance particles
    void advance(const real& dt);

  private:
    //! Don't permit copy constructor
    SimplifiedLangevin(const SimplifiedLangevin&) = delete;
    //! Don't permit copy assigment
    SimplifiedLangevin& operator=(const SimplifiedLangevin&) = delete;
    //! Don't permit move constructor
    SimplifiedLangevin(SimplifiedLangevin&&) = delete;
    //! Don't permit move assigment
    SimplifiedLangevin& operator=(SimplifiedLangevin&&) = delete;

    const real m_C0;                //!< Parameter C0 in SLM
};

} // namespace Quinoa

#endif // SimplifiedLangevin_h
