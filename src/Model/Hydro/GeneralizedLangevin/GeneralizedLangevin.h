//******************************************************************************
/*!
  \file      src/Model/Hydro/GeneralizedLangevin/GeneralizedLangevin.h
  \author    J. Bakosi
  \date      Fri May 10 17:48:14 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Generalized Langevin hydrodynamics model
  \details   Generalized Langevin hydrodynamics model
*/
//******************************************************************************
#ifndef GeneralizedLangevin_h
#define GeneralizedLangevin_h

#include <Memory.h>
#include <Hydro.h>

namespace Quinoa {

class Memory;
class Paradigm;
class MKLRandom;
class MKLRndStream;
class MemoryEntry;
class JPDF;

//! GeneralizedLangevin : Hydro<GeneralizedLangevin> child for CRTP
//! See: http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
class GeneralizedLangevin : public Hydro<GeneralizedLangevin> {

  public:
    //! Constructor
    explicit GeneralizedLangevin(Memory* const memory,
                                 Paradigm* const paradigm,
                                 Control* const control,
                                 real* const velocities);

    //! Destructor
    virtual ~GeneralizedLangevin() noexcept = default;

    //! Initialize particles
    void init();

    //! Advance particles
    void advance(const real& dt);

    //! Echo information on the generalized Langevin model
    void echo() const;

  private:
    //! Don't permit copy constructor
    GeneralizedLangevin(const GeneralizedLangevin&) = delete;
    //! Don't permit copy assigment
    GeneralizedLangevin& operator=(const GeneralizedLangevin&) = delete;
    //! Don't permit move constructor
    GeneralizedLangevin(GeneralizedLangevin&&) = delete;
    //! Don't permit move assigment
    GeneralizedLangevin& operator=(GeneralizedLangevin&&) = delete;
};

} // namespace Quinoa

#endif // GeneralizedLangevin_h
