//******************************************************************************
/*!
  \file      src/Model/Hydro/GeneralizedLangevin/GeneralizedLangevin.h
  \author    J. Bakosi
  \date      Thu May  9 19:05:22 2013
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

//! GeneralizedLangevin : Hydro
class GeneralizedLangevin : public Hydro {

  public:
    //! Constructor
    explicit GeneralizedLangevin(Memory* const memory,
                                 Paradigm* const paradigm,
                                 Control* const control,
                                 real* const velocities);

    //! Destructor
    virtual ~GeneralizedLangevin() noexcept = default;

    //! Initialize particles
    virtual void init();

    //! Advance particles
    virtual void advance(const real dt);

    //! Echo information on the generalized Langevin model
    virtual void echo() const;

  private:
    //! Don't permit copy constructor
    GeneralizedLangevin(const GeneralizedLangevin&) = delete;
    //! Don't permit copy assigment
    GeneralizedLangevin& operator=(const GeneralizedLangevin&) = delete;
    //! Don't permit move constructor
    GeneralizedLangevin(GeneralizedLangevin&&) = delete;
    //! Don't permit move assigment
    GeneralizedLangevin& operator=(GeneralizedLangevin&&) = delete;

    real* const m_velocities;        //!< Raw pointer to particle velocities
};

} // namespace Quinoa

#endif // GeneralizedLangevin_h
