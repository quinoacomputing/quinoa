//******************************************************************************
/*!
  \file      src/Model/Hydro/SimplifiedLangevin/SimplifiedLangevin.h
  \author    J. Bakosi
  \date      Fri May 10 17:52:48 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Simplified Langevin hydrodynamics model
  \details   Simplified Langevin hydrodynamics model
*/
//******************************************************************************
#ifndef SimplifiedLangevin_h
#define SimplifiedLangevin_h

#include <mkl_vsl.h>

#include <Memory.h>
#include <Hydro.h>

namespace Quinoa {

class Memory;
class Paradigm;
class MKLRandom;
class MKLRndStream;
class MemoryEntry;
class JPDF;

//! SimplifiedLangevin : Hydro<SimplifiedLangevin> child for CRTP
//! See: http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
class SimplifiedLangevin : public Hydro<SimplifiedLangevin> {

  public:
    //! Constructor
    explicit SimplifiedLangevin(Memory* const memory,
                                Paradigm* const paradigm,
                                Control* const control,
                                real* const velocities);

    //! Destructor
    virtual ~SimplifiedLangevin() noexcept;

    //! Initialize particles
    void init();

    //! Advance particles
    void advance(const real& dt);

    //! Echo information on the simplified Langevin model
    void echo() const;

  private:
    //! Don't permit copy constructor
    SimplifiedLangevin(const SimplifiedLangevin&) = delete;
    //! Don't permit copy assigment
    SimplifiedLangevin& operator=(const SimplifiedLangevin&) = delete;
    //! Don't permit move constructor
    SimplifiedLangevin(SimplifiedLangevin&&) = delete;
    //! Don't permit move assigment
    SimplifiedLangevin& operator=(SimplifiedLangevin&&) = delete;

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    void finalize() noexcept;

    const real m_C0;                //!< Parameter C0 in SLM

    const VSLStreamStatePtr* m_str; //!< Array of MKL VSL stream state pointers
    MKLRandom* m_random;            //!< Random number generator object
    MKLRndStream* m_rndStr;         //!< Random number stream object
};

} // namespace Quinoa

#endif // SimplifiedLangevin_h
