//******************************************************************************
/*!
  \file      src/Model/Hydro/SimplifiedLangevin/SimplifiedLangevin.h
  \author    J. Bakosi
  \date      Sat 30 Mar 2013 01:20:51 PM MDT
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
class MKLRandom;
class MKLRndStream;
class MemoryEntry;
class JPDF;

//! SimplifiedLangevin : Hydro
class SimplifiedLangevin : public Hydro {

  public:
    //! Constructor
    SimplifiedLangevin(Memory* const memory,
                       Paradigm* const paradigm,
                       Control* const control);

    //! Destructor
    virtual ~SimplifiedLangevin() {}

    //! Initialize particles
    virtual void init();

    //! Advance particles
    virtual void advance(const real dt);

    //! Echo information on the simplified Langevin model
    virtual void echo();

    //! Constant accessor to particle properties pointer
    virtual const real* particles() const { return m_particles.ptr; }    

  private:
    //! Don't permit copy constructor
    SimplifiedLangevin(const SimplifiedLangevin&) = delete;
    //! Don't permit copy assigment
    SimplifiedLangevin& operator=(const SimplifiedLangevin&) = delete;
    //! Don't permit move constructor
    SimplifiedLangevin(SimplifiedLangevin&&) = delete;
    //! Don't permit move assigment
    SimplifiedLangevin& operator=(SimplifiedLangevin&&) = delete;

    Data<real> m_particles;        //!< Particle properties
};

} // namespace Quinoa

#endif // SimplifiedLangevin_h
