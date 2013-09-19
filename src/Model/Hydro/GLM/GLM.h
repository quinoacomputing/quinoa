//******************************************************************************
/*!
  \file      src/Model/Hydro/GLM/GLM.h
  \author    J. Bakosi
  \date      Wed Sep 18 14:02:07 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Generalized Langevin hydrodynamics model
  \details   Generalized Langevin hydrodynamics model
*/
//******************************************************************************
#ifndef GeneralizedLangevin_h
#define GeneralizedLangevin_h

#include <Memory.h>
#include <Hydro.h>

namespace quinoa {

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
    explicit GeneralizedLangevin(const Base& base, real* const particles) :
      Hydro(base, particles) {
      // Error out if hydro model selected at compile time does not match that
      // whose options are given in control file
      //control->matchModels<select::Hydro, select::HydroType, control::HYDRO>(
      //  select::HydroType::GLM);
    }

    //! Destructor
    ~GeneralizedLangevin() noexcept override = default;

    //! Initialize particles
    void init() override;

    //! Advance particles
    void advance(const real& dt) override;

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

} // namespace quinoa

#endif // GLM_h
