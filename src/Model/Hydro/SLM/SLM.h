//******************************************************************************
/*!
  \file      src/Model/Hydro/SLM/SLM.h
  \author    J. Bakosi
  \date      Thu Sep 19 09:15:59 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Simplified Langevin hydrodynamics model
  \details   Simplified Langevin hydrodynamics model
*/
//******************************************************************************
#ifndef SLM_h
#define SLM_h

#include <Memory.h>
#include <Hydro.h>

namespace quinoa {

class Memory;
class Paradigm;
class JPDF;

//! SimplifiedLangevin : Hydro
class SimplifiedLangevin : public Hydro {

  public:
    //! Constructor
    explicit SimplifiedLangevin(const Base& base, real* const particles) :
      Hydro(base, particles),
      m_C0(base.control.get<ctr::param, ctr::slm, ctr::c0>()) {
      //ErrChk on m_C0
    }

    //! Destructor
    ~SimplifiedLangevin() noexcept override = default;

    //! Initialize particles
    void init() override;

    //! Advance particles
    void advance(const real& dt) override;

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

} // namespace quinoa

#endif // SLM_h
