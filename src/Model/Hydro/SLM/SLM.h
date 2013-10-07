//******************************************************************************
/*!
  \file      src/Model/Hydro/SLM/SLM.h
  \author    J. Bakosi
  \date      Mon Oct  7 10:32:09 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Simplified Langevin hydrodynamics model
  \details   Simplified Langevin hydrodynamics model
*/
//******************************************************************************
#ifndef SLM_h
#define SLM_h

#include <Hydro/Hydro.h>

namespace quinoa {

//! SimplifiedLangevin : Hydro
class SimplifiedLangevin : public Hydro {

  public:
    //! Constructor
    explicit SimplifiedLangevin(const Base& base, tk::real* const particles) :
      Hydro(base, particles),
      m_C0(base.control.get<ctr::param, ctr::slm, ctr::c0>()) {
      //ErrChk on m_C0
    }

    //! Destructor
    ~SimplifiedLangevin() noexcept override = default;

    //! Initialize particles
    void init() override;

    //! Advance particles
    void advance(int p, int tid, tk::real dt) override;

  private:
    //! Don't permit copy constructor
    SimplifiedLangevin(const SimplifiedLangevin&) = delete;
    //! Don't permit copy assigment
    SimplifiedLangevin& operator=(const SimplifiedLangevin&) = delete;
    //! Don't permit move constructor
    SimplifiedLangevin(SimplifiedLangevin&&) = delete;
    //! Don't permit move assigment
    SimplifiedLangevin& operator=(SimplifiedLangevin&&) = delete;

    const tk::real m_C0;                //!< Parameter C0 in SLM
};

} // quinoa::

#endif // SLM_h
