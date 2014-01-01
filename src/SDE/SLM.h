//******************************************************************************
/*!
  \file      src/SDE/SLM.h
  \author    J. Bakosi
  \date      Wed 01 Jan 2014 01:19:39 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Simplified Langevin hydrodynamics model
  \details   Simplified Langevin hydrodynamics model
*/
//******************************************************************************
#ifndef SLM_h
#define SLM_h

#include <Hydro.h>

namespace quinoa {

//! SLM : Hydro
class SLM : public Hydro {

  public:
    //! Constructor
    explicit SLM() = default;
//     explicit SLM(const Base& base, tk::real* const particles) :
//       Hydro(base, particles),
//       m_C0(base.control.get<ctr::param, ctr::slm, ctr::c0>()) {
//       //ErrChk on m_C0
//     }

    //! Destructor
    ~SLM() noexcept override = default;

//     //! Initialize particles
//     void init() override;
// 
//     //! Advance particles
//     void advance(int p, int tid, tk::real dt) override;

  private:
    //! Don't permit copy constructor
    SLM(const SLM&) = delete;
    //! Don't permit copy assigment
    SLM& operator=(const SLM&) = delete;
    //! Don't permit move constructor
    SLM(SLM&&) = delete;
    //! Don't permit move assigment
    SLM& operator=(SLM&&) = delete;

//     const tk::real m_C0;                //!< Parameter C0 in SLM
};

} // quinoa::

#endif // SLM_h
