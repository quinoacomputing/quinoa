//******************************************************************************
/*!
  \file      src/Model/Mass/Beta/Beta.h
  \author    J. Bakosi
  \date      Mon 30 Sep 2013 08:35:19 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Beta mass model
  \details   Beta mass model
*/
//******************************************************************************
#ifndef Beta_h
#define Beta_h

#include <Mass/Mass.h>

namespace quinoa {

class Memory;
class Paradigm;
class JPDF;

//! Beta : Mass
class Beta : public Mass {

  public:
    //! Constructor
    explicit Beta(const Base& base, real* const particles) :
      Mass(base, particles),
      m_At(base.control.get<ctr::param, ctr::beta, ctr::atwood>()) {
      // ErrChk on m_At
    }

    //! Destructor
    ~Beta() noexcept override = default;

    //! Initialize particles
    void init() override;

    //! Advance particles
    void advance(int p, int tid, real dt) override;

    //! Estimate joint scalar PDF
    void jpdf(JPDF& jpdf);

  private:
    //! Don't permit copy constructor
    Beta(const Beta&) = delete;
    //! Don't permit copy assigment
    Beta& operator=(const Beta&) = delete;
    //! Don't permit move constructor
    Beta(Beta&&) = delete;
    //! Don't permit move assigment
    Beta& operator=(Beta&&) = delete;

    const real m_At;            //!< Atwood-number
};

} // namespace quinoa

#endif // Beta_h
