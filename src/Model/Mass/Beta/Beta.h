//******************************************************************************
/*!
  \file      src/Model/Mass/Beta/Beta.h
  \author    J. Bakosi
  \date      Wed 04 Sep 2013 07:25:13 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Beta mass model
  \details   Beta mass model
*/
//******************************************************************************
#ifndef Beta_h
#define Beta_h

#include <Mass.h>

namespace quinoa {

class Memory;
class Paradigm;
class JPDF;

//! Beta : Mass<Beta> child for CRTP
//! See: http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
class Beta : public Mass<Beta> {

  public:
    //! Constructor
    explicit Beta(Memory* const memory,
                  Paradigm* const paradigm,
                  const QuinoaControl& control,
                  real* const particles) :
      Mass<Beta>(memory, paradigm, control, particles),
      m_At(control.get<control::parameter, control::beta, control::atwood>()) {
      // Error out if mass model selected at compile time does not match that
      // whose options are given in control file
      //control->matchModels<select::Mass, select::MassType, control::MASS>(
      //  select::MassType::BETA);
      // ErrChk on m_At
    }

    //! Destructor
    ~Beta() noexcept override = default;

    //! Initialize particles
    void init();

    //! Advance particles
    void advance(const real& dt);

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
