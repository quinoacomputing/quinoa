//******************************************************************************
/*!
  \file      src/Model/Mass/Beta/Beta.h
  \author    J. Bakosi
  \date      Mon Oct 28 07:27:14 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Beta mass model
  \details   Beta mass model
*/
//******************************************************************************
#ifndef Beta_h
#define Beta_h

#include <Mass/Mass.h>

namespace quinoa {

//! Beta : Mass
class Beta : public Mass {

  public:
    //! Constructor
    explicit Beta() {}
//     explicit Beta(const Base& base, tk::real* const particles) :
//       Mass(base, particles),
//       m_At(base.control.get<ctr::param, ctr::beta, ctr::atwood>()) {
//       // ErrChk on m_At
//     }

    //! Destructor
    ~Beta() noexcept override = default;

//     //! Initialize particles
//     void init() override;
// 
//     //! Advance particles
//     void advance(int p, int tid, tk::real dt) override;
// 
//     //! Estimate joint scalar PDF
//     void jpdf(tk::JPDF& jpdf);

  private:
    //! Don't permit copy constructor
    Beta(const Beta&) = delete;
    //! Don't permit copy assigment
    Beta& operator=(const Beta&) = delete;
    //! Don't permit move constructor
    Beta(Beta&&) = delete;
    //! Don't permit move assigment
    Beta& operator=(Beta&&) = delete;

//    const tk::real m_At;            //!< Atwood-number
};

} // quinoa::

#endif // Beta_h
