//******************************************************************************
/*!
  \file      src/Model/Mix/Dirichlet/Dirichlet.h
  \author    J. Bakosi
  \date      Sun 12 May 2013 03:40:52 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet mix model
  \details   Dirichlet mix model
*/
//******************************************************************************
#ifndef Dirichlet_h
#define Dirichlet_h

#include <vector>

#include <Mix.h>

namespace Quinoa {

class Memory;
class Paradigm;
class JPDF;

//! Dirichlet : Mix<Dirichlet> child for CRTP
//! See: http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
class Dirichlet : public Mix<Dirichlet> {

  public:
    //! Constructor
    explicit Dirichlet(Memory* const memory,
                       Paradigm* const paradigm,
                       Control* const control,
                       real* const scalars);

    //! Destructor
    virtual ~Dirichlet() noexcept = default;

    //! Initialize particles
    void init();

    //! Advance particles
    void advance(const real& dt);

    //! Estimate joint scalar PDF
    void jpdf(JPDF& jpdf);

  private:
    //! Don't permit copy constructor
    Dirichlet(const Dirichlet&) = delete;
    //! Don't permit copy assigment
    Dirichlet& operator=(const Dirichlet&) = delete;
    //! Don't permit move constructor
    Dirichlet(Dirichlet&&) = delete;
    //! Don't permit move assigment
    Dirichlet& operator=(Dirichlet&&) = delete;

    const vector<real> m_b;         //!< SDE coefficients
    const vector<real> m_S;
    const vector<real> m_k;
};

} // namespace Quinoa

#endif // Dirichlet_h
