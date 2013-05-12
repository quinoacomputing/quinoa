//******************************************************************************
/*!
  \file      src/Model/Mix/GeneralizedDirichlet/GeneralizedDirichlet.h
  \author    J. Bakosi
  \date      Sun 12 May 2013 03:40:49 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     The generalized Dirichlet mix model
  \details   The generalized Dirichlet mix model
*/
//******************************************************************************
#ifndef GeneralizedDirichlet_h
#define GeneralizedDirichlet_h

#include <vector>

#include <Mix.h>

namespace Quinoa {

class Memory;
class Paradigm;
class JPDF;

//! GeneralizedDirichlet : Mix<GeneralizedDirichlet> child for CRTP
//! See: http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
class GeneralizedDirichlet : public Mix<GeneralizedDirichlet> {

  public:
    //! Constructor
    explicit GeneralizedDirichlet(Memory* const memory,
                                  Paradigm* const paradigm,
                                  Control* const control,
                                  real* const scalars);

    //! Destructor
    virtual ~GeneralizedDirichlet() noexcept = default;

    //! Initialize particles
    void init();

    //! Advance particles
    void advance(const real& dt);

    //! Estimate joint scalar PDF
    void jpdf(JPDF& jpdf);

  private:
    //! Don't permit copy constructor
    GeneralizedDirichlet(const GeneralizedDirichlet&) = delete;
    //! Don't permit copy assigment
    GeneralizedDirichlet& operator=(const GeneralizedDirichlet&) = delete;
    //! Don't permit move constructor
    GeneralizedDirichlet(GeneralizedDirichlet&&) = delete;
    //! Don't permit move assigment
    GeneralizedDirichlet& operator=(GeneralizedDirichlet&&) = delete;

    const vector<real> m_b;         //!< SDE coefficients
    const vector<real> m_S;
    const vector<real> m_k;
    const vector<real> m_c;
};

} // namespace Quinoa

#endif // GeneralizedDirichlet_h
