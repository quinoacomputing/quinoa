//******************************************************************************
/*!
  \file      src/Model/Mix/GeneralizedDirichlet/GeneralizedDirichlet.h
  \author    J. Bakosi
  \date      Tue Jul  2 16:06:07 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     The generalized Dirichlet mix model
  \details   The generalized Dirichlet mix model
*/
//******************************************************************************
#ifndef GeneralizedDirichlet_h
#define GeneralizedDirichlet_h

#include <vector>

#include <Macro.h>
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
                                  real* const particles) :
      Mix<GeneralizedDirichlet>(memory, paradigm, control, particles),
      m_b(control->get<control::B>()),
      m_S(control->get<control::S>()),
      m_k(control->get<control::KAPPA>()),
      m_c(control->get<control::C>()) {
      ErrChk(m_b.size() == static_cast<unsigned int>(m_nscalar),
             ExceptType::FATAL,
             "Wrong number of generalized Dirichlet model parameters 'b'");
      ErrChk(m_S.size() == static_cast<unsigned int>(m_nscalar),
             ExceptType::FATAL, 
             "Wrong number of generalized Dirichlet model parameters 'S'");
      ErrChk(m_k.size() == static_cast<unsigned int>(m_nscalar),
             ExceptType::FATAL,
             "Wrong number of generalized Dirichlet model parameters 'k'");
      ErrChk(m_c.size() == static_cast<unsigned int>(m_nscalar*(m_nscalar-1)/2),
             ExceptType::FATAL,
             "Wrong number of generalized Dirichlet model parameters 'c'");
    }


    //! Destructor
    virtual ~GeneralizedDirichlet() noexcept = default;

    //! Return mix model identification
    select::MixTypes id() const noexcept {
      return select::MixTypes::GENERALIZED_DIRICHLET;
    }

    //! Initialize particles
    void init(int p, int tid) { initZero(p); IGNORE(tid); }

    //! Advance particles
    void advance(int p, int tid, real dt);

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

    const std::vector<real> m_b;         //!< SDE coefficients
    const std::vector<real> m_S;
    const std::vector<real> m_k;
    const std::vector<real> m_c;
};

} // namespace Quinoa

#endif // GeneralizedDirichlet_h
