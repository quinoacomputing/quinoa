//******************************************************************************
/*!
  \file      src/Model/Mix/GenDirichlet/GenDirichlet.h
  \author    J. Bakosi
  \date      Wed Sep  4 07:57:12 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     The generalized Dirichlet mix model
  \details   The generalized Dirichlet mix model
*/
//******************************************************************************
#ifndef GenDirichlet_h
#define GenDirichlet_h

#include <vector>

#include <Macro.h>
#include <Mix.h>

namespace quinoa {

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
                                  const QuinoaControl& control,
                                  real* const particles) :
      Mix<GeneralizedDirichlet>(memory, paradigm, control, particles),
      m_b(control.get<control::parameter>().get<control::gendirichlet>().b),
      m_S(control.get<control::parameter>().get<control::gendirichlet>().S),
      m_k(control.get<control::parameter>().get<control::gendirichlet>().kappa),
      m_c(control.get<control::parameter>().get<control::gendirichlet>().c) {
      // Error out if mix model selected at compile time does not match that
      // whose options are given in control file
      //control->matchModels<select::Mix, select::MixType, control::MIX>(
      //  select::MixType::GENERALIZED_DIRICHLET);
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
    ~GeneralizedDirichlet() noexcept override = default;

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

} // namespace quinoa

#endif // GenDirichlet_h
