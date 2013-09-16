//******************************************************************************
/*!
  \file      src/Model/Mix/GenDirichlet/GenDirichlet.h
  \author    J. Bakosi
  \date      Sun 15 Sep 2013 05:28:21 PM MDT
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
    explicit GeneralizedDirichlet(const Base& base, real* const particles) :
      Mix<GeneralizedDirichlet>(base, particles),
      m_b(base.control.get<control::param, control::gendirichlet, control::b>()),
      m_S(base.control.get<control::param, control::gendirichlet, control::S>()),
      m_k(base.control.get<control::param, control::gendirichlet, control::kappa>()),
      m_c(base.control.get<control::param, control::gendirichlet, control::c>()) {
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
