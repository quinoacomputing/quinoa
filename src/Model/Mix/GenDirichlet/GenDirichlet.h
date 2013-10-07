//******************************************************************************
/*!
  \file      src/Model/Mix/GenDirichlet/GenDirichlet.h
  \author    J. Bakosi
  \date      Mon Oct  7 10:34:18 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     The generalized Dirichlet mix model
  \details   The generalized Dirichlet mix model
*/
//******************************************************************************
#ifndef GenDirichlet_h
#define GenDirichlet_h

#include <vector>

#include <Macro.h>
#include <Mix/Mix.h>

namespace quinoa {

class Memory;
class Paradigm;
class JPDF;

//! GeneralizedDirichlet : Mix
class GeneralizedDirichlet : public Mix {

  public:
    //! Constructor
    explicit GeneralizedDirichlet(const Base& base, tk::real* const particles) :
      Mix(base, particles),
      m_b(base.control.get<ctr::param, ctr::gendirichlet, ctr::b>()),
      m_S(base.control.get<ctr::param, ctr::gendirichlet, ctr::S>()),
      m_k(base.control.get<ctr::param, ctr::gendirichlet, ctr::kappa>()),
      m_c(base.control.get<ctr::param, ctr::gendirichlet, ctr::c>()) {
      ErrChk(m_b.size() == static_cast<unsigned int>(m_nscalar),
             tk::ExceptType::FATAL,
             "Wrong number of generalized Dirichlet model parameters 'b'");
      ErrChk(m_S.size() == static_cast<unsigned int>(m_nscalar),
             tk::ExceptType::FATAL, 
             "Wrong number of generalized Dirichlet model parameters 'S'");
      ErrChk(m_k.size() == static_cast<unsigned int>(m_nscalar),
             tk::ExceptType::FATAL,
             "Wrong number of generalized Dirichlet model parameters 'k'");
      ErrChk(m_c.size() == static_cast<unsigned int>(m_nscalar*(m_nscalar-1)/2),
             tk::ExceptType::FATAL,
             "Wrong number of generalized Dirichlet model parameters 'c'");
    }


    //! Destructor
    ~GeneralizedDirichlet() noexcept override = default;

    //! Initialize particles
    void init(int p, int tid) override { initZero(p); IGNORE(tid); }

    //! Advance particles
    void advance(int p, int tid, tk::real dt) override;

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

    const std::vector<tk::real> m_b;         //!< SDE coefficients
    const std::vector<tk::real> m_S;
    const std::vector<tk::real> m_k;
    const std::vector<tk::real> m_c;
};

} // quinoa::

#endif // GenDirichlet_h
