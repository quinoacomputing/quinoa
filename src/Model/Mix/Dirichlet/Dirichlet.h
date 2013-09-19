//******************************************************************************
/*!
  \file      src/Model/Mix/Dirichlet/Dirichlet.h
  \author    J. Bakosi
  \date      Thu Sep 19 09:14:59 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet mix model
  \details   Dirichlet mix model
*/
//******************************************************************************
#ifndef Dirichlet_h
#define Dirichlet_h

#include <vector>

#include <Macro.h>
#include <Mix.h>

namespace quinoa {

class Memory;
class Paradigm;
class JPDF;

//! Dirichlet : Mix
class Dirichlet : public Mix {

  public:
    //! Constructor
    explicit Dirichlet(const Base& base, real* const particles) :
      Mix(base, particles),
      m_b(base.control.get<ctr::param, ctr::dirichlet, ctr::b>()),
      m_S(base.control.get<ctr::param, ctr::dirichlet, ctr::S>()),
      m_k(base.control.get<ctr::param, ctr::dirichlet, ctr::kappa>()) {
      ErrChk(m_b.size() == static_cast<unsigned int>(m_nscalar),
             ExceptType::FATAL,
             "Wrong number of Dirichlet model parameters 'b'");
      ErrChk(m_S.size() == static_cast<unsigned int>(m_nscalar),
             ExceptType::FATAL,
             "Wrong number of Dirichlet model parameters 'S'");
      ErrChk(m_k.size() == static_cast<unsigned int>(m_nscalar),
             ExceptType::FATAL,
             "Wrong number of Dirichlet model parameters 'k'");
    }

    //! Destructor
    ~Dirichlet() noexcept override = default;

    //! Initialize particles
    void init(int p, int tid) override { initZero(p); IGNORE(tid); }

    //! Advance particles
    void advance(int p, int tid, real dt) override;

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

    const std::vector<real> m_b;         //!< SDE coefficients
    const std::vector<real> m_S;
    const std::vector<real> m_k;
};

} // namespace quinoa

#endif // Dirichlet_h
