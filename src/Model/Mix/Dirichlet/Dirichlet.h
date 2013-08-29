//******************************************************************************
/*!
  \file      src/Model/Mix/Dirichlet/Dirichlet.h
  \author    J. Bakosi
  \date      Wed 28 Aug 2013 08:58:56 PM MDT
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
                       QuinoaControl* const control,
                       real* const particles) :
      Mix<Dirichlet>(memory, paradigm, control, particles),
      m_b(control->get<control::B>()),
      m_S(control->get<control::S>()),
      m_k(control->get<control::KAPPA>()) {
      // Error out if mix model selected at compile time does not match that
      // whose options are given in control file
      control->matchModels<select::Mix, select::MixType, control::MIX>(
        select::MixType::DIRICHLET);
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
    virtual ~Dirichlet() noexcept = default;

    //! Initialize particles
    void init(int p, int tid) { initZero(p); IGNORE(tid); }

    //! Advance particles
    void advance(int p, int tid, real dt);

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

} // namespace Quinoa

#endif // Dirichlet_h
