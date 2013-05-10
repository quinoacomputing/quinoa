//******************************************************************************
/*!
  \file      src/Model/Mix/GeneralizedDirichlet/GeneralizedDirichlet.h
  \author    J. Bakosi
  \date      Fri May 10 16:59:47 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     The generalized Dirichlet mix model
  \details   The generalized Dirichlet mix model
*/
//******************************************************************************
#ifndef GeneralizedDirichlet_h
#define GeneralizedDirichlet_h

#include <vector>

#include <mkl_vsl.h>

#include <Memory.h>
#include <Mix.h>

namespace Quinoa {

class Memory;
class Paradigm;
class MKLRandom;
class MKLRndStream;
class MemoryEntry;
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
    virtual ~GeneralizedDirichlet() noexcept;

    //! Initialize particles
    void init();

    //! Advance particles
    void advance(const real& dt);

    //! Echo information on Dirichlet model
    void echo() const;

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

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    void finalize() noexcept;

    //! Initialize scalars with unirom PDF with the last constrained
    void initUniform();

    const vector<real> m_b;         //!< SDE coefficients
    const vector<real> m_S;
    const vector<real> m_k;
    const vector<real> m_c;
    const VSLStreamStatePtr* m_str; //!< Array of MKL VSL stream state pointers

    MKLRandom* m_random;            //!< Random number generator object
    MKLRndStream* m_rndStr;         //!< Random number stream object
};

} // namespace Quinoa

#endif // GeneralizedDirichlet_h
