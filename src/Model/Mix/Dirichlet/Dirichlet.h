//******************************************************************************
/*!
  \file      src/Model/Mix/Dirichlet/Dirichlet.h
  \author    J. Bakosi
  \date      Fri May 10 17:11:58 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet mix model
  \details   Dirichlet mix model
*/
//******************************************************************************
#ifndef Dirichlet_h
#define Dirichlet_h

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
    virtual ~Dirichlet() noexcept;

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
    Dirichlet(const Dirichlet&) = delete;
    //! Don't permit copy assigment
    Dirichlet& operator=(const Dirichlet&) = delete;
    //! Don't permit move constructor
    Dirichlet(Dirichlet&&) = delete;
    //! Don't permit move assigment
    Dirichlet& operator=(Dirichlet&&) = delete;

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    void finalize() noexcept;

    //! Initialize scalars with unirom PDF with the last constrained
    void initUniform();

    //! Initialize scalars with Gaussian PDF
    void initGaussian();

    const vector<real> m_b;         //!< SDE coefficients
    const vector<real> m_S;
    const vector<real> m_k;
    const VSLStreamStatePtr* m_str; //!< Array of MKL VSL stream state pointers

    MKLRandom* m_random;            //!< Random number generator object
    MKLRndStream* m_rndStr;         //!< Random number stream object
};

} // namespace Quinoa

#endif // Dirichlet_h
