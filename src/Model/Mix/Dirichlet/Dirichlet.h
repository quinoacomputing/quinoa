//******************************************************************************
/*!
  \file      src/Model/Mix/Dirichlet/Dirichlet.h
  \author    J. Bakosi
  \date      Fri Apr 26 16:30:57 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet mix model
  \details   Dirichlet mix model
*/
//******************************************************************************
#ifndef Dirichlet_h
#define Dirichlet_h

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

//! Dirichlet : Mix
class Dirichlet : public Mix {

  public:
    //! Constructor
    explicit Dirichlet(Memory* const memory,
                       Paradigm* const paradigm,
                       Control* const control);

    //! Destructor
    virtual ~Dirichlet() noexcept;

    //! Initialize particles
    virtual void init();

    //! Advance particles
    virtual void advance(const real dt);

    //! Echo information on Dirichlet model
    virtual void echo() const;

    //! Estimate joint scalar PDF
    virtual void jpdf(JPDF& jpdf);

    //! Constant accessor to particle properties (scalars) pointer
    //! \return Number of particle scalars
    virtual const real* particles() const { return m_allScalars.ptr; }

  private:
    //! Don't permit copy constructor
    Dirichlet(const Dirichlet&) = delete;
    //! Don't permit copy assigment
    Dirichlet& operator=(const Dirichlet&) = delete;
    //! Don't permit move constructor
    Dirichlet(Dirichlet&&) = delete;
    //! Don't permit move assigment
    Dirichlet& operator=(Dirichlet&&) = delete;

    //! Initialize scalars with unirom PDF with the last constrained
    void initUniform();

    //! Initialize scalars with Gaussian PDF
    void initGaussian();

    const VSLStreamStatePtr* m_str; //!< Array of MKL VSL stream state pointers
    MKLRandom* m_random;            //!< Random number generator object
    MKLRndStream* m_rndStr;         //!< Random number stream object
    Data<real> m_allScalars;        //!< All particle scalars
    const vector<real> m_b;         //!< SDE coefficients
    const vector<real> m_S;
    const vector<real> m_k;
};

} // namespace Quinoa

#endif // Dirichlet_h
