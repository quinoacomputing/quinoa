//******************************************************************************
/*!
  \file      src/Model/Mix/Dirichlet/Dirichlet.h
  \author    J. Bakosi
  \date      Sat 02 Mar 2013 10:03:54 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet mix model
  \details   Dirichlet mix model
*/
//******************************************************************************
#ifndef Dirichlet_h
#define Dirichlet_h

#include <mkl_vsl.h>

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
    Dirichlet(Memory* const memory,
              Paradigm* const paradigm,
              Control* const control);

    //! Destructor
    virtual ~Dirichlet();

    //! Echo information on Dirichlet model
    virtual void echo();

    //! Initialize particles
    virtual void init();

    //! Advance particles
    virtual void advance(const real dt);

    //! Estimate joint scalar PDF
    virtual void jpdf(JPDF& jpdf);

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

    MemoryEntry* m_allScalars;      //!< Memory entry storing all the scalars
    MemoryEntry* m_parScalars;      //!< Memory entry storing particle scalars
    real* m_rawAllScalars;          //!< Raw pointer to all scalars
    real* m_rawParScalars;          //!< Raw pointer to particle scalars

    const vector<real> m_b;         //!< SDE coefficients
    const vector<real> m_S;
    const vector<real> m_k;
};

} // namespace Quinoa

#endif // Dirichlet_h
