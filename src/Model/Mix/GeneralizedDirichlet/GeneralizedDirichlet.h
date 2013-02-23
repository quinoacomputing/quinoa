//******************************************************************************
/*!
  \file      src/Model/Mix/GeneralizedDirichlet/GeneralizedDirichlet.h
  \author    J. Bakosi
  \date      Sat 23 Feb 2013 12:10:11 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     The generalized Dirichlet mix model
  \details   The generalized Dirichlet mix model
*/
//******************************************************************************
#ifndef GeneralizedDirichlet_h
#define GeneralizedDirichlet_h

#include <mkl_vsl.h>

#include <Mix.h>

namespace Quinoa {

class Memory;
class Paradigm;
class MKLRandom;
class MKLRndStream;
class MemoryEntry;
class JPDF;

//! GeneralizedDirichlet : Mix
class GeneralizedDirichlet : public Mix {

  public:
    //! Constructor
    GeneralizedDirichlet(Memory* const memory,
                         Paradigm* const paradigm,
                         Control* const control);

    //! Destructor
    virtual ~GeneralizedDirichlet();

    //! Echo information on Generalized Dirichlet model
    virtual void echo();

    //! Initialize particles
    virtual void init();

    //! Advance particles
    virtual void advance(const real dt);

    //! Estimate joint scalar PDF
    virtual void jpdf(JPDF& jpdf);

  private:
    //! Don't permit copy constructor
    GeneralizedDirichlet(const GeneralizedDirichlet&) = delete;
    //! Don't permit copy assigment
    GeneralizedDirichlet& operator=(const GeneralizedDirichlet&) = delete;
    //! Don't permit move constructor
    GeneralizedDirichlet(GeneralizedDirichlet&&) = delete;
    //! Don't permit move assigment
    GeneralizedDirichlet& operator=(GeneralizedDirichlet&&) = delete;

    //! Initialize scalars with unirom PDF with the last constrained
    void initUniform();

    const VSLStreamStatePtr* m_str; //!< Array of MKL VSL stream state pointers
    MKLRandom* m_random;            //!< Random number generator object
    MKLRndStream* m_rndStr;         //!< Random number stream object
    MemoryEntry* m_MEscalar;        //!< Memory entry storing the scalars
    real* m_scalar;                 //!< Raw pointer to scalars
    const vector<real> m_b;         //!< SDE coefficients
    const vector<real> m_S;
    const vector<real> m_k;
    const vector<real> m_c;
};

} // namespace Quinoa

#endif // GeneralizedDirichlet_h
