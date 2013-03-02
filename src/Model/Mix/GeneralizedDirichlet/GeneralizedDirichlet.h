//******************************************************************************
/*!
  \file      src/Model/Mix/GeneralizedDirichlet/GeneralizedDirichlet.h
  \author    J. Bakosi
  \date      Sat 02 Mar 2013 10:21:20 AM MST
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

    MemoryEntry* m_allScalars;      //!< Memory entry storing all the scalars
    MemoryEntry* m_parScalars;      //!< Memory entry storing particle scalars
    MemoryEntry* m_Y;               //!< Memory entry storing particle scalars
    MemoryEntry* m_U;               //!< Memory entry storing particle scalars
    real* m_rawAllScalars;          //!< Raw pointer to all scalars
    real* m_rawParScalars;          //!< Raw pointer to particle scalars
    real* m_rawY;                   //!< Raw pointer to particle scalars
    real* m_rawU;                   //!< Raw pointer to particle scalars

    const vector<real> m_b;         //!< SDE coefficients
    const vector<real> m_S;
    const vector<real> m_k;
    const vector<real> m_c;
};

} // namespace Quinoa

#endif // GeneralizedDirichlet_h
