//******************************************************************************
/*!
  \file      src/Physics/HomMix/HomMix.h
  \author    J. Bakosi
  \date      Mon 04 Feb 2013 09:30:20 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous material mix model
  \details   Homogeneous material mix model
*/
//******************************************************************************
#ifndef HomMix_h
#define HomMix_h

#include <limits>

#include <mkl_vsl.h>

#include <Physics.h>

namespace Quinoa {

class Memory;
class Paradigm;
class MKLRandom;
class MKLRndStream;
class Dirichlet;

//! HomMix : Physics
class HomMix : public Physics {

  public:
    //! Constructor
    HomMix(Memory* memory,
           Paradigm* paradigm,
           const int& nscalar,
           const int& npar,
           const real time,
           const int echo = 1,
           const int nstep = numeric_limits<int>::max());

    //! Destructor
    virtual ~HomMix();

    //! Echo informaion on model
    virtual void echo();

    //! Initialize model
    virtual void init();

    //! Solve model
    virtual void solve();

  private:
    //! Don't permit copy constructor
    HomMix(const HomMix&) = delete;
    //! Don't permit copy assigment
    HomMix& operator=(const HomMix&) = delete;
    //! Don't permit move constructor
    HomMix(HomMix&&) = delete;
    //! Don't permit move assigment
    HomMix& operator=(HomMix&&) = delete;

    //! Initialize scalars with unirom PDF with the last constrained
    void initUniform();

    //! Initialize scalars with Gaussian PDF
    void initGaussian();

    //! Advance particles
    void advance(const real dt);

    //! Output joint scalar PDF
    void outJPDF();

    MKLRandom* m_random;            //!< Random number generator object
    MKLRndStream* m_rndStr;         //!< Random number stream object
    Dirichlet* m_dir;               //!< Dirichlet mix model object
    const int m_N;                  //!< Number of mixing scalars
    const int m_npar;               //!< Number of particles
    MemoryEntry* m_MEscalar;        //!< Memory entry storing the scalars
    real* m_scalar;                 //!< Raw pointer to scalars
    const VSLStreamStatePtr* m_str; //!< Array of MKL VSL stream state pointers
};

} // namespace Quinoa

#endif // HomMix_h
