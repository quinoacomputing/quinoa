//******************************************************************************
/*!
  \file      src/Physics/HomogeneousDirichlet/HomogeneousDirichlet.h
  \author    J. Bakosi
  \date      Mon 21 Jan 2013 10:53:11 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous Dirichlet model
  \details   Homogeneous Dirichlet model
*/
//******************************************************************************
#ifndef HomogeneousDirichlet_h
#define HomogeneousDirichlet_h

#include <limits>

#include <mkl_vsl.h>

#include <Physics.h>

namespace Quinoa {

class Memory;
class Paradigm;
class MKLRandom;
class MKLRndStream;
class Dirichlet;

//! HomogeneousDirichlet : Physics
class HomogeneousDirichlet : public Physics {

  public:
    //! Constructor
    HomogeneousDirichlet(Memory* memory,
                         Paradigm* paradigm,
                         const int& nscalar,
                         const int& npar,
                         const real time,
                         const int echo = 1,
                         const int nstep = numeric_limits<int>::max());

    //! Destructor
    virtual ~HomogeneousDirichlet();

    //! Echo informaion on model
    virtual void echo();

    //! Initialize model
    virtual void init();

    //! Solve model
    virtual void solve();

  private:
    //! Don't permit copy constructor
    HomogeneousDirichlet(const HomogeneousDirichlet&) = delete;
    //! Don't permit copy assigment
    HomogeneousDirichlet& operator=(const HomogeneousDirichlet&) = delete;
    //! Don't permit move constructor
    HomogeneousDirichlet(HomogeneousDirichlet&&) = delete;
    //! Don't permit move assigment
    HomogeneousDirichlet& operator=(HomogeneousDirichlet&&) = delete;

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

#endif // HomogeneousDirichlet_h
