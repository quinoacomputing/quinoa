//******************************************************************************
/*!
  \file      src/Physics/SPINSFlow/SPINSFlow.h
  \author    J. Bakosi
  \date      Sun 20 Jan 2013 05:47:26 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Standalone-Particle Incompressible Navier-Stokes Flow
  \details   Standalone-Particle Incompressible Navier-Stokes Flow
*/
//******************************************************************************
#ifndef SPINSFlow_h
#define SPINSFlow_h

#include <limits>

#include <mkl_vsl.h>

#include <Physics.h>

namespace Quinoa {

class Memory;
class Paradigm;
class MKLRandom;
class MKLRndStream;

//! SPINSFlow : Physics
class SPINSFlow : public Physics {

  public:
    //! Constructor
    SPINSFlow(Memory* memory,
              Paradigm* paradigm,
              const int& npar,
              const real time,
              const int echo = 1,
              const int nstep = numeric_limits<int>::max());

    //! Destructor
    virtual ~SPINSFlow();

    //! Echo informaion on model
    virtual void echo();

    //! Initialize model
    virtual void init();

    //! Solve model
    virtual void solve();

  private:
    //! Don't permit copy constructor
    SPINSFlow(const SPINSFlow&) = delete;
    //! Don't permit copy assigment
    SPINSFlow& operator=(const SPINSFlow&) = delete;
    //! Don't permit move constructor
    SPINSFlow(SPINSFlow&&) = delete;
    //! Don't permit move assigment
    SPINSFlow& operator=(SPINSFlow&&) = delete;

    //! Advance particles
    void advance(const real dt);

    MKLRandom* m_random;            //!< Ptr to random number generator object
    MKLRndStream* m_rndStr;         //!< Ptr to random number stream object
    const int m_npar;               //!< Number of particles
    const VSLStreamStatePtr* m_str; //!< Array of MKL VSL stream state pointers
};

} // namespace Quinoa

#endif // SPINSFlow_h
