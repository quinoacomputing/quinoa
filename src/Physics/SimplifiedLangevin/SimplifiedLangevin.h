//******************************************************************************
/*!
  \file      src/Physics/SimplifiedLangevin/SimplifiedLangevin.h
  \author    J. Bakosi
  \date      Sun 20 Jan 2013 01:19:34 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     The simplified Langevin model
  \details   The simplified Langevin model
*/
//******************************************************************************
#ifndef SimplifiedLangevin_h
#define SimplifiedLangevin_h

#include <limits>

#include <mkl_vsl.h>

#include <Physics.h>

namespace Quinoa {

class Memory;
class Paradigm;
class MKLRandom;
class MKLRndStream;

//! SimplifiedLangevin : Physics
class SimplifiedLangevin : public Physics {

  public:
    //! Constructor
    SimplifiedLangevin(Memory* memory,
                       Paradigm* paradigm,
                       const int& npar,
                       const real time,
                       const int echo = 1,
                       const int nstep = numeric_limits<int>::max());

    //! Destructor
    virtual ~SimplifiedLangevin();

    //! Echo informaion on model
    virtual void echo();

    //! Initialize model
    virtual void init();

    //! Solve model
    virtual void solve();

  private:
    //! Don't permit copy constructor
    SimplifiedLangevin(const SimplifiedLangevin&) = delete;
    //! Don't permit copy assigment
    SimplifiedLangevin& operator=(const SimplifiedLangevin&) = delete;
    //! Don't permit move constructor
    SimplifiedLangevin(SimplifiedLangevin&&) = delete;
    //! Don't permit move assigment
    SimplifiedLangevin& operator=(SimplifiedLangevin&&) = delete;

    //! Advance particles
    void advance(const real dt);

    MKLRandom* m_random;            //!< Ptr to random number generator object
    MKLRndStream* m_rndStr;         //!< Ptr to random number stream object
    const int m_npar;               //!< Number of particles
    const VSLStreamStatePtr* m_str; //!< Array of MKL VSL stream state pointers
};

} // namespace Quinoa

#endif // SimplifiedLangevin_h
