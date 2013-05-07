//******************************************************************************
/*!
  \file      src/Physics/SPINSFlow/SPINSFlow.h
  \author    J. Bakosi
  \date      Tue May  7 10:49:49 2013
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
class Control;
class MKLRandom;
class MKLRndStream;
class UnsMesh;
class Hydro;
class Timer;

//! SPINSFlow : Physics
class SPINSFlow : public Physics {

  public:
    //! Constructor
    explicit SPINSFlow(Memory* const memory,
                       Paradigm* const paradigm,
                       Control* const control,
                       Timer* const timer,
                       const string& filename);

    //! Destructor
    virtual ~SPINSFlow() noexcept;

    //! Echo informaion on model
    virtual void echo() const;

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

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    void finalize() noexcept;

    //! Advance particles
    void advance(const real dt);

    MKLRandom* m_random;            //!< Random number generator object
    MKLRndStream* m_rndStr;         //!< Random number stream object
    Hydro* m_hydro;                 //!< Hydro model object
    UnsMesh* m_mesh;                //!< Unstructured mesh object
    const int m_npar;               //!< Number of particles
    const VSLStreamStatePtr* m_str; //!< Array of MKL VSL stream state pointers
    const string m_filename;        //!< Unstructured mesh object
};

} // namespace Quinoa

#endif // SPINSFlow_h
