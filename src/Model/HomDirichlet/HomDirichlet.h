//******************************************************************************
/*!
  \file      src/Model/HomDirichlet/HomDirichlet.h
  \author    J. Bakosi
  \date      Fri 16 Nov 2012 09:49:32 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous Dirichlet model
  \details   Homogeneous Dirichlet model
*/
//******************************************************************************
#ifndef HomDirichlet_h
#define HomDirichlet_h

#include <Model.h>

namespace Quinoa {

class Memory;
class Paradigm;
class MKLRandom;
class MKLRndStream;
class Dirichlet;

//! HomDirichlet : Model
class HomDirichlet : public Model {

  public:
    //! Constructor
    HomDirichlet(Memory* memory,
                 Paradigm* paradigm,
                 const int& nscalar,
                 const int& npar,
                 const real time,
                 const int nstep = numeric_limits<int>::max());

    //! Destructor
    virtual ~HomDirichlet();

    //! Echo informaion on model
    virtual void echo();

    //! Initialize model
    virtual void init();

    //! Solve model
    virtual void solve();

  private:
    //! Don't permit copy constructor
    HomDirichlet(const HomDirichlet&) = delete;
    //! Don't permit copy assigment
    HomDirichlet& operator=(const HomDirichlet&) = delete;
    //! Don't permit move constructor
    HomDirichlet(HomDirichlet&&) = delete;
    //! Don't permit move assigment
    HomDirichlet& operator=(HomDirichlet&&) = delete;

    //! Initialize scalars with unirom PDF with the last constrained
    void initUniform();

    //! One-liner report
    void report(const int it, const real t, const real dt,
                long int& hrs2beg, long int& mins2beg, long int& secs2beg,
                long int& hrs2end, long int& mins2end, long int& secs2end);

    MKLRandom* m_random;          //!< Pointer to random number generator object
    MKLRndStream* m_rndStr;       //!< Pointer to random number stream object
    Dirichlet* m_dir;             //!< Pointer to Dirichlet mix model object
    const int m_nscalar;          //!< Number of mixing scalars
    const int m_npar;             //!< Number of particles
    MemoryEntry* m_scalar;        //!< Memory entry storing the scalars
    struct timeval m_startTime;   //!< Date/time loop started
};

} // namespace Quinoa

#endif // HomDirichlet_h
