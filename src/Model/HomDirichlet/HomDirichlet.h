//******************************************************************************
/*!
  \file      src/Model/HomDirichlet/HomDirichlet.h
  \author    J. Bakosi
  \date      Fri Nov 16 08:23:16 2012
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
                 const int& npar);

    //! Destructor
    virtual ~HomDirichlet();

    //! Echo informaion on model
    virtual void echo();

    //! Initialize model
    virtual void init();

  private:
    //! Don't permit copy constructor
    HomDirichlet(const HomDirichlet&) = delete;
    //! Don't permit copy assigment
    HomDirichlet& operator=(const HomDirichlet&) = delete;
    //! Don't permit move constructor
    HomDirichlet(HomDirichlet&&) = delete;
    //! Don't permit move assigment
    HomDirichlet& operator=(HomDirichlet&&) = delete;

    //! Initialize with N-leak delta
    void initNpeakDelta();

    MKLRandom* m_random;          //!< Pointer to random number generator object
    MKLRndStream* m_rndStr;       //!< Pointer to random number stream object
    Dirichlet* m_dir;             //!< Pointer to Dirichlet mix model object
    const int m_nscalar;          //!< Number of mixing scalars
    const int m_npar;             //!< Number of particles
    MemoryEntry* m_scalar;        //!< Memory entry storing the scalars
};

} // namespace Quinoa

#endif // HomDirichlet_h
