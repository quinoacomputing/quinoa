//******************************************************************************
/*!
  \file      HomogeneousGeneralizedDirichlet.h
  \author    J. Bakosi
  \date      Sun 20 Jan 2013 01:54:24 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous generalized Dirichlet model
  \details   Homogeneous generalized Dirichlet model
*/
//******************************************************************************
#ifndef HomogeneousGeneralizedDirichlet_h
#define HomogeneousGeneralizedDirichlet_h

#include <limits>

#include <Physics.h>

namespace Quinoa {

class Memory;
class MixModel;
class Paradigm;
class MKLRandom;

//! HomogeneousGeneralizedDirichlet : Physics
class HomogeneousGeneralizedDirichlet : public Physics {

  public:
    //! Constructor
    HomogeneousGeneralizedDirichlet(Memory* memory,
                                    Paradigm* paradigm,
                                    const int nscalar,
                                    const real time,
                                    const int echo = 1,
                                    const int nstep =
                                      numeric_limits<int>::max());

    //! Destructor
    virtual ~HomogeneousGeneralizedDirichlet();

    //! Echo informaion on model
    virtual void echo();

    //! Initialize model
    virtual void init();

    //! Solve model
    virtual void solve();

  private:
    //! Don't permit copy constructor
    HomogeneousGeneralizedDirichlet(const HomogeneousGeneralizedDirichlet&) =
      delete;
    //! Don't permit copy assigment
    HomogeneousGeneralizedDirichlet& operator=
      (const HomogeneousGeneralizedDirichlet&) = delete;
    //! Don't permit move constructor
    HomogeneousGeneralizedDirichlet(HomogeneousGeneralizedDirichlet&&) = delete;
    //! Don't permit move assigment
    HomogeneousGeneralizedDirichlet& operator=
      (HomogeneousGeneralizedDirichlet&&) = delete;

    const int m_nscalar;          //!< Number of mixing scalars
    MKLRandom* m_random;          //!< Pointer to random number generator object
    MixModel* m_mixModel;         //!< Pointer to MixModel object
};

} // namespace Quinoa

#endif // HomogeneousGeneralizedDirichlet_h
