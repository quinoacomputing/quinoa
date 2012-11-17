//******************************************************************************
/*!
  \file      src/Model/HomGenDirichlet/HomGenDirichlet.h
  \author    J. Bakosi
  \date      Sat 17 Nov 2012 08:15:15 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous generalized Dirichlet model
  \details   Homogeneous generalized Dirichlet model
*/
//******************************************************************************
#ifndef HomGenDirichlet_h
#define HomGenDirichlet_h

#include <limits>

#include <Model.h>

namespace Quinoa {

class Memory;
class MixModel;
class Paradigm;
class MKLRandom;

//! HomGenDirichlet : Model
class HomGenDirichlet : public Model {

  public:
    //! Constructor
    HomGenDirichlet(Memory* memory,
                    Paradigm* paradigm,
                    const int nscalar,
                    const real time,
                    const int echo = 1,
                    const int nstep = numeric_limits<int>::max());

    //! Destructor
    virtual ~HomGenDirichlet();

    //! Echo informaion on model
    virtual void echo();

    //! Initialize model
    virtual void init();

    //! Solve model
    virtual void solve();

  private:
    //! Don't permit copy constructor
    HomGenDirichlet(const HomGenDirichlet&) = delete;
    //! Don't permit copy assigment
    HomGenDirichlet& operator=(const HomGenDirichlet&) = delete;
    //! Don't permit move constructor
    HomGenDirichlet(HomGenDirichlet&&) = delete;
    //! Don't permit move assigment
    HomGenDirichlet& operator=(HomGenDirichlet&&) = delete;

    const int m_nscalar;          //!< Number of mixing scalars
    MKLRandom* m_random;          //!< Pointer to random number generator object
    MixModel* m_mixModel;         //!< Pointer to MixModel object
};

} // namespace Quinoa

#endif // HomGenDirichlet_h
