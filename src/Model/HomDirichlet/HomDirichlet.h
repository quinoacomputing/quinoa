//******************************************************************************
/*!
  \file      src/Model/HomDirichlet/HomDirichlet.h
  \author    J. Bakosi
  \date      Thu Nov 15 15:07:42 2012
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
class MixModel;
class Paradigm;
class MKLRandom;

//! HomDirichlet : Model
class HomDirichlet : public Model {

  public:
    //! Constructor
    HomDirichlet(Memory* memory, Paradigm* paradigm, const int nscalar);

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

    const int m_nscalar;          //!< Number of mixing scalars
    MKLRandom* m_random;          //!< Pointer to random number generator object
    MixModel* m_mixModel;         //!< Pointer to MixModel object
};

} // namespace Quinoa

#endif // HomDirichlet_h
