//******************************************************************************
/*!
  \file      src/Model/HomDirichlet/HomDirichlet.h
  \author    J. Bakosi
  \date      Thu Nov 15 15:30:38 2012
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
class Dirichlet;

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

    MKLRandom* m_random;          //!< Pointer to random number generator object
    Dirichlet* m_dir;             //!< Pointer to Dirichlet mix model object
};

} // namespace Quinoa

#endif // HomDirichlet_h
