//******************************************************************************
/*!
  \file      src/SDE/Model.h
  \author    J. Bakosi
  \date      Mon 20 Jan 2014 05:04:39 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Model
  \details   Model
*/
//******************************************************************************
#ifndef Model_h
#define Model_h

#include <Types.h>
#include <Options/RNG.h>

namespace quinoa {

//! Model
class Model {

  public:
    //! Constructor
    explicit Model() = default;

    //! Advance particles
    virtual void advance(int p, int tid, tk::real dt) = 0;

    //! Destructor
    virtual ~Model() noexcept = default;

    //! Return true if model is stochastic (false if deterministic)
    virtual bool stochastic() const noexcept = 0;

    //! Return RNG type if stochastic, NO_RNG if deterministic
    virtual tk::ctr::RNGType rng() const noexcept = 0;

  private:
    //! Don't permit copy constructor
    Model(const Model&) = delete;
    //! Don't permit copy assigment
    Model& operator=(const Model&) = delete;
    //! Don't permit move constructor
    Model(Model&&) = delete;
    //! Don't permit move assigment
    Model& operator=(Model&&) = delete;
};

} // quinoa::

#endif // Model_h
