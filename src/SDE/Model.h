//******************************************************************************
/*!
  \file      src/SDE/Model.h
  \author    J. Bakosi
  \date      Tue Jan 14 08:54:50 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Model
  \details   Model
*/
//******************************************************************************
#ifndef Model_h
#define Model_h

#include <Types.h>

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
