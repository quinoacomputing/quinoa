//******************************************************************************
/*!
  \file      src/Model/Model.h
  \author    J. Bakosi
  \date      Thu 08 Nov 2012 09:14:27 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Model base
  \details   Model base
*/
//******************************************************************************
#ifndef Model_h
#define Model_h

namespace Quinoa {

//! Model base
class Model {

  public:
    //! Constructor
    Model() = default;

    //! Destructor
    virtual ~Model() {};

    //! Interface for setting initial conditions
    virtual void setIC() = 0;

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

} // namespace Quinoa

#endif // Model_h
