//******************************************************************************
/*!
  \file      src/Physics/PhysicsException.h
  \author    J. Bakosi
  \date      Mon Apr 29 15:55:45 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics Exception handling
  \details   Physics Exception handling
*/
//******************************************************************************
#ifndef PhysicsException_h
#define PhysicsException_h

#include <string>

using namespace std;

#include <Exception.h>

namespace Quinoa {

//! Physics exception types
enum PhysicsExceptType { PHYSICS_UNIMPLEMENTED=0,    //!< Physics unimplemented
                         NO_PHYSICS,                 //!< No physics selected
                         NUM_PHYSICS_EXCEPT
};

//! Physics exception error messages
const string PhysicsMsg[NUM_PHYSICS_EXCEPT] = {
  "Selected physics not implemented",
  "No physics selected"
};

//! PhysicsException : Exception
class PhysicsException : public Exception {

  public:
    //! Constructor
    explicit PhysicsException(const ExceptType except,
                              const PhysicsExceptType physicsExcept,
                              const string& file,
                              const string& func,
                              unsigned int line) noexcept :
      Exception(except,
                PhysicsMsg[static_cast<int>(physicsExcept)],
                file,
                func,
                line) {}

    //! Move constructor for throws, default compiler generated
    PhysicsException(PhysicsException&&) = default;

    //! Destructor
    virtual ~PhysicsException() noexcept = default;

  private:
    //! Don't permit copy constructor
    PhysicsException(const PhysicsException&) = default;
    //! Don't permit copy assignment
    PhysicsException& operator=(const PhysicsException&) = delete;
    //! Don't permit move assignment
    PhysicsException& operator=(PhysicsException&&) = delete;
};

} // namespace Quinoa

#endif // PhysicsException_h
