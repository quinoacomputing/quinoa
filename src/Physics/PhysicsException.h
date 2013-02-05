//******************************************************************************
/*!
  \file      src/Physics/PhysicsException.h
  \author    J. Bakosi
  \date      Mon 04 Feb 2013 09:54:05 PM MST
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
    PhysicsException(ExceptType except,
                    PhysicsExceptType physicsExcept,
                    const string& file,
                    const string& func,
                    const unsigned int& line) :
      Exception(except, file, func, line), m_except(physicsExcept) {}

    //! Move constructor, necessary for throws, default compiler generated
    PhysicsException(PhysicsException&&) = default;

    //! Don't permit copy constructor
    // ICC: should be deleted and private
    PhysicsException(const PhysicsException&);

    //! Destructor
    //virtual ~PhysicsException() {}

    //! Handle PhysicsException
    virtual ErrCode handleException(Driver* driver);

  private:
    //! Don't permit copy assignment
    PhysicsException& operator=(const PhysicsException&) = delete;
    //! Don't permit move assignment
    PhysicsException& operator=(PhysicsException&&) = delete;

    //! Physics exception type
    PhysicsExceptType m_except;
};

} // namespace Quinoa

#endif // PhysicsException_h
