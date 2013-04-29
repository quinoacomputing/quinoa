//******************************************************************************
/*!
  \file      src/Random/RandomException.h
  \author    J. Bakosi
  \date      Mon Apr 29 15:56:03 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RandomException class declaration
  \details   RandomException class declaration
*/
//******************************************************************************
#ifndef RandomException_h
#define RandomException_h

#include <string>

using namespace std;

#include <Exception.h>

namespace Quinoa {

//! RandomException types
enum RndExceptType { RND_MKL=0,             //!< MKL error
                     RND_UNIMPLEMENTED,     //!< Unimplemented exception
                     NUM_RND_EXCEPT
};

//! Random exception error messages
const string RndMsg[NUM_RND_EXCEPT] = {
  "MKL exception: ",
  "Unimplemented random number generator exception"
};

//! RandomException : Exception
class RandomException : public Exception {

  public:
    //! Constructor
    explicit RandomException(const ExceptType except,
                             const RndExceptType rndExcept,
                             const string& file,
                             const string& func,
                             unsigned int line,
                             const string& message = "") noexcept :
      Exception(except,
                RndMsg[static_cast<int>(rndExcept)] + message,
                file,
                func,
                line) {}

    //! Destructor
    virtual ~RandomException() noexcept = default;

  protected:
    //! Permit copy constructor only for children
    RandomException(const RandomException&) = default;
    //! Permit move constructor only for children
    RandomException(RandomException&&) = default;

  private:
    //! Don't permit copy assignment
    RandomException& operator=(const RandomException&) = delete;
    //! Don't permit move assignment
    RandomException& operator=(RandomException&&) = delete;
};

} // namespace Quinoa

#endif // RandomException_h
