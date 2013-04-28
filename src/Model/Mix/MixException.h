//******************************************************************************
/*!
  \file      src/Model/Mix/MixException.h
  \author    J. Bakosi
  \date      Sat 27 Apr 2013 08:27:26 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mix model exception handler
  \details   Mix model exception handler
*/
//******************************************************************************
#ifndef MixException_h
#define MixException_h

#include <string>

using namespace std;

#include <ModelException.h>

namespace Quinoa {

//! Mix exception types
enum MixExceptType { BAD_NSCALAR=0,          //!< Wrong number of scalars
                     BAD_MODEL_PARAMETERS,   //!< Wrong number of model params
                     NO_SUCH_MIX,            //!< No mix model selected
                     MIX_UNIMPLEMENTED,      //!< Mix model unimplemented
                     NUM_MIX_EXCEPT
};

//! Mix exception error messages
const string MixMsg[NUM_MIX_EXCEPT] = {
  "Wrong number of scalars",
  "Wrong number of model parameters",
  "No mix model selected",
  "Mix model not implemented"
};

//! MixException : ModelException
class MixException : public ModelException {

  public:
    //! Constructor
    explicit MixException(const ExceptType except,
                          const MixExceptType mixExcept,
                          const string& file,
                          const string& func,
                          const unsigned int& line) noexcept :
      ModelException(except,
                     MIX_EXCEPT,
                     file,
                     func,
                     line,
                     MixMsg[static_cast<int>(mixExcept)]) {}

    //! Move constructor for throws, default compiler generated
    MixException(MixException&&) = default;

    //! Destructor
    virtual ~MixException() noexcept = default;

  private:
    //! Don't permit copy constructor
    MixException(const MixException&) = delete;
    //! Don't permit copy assignment
    MixException& operator=(const MixException&) = delete;
    //! Don't permit move assignment
    MixException& operator=(MixException&&) = delete;
};

} // namespace Quinoa

#endif // MixException_h
