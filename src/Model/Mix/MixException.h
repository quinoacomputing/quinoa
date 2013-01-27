//******************************************************************************
/*!
  \file      src/Model/Mix/MixException.h
  \author    J. Bakosi
  \date      Sun 27 Jan 2013 12:18:31 PM MST
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
enum MixExceptType { BAD_NSCALARS=0,          //!< Wrong number of scalars
                     NUM_MIX_EXCEPT
};

//! Mix exception error messages
const string MixMsg[NUM_MIX_EXCEPT] = {
  "Wrong number of scalars"
};

//! MixException : ModelException
class MixException : public ModelException {

  public:
    //! Constructor
    MixException(ExceptType except,
                      MixExceptType mixExcept,
                      const string& file,
                      const string& func,
                      const unsigned int& line) :
      ModelException(except, MIX_EXCEPT, file, func, line),
      m_except(mixExcept) {}

    //! Move constructor, necessary for throws, default compiler generated
    MixException(MixException&&) = default;

    //! Don't permit copy constructor
    // ICC: should be deleted and private
    MixException(const MixException&);

    //! Destructor
    //virtual ~MixException() {}

    //! Handle MixException
    virtual ErrCode handleException(Driver* driver);

  private:
    //! Don't permit copy assignment
    MixException& operator=(const MixException&) = delete;
    //! Don't permit move assignment
    MixException& operator=(MixException&&) = delete;

    //! Mix exception type (BAD_SCALARS, etc.)
    MixExceptType m_except;
};

} // namespace Quinoa

#endif // MixException_h
