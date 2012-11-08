//******************************************************************************
/*!
  \file      src/Random/VSLException.h
  \author    J. Bakosi
  \date      Wed Nov  7 18:03:42 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     VSLException class declaration
  \details   VSLException class declaration
*/
//******************************************************************************
#ifndef VSLException_h
#define VSLException_h

#include <string>
#include <map>

using namespace std;

#include <MKLException.h>

namespace Quinoa {

//! VSL exception types
enum VSLExceptType { VSL_UNIMPLEMENTED=0,
                     VSL_UNKNOWN,
                     VSL_BADARGS,
                     VSL_MEM_FAILURE,
                     VSL_NULL_PTR,
                     VSL_INVALID_BRNG_INDEX,
                     VSL_LEAPFROG_UNSUPPORTED,
                     VSL_SKIPAHEAD_UNSUPPORTED,
                     VSL_BRNGS_INCOMPATIBLE,
                     VSL_BAD_STREAM,
                     VSL_BRNG_TABLE_FULL,
                     VSL_BAD_STREAM_STATE_SIZE,
                     VSL_BAD_WORD_SIZE,
                     VSL_BAD_NSEEDS,
                     VSL_BAD_NBITS,
                     VSL_BAD_UPDATE,
                     VSL_NO_NUMBERS,
                     VSL_INVALID_ABSTRACT_STREAM,
                     VSL_FILE_CLOSE,
                     VSL_FILE_OPEN,
                     VSL_FILE_WRITE,
                     VSL_FILE_READ,
                     VSL_BAD_FILE_FORMAT,
                     VSL_UNSUPPORTED_FILE_VER,
                     NUM_VSL_EXCEPT
};

//! VSL exception error messages
const string VSLMsg[NUM_VSL_EXCEPT] = {
  "feature not yet implemented",
  "unknown error",
  "bad arguments",
  "memory allocation problem",
  "null pointer",
  "invalid BRNG index",
  "LeapFrog initialization unsupported",
  "SkipAhead initialization unsupported",
  "BRNGs are not compatible for the operation",
  "Random stream is invalid",
  "table of registered BRNGs is full",
  "value in StreamStateSize field is bad",
  "value in WordSize field is bad",
  "value in NSeeds field is bad",
  "value in NBits field is bad",
  "number of updated entries in buffer is invalid",
  "zero number of updated entries in buffer",
  "abstract random stream is invalid",
  "cannot close file",
  "cannot open file",
  "cannot write to file",
  "cannot read from file",
  "file format is unknown",
  "unsupported file version"
};

//! VSLException : RandomException
class VSLException : public MKLException {

  public:
    //! Constructor: fill VSLErrMap
    VSLException(ExceptType except, int vslerr);

    //! Move constructor, necessary for throws, default compiler generated
    VSLException(VSLException&&) = default;

    //! Destructor
    ~VSLException() = default;

    //! Handle VSLException
    virtual ErrCode handleException(Driver* driver);

    //! Get VSLException based on VSLError
    VSLExceptType getException(int vslerr);

  private:
    //! Don't permit copy constructor
    VSLException(const VSLException&) = delete;
    //! Don't permit copy assignment
    VSLException& operator=(const VSLException&) = delete;
    //! Don't permit move assignment
    VSLException& operator=(VSLException&&) = delete;

    //! VSL exception type (VSL_UNIMPLEMENTED, VSL_UNKNOWN, etc.)
    VSLExceptType m_except;

    //! VSLError -> VSLExceptType map
    map<int,VSLExceptType> m_VSLErrMap;
};

} // namespace Quinoa

#endif // VSLException_h
