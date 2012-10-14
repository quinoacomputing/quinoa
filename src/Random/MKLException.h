//******************************************************************************
/*!
  \file      src/Random/MKLException.h
  \author    J. Bakosi
  \date      Sat 13 Oct 2012 10:54:53 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKLException class declaration
  \details   MKLException class declaration
*/
//******************************************************************************
#ifndef MKLException_h
#define MKLException_h

#include <string>
#include <map>

using namespace std;

#include <RandomException.h>

namespace Quinoa {

//! MKL exception types
enum MKLExceptType { MKLEXCEPTION_UNIMPLEMENTED=0,
                     MKL_UNIMPLEMENTED,
                     MKL_UNKNOWN,
                     MKL_BADARGS,
                     MKL_MEM_FAILURE,
                     MKL_NULL_PTR,
                     MKL_INVALID_BRNG_INDEX,
                     MKL_LEAPFROG_UNSUPPORTED,
                     MKL_SKIPAHEAD_UNSUPPORTED,
                     MKL_BRNGS_INCOMPATIBLE,
                     MKL_BAD_STREAM,
                     MKL_BRNG_TABLE_FULL,
                     MKL_BAD_STREAM_STATE_SIZE,
                     MKL_BAD_WORD_SIZE,
                     MKL_BAD_NSEEDS,
                     MKL_BAD_NBITS,
                     MKL_BAD_UPDATE,
                     MKL_NO_NUMBERS,
                     MKL_INVALID_ABSTRACT_STREAM,
                     MKL_ERROR_FILE_CLOSE,
                     MKL_ERROR_FILE_OPEN,
                     MKL_ERROR_FILE_WRITE,
                     MKL_ERROR_FILE_READ,
                     MKL_BAD_FILE_FORMAT,
                     MKL_UNSUPPORTED_FILE_VER,
                     NUM_MKL_EXCEPT
};

//! MKL exception error messages
const string MKLMsg[NUM_MKL_EXCEPT] = {
  "MKL exception type unimplemented",
  "VSL feature not yet implemented",
  "VSL unknown error",
  "VSL bad arguments",
  "VSL memory allocation problem",
  "VSL null pointer",
  "VSL invalid BRNG index",
  "VSL LeapFrog initialization unsupported",
  "VSL SkipAhead initialization unsupported",
  "VSL BRNGs are not compatible for the operation",
  "VSL random stream is invalid",
  "VSL table of registered BRNGs is full",
  "VSL value in StreamStateSize field is bad",
  "VSL value in WordSize field is bad",
  "VSL value in NSeeds field is bad",
  "VSL value in NBits field is bad",
  "VSL number of updated entries in buffer is invalid",
  "VSL zero number of updated entries in buffer",
  "VSL abstract random stream is invalid",
  "VSL can`t close file",
  "VSL can`t open file",
  "VSL can`t write to file",
  "VSL can`t read from file",
  "VSL file format is unknown",
  "VSL unsupported file version"
};

//! MKLException : RandomException
class MKLException : protected RandomException {

  public:
    //! Constructor: fill VSLErrMap
    MKLException(ExceptType except, Int vslerr);

    //! Move constructor, necessary for throws, default compiler generated
    MKLException(MKLException&&) = default;

    //! Destructor
    ~MKLException() = default;

    //! Handle MKLException
    ErrCode handleException(Driver* driver);

    //! Get MKLException based on VSLError
    MKLExceptType getException(Int vslerr);

  private:
    //! Don't permit copy constructor
    MKLException(const MKLException&) = delete;
    //! Don't permit copy assignment
    MKLException& operator=(const MKLException&) = delete;
    //! Don't permit move assignment
    MKLException& operator=(MKLException&&) = delete;

    //! MKL exception type (MKL_UNIMPLEMENTED, MKL_UNKNOWN, etc.)
    MKLExceptType m_except;

    //! VSLError -> MKLExceptType map
    map<Int,MKLExceptType> m_VSLErrMap;
};

} // namespace Quinoa

#endif // MKLException_h
