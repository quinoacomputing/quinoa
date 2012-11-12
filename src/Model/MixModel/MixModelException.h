//******************************************************************************
/*!
  \file      src/Model/MixModel/MixModelException.h
  \author    J. Bakosi
  \date      Mon 12 Nov 2012 10:06:17 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mix model exception handler
  \details   Mix model exception handler
*/
//******************************************************************************
#ifndef MixModelException_h
#define MixModelException_h

#include <string>

using namespace std;

#include <ModelException.h>

namespace Quinoa {

//! MixModel exception types
enum MixModelExceptType { BAD_SCALARS=0,          //!< Wrong number of scalars
                          NUM_MIXMODEL_EXCEPT
};

//! MixModel exception error messages
const string MixModelMsg[NUM_MIXMODEL_EXCEPT] = {
  "Wrong number of scalars"
};

//! MixModelException : ModelException
class MixModelException : public ModelException {

  public:
    //! Constructor
    MixModelException(ExceptType except,
                      MixModelExceptType mixModelExcept,
                      const string& file,
                      const string& func,
                      const unsigned int& line) :
      ModelException(except, MIXMODEL_EXCEPT, file, func, line),
      m_except(mixModelExcept) {}

    //! Move constructor, necessary for throws, default compiler generated
    MixModelException(MixModelException&&) = default;

    //! Don't permit copy constructor
    // ICC: should be deleted and private
    MixModelException(const MixModelException&);

    //! Destructor
    virtual ~MixModelException() {}

    //! Handle MixModelException
    virtual ErrCode handleException(Driver* driver);

  private:
    //! Don't permit copy assignment
    MixModelException& operator=(const MixModelException&) = delete;
    //! Don't permit move assignment
    MixModelException& operator=(MixModelException&&) = delete;

    //! MixModel exception type (BAD_SCALARS, etc.)
    MixModelExceptType m_except;
};

} // namespace Quinoa

#endif // MixModelException_h
