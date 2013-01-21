//******************************************************************************
/*!
  \file      src/Model/HydroModel/HydroModelException.h
  \author    J. Bakosi
  \date      Mon 21 Jan 2013 10:47:45 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Hydro model exception handler
  \details   Hydro model exception handler
*/
//******************************************************************************
#ifndef HydroModelException_h
#define HydroModelException_h

#include <string>

using namespace std;

#include <ModelException.h>

namespace Quinoa {

//! HydroModel exception types
enum HydroModelExceptType { BAD_HYDRO=0,          //!< Placeholder for now
                            NUM_HYDROMODEL_EXCEPT
};

//! HydroModel exception error messages
const string HydroModelMsg[NUM_HYDROMODEL_EXCEPT] = {
  "Bad hydro"
};

//! HydroModelException : ModelException
class HydroModelException : public ModelException {

  public:
    //! Constructor
    HydroModelException(ExceptType except,
                      HydroModelExceptType mixModelExcept,
                      const string& file,
                      const string& func,
                      const unsigned int& line) :
      ModelException(except, HYDROMODEL_EXCEPT, file, func, line),
      m_except(mixModelExcept) {}

    //! Move constructor, necessary for throws, default compiler generated
    HydroModelException(HydroModelException&&) = default;

    //! Don't permit copy constructor
    // ICC: should be deleted and private
    HydroModelException(const HydroModelException&);

    //! Destructor
    virtual ~HydroModelException() {}

    //! Handle HydroModelException
    virtual ErrCode handleException(Driver* driver);

  private:
    //! Don't permit copy assignment
    HydroModelException& operator=(const HydroModelException&) = delete;
    //! Don't permit move assignment
    HydroModelException& operator=(HydroModelException&&) = delete;

    //! HydroModel exception type
    HydroModelExceptType m_except;
};

} // namespace Quinoa

#endif // HydroModelException_h
