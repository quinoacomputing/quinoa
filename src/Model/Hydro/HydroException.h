//******************************************************************************
/*!
  \file      src/Model/Hydro/HydroException.h
  \author    J. Bakosi
  \date      Sat 30 Mar 2013 11:48:15 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Hydro model exception handler
  \details   Hydro model exception handler
*/
//******************************************************************************
#ifndef HydroException_h
#define HydroException_h

#include <string>

using namespace std;

#include <ModelException.h>

namespace Quinoa {

//! Hydro exception types
enum HydroExceptType { NO_SUCH_HYDRO=0,          //!< No such hydro model
                       HYDRO_UNIMPLEMENTED,      //!< Hydro model unimplemented
                       NUM_HYDRO_EXCEPT
};

//! Hydro exception error messages
const string HydroMsg[NUM_HYDRO_EXCEPT] = {
  "No such hydrodynamics model",
  "Hydro model not implemented"
};

//! HydroException : ModelException
class HydroException : public ModelException {

  public:
    //! Constructor
    HydroException(ExceptType except,
                   HydroExceptType mixExcept,
                   const string& file,
                   const string& func,
                   const unsigned int& line) :
      ModelException(except, HYDROMODEL_EXCEPT, file, func, line),
      m_except(mixExcept) {}

    //! Move constructor, necessary for throws, default compiler generated
    HydroException(HydroException&&) = default;

    //! Don't permit copy constructor
    // ICC: should be deleted and private
    HydroException(const HydroException&);

    //! Destructor
    //virtual ~HydroException() {}

    //! Handle HydroException
    virtual ErrCode handleException(Driver* driver);

  private:
    //! Don't permit copy assignment
    HydroException& operator=(const HydroException&) = delete;
    //! Don't permit move assignment
    HydroException& operator=(HydroException&&) = delete;

    //! Hydro exception type
    HydroExceptType m_except;
};

} // namespace Quinoa

#endif // HydroException_h
