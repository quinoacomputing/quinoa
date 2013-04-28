//******************************************************************************
/*!
  \file      src/Model/Hydro/HydroException.h
  \author    J. Bakosi
  \date      Sat 27 Apr 2013 08:27:46 PM MDT
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
enum HydroExceptType { NO_SUCH_HYDRO=0,        //!< No such hydro model
                       BAD_NPROP,              //!< Wrong number of properties
                       HYDRO_UNIMPLEMENTED,    //!< Hydro model unimplemented
                       NUM_HYDRO_EXCEPT
};

//! Hydro exception error messages
const string HydroMsg[NUM_HYDRO_EXCEPT] = {
  "No such hydrodynamics model",
  "Wrong number of particle properties",
  "Hydro model not implemented"
};

//! HydroException : ModelException
class HydroException : public ModelException {

  public:
    //! Constructor
    explicit HydroException(const ExceptType except,
                            const HydroExceptType hydroExcept,
                            const string& file,
                            const string& func,
                            const unsigned int& line) noexcept :
      ModelException(except,
                     HYDROMODEL_EXCEPT,
                     file,
                     func,
                     line,
                     HydroMsg[static_cast<int>(hydroExcept)]) {}

    //! Move constructor for throws, default compiler generated
    HydroException(HydroException&&) = default;

    //! Destructor
    virtual ~HydroException() noexcept = default;

  private:
    //! Don't permit copy constructor
    HydroException(const HydroException&) = delete;
    //! Don't permit copy assignment
    HydroException& operator=(const HydroException&) = delete;
    //! Don't permit move assignment
    HydroException& operator=(HydroException&&) = delete;
};

} // namespace Quinoa

#endif // HydroException_h
