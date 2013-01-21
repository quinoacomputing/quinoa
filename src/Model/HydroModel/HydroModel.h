//******************************************************************************
/*!
  \file      src/Model/HydroModel/HydroModel.h
  \author    J. Bakosi
  \date      Mon 21 Jan 2013 10:31:16 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     HydroModel base
  \details   HydroModel base
*/
//******************************************************************************
#ifndef HydroModel_h
#define HydroModel_h

#include <string>

#include <QuinoaTypes.h>

namespace Quinoa {

using namespace std;

//! HydroModel base
class HydroModel {

  public:
    //! Constructor
    HydroModel(const string& name);

    //! Destructor
    virtual ~HydroModel() {}

    //! Interface for echo information on mix model
    virtual void echo() = 0;

  protected:
    const string m_name;           //!< Name of mix model

  private:
    //! Don't permit copy constructor
    HydroModel(const HydroModel&) = delete;
    //! Don't permit copy assigment
    HydroModel& operator=(const HydroModel&) = delete;
    //! Don't permit move constructor
    HydroModel(HydroModel&&) = delete;
    //! Don't permit move assigment
    HydroModel& operator=(HydroModel&&) = delete;
};

} // namespace Quinoa

#endif // HydroModel_h
