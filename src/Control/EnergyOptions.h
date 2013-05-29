//******************************************************************************
/*!
  \file      src/Control/EnergyOptions.h
  \author    J. Bakosi
  \date      Wed May 29 07:10:32 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Energy model options and associations
  \details   Energy model options and associations
*/
//******************************************************************************
#ifndef EnergyOptions_h
#define EnergyOptions_h

#include <map>

#include <Exception.h>
#include <Toggle.h>

namespace Quinoa {

namespace select {

//! Energy model types
enum class EnergyTypes { NO_ENERGY=0,
                       };

//! Class with base templated on the above enum class with associations
class Energy : public Toggle<EnergyTypes> {

  public:
    //! Constructor initializing associations
    // ICC: use initializer lists
    Energy() : Toggle<EnergyTypes>(names, values) {
      //! Enums -> names
      names[EnergyTypes::NO_ENERGY] = "No energy";
      //! keywords -> Enums
      values["no_energy"] = EnergyTypes::NO_ENERGY;
    }

  private:
    //! Don't permit copy constructor
    Energy(const Energy&) = delete;
    //! Don't permit copy assigment
    Energy& operator=(const Energy&) = delete;
    //! Don't permit move constructor
    Energy(Energy&&) = delete;
    //! Don't permit move assigment
    Energy& operator=(Energy&&) = delete;

    std::map<EnergyTypes, std::string> names;
    std::map<std::string, EnergyTypes> values;
};

} // namespace select

} // namespace Quinoa

#endif // EnergyOptions_h
