//******************************************************************************
/*!
  \file      src/Control/HydroOptions.h
  \author    J. Bakosi
  \date      Mon 27 May 2013 02:45:30 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Hydro model options and associations
  \details   Hydro model options and associations
*/
//******************************************************************************
#ifndef HydroOptions_h
#define HydroOptions_h

#include <map>

#include <Exception.h>

namespace Quinoa {

namespace select {

class Hydro {

  public:
    //! Hydro model types
    enum class Enum { NO_HYDRO=0,
                      SLM,
                      GLM };

    //! Operator << for writing Enum to output streams
    friend std::ostream& operator<< (std::ostream& os, const Enum& e) {
      os << static_cast<std::underlying_type<Enum>::type>(e);
      return os;
    }

    //! Operator + for adding Enum to a std::string
    friend std::string operator+ (const std::string& lhs, Enum e) {
      std::stringstream ss;
      ss << lhs << e ;
      std::string rhs = ss.str();
      return rhs;
    }

    //! Constructor initializing associations
    Hydro() :
      //! Enums -> names
      names{ { Enum::NO_HYDRO, "No hydro" },
             { Enum::SLM, "Simplified Langevin"},
             { Enum::GLM, "Generalized Langevin" } },
      //! keywords -> Enums
      values{ { "no_hydro",  Enum::NO_HYDRO},
              { "slm",       Enum::SLM},
              { "glm",       Enum::GLM} } {}

    //! Lookup Enum value based on keyword
    Enum value(const std::string keyword) const {
      auto it = values.find(keyword);
      Assert(it != values.end(), FATAL,
             "Cannot find value for keyword \"" + keyword + "\"");
      return it->second;
    }

    //! Lookup option name based on Enum
    const std::string& name(Enum value) const {
      auto it = names.find(value);
      Assert(it != names.end(), FATAL,
             std::string("Cannot find name for value \"") + value + "\"");
      return it->second;
    }

  private:
    //! Don't permit copy constructor
    Hydro(const Hydro&) = delete;
    //! Don't permit copy assigment
    Hydro& operator=(const Hydro&) = delete;
    //! Don't permit move constructor
    Hydro(Hydro&&) = delete;
    //! Don't permit move assigment
    Hydro& operator=(Hydro&&) = delete;

    std::map<Enum, std::string> names;
    std::map<std::string, Enum> values;
};

} // namespace select

} // namespace Quinoa

#endif // HydroOptions_h
