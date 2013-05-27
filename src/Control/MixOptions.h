//******************************************************************************
/*!
  \file      src/Control/MixOptions.h
  \author    J. Bakosi
  \date      Mon 27 May 2013 02:41:33 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Mix model options and associations
  \details   Mix model options and associations
*/
//******************************************************************************
#ifndef MixOptions_h
#define MixOptions_h

#include <map>
#include <sstream>

#include <Exception.h>

namespace Quinoa {

namespace select {

class Mix {

  public:
    //! Mix model types
    enum class Enum { NO_MIX=0,
                      IEM,
                      IECM,
                      DIRICHLET,
                      GENERALIZED_DIRICHLET };

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
    Mix() :
      //! Enums -> names
      names{ { Enum::NO_MIX, "No mix" },
             { Enum::IEM, "Interaction by exchange with the mean"},
             { Enum::IECM,
               "Interaction by exchange with the conditional mean" },
             { Enum::DIRICHLET, "Dirichlet"},
             { Enum::GENERALIZED_DIRICHLET, "Dirichlet"} },
      //! keywords -> Enums
      values{ { "no_mix",    Enum::NO_MIX},
              { "iem",       Enum::IEM},
              { "iecm",      Enum::IECM},
              { "dir",       Enum::DIRICHLET},
              { "gendir",    Enum::GENERALIZED_DIRICHLET} } {}

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
    Mix(const Mix&) = delete;
    //! Don't permit copy assigment
    Mix& operator=(const Mix&) = delete;
    //! Don't permit move constructor
    Mix(Mix&&) = delete;
    //! Don't permit move assigment
    Mix& operator=(Mix&&) = delete;

    std::map<Enum, std::string> names;
    std::map<std::string, Enum> values;
};

} // namespace select

} // namespace Quinoa

#endif // MixOptions_h
