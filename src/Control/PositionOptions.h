//******************************************************************************
/*!
  \file      src/Control/PositionOptions.h
  \author    J. Bakosi
  \date      Mon 27 May 2013 02:44:02 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Position model options and associations
  \details   Position model options and associations
*/
//******************************************************************************
#ifndef PositionOptions_h
#define PositionOptions_h

#include <map>

#include <Exception.h>

namespace Quinoa {

namespace select {

class Position {

  public:
    //! Position model types
    enum class Enum { NO_POSITION=0,
                      INVISCID,
                      VISCOUS };

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
    Position() :
      //! Enums -> names
      names{ { Enum::NO_POSITION, "No position" },
             { Enum::INVISCID, "Inviscid"},
             { Enum::VISCOUS, "Viscous"} },
      //! keywords -> Enums
      values{ { "no_position", Enum::NO_POSITION},
              { "invpos",      Enum::INVISCID},
              { "vispos",      Enum::INVISCID } } {}

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
    Position(const Position&) = delete;
    //! Don't permit copy assigment
    Position& operator=(const Position&) = delete;
    //! Don't permit move constructor
    Position(Position&&) = delete;
    //! Don't permit move assigment
    Position& operator=(Position&&) = delete;

    std::map<Enum, std::string> names;
    std::map<std::string, Enum> values;
};

} // namespace select

} // namespace Quinoa

#endif // PositionOptions_h
