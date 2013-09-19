//******************************************************************************
/*!
  \file      src/Control/GeometryOptions.h
  \author    J. Bakosi
  \date      Thu Sep 19 09:42:05 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Geometry options and associations
  \details   Geometry options and associations
*/
//******************************************************************************
#ifndef GeometryOptions_h
#define GeometryOptions_h

#include <map>

#include <Exception.h>
#include <Toggle.h>

namespace quinoa {
namespace sel {

//! Geometry definition types
enum class GeometryType : uint8_t { NO_GEOMETRY=0,
                                    ANALYTIC,
                                    DISCRETE };

//! Class with base templated on the above enum class with associations
class Geometry : public Toggle<GeometryType> {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit Geometry() : Toggle<GeometryType>(names, values) {}

  private:
    //! Don't permit copy constructor
    Geometry(const Geometry&) = delete;
    //! Don't permit copy assigment
    Geometry& operator=(const Geometry&) = delete;
    //! Don't permit move constructor
    Geometry(Geometry&&) = delete;
    //! Don't permit move assigment
    Geometry& operator=(Geometry&&) = delete;

    //! Enums -> names
    const std::map<GeometryType, std::string> names {
      { GeometryType::NO_GEOMETRY, "" },
      { GeometryType::ANALYTIC, "Analytic" },
      { GeometryType::DISCRETE, "Discrete" }
    };

    //! keywords -> Enums
    const std::map<std::string, GeometryType> values {
      { "no_geometry", GeometryType::NO_GEOMETRY },
      { "analytic_geometry", GeometryType::ANALYTIC },
      { "discrete_geometry", GeometryType::DISCRETE }
    };
};

} // sel::
} // quinoa::

#endif // GeometryOptions_h
