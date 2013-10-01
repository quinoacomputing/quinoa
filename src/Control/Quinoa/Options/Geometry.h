//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Geometry.h
  \author    J. Bakosi
  \date      Mon 30 Sep 2013 09:30:43 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Geometry options and associations
  \details   Geometry options and associations
*/
//******************************************************************************
#ifndef QuinoaGeometryOptions_h
#define QuinoaGeometryOptions_h

#include <map>

#include <Exception.h>
#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

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
    explicit Geometry() : Toggle<GeometryType>("Geometry", names, values) {}

  private:
    //! Don't permit copy constructor
    Geometry(const Geometry&) = delete;
    //! Don't permit copy assigment
    Geometry& operator=(const Geometry&) = delete;
    //! Don't permit move constructor
    Geometry(Geometry&&) = delete;
    //! Don't permit move assigment
    Geometry& operator=(Geometry&&) = delete;

    //! Get access to geometry keywords
    const grm::kw::analytic_geometry ag {};
    const grm::kw::discrete_geometry dg {};

    //! Enums -> names
    const std::map<GeometryType, std::string> names {
      { GeometryType::NO_GEOMETRY, "n/a" },
      { GeometryType::ANALYTIC, ag.name() },
      { GeometryType::DISCRETE, dg.name() }
    };

    //! keywords -> Enums
    const std::map<std::string, GeometryType> values {
      { "no_geometry", GeometryType::NO_GEOMETRY },
      { ag.string(), GeometryType::ANALYTIC },
      { dg.string(), GeometryType::DISCRETE }
    };
};

} // sel::
} // quinoa::

#endif // QuinoaGeometryOptions_h
