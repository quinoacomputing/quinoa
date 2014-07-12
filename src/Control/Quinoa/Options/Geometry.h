//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/Geometry.h
  \author    J. Bakosi
  \date      Mon Oct  7 09:16:38 2013
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Geometry options and associations
  \details   Geometry options and associations
*/
//******************************************************************************
#ifndef QuinoaGeometryOptions_h
#define QuinoaGeometryOptions_h

#include <map>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace ctr {

//! Geometry definition types
enum class GeometryType : uint8_t { NO_GEOMETRY=0,
                                    ANALYTIC,
                                    DISCRETE };

//! Class with base templated on the above enum class with associations
class Geometry : public tk::Toggle<GeometryType> {

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
    const kw::analytic_geometry ag {};
    const kw::discrete_geometry dg {};

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

} // ctr::
} // quinoa::

#endif // QuinoaGeometryOptions_h
