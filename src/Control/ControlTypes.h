//******************************************************************************
/*!
  \file      src/Control/ControlTypes.h
  \author    J. Bakosi
  \date      Fri Aug  2 17:57:21 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Types for control and parsing
  \details   Types for control and parsing
*/
//******************************************************************************
#ifndef ControlTypes_h
#define ControlTypes_h

#include <string>
#include <vector>
#include <tuple>
#include <ostream>
#include <sstream>

#include <QuinoaTypes.h>
#include <GeometryOptions.h>
#include <PhysicsOptions.h>
#include <PositionOptions.h>
#include <MassOptions.h>
#include <HydroOptions.h>
#include <EnergyOptions.h>
#include <MixOptions.h>
#include <FrequencyOptions.h>
#include <MixRateOptions.h>
#include <RNGTestOptions.h>
#include <RNGOptions.h>

namespace Quinoa {

namespace control {

const int NCOMP_POS = 3;        //!< Number of position components for a field

//! Quantities whose statistics can be estimated. If you change this, make sure
//! you change Control::termOffset() as well.
enum class Quantity : uint8_t { POSITION=0,
                                DENSITY,
                                VELOCITY_X,
                                VELOCITY_Y,
                                VELOCITY_Z,
                                SCALAR
};

//! Moment specifies which type of moment is computed for a Quantity in a Term
enum class Moment : uint8_t { ORDINARY=0,      //!< Full variable
                              CENTRAL          //!< Fluctuation
};

//! Term is a Moment of a Quantity with a field ID to be ensemble averaged.
//! Internally the Numbering of field IDs starts from 0, but presented to the
//! user as starting from 1. Examples: 1st pressure fluctuation: {0, PRESSURE,
//! CENTRAL, ...}, mean of 2nd scalar: {1, TRANSPORTED_SCALAR, ORDINARY, ...}.
struct Term {
  int field;         //!< Field ID
  Quantity quantity; //!< Physical quantity
  Moment moment;     //!< Moment type: ordinary, central
  int name;          //!< Character code as name, converted only for output
  bool plot;   //!< Shows whether the variable will be plotted
  // Conceptually, plot should be in Product, since plot will only be false for
  // a mean that was triggered by a central moment by one of the Terms of a
  // Product requesting the mean or a model. However, that would require
  // Product to be a vector<struct>, which then would need custom comparitors
  // for std::sort() and std::unique() in Parser::unique(). Since this is not a
  // performance issue, plot is here in Term.

  //! Constructor
  explicit Term(int f,
                Quantity q,
                Moment m,
                char n,
                bool p) : field(f),
                          quantity(q),
                          moment(m),
                          name(n),
                          plot(p) {}

  //! Equal operator for finding unique elements, used by e.g., std::unique()
  //! Test only on field, quantity, and moment
  bool operator== (const Term& term) const {
    if (field == term.field &&
        quantity == term.quantity &&
        moment == term.moment) {
      return true;
    } else {
      return false;
    }
  }

  //! Less than operator for ordering, used by e.g., std::sort().
  //! Test on field, quantity, term, moment, and !plot.
  //! Using operator >, instead of operator <, on plot ensures that if a Term is
  //! user-requested, i.e., plotted, and also triggered by e.g., a model, the
  //! user-requested Term will take precendence.
  bool operator< (const Term& term) const {
    // test on everything except name
    if (field < term.field) {
      return true;
    } else if (field == term.field && quantity < term.quantity) {
      return true;
    } else if (field == term.field && quantity == term.quantity &&
               moment < term.moment) {
      return true;
    } else if (field == term.field && quantity == term.quantity &&
               moment == term.moment && plot > term.plot) {
      return true;
    } else {
      return false;
    }
  }

  //! Operator + for adding Term (name+field ID) to a std::string
  friend std::string operator+ (const std::string& lhs,
                                const control::Term& term) {
    std::stringstream ss;
    ss << lhs << char(term.name) << term.field+1;
    std::string rhs = ss.str();
    return rhs;
  }

  //! Operator << for writing Term to output streams
  friend std::ostream& operator<< (std::ostream& os, const Term& term) {
    os << char(term.name) << term.field+1;
    return os;
  }

  //! Operator << for writing vector<Term> to output streams
  friend std::ostream& operator<< (std::ostream& os,
                                   const std::vector<Term>& vec) {
    os << " <";
    for (auto& w : vec) os << w;
    os << ">";
    return os;
  }

  //! Operator <<= for writing requested vector<Term> to output streams
  friend std::ostream& operator<<= (std::ostream& os,
                                    const std::vector<Term>& vec) {
    if (vec[0].plot) {
      os << " <";
      for (auto& w : vec) os << w;
      os << ">";
    }
    return os;
  }
};

//! Lighter-weight structure for field names
struct FieldName {
  int name;
  int field;

  //! Constructor
  explicit FieldName(const int n = 0, const int f = 0) :
    name(n), field(f) {}

  //! Operator << for writing FieldName to output streams
  friend std::ostream& operator<< (std::ostream& os, const FieldName& fn) {
     os << char(fn.name) << fn.field+1;
     return os;
  }

  //! Operator += for adding FieldName to std::string
  friend std::string& operator+= (std::string& os, const FieldName& fn) {
     std::stringstream ss;
     ss << os << char(fn.name) << fn.field+1;
     os = ss.str();
     return os;
  }
};

//! Products are N Terms to be multiplied and ensemble averaged
//! E.g. the scalar flux in x direction needs two terms for ensemble averaging:
//! (Y-<Y>) and (U-<U>), then the moment is <yu> = <(Y-<Y>)(U-<U>)>
//! E.g the third mixed central moment of three scalars needs three terms for
//! ensemble averaging: (Y1-<Y1>), (Y2-<Y2>), and (Y3-<Y3>), then the moment is
//! <y1y2y3> = <(Y1-<Y1>)(Y2-<Y2>)(Y3-<Y3>)>
using Product = std::vector<Term>;

//! Position enum for accessing fields of tuple Bundle using names as in struct
enum BundlePosition { TITLE=0,
                      GEOMETRY,
                      PHYSICS,
                      POSITION,
                      MASS,
                      HYDRO,
                      ENERGY,
                      MIX,
                      FREQUENCY,
                      MIXRATE,
                      RNGTEST,
                      NSTEP,
                      TERM,
                      DT,
                      NPOSITION,
                      NDENSITY,
                      NVELOCITY,
                      NSCALAR,
                      NFREQUENCY,
                      NPAR,
                      TTYI,
                      DMPI,
                      STAI,
                      PDFI,
                      GLBI,
                      INPUT,
                      OUTPUT,
                      PDFNAME,
                      GLOBNAME,
                      STATNAME,
                      B,
                      S,
                      KAPPA,
                      C,
                      RNGS,
                      C0,
                      AT,
                      FREQ_GAMMA_C1,
                      FREQ_GAMMA_C2,
                      FREQ_GAMMA_C3,
                      FREQ_GAMMA_C4,
                      BOXES,
                      STATISTICS
};

//! Storage bundle for parsed data
using Bundle = std::tuple<
  std::string,             //!< Problem Title
  select::GeometryType,    //!< Selected geometry definition
  select::PhysicsType,     //!< Selected physics
  select::PositionType,    //!< Selected position model
  select::MassType,        //!< Selected mass model
  select::HydroType,       //!< Selected hydrodynamics model
  select::EnergyType,      //!< Selected internal energy model
  select::MixType,         //!< Selected material mix model
  select::FrequencyType,   //!< Selected turbulence frequency model
  select::MixRateType,     //!< Selected material mix rate model
  select::RNGTestType,     //!< Selected RNG test suite
  uint64_t,                //!< Number of time steps to take
  real,                    //!< Time to terminate time stepping
  real,                    //!< Size of time step
  int,                     //!< Number of position components in position model
  int,                     //!< Number of density components in mass model
  int,                     //!< Number of velocity components in hydro model
  int,                     //!< Number of mixing scalars in material mix model
  int,                     //!< Number of frequencies in turb. frequency model
  uint64_t,                //!< Total number of particles
  int,                     //!< TTY output interval
  int,                     //!< Dump output interval
  int,                     //!< Plot output interval
  int,                     //!< PDF output interval
  int,                     //!< Glob output interval
  std::string,             //!< Input filename
  std::string,             //!< Output filename
  std::string,             //!< PDF filename
  std::string,             //!< Glob filename
  std::string,             //!< Statistics filename
  std::vector<real>,       //!< Parameters 'b' in Dirichlet mix models
  std::vector<real>,       //!< Parameters 'S' in Dirichlet mix models
  std::vector<real>,       //!< Parameters 'kappa' in Dirichlet mix models
  std::vector<real>,       //!< Parameters 'c_ij' in GenDirichlet mix models
  std::vector<select::RNGType>,  //!< Random number generators
  real,                    //!< Parameter C0 in the simplified Langevin model
  real,                    //!< Atwood number in beta model
  real,                    //!< C1 in gamma frequency model
  real,                    //!< C2 in gamma frequency model
  real,                    //!< C3 in gamma frequency model
  real,                    //!< C4 in gamma frequency model
  std::vector<real>,       //!< Sextet of box coordinates for anal. geometry def
  std::vector<Product>     //!< Requested (and triggered) statistics
>;

//! Vector of bools indicating whether data is set in Bundle during parsing
using BoolBundle = std::vector<bool>;

} // namespace control

} // namespace Quinoa

#endif // ControlTypes_h
