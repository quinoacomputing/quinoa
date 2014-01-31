//******************************************************************************
/*!
  \file      src/Control/Quinoa/Types.h
  \author    J. Bakosi
  \date      Fri Jan 31 10:44:09 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Types for Quinoa's parsers
  \details   Types for Quinoa's parsers
*/
//******************************************************************************
#ifndef QuinoaTypes_h
#define QuinoaTypes_h

#include <Types.h>
#include <Quinoa/Tags.h>
#include <Quinoa/Options/Geometry.h>
#include <Quinoa/Options/MonteCarlo.h>
#include <Quinoa/Options/Position.h>
#include <Quinoa/Options/Mass.h>
#include <Quinoa/Options/Hydro.h>
#include <Quinoa/Options/Energy.h>
#include <Quinoa/Options/Mix.h>
#include <Quinoa/Options/Frequency.h>
#include <Quinoa/Options/MixRate.h>
#include <Quinoa/Options/SDE.h>
#include <Options/RNG.h>

namespace quinoa {
//! control and parsing
namespace ctr {

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
//! Internally the numbering of field IDs starts from 0, but presented to the
//! user as starting from 1. Examples: 1st pressure fluctuation: {0, PRESSURE,
//! CENTRAL, ...}, mean of 2nd scalar: {1, SCALAR, ORDINARY, ...}.
struct Term {
  int field;         //!< Field ID
  Quantity quantity; //!< Physical quantity
  Moment moment;     //!< Moment type: ordinary, central
  int name;          //!< Character code as name, converted only for output
  bool plot;         //!< Indicates whether the variable will be plotted
  // Conceptually, plot should be in Product, since plot will only be false for
  // a mean that was triggered by a central moment by one of the Terms of a
  // Product requesting the mean or a model. However, that would require
  // Product to be a vector<struct>, which then would need custom comparitors
  // for std::sort() and std::unique() in Parser::unique(). Since this is not a
  // performance issue, plot is here in Term.

  //! Constructor
  explicit Term( int f, Quantity q, Moment m, char n, bool p ) :
    field(f), quantity(q), moment(m), name(n), plot(p) {}

  //! Equal operator for finding unique elements, used by e.g., std::unique()
  //! Test only on field, quantity, and moment
  bool operator== ( const Term& term ) const {
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
  bool operator< ( const Term& term ) const {
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
  friend std::string operator+ ( const std::string& lhs, const Term& term ) {
    std::stringstream ss;
    ss << lhs << char(term.name) << term.field+1;
    std::string rhs = ss.str();
    return rhs;
  }

  //! Operator << for writing Term to output streams
  friend std::ostream& operator<< ( std::ostream& os, const Term& term ) {
    os << char(term.name) << term.field+1;
    return os;
  }

  //! Operator << for writing vector<Term> to output streams
  friend std::ostream&
  operator<< ( std::ostream& os, const std::vector< Term >& vec ) {
    os << "<";
    for (auto& w : vec) os << w;
    os << ">";
    return os;
  }

  //! Operator <<= for writing requested vector<Term> to output streams
  friend std::ostream&
  operator<<= ( std::ostream& os, const std::vector< Term >& vec ) {
    if (vec[0].plot) {
      os << "<";
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
  explicit FieldName( const int n = 0, const int f = 0 ) :
    name(n), field(f) {}

  //! Operator << for writing FieldName to output streams
  friend std::ostream& operator<< ( std::ostream& os, const FieldName& fn ) {
     os << char(fn.name) << fn.field+1;
     return os;
  }

  //! Operator += for adding FieldName to std::string
  friend std::string& operator+= ( std::string& os, const FieldName& fn ) {
     std::stringstream ss;
     ss << os << char(fn.name) << fn.field+1;
     os = ss.str();
     return os;
  }
};

//! Products are N Terms to be multiplied and ensemble averaged
//! E.g. the scalar flux in x direction needs two terms for ensemble averaging:
//! (Y-\<Y\>) and (U-\<U\>), then the moment is \<yu\> = <(Y-\<Y\>)(U-\<U\>)>
//! E.g the third mixed central moment of three scalars needs three terms for
//! ensemble averaging: (Y1-\<Y1\>), (Y2-\<Y2\>), and (Y3-\<Y3\>), then the
//! moment is \<y1y2y3\> = \<(Y1-\<Y1\>)(Y2-\<Y2\>)(Y3-\<Y3\>)\>
using Product = std::vector<Term>;

//! Storage of selected options
using selects = tk::tuple::tagged_tuple<
  tag::geometry,   ctr::GeometryType,   //!< Selected geometry definition
  tag::montecarlo, ctr::MonteCarloType, //!< Selected physics
  tag::position,   ctr::PositionType,   //!< Selected position model
  tag::mass,       ctr::MassType,       //!< Selected mass model
  tag::hydro,      ctr::HydroType,      //!< Selected hydrodynamics model
  tag::energy,     ctr::EnergyType,     //!< Selected internal energy model
  tag::mix,        ctr::MixType,        //!< Selected material mix model
  tag::frequency,  ctr::FrequencyType,  //!< Selected turbulence frequency model
  tag::mixrate,    ctr::MixRateType,    //!< Selected material mix rate model
  tk::tag::rng,    std::vector< tk::ctr::RNGType >, //!< Selected RNGs
  tag::sde,        std::vector< ctr::SDEType >      //!< Selected SDEs
>;

//! Time incrementation parameters storage
using incpars = tk::tuple::tagged_tuple<
  tag::nstep, uint64_t,  //!< Number of time steps to take
  tag::term,  tk::real,  //!< Time to terminate time stepping
  tag::dt,    tk::real   //!< Size of time step
>;

//! Components storage
using components = tk::tuple::tagged_tuple<
  tag::nposition,  uint32_t, //!< Number of position components in position model
  tag::ndensity,   uint32_t, //!< Number of density components in mass model
  tag::nvelocity,  uint32_t, //!< Number of velocity components in hydro model
  tag::nscalar,    uint32_t, //!< Number of mixing scalars in material mix model
  tag::nfrequency, uint32_t, //!< Number of frequencies in turb. frequency model
  tag::npar,       uint64_t, //!< Total number of particles
  tag::ndirichlet, uint32_t, //!< Number of components in the Dirichlet SDE
  tag::ngendir,    uint32_t  //!< Number of components in the gen. Dirichlet SDE
>;

//! Output intervals storage
using intervals = tk::tuple::tagged_tuple<
  tag::tty,  uint32_t,  //!< TTY output interval
  tag::dump, uint32_t,  //!< Dump output interval
  tag::plot, uint32_t,  //!< Plot output interval
  tag::pdf,  uint32_t,  //!< PDF output interval
  tag::glob, uint32_t   //!< Glob output interval
>;

//! IO parameters storage
using ios = tk::tuple::tagged_tuple<
  tag::control,     std::string,  //!< Control filename
  tag::input,       std::string,  //!< Input filename
  tag::output,      std::string,  //!< Output filename
  tag::pdf,         std::string,  //!< PDF filename
  tag::glob,        std::string,  //!< Glob filename
  tag::stat,        std::string   //!< Statistics filename
>;

//! Beta mass model parameters storage
using BetaParameters = tk::tuple::tagged_tuple<
  tag::atwood,     tk::real,
  tk::tag::rng,    tk::ctr::RNGType            //!< RNG used
>;

//! Dirichlet mix model parameters storage
using DirichletParameters = tk::tuple::tagged_tuple<
  tag::b,         std::vector<tk::real>,
  tag::S,         std::vector<tk::real>,
  tag::kappa,     std::vector<tk::real>,
  tk::tag::rng,   tk::ctr::RNGType            //!< RNG used
>;

//! Generalized Dirichlet mix model parameters storage
using GenDirichletParameters = tk::tuple::tagged_tuple<
  tag::b,         std::vector<tk::real>,
  tag::S,         std::vector<tk::real>,
  tag::kappa,     std::vector<tk::real>,
  tag::c,         std::vector<tk::real>,
  tk::tag::rng,   tk::ctr::RNGType            //!< RNG used
>;

//! Gamma mix model parameters storage
using GammaParameters = tk::tuple::tagged_tuple<
  tag::c1,        tk::real,
  tag::c2,        tk::real,
  tag::c3,        tk::real,
  tag::c4,        tk::real,
  tk::tag::rng,   tk::ctr::RNGType            //!< RNG used
>;

//! Simplified Langevin hydro model parameters storage
using SLMParameters = tk::tuple::tagged_tuple<
  tag::c0, tk::real
>;

//! Generalized Langevin hydro model parameters storage
using GLMParameters = tk::tuple::tagged_tuple<
  tag::c0, tk::real
>;

//! Skew-normal parameters storage
using SkewNormalParameters = tk::tuple::tagged_tuple<
  tag::sigma,     tk::real,
  tag::timescale, tk::real,
  tag::lambda,    tk::real,
  tk::tag::rng,   tk::ctr::RNGType            //!< RNG used
>;

//! Ornstein-Uhlenbeck parameters storage
using OUParameters = tk::tuple::tagged_tuple<
  tag::sigma,     tk::real,
  tag::timescale, tk::real,
  tk::tag::rng,   tk::ctr::RNGType            //!< RNG used
>;

//! Parameters storage
using parameters = tk::tuple::tagged_tuple<
  tk::tag::mklrng,   tk::ctr::MKLRNGParameters,         // MKL RNG parameters
  tk::tag::rngsse,   tk::ctr::RNGSSEParameters,         // RNGSSE RNG parameters
  tag::beta,         BetaParameters,
  tag::dirichlet,    DirichletParameters,
  tag::gendir,       GenDirichletParameters,
  tag::gamma,        GammaParameters,
  tag::slm,          SLMParameters,
  tag::glm,          GLMParameters,
  tag::skewnormal,   SkewNormalParameters
>;

//! PEGTL location type to use throughout all of Quinoa's parsers
using Location = pegtl::ascii_location;

} // ctr::
} // quinoa::

#endif // QuinoaTypes_h
