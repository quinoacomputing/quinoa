//******************************************************************************
/*!
  \file      src/Control/Quinoa/Types.h
  \author    J. Bakosi
  \date      Sun 01 Jun 2014 11:46:18 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
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

//! Moment specifies which type of moment is computed for a quantity in a Term
enum class Moment : uint8_t { ORDINARY=0,      //!< Full variable
                              CENTRAL          //!< Fluctuation
};

//! Term is a Moment of a quantity with a field ID to be ensemble averaged.
//! Internally the numbering of field IDs starts from 0, but presented to the
//! user as starting from 1.
struct Term {
  int field;         //!< Field ID
  Moment moment;     //!< Moment type: ordinary, central
  char var;          //!< Dependent variable
  bool plot;         //!< Indicates whether the variable will be plotted
  // Conceptually, plot should be in Product, since plot will only be false for
  // a mean that was triggered by a central moment by one of the Terms of a
  // Product requesting the mean or a model. However, that would require
  // Product to be a vector<struct>, which then would need custom comparitors
  // for std::sort() and std::unique() in Parser::unique(). Since this is not a
  // performance issue, plot is here in Term.

  //! Constructor
  explicit Term( int f, Moment m, char v, bool p ) :
    field(f), moment(m), var(v), plot(p) {}

  //! Equal operator for finding unique elements, used by e.g., std::unique()
  //! Test only on field and moment
  bool operator== ( const Term& term ) const {
    if (field == term.field && moment == term.moment && var == term.var) {
      return true;
    } else {
      return false;
    }
  }

  //! Less than operator for ordering, used by e.g., std::sort().
  //! Test on field, term, moment, and !plot.
  //! Using operator >, instead of operator <, on plot ensures that if a Term is
  //! user-requested, i.e., plotted, and also triggered by e.g., a model, the
  //! user-requested Term will take precendence.
  bool operator< ( const Term& term ) const {
    // test on everything except var
    if (field < term.field) {
      return true;
    } else if (field == term.field && moment < term.moment) {
      return true;
    } else if (field == term.field && moment == term.moment && var < term.var) {
      return true;
    } else if (field == term.field && moment == term.moment &&
               var == term.var && plot > term.plot) {
      return true;
    } else {
      return false;
    }
  }

  //! Operator + for adding Term (var+field ID) to a std::string
  friend std::string operator+ ( const std::string& lhs, const Term& term ) {
    std::stringstream ss;
    ss << lhs << char(term.var) << term.field+1;
    std::string rhs = ss.str();
    return rhs;
  }

  //! Operator << for writing Term to output streams
  friend std::ostream& operator<< ( std::ostream& os, const Term& term ) {
    os << char(term.var) << term.field+1;
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

//! Lighter-weight (lighter than Term) structure for var+field
struct FieldVar {
  char var;
  int field;

  //! Constructor
  explicit FieldVar( const char n = '\0', const int f = 0 ) :
    var(n), field(f) {}

  //! Operator << for writing FieldVar to output streams
  friend std::ostream& operator<< ( std::ostream& os, const FieldVar& fn ) {
     os << char(fn.var) << fn.field+1;
     return os;
  }

  //! Operator += for adding FieldVar to std::string
  friend std::string& operator+= ( std::string& os, const FieldVar& fn ) {
     std::stringstream ss;
     ss << os << fn.var << fn.field+1;
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
using Product = std::vector< Term >;

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
  tag::sde,        std::vector< ctr::SDEType >,    //!< Selected SDEs
  tk::tag::rng,    std::vector< tk::ctr::RNGType > //!< Selected RNGs
>;

//! Time incrementation parameters storage
using incpars = tk::tuple::tagged_tuple<
  tag::npar,  uint64_t,  //!< Total number of particles
  tag::nstep, uint64_t,  //!< Number of time steps to take
  tag::term,  tk::real,  //!< Time to terminate time stepping
  tag::dt,    tk::real   //!< Size of time step
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

//! Position parameters storage
using PositionParameters = tk::tuple::tagged_tuple<
  tag::depvar,      char
>;

//! Mass parameters storage
using MassParameters = tk::tuple::tagged_tuple<
  tag::depvar,      char
>;

//! Hydro parameters storage
using HydroParameters = tk::tuple::tagged_tuple<
  tag::depvar,      char
>;

//! Mix parameters storage
using MixParameters = tk::tuple::tagged_tuple<
  tag::depvar,      char
>;

//! Frequency parameters storage
using FrequencyParameters = tk::tuple::tagged_tuple<
  tag::depvar,      char
>;

//! Dirichlet mix model parameters storage
using DirichletParameters = tk::tuple::tagged_tuple<
  tag::depvar,      char,
  tag::b,           std::vector< tk::real >,
  tag::S,           std::vector< tk::real >,
  tag::kappa,       std::vector< tk::real >,
  tk::tag::rng,     tk::ctr::RNGType,
  tag::initpolicy,  ctr::InitPolicyType,
  tag::coeffpolicy, ctr::CoeffPolicyType
>;

//! Generalized Dirichlet parameters storage
using GenDirichletParameters = tk::tuple::tagged_tuple<
  tag::depvar,      char,
  tag::b,           std::vector< tk::real >,
  tag::S,           std::vector< tk::real >,
  tag::kappa,       std::vector< tk::real >,
  tag::c,           std::vector< tk::real >,
  tk::tag::rng,     tk::ctr::RNGType,
  tag::initpolicy,  ctr::InitPolicyType,
  tag::coeffpolicy, ctr::CoeffPolicyType
>;

//! Ornstein-Uhlenbeck parameters storage
using OrnsteinUhlenbeckParameters = tk::tuple::tagged_tuple<
  tag::depvar,      char,
  tag::sigma,       tk::real,
  tag::timescale,   tk::real,
  tk::tag::rng,     tk::ctr::RNGType
>;

//! Log-normal parameters storage
using LogNormalParameters = tk::tuple::tagged_tuple<
  tag::depvar,      char,
  tag::sigma,       tk::real,
  tag::timescale,   tk::real,
  tk::tag::rng,     tk::ctr::RNGType
>;

//! Skew-normal parameters storage
using SkewNormalParameters = tk::tuple::tagged_tuple<
  tag::depvar,      char,
  tag::sigma,       tk::real,
  tag::timescale,   tk::real,
  tag::lambda,      tk::real,
  tk::tag::rng,     tk::ctr::RNGType
>;

//! Gamma parameters storage
using GammaParameters = tk::tuple::tagged_tuple<
  tag::depvar,      char,
  tag::c1,          tk::real,
  tag::c2,          tk::real,
  tag::c3,          tk::real,
  tag::c4,          tk::real,
  tk::tag::rng,     tk::ctr::RNGType
>;

//! Beta parameters storage
using BetaParameters = tk::tuple::tagged_tuple<
  tag::depvar,      char,
  tag::atwood,      tk::real,
  tag::b,           tk::real,
  tag::S,           tk::real,
  tag::kappa,       tk::real,
  tk::tag::rng,     tk::ctr::RNGType
>;

//! Simplified Langevin hydro model parameters storage
using SLMParameters = tk::tuple::tagged_tuple<
  tag::c0, tk::real
>;

//! Generalized Langevin hydro model parameters storage
using GLMParameters = tk::tuple::tagged_tuple<
  tag::c0, tk::real
>;

//! Parameters storage
using parameters = tk::tuple::tagged_tuple<
  #ifdef HAS_MKL
  tk::tag::rngmkl,   tk::ctr::RNGMKLParameters,   //!< MKL RNG parameters
  #endif
  tk::tag::rngsse,   tk::ctr::RNGSSEParameters,   //!< RNGSSE RNG parameters
  tag::position,     PositionParameters,
  tag::mass,         MassParameters,
  tag::hydro,        HydroParameters,
  tag::mix,          MixParameters,
  tag::frequency,    FrequencyParameters,
  tag::slm,          SLMParameters,
  tag::glm,          GLMParameters,
  tag::dirichlet,    DirichletParameters,
  tag::gendir,       GenDirichletParameters,
  tag::ou,           OrnsteinUhlenbeckParameters,
  tag::lognormal,    LogNormalParameters,
  tag::skewnormal,   SkewNormalParameters,
  tag::gamma,        GammaParameters,
  tag::beta,         BetaParameters
>;

//! PEGTL location type to use throughout all of Quinoa's parsers
using Location = pegtl::ascii_location;

} // ctr::
} // quinoa::

#endif // QuinoaTypes_h
