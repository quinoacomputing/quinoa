//******************************************************************************
/*!
  \file      src/Control/Quinoa/Types.h
  \author    J. Bakosi
  \date      Fri 05 Dec 2014 02:34:47 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Types for Quinoa's parsers
  \details   Types for Quinoa's parsers
*/
//******************************************************************************
#ifndef QuinoaTypes_h
#define QuinoaTypes_h

#include <tkTags.h>
#include <tkTypes.h>
#include <Types.h>
#include <Quinoa/Tags.h>
#include <Quinoa/Options/MonteCarlo.h>
#include <Quinoa/Options/Position.h>
#include <Quinoa/Options/Mass.h>
#include <Quinoa/Options/Hydro.h>
#include <Quinoa/Options/Energy.h>
#include <Quinoa/Options/Mix.h>
#include <Quinoa/Options/Frequency.h>
#include <Quinoa/Options/MixRate.h>
#include <Quinoa/Options/DiffEq.h>
#include <Quinoa/Options/InitPolicy.h>
#include <Quinoa/Options/CoeffPolicy.h>
#include <Quinoa/Options/PDFFile.h>
#include <Quinoa/Options/PDFPolicy.h>
#include <Quinoa/Options/PDFCentering.h>
#include <Quinoa/Options/TxtFloatFormat.h>
#include <Options/RNG.h>
#include <PUPUtil.h>

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
  // for std::sort() and std::unique() in, e.g, Parser::unique(). Since this is
  // not a performance issue, plot is here, redundantly, in Term.

  //! Pack/Unpack
  void pup( PUP::er& p ) {
    p | field;
    PUP::pup( p, moment );
    p | var;
    p | plot;
  }
  friend void operator|( PUP::er& p, Term& t ) { t.pup(p); } 

  //! Empty constructor for Charm++
  explicit Term() : field(0), moment(Moment::ORDINARY), var(0), plot(false) {}

  //! Constructor
  explicit Term( int f, Moment m, char v, bool p ) :
    field(f), moment(m), var(v), plot(p) {}

  //! Equal operator for, e.g., finding unique elements, used by, e.g.,
  //! std::unique(). Test on field, moment, and var, ignore plot.
  bool operator== ( const Term& term ) const {
    if (field == term.field && moment == term.moment && var == term.var)
      return true;
    else
      return false;
  }

  //! Less than operator for ordering, used by e.g., std::sort().
  //! Test on field, term, moment, and !plot.
  //! Using operator >, instead of operator <, on plot ensures that if a Term is
  //! user-requested, i.e., plotted, and also triggered by e.g., a model, the
  //! user-requested Term will take precendence.
  bool operator< ( const Term& term ) const {
    // test on everything except var
    if (field < term.field)
      return true;
    else if (field == term.field && moment < term.moment)
      return true;
    else if (field == term.field && moment == term.moment && var < term.var)
      return true;
    else if (field == term.field && moment == term.moment &&
             var == term.var && plot > term.plot)
      return true;
    else
      return false;
  }
};

//! Pack/Unpack Term
inline void pup( PUP::er& p, Term& t ) { t.pup(p); }

//! Operator + for adding Term (var+field ID) to a std::string
static std::string operator+ ( const std::string& lhs, const Term& term ) {
  std::stringstream ss;
  ss << lhs << char(term.var) << term.field+1;
  std::string rhs = ss.str();
  return rhs;
}

//! Operator << for writing Term to output streams
static std::ostream& operator<< ( std::ostream& os, const Term& term ) {
  os << char(term.var) << term.field+1;
  return os;
}

//! Lighter-weight (lighter than Term) structure for var+field. Used for
//! representing the variable + field ID in e.g., statistics or sample space.
struct FieldVar {
  char var;
  int field;

  //! Constructor
  explicit FieldVar( const char v='\0', const int f=0 ) : var(v), field(f) {}

  //! Equal operator for, e.g., testing on equality of containers containing
  //! FieldVars in any way finding, e.g., InputDeck< tag::pdf >. Test on both
  //! var and field.
  bool operator== ( const FieldVar& f ) const {
    if (field == f.field && var == f.var)
      return true;
    else
      return false;
  }

  //! Operator += for adding FieldVar to std::string
  friend std::string& operator+= ( std::string& os, const FieldVar& f ) {
     std::stringstream ss;
     ss << os << f.var << f.field+1;
     os = ss.str();
     return os;
  }
};

//! Function for writing std::vector< Term > to output streams
static
std::ostream& estimated( std::ostream& os, const std::vector< Term >& vec ) {
  os << "<";
  for (const auto& w : vec) os << w;
  os << "> ";
  return os;
}

//! Function for writing requested statistics terms to output streams
static
std::ostream& requested( std::ostream& os, const std::vector< Term >& vec ) {
  if (!vec.empty() && vec[0].plot) {
    os << "<";
    for (const auto& w : vec) os << w;
    os << "> ";
  }
  return os;
}

//! Function for writing triggered statistics terms to output streams
static
std::ostream& triggered( std::ostream& os, const std::vector< Term >& vec ) {
  if (!vec.empty() && !vec[0].plot) {
    os << "<";
    for (const auto& w : vec) os << w;
    os << "> ";
  }
  return os;
}

//! Function for writing pdf sample space variables to output streams
static
std::ostream& pdf( std::ostream& os,
                   const std::vector< Term >& var,
                   const std::vector< tk::real >& bin,
                   const std::string& name,
                   const std::vector< tk::real >& ext )
{
  Assert( !var.empty(), "var is empty in sample_space()" );
  Assert( !bin.empty(), "bin is empty in sample_space()" );
  Assert( var.size() == bin.size(),
          "var.size and bin.size() must equal in ctr::pdf()" );

  os << name << '(';
  std::size_t i;
  // sample space variables
  for (i=0; i<var.size()-1; ++i) os << var[i] << ',';
  os << var[i] << ':';
  // sample space bin sizes
  for (i=0; i<bin.size()-1; ++i) os << bin[i] << ',';
  os << bin[i];
  // sample space extents
  if (!ext.empty()) {
    os << ';';
    for (i=0; i<ext.size()-1; ++i) os << ext[i] << ',';
    os << ext[i];
  }
  os << ") ";
  return os;
}

//! Case-insensitive character comparison functor
struct CaseInsensitiveCharLess {
  bool operator() ( char lhs, char rhs ) const {
    return std::tolower( lhs ) < std::tolower( rhs );
  }
};

//! Products are arbitrary number of Terms to be multiplied and ensemble
//! averaged, an example is the scalar flux in x direction which needs two terms
//! for ensemble averaging: (Y-\<Y\>) and (U-\<U\>), then the moment is \<yu\> =
//! <(Y-\<Y\>)(U-\<U\>)>, another example is the third mixed central moment of
//! three scalars which needs three terms for ensemble averaging: (Y1-\<Y1\>),
//! (Y2-\<Y2\>), and (Y3-\<Y3\>), then the moment is \<y1y2y3\> =
//! \<(Y1-\<Y1\>)(Y2-\<Y2\>)(Y3-\<Y3\>)\>
using Product = std::vector< Term >;

//! Find out if a vector of Terms only contains ordinary moment terms
//! \details Iff all terms are ordinary, the vector of Terms is ordinary.
static inline bool ordinary( const std::vector< ctr::Term >& vec ) {
  bool ord = true;
  for (auto& term : vec) if (term.moment == ctr::Moment::CENTRAL) ord = false;
  return ord;
}

//! Find out if a vector of Terms  contains any central moment terms
//! \details If any term is central, the vector of Terms is central.
static inline bool central( const std::vector< ctr::Term >& vec )
{ return !ordinary( vec ); }

//! Probability density function (sample space variables)
using Probability = std::vector< Term >;

//! Storage of selected options
using selects = tk::tuple::tagged_tuple<
  tag::montecarlo,   ctr::MonteCarloType,   //!< Physics
  tag::position,     ctr::PositionType,     //!< Position model
  tag::mass,         ctr::MassType,         //!< Mass model
  tag::hydro,        ctr::HydroType,        //!< Hydrodynamics model
  tag::energy,       ctr::EnergyType,       //!< Internal energy model
  tag::mix,          ctr::MixType,          //!< Material mix model
  tag::frequency,    ctr::FrequencyType,    //!< Turbulence frequency model
  tag::mixrate,      ctr::MixRateType,      //!< Material mix rate model
  tag::diffeq,       std::vector< ctr::DiffEqType >,  //!< Differential eqs
  tk::tag::rng,      std::vector< tk::ctr::RNGType >, //!< RNGs
  tag::pdffiletype,  ctr::PDFFileType,      //!< PDF output file type
  tag::pdfpolicy,    ctr::PDFPolicyType,    //!< PDF output file policy
  tag::pdfctr,       ctr::PDFCenteringType, //!< PDF output file centering
  tag::float_format, ctr::TxtFloatFormatType//!< Text floating-point format
>;

//! Discretization parameters storage
using discretization = tk::tuple::tagged_tuple<
  tag::npar,      uint64_t,  //!< Total number of particles
  tag::nstep,     uint64_t,  //!< Number of time steps to take
  tag::term,      tk::real,  //!< Time to terminate time stepping
  tag::dt,        tk::real,  //!< Size of time step
  tag::binsize,   std::vector< std::vector< tk::real > >, //!< PDF binsizes
  tag::extent,    std::vector< std::vector< tk::real > >, //!< PDF extents
  tag::precision, std::streamsize  //!< Precision in digits
>;

//! Output intervals storage
using intervals = tk::tuple::tagged_tuple<
  tag::tty,  uint32_t,  //!< TTY output interval
  tag::dump, uint32_t,  //!< Dump output interval
  tag::stat, uint32_t,  //!< Statistics output interval
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
  tag::stat,        std::string,  //!< Statistics filename
  tag::pdfnames,    std::vector< std::string >  //!< PDF identifiers
>;

//! Position parameters storage
using PositionParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >
>;

//! Mass parameters storage
using MassParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >
>;

//! Hydro parameters storage
using HydroParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >
>;

//! Mix parameters storage
using MixParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >
>;

//! Frequency parameters storage
using FrequencyParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >
>;

//! Dirichlet mix model parameters storage
using DirichletParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::b,           std::vector< std::vector< tk::real > >,
  tag::S,           std::vector< std::vector< tk::real > >,
  tag::kappa,       std::vector< std::vector< tk::real > >,
  tk::tag::rng,     std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< ctr::CoeffPolicyType >
>;

//! Generalized Dirichlet parameters storage
using GenDirichletParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::b,           std::vector< std::vector< tk::real > >,
  tag::S,           std::vector< std::vector< tk::real > >,
  tag::kappa,       std::vector< std::vector< tk::real > >,
  tag::c,           std::vector< std::vector< tk::real > >,
  tk::tag::rng,     std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< ctr::CoeffPolicyType >
>;

//! Wright-Fisher parameters storage
using WrightFisherParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::omega,       std::vector< std::vector< tk::real > >,
  tk::tag::rng,     std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< ctr::CoeffPolicyType >
>;

//! Ornstein-Uhlenbeck parameters storage
using OrnsteinUhlenbeckParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::sigma,       std::vector< std::vector< tk::real > >,
  tag::theta,       std::vector< std::vector< tk::real > >,
  tag::mu,          std::vector< std::vector< tk::real > >,
  tk::tag::rng,     std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< ctr::CoeffPolicyType >
>;

//! Diagonal Ornstein-Uhlenbeck parameters storage
using DiagOrnsteinUhlenbeckParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::sigma,       std::vector< std::vector< tk::real > >,
  tag::theta,       std::vector< std::vector< tk::real > >,
  tag::mu,          std::vector< std::vector< tk::real > >,
  tk::tag::rng,     std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< ctr::CoeffPolicyType >
>;

//! Log-normal parameters storage
using LogNormalParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::sigma,       tk::real,
  tag::timescale,   tk::real,
  tk::tag::rng,     std::vector< tk::ctr::RNGType >
>;

//! Skew-normal parameters storage
using SkewNormalParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::timescale,   std::vector< std::vector< tk::real > >,
  tag::sigma,       std::vector< std::vector< tk::real > >,
  tag::lambda,      std::vector< std::vector< tk::real > >,
  tk::tag::rng,     std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< ctr::CoeffPolicyType >
>;

//! Gamma parameters storage
using GammaParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::b,           std::vector< std::vector< tk::real > >,
  tag::S,           std::vector< std::vector< tk::real > >,
  tag::kappa,       std::vector< std::vector< tk::real > >,
  tk::tag::rng,     std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< ctr::CoeffPolicyType >
>;

//! Beta parameters storage
using BetaParameters = tk::tuple::tagged_tuple<
  tag::depvar,      std::vector< char >,
  tag::b,           std::vector< std::vector< tk::real > >,
  tag::S,           std::vector< std::vector< tk::real > >,
  tag::kappa,       std::vector< std::vector< tk::real > >,
  tk::tag::rng,     std::vector< tk::ctr::RNGType >,
  tag::initpolicy,  std::vector< ctr::InitPolicyType >,
  tag::coeffpolicy, std::vector< ctr::CoeffPolicyType >
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
  tag::wrightfisher, WrightFisherParameters,
  tag::ou,           OrnsteinUhlenbeckParameters,
  tag::diagou,       DiagOrnsteinUhlenbeckParameters,
  tag::lognormal,    LogNormalParameters,
  tag::skewnormal,   SkewNormalParameters,
  tag::gamma,        GammaParameters,
  tag::beta,         BetaParameters
>;

//! PEGTL location type to use throughout Quinoa's parsers
using Location = pegtl::ascii_location;

} // ctr::
} // quinoa::

#endif // QuinoaTypes_h
