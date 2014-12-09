//******************************************************************************
/*!
  \file      src/Control/ControlTypes.h
  \author    J. Bakosi
  \date      Tue 09 Dec 2014 09:36:58 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Types for tk control
  \details   Types for tk control
*/
//******************************************************************************
#ifndef ControlTypes_h
#define ControlTypes_h

#include <TaggedTuple.h>
#include <Tags.h>
#include <Types.h>
#include <Options/RNG.h>
#include <Options/RNGSSESeqLen.h>

#ifdef HAS_MKL
#include <Options/MKLUniformMethod.h>
#include <Options/MKLGaussianMethod.h>
#endif

namespace tk {
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

//! RNGSSE random number generator parameters storage
using RNGSSEParam = tk::tuple::tagged_tuple<
  tag::seed,          unsigned int,              //!< seed
  tag::seqlen,        RNGSSESeqLenType           //!< sequence length type
>;
//! RNGSSE parameters bundle
using RNGSSEParameters = std::map< RNGType, RNGSSEParam >;

#ifdef HAS_MKL
//! MKL random number generator parameters storage
using RNGMKLParam = tk::tuple::tagged_tuple<
  tag::seed,             unsigned int,              //!< seed
  tag::uniform_method,   MKLUniformMethodType,      //!< uniform method type
  tag::gaussian_method,  MKLGaussianMethodType      //!< Gaussian method type
>;
//! MKL RNG parameters bundle
using RNGMKLParameters = std::map< RNGType, RNGMKLParam >;
#endif

} // ctr::
} // tk::

#endif // ControlTypes_h
