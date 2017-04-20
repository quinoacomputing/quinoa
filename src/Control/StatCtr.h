// *****************************************************************************
/*!
  \file      src/Control/StatCtr.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Types and associated functions to deal with moments and PDFs
  \details   Types and associated functions to deal with statistical moments and
    probability density functions.
*/
// *****************************************************************************
#ifndef StatControl_h
#define StatControl_h

#include "Types.h"
#include "Exception.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace tk {
namespace ctr {

//! \brief Moment specifies which type of moment is computed for a quantity in
//!    a Term
//! \author J. Bakosi
enum class Moment : uint8_t { ORDINARY=0,      //!< Full variable
                              CENTRAL          //!< Fluctuation
};

//! \brief Term is a Moment of a quantity with a field ID to be ensemble
//!    averaged
//! \details Internally the numbering of field IDs starts from 0, but presented
//!    to the user, e.g., in screen-output, as starting from 1.
//! \author J. Bakosi
struct Term {
  using ncomp_t = kw::ncomp::info::expect::type;

  char var;             //!< Variable name
  ncomp_t field;        //!< Field ID
  Moment moment;        //!< Moment type: ordinary, central

  /** @name Pack/Unpack: Serialize Term object for Charm++ */
  ///@{
  //! Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \author J. Bakosi
  void pup( PUP::er& p ) {
    p | var;
    p | field;
    PUP::pup( p, moment );
  }
  //! \brief Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] t Term object reference
  //! \author J. Bakosi
  friend void operator|( PUP::er& p, Term& t ) { t.pup(p); } 
  ///@}

  //! \brief Constructor: initialize all state data
  //! \param[in] v Variable name
  //! \param[in] f Field ID
  //! \param[in] m Moment type enum: Moment::ORDINARY or Moment::CENTRAL
  //! \author J. Bakosi
  explicit Term( char v = 0, ncomp_t f = 0, Moment m = Moment::ORDINARY ) :
    var( v ), field( f ), moment( m ) {}

  //! \brief Equal operator for, e.g., finding unique elements, used by, e.g.,
  //!    std::unique().
  //! \details Test on field, moment, and var
  //! \param[in] term Term to compare
  //! \return Boolean indicating if term equals 'this'
  //! \author J. Bakosi
  bool operator== ( const Term& term ) const {
    if (var == term.var && field == term.field && moment == term.moment)
      return true;
    else
      return false;
  }

  //! \brief Less-than operator for ordering, used by, e.g., std::sort().
  //! \details Test on var, field, and moment.
  //! \param[in] term Term to compare
  //! \return Boolean indicating if term is less than 'this'
  //! \author J. Bakosi
  bool operator< ( const Term& term ) const {
    if (var < term.var)
      return true;
    else if (var == term.var && field < term.field)
      return true;
    else if (var == term.var && field == term.field && moment < term.moment)
      return true;
    else
      return false;
  }
};

//! \brief Pack/Unpack: Namespace-scope serialize Term object for Charm++
//! \param[in,out] p Charm++'s PUP::er serializer object reference
//! \param[in,out] t Term object reference
//! \author J. Bakosi
inline void pup( PUP::er& p, Term& t ) { t.pup(p); }

//! \brief Products are arbitrary number of Terms to be multiplied and ensemble
//!   averaged.
//! \details An example is the scalar flux in x direction which needs two terms
//! for ensemble averaging: (Y-\<Y\>) and (U-\<U\>), then the central moment is
//! \<yu\> = <(Y-\<Y\>)(U-\<U\>)>, another example is the third mixed central
//! moment of three scalars which needs three terms for ensemble averaging:
//! (Y1-\<Y1\>), (Y2-\<Y2\>), and (Y3-\<Y3\>), then the central moment is
//! \<y1y2y3\> = \<(Y1-\<Y1\>)(Y2-\<Y2\>)(Y3-\<Y3\>)\>.
//! \author J. Bakosi
using Product = std::vector< Term >;


// The following functions are useful for debugging, and unused.
#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-function"
#endif

//! \brief Operator + for adding Term (var+field) to a std::string
//! \param[in] lhs std::string to add to
//! \param[in] term Term to add
//! \return Updated std::string
//! \author J. Bakosi
static std::string operator+ ( const std::string& lhs, const Term& term ) {
  std::stringstream ss;
  ss << lhs << term.var << term.field+1;
  std::string rhs = ss.str();
  return rhs;
}

//! \brief Operator += for adding Term (var+field) to a std::string
//! \param[in,out] os std::string to add to
//! \param[in] term Term to add
//! \return Updated std::string
//! \author J. Bakosi
static std::string& operator+= ( std::string& os, const Term& term ) {
  std::stringstream ss;
  ss << os << term.var << term.field+1;
  os = ss.str();
  return os;
}

//! \brief Operator << for writing Term to output streams
//! \param[in,out] os Output stream to write to
//! \param[in] term Term to write
//! \return Updated output stream
//! \author J. Bakosi
static std::ostream& operator<< ( std::ostream& os, const Term& term ) {
  os << term.var << term.field+1;
  return os;
}

//! \brief Operator + for adding products (var+field) to a std::string
//! \param[in] lhs std::string to add to
//! \param[in] p Product to add
//! \return Updated std::string
//! \author J. Bakosi
static std::string operator+ ( const std::string& lhs, const Product& p ) {
  std::stringstream ss;
  ss << lhs;
  if (!p.empty()) {
    ss << "<";
    for (const auto& w : p) ss << w;
    ss << ">";
  }
  std::string rhs = ss.str();
  return rhs;
}

//! \brief Operator << for writing products to output streams
//! \param[in,out] os Output stream to write to
//! \param[in] p Product, std::vector< Term >, to write
//! \return Updated output stream
//! \author J. Bakosi
static
std::ostream& operator<< ( std::ostream& os, const Product& p ) {
  if (!p.empty()) {
    os << "<";
    for (const auto& w : p) os << w;
    os << "> ";
  }
  return os;
}

//! \brief Function for writing PDF sample space variables to output streams
//! \param[in,out] os Output stream to write to
//! \param[in] var Vector of Terms to write
//! \param[in] bin Vector of PDF bin sizes
//! \param[in] name Name of PDF
//! \param[in] ext Vector of sample space extents
//! \return Updated output stream
//! \author J. Bakosi
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

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

//! \brief Case-insensitive character comparison functor
//! \author J. Bakosi
struct CaseInsensitiveCharLess {
  //! Function call operator
  //! \param[in] lhs Left character of the comparitor operand
  //! \param[in] rhs Right character of the comparitor operand
  //! \return Boolean indicating the result of the comparison
  //! \author J. Bakosi
  bool operator() ( char lhs, char rhs ) const {
    return std::tolower( lhs ) < std::tolower( rhs );
  }
};

//! \brief Find out if a vector of Terms only contains ordinary moment terms
//! \details If and only if all terms are ordinary, the vector of Terms is
//!    ordinary.
//! \param[in] vec Vector of Terms to check
//! \return Boolean indicating if all terms are ordinary
//! \author J. Bakosi
static inline bool
ordinary( const std::vector< ctr::Term >& vec ) {
  bool ord = true;
  for (auto& term : vec) if (term.moment == ctr::Moment::CENTRAL) ord = false;
  return ord;
}

//! \brief Find out if a vector of Terms contains any central moment terms
//! \details If any term is central, the vector of Terms is central.
//! \param[in] vec Vector of Terms to check
//! \return Boolean indicating of any term is central
//! \author J. Bakosi
static inline bool
central( const std::vector< ctr::Term >& vec )
{ return !ordinary( vec ); }

//! \brief Probability density function (vector of sample space variables)
//! \author J. Bakosi
using Probability = std::vector< Term >;

//! \brief PDF information bundle
//! \note If the user did not specify extents for a PDF, the corresponding
//!   extents vector still exists but it is empty.
//! \author J. Bakosi
struct PDFInfo {
  const std::string& name;                  //!< PDF identifier, i.e., name
  const std::vector< tk::real >& exts;      //!< Sample space extents
  std::vector< std::string > vars;          //!< List of sample space ariables
};

//! \brief Find PDF information, see tk::ctr::PDFInfo
//! \note Size of binsizes, names, pdfs, and exts must all be equal
//! \note idx must be less than the length of binsizes, names, and pdfs
//! \param[in] binsizes Vector of sample space bin sizes for multiple PDFs
//! \param[in] names Vector of PDF names
//! \param[in] exts Vector of sample space extents. Note: if the user did not
//!   specify extents for a PDF, the corresponding extents vector still exists
//!   but it is empty.
//! \param[in] pdfs Vector of PDFs
//! \param[in] m ORDINARY or CENTRAL PDF we are looking for
//! \param[in] idx Index of the PDF within the list of matching (based on D and
//!   m) PDFs requested
//! \return The PDF metadata requested
//! \details Find PDF information given the sample space dimension (template
//!   argument D), ordinary or central PDF (m), and the index within the list of
//!   matching (based on D and m) PDFs requested (idx). This function must find
//!   the PDF, if it does not, it throws an exception.
//! \see walker::Distributor
//! \author J. Bakosi
template< std::size_t D >
PDFInfo pdfInfo( const std::vector< std::vector< tk::real > >& binsizes,
                 const std::vector< std::string >& names,
                 const std::vector< std::vector< tk::real > >& exts,
                 const std::vector< Probability >& pdfs,
                 tk::ctr::Moment m,
                 std::size_t idx )
{
  Assert( binsizes.size() == names.size(),
          "Length of binsizes vector and that of PDF names must equal" );
  Assert( binsizes.size() == pdfs.size(),
          "Length of binsizes vector and that of PDFs must equal" );
  Assert( binsizes.size() == exts.size(),
          "Length of binsizes vector and that of PDF extents must equal" );
  Assert( binsizes.size() > idx, "Indexing out of bounds" );

  std::size_t i = 0;  // will count all PDFs queried
  std::size_t n = 0;  // will count PDFs with sample space dimensions and type
                      // (orindary or central) we are looking for
  for (const auto& bs : binsizes) {
    if ( bs.size() == D &&
         ((m == Moment::ORDINARY && ordinary(pdfs[i])) ||
          (m == Moment::CENTRAL && central(pdfs[i]))) ) ++n;
    if (n == idx+1) {
      std::vector< std::string > vars;
      for (const auto& term : pdfs[i])
        vars.push_back( term.var + std::to_string(term.field+1) );
      return { names[i], exts[i], std::move(vars) };
    }
    ++i;
  }
  Throw( "Cannot find PDF." );
}

//! Extract number of PDFs given sample space dimension
//! \details Count number of PDFs given the sample space dimension (template
//!   argument D) and whether the PDF is ordinary or central (m)
//! \note Size of binsizes, names, pdfs, and exts must all be equal
//! \param[in] pdfs Vector of PDFs
//! \param[in] m ORDINARY or CENTRAL PDF we are looking for
//! \return The number of PDFs matchin the criteria discussed above
//! \author J. Bakosi
template< std::size_t D >
std::size_t numPDF( const std::vector< std::vector< tk::real > >& binsizes,
                    const std::vector< Probability >& pdfs,
                    ctr::Moment m )
{
  Assert( binsizes.size() == pdfs.size(),
          "Length of binsizes vector and that of PDFs must equal" );
  auto& kind = (m == Moment::ORDINARY ? ordinary : central);
  std::size_t i=0, n=0;
  for (const auto& p : pdfs) {
    const auto& bs = binsizes[i++];
    if (kind(p) && bs.size() == D) ++n;
  }
  return n;
}

//! Lookup moment in moments map based on product key
static inline tk::real
lookup( const Product& p, const std::map< Product, tk::real >& moments ) {
  const auto& it = moments.find( p );
  if (it != end(moments))
    return it->second;
  else
    Throw( "Cannot find moment " + p + " in moments map" );
}

//! Construct mean
//! \param[in] var Variable
//! \param[in] c Component number
//! \return Constructed vector< Term > identifying the first ordinary moment
//!   (mean) of field (component) c of variable var
//! \author J. Bakosi
static inline Product
mean( char var, kw::ncomp::info::expect::type c ) {
  tk::ctr::Term m( static_cast<char>(std::toupper(var)), c, Moment::ORDINARY );
  return tk::ctr::Product( { m } );
}
//! Construct variance
//! \param[in] var Variable
//! \param[in] c Component number
//! \return Constructed vector< Term > identifying the second central moment
//!   (variance) of field (component) c of variable var
//! \author J. Bakosi
static inline Product
variance( char var, kw::ncomp::info::expect::type c ) {
  tk::ctr::Term f( static_cast<char>(std::tolower(var)), c, Moment::CENTRAL );
  return tk::ctr::Product( { f, f } );
}
//! Construct third central moment
//! \param[in] var Variable
//! \param[in] c Component number
//! \return Constructed vector< Term > identifying the third central moment
//!   of field (component) c of variable var
//! \author J. Bakosi
static inline Product
cen3( char var, kw::ncomp::info::expect::type c ) {
  tk::ctr::Term f( static_cast<char>(std::tolower(var)), c, Moment::CENTRAL );
  return tk::ctr::Product( { f, f, f } );
}

//! Construct second ordinary moment
//! \param[in] var Variable
//! \param[in] c Component number
//! \return Constructed vector< Term > identifying the second ordinary moment
//!   of field (component) c of variable var
//! \author J. Bakosi
static inline Product
ord2( char var, kw::ncomp::info::expect::type c ) {
  tk::ctr::Term f( static_cast<char>(std::toupper(var)), c, Moment::ORDINARY );
  return tk::ctr::Product( { f, f } );
}

} // ctr::
} // tk::

#endif // StatControl_h
