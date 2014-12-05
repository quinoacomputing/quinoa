//******************************************************************************
/*!
  \file      src/Main/QuinoaPrint.h
  \author    J. Bakosi
  \date      Fri 05 Dec 2014 11:36:29 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Quinoa's printer
  \details   Quinoa's printer
*/
//******************************************************************************
#ifndef QuinoaPrint_h
#define QuinoaPrint_h

#include <algorithm>

#include <boost/algorithm/string/replace.hpp>

#include <RNGPrint.h>
#include <DiffEq.h>
#include <Quinoa/Types.h>
#include <Quinoa/Options/DiffEq.h>
#include <Quinoa/InputDeck/InputDeck.h>

namespace quinoa {

extern ctr::InputDeck g_inputdeck_defaults;
extern ctr::InputDeck g_inputdeck;

//! QuinoaPrint : RNGPrint
class QuinoaPrint : public tk::RNGPrint {

  public:
    //! Constructor
    explicit QuinoaPrint( std::ostream& str = std::clog,
                          std::ostream& qstr = std::cout ) :
      RNGPrint( str, qstr ) {}

    //! Print control option: 'group : option' only if differs from its default
    template< typename Option, typename... tags >
    void Section() const {
      if (g_inputdeck.get< tags... >() !=
            g_inputdeck_defaults.get< tags... >() ) {
        Option opt;
        auto& group = opt.group();
        auto& value = opt.name( g_inputdeck.get< tags... >() );
        m_stream << m_section_title_value_fmt % m_section_indent
                                              % m_section_bullet
                                              % group
                                              % value;
        m_stream << m_section_underline_fmt
                    % m_section_indent
                    % std::string( m_section_indent.size() + 3 +
                                   group.size() + value.size(), '-' );
      }
    }

    //! Print item: 'name : value' only if differs from its default
    template<typename... tags>
    void Item(const std::string& name) const {
      if (g_inputdeck.get< tags... >() !=
            g_inputdeck_defaults.get< tags... >() )
        m_stream << m_item_name_value_fmt
                    % m_item_indent % name % g_inputdeck.get<tags...>();
    }

    //! Print control option: 'group : option' only if differs from its default
    template<typename Option, typename... tags>
    void Item() const {
      if (g_inputdeck.get< tags... >() !=
            g_inputdeck_defaults.get< tags... >() ) {
        Option opt;
        m_stream << m_item_name_value_fmt
                    % m_item_indent % opt.group()
                    % opt.name( g_inputdeck.get< tags... >() );
      }
    }

    // Helper class for compact output of diff eq policies
    class Policies {
      public:
        // Default constructor
        explicit Policies() {}
        // Initializer constructor
        explicit Policies( const std::string& i, const std::string& c ) :
          init(i), coef(c) {}
        // Operator += for adding up two Policies structs
        Policies& operator+= ( const Policies& p ) {
          init += p.init;
          coef += p.coef;
          return *this;
        }
        // Output unique policies to output stream
        friend std::ostream& operator<< ( std::ostream& os, const Policies& p ) {
          Policies copy( p );     // copy policies
          copy.unique();          // get rid of duplicate policies
          os << "i:" << copy.init << ", c:" << copy.coef;
          return os;
        }

      private:
        // Make list of policies unique
        void unique( std::string& list ) {
          std::sort( begin(list), end(list) );
          auto it = std::unique( begin(list), end(list) );
          list.resize( std::distance( begin(list), it ) );
        }
        // Make all policies unique
        void unique() {
          unique( init );
          unique( coef );
        }

        std::string init;
        std::string coef;
    };

    //! Print list: name: option names with policies
    template< class Factory >
    void eqlist( const std::string& title,
                 const Factory& factory,
                 std::size_t ntypes ) const {
      if (!factory.empty()) {
        section( title );
        item( "Unique equation types", ntypes );
        item( "With Cartesian product of policy combinations", factory.size() );
        raw( '\n' );
        raw( m_item_indent + "Legend: equation name : supported policies\n" );
        raw( '\n' );
        raw( m_item_indent + "Policy codes:\n" +
             m_item_indent + " * i: initialization policy: "
                                    "R-raw, "
                                    "Z-zero\n" +
             m_item_indent + " * c: coefficients policy: "
                                    "C-constant\n\n" );
        // extract eqname and supported policies
        const auto ip = ctr::InitPolicy();
        const auto cp = ctr::CoeffPolicy();
        std::map< std::string, Policies > eqs;      // eqname : policies
        for (const auto& f : factory)
          eqs[ DiffEqName( f.first ) ] +=
            Policies( ip.name( f.first.template get< tag::initpolicy >() ),
                      cp.name( f.first.template get< tag::coeffpolicy >() ) );
        // output eqname and supported policies
        for (const auto& e : eqs)
          m_stream << m_item_name_value_fmt % m_item_indent % e.first % e.second;
      }
    }

    //! Print time integration header
    void inthead( const std::string& title, const std::string& name,
                  const std::string& legend, const std::string& head ) const;

    //! Print statistics and PDFs
    void statistics( const std::string& title ) const;

    //! Print configuration of a stack of differential equations
    void diffeqs( const std::string& title,
      const std::vector< std::vector< std::pair< std::string, std::string > > >&
        info ) const;

  private:
    //! Return differential equation name
    template< class Key >
    std::string DiffEqName ( const Key& key ) const
    { return ctr::DiffEq().name( key.template get< tag::diffeq >() ); }

    //! Echo statistics container contents if differs from default applying op
    void stats( const std::string& msg, std::function< std::ostream& (
      std::ostream&, const std::vector< ctr::Term >& ) > op ) const;

    //! Echo pdfs container contents if differs from default applying op
    void pdfs( const std::string& msg,
               std::function<
                 std::ostream& ( std::ostream&,
                                 const std::vector< ctr::Term >&,
                                 const std::vector< tk::real >&,
                                 const std::string&,
                                 const std::vector< tk::real >& ) > op ) const;
};

} // quinoa::

#endif // QuinoaPrint_h
