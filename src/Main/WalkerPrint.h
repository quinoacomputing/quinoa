// *****************************************************************************
/*!
  \file      src/Main/WalkerPrint.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Walker-specific pretty printer functionality
  \details   Walker-specific pretty printer functionality.
*/
// *****************************************************************************
#ifndef WalkerPrint_h
#define WalkerPrint_h

#include <functional>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <utility>
#include <cstddef>

#include "NoWarning/format.h"
#include "NoWarning/for_each.h"

#include "Keywords.h"
#include "Print.h"
#include "Types.h"
#include "Tags.h"
#include "RNGPrint.h"
#include "ContainerUtil.h"
#include "Walker/Options/DiffEq.h"
#include "Walker/Options/InitPolicy.h"
#include "Walker/Options/CoeffPolicy.h"
#include "Walker/InputDeck/InputDeck.h"

namespace tk { namespace ctr { struct Term; } }

namespace walker {

extern ctr::InputDeck g_inputdeck_defaults;
extern ctr::InputDeck g_inputdeck;

//! WalkerPrint : tk::RNGPrint
class WalkerPrint : public tk::RNGPrint {

  public:
    //! Constructor
    //! \param[in,out] str Verbose stream
    //! \param[in,out] qstr Quiet stream
    //! \see tk::RNGPrint::RNGPrint and tk::Print::Print
    //! \author J. Bakosi
    explicit WalkerPrint( std::ostream& str = std::clog,
                          std::ostream& qstr = std::cout ) :
      RNGPrint( str, qstr ) {}

    //! Print section only if differs from its default
    //! \author J. Bakosi
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
    //! \param[in] name Name of item
    //! \author J. Bakosi
    template< typename... tags >
    void Item( const std::string& name ) const {
      if (g_inputdeck.get< tags... >() !=
            g_inputdeck_defaults.get< tags... >() )
        m_stream << m_item_name_value_fmt
                    % m_item_indent % name % g_inputdeck.get<tags...>();
    }

    //! Print control option: 'group : option' only if differs from its default
    //! \author J. Bakosi
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
    //! \author J. Bakosi
    class Policies {
      public:
        // Default constructor
        explicit Policies() : init(), coef() {}
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
        friend std::ostream& operator<< ( std::ostream& os, const Policies& p )
        {
          Policies copy( p );     // copy policies
          copy.unique();          // get rid of duplicate policies
          os << "i:" << copy.init << ", c:" << copy.coef;
          return os;
        }

      private:
        // Make all policies unique
        void unique() {
          tk::unique( init );
          tk::unique( coef );
        }

        std::string init;
        std::string coef;
    };

    //! Print equation list with policies
    //! \param[in] t Section title
    //! \param[in] factory Factory to get equation data from
    //! \param[in] ntypes Unique equation types
    //! \author J. Bakosi
    template< class Factory >
    void eqlist( const std::string& t,
                 const Factory& factory,
                 std::size_t ntypes ) const
    {
      if (!factory.empty()) {
        section( t );
        item( "Unique equation types", ntypes );
        item( "With all policy combinations", factory.size() );
        raw( '\n' );
        raw( m_item_indent + "Legend: equation name : supported policies\n" );
        raw( '\n' );
        raw( m_item_indent + "Policy codes:\n" );
        static_assert( tk::HasTypedefCode< kw::init::info >::value,
                       "Policy code undefined for keyword" );
        static_assert( tk::HasTypedefCode< kw::coeff::info >::value,
                       "Policy code undefined for keyword" );
        raw( m_item_indent + " * " + *kw::init::code() + ": "
                           + kw::init::name() + ":\n" );
        boost::mpl::for_each< ctr::InitPolicy::keywords >( echoPolicies(this) );
        raw( m_item_indent + " * " + *kw::coeff::code() + ": "
                           + kw::coeff::name() + ":\n" );
        boost::mpl::for_each< ctr::CoeffPolicy::keywords >( echoPolicies(this) );
        raw( '\n' );
        // extract eqname and supported policies
        const auto ip = ctr::InitPolicy();
        const auto cp = ctr::CoeffPolicy();
        std::map< std::string, Policies > eqs;      // eqname : policies
        for (const auto& f : factory)
          eqs[ DiffEqName( f.first ) ] +=
            Policies( ip.code( f.first.template get< tag::initpolicy >() ),
                      cp.code( f.first.template get< tag::coeffpolicy >() ) );
        // output eqname and supported policies
        for (const auto& e : eqs)
          m_stream << m_item_name_value_fmt % m_item_indent
                                            % e.first % e.second;
      }
    }

    //! Print time integration header
    void inthead( const std::string& t, const std::string& name,
                  const std::string& legend, const std::string& head ) const;

    //! Print statistics and PDFs
    void statistics( const std::string& t ) const;

    //! Print configuration of a stack of differential equations
    void diffeqs( const std::string& t,
      const std::vector< std::vector< std::pair< std::string, std::string > > >&
        info ) const;

  private:
    //! Return differential equation name
    //! \param[in] key Equation key
    //! \return Differential equation name based on key
    template< class Key >
    std::string DiffEqName ( const Key& key ) const
    { return ctr::DiffEq().name( key.template get< tag::diffeq >() ); }

    //! Echo statistics container contents if differs from default
    void stats( const std::string& msg ) const;

    //! Echo pdfs container contents if differs from default applying op
    void pdfs( const std::string& msg,
               std::function<
                 std::ostream& ( std::ostream&,
                                 const std::vector< tk::ctr::Term >&,
                                 const std::vector< tk::real >&,
                                 const std::string&,
                                 const std::vector< tk::real >& ) > op ) const;
};

} // walker::

#endif // WalkerPrint_h
