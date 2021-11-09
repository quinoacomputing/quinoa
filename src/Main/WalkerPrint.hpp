// *****************************************************************************
/*!
  \file      src/Main/WalkerPrint.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
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

#include <brigand/algorithms/for_each.hpp>

#include "NoWarning/format.hpp"

#include "Keywords.hpp"
#include "Print.hpp"
#include "Types.hpp"
#include "Tags.hpp"
#include "RNGPrint.hpp"
#include "ContainerUtil.hpp"
#include "Walker/Options/DiffEq.hpp"
#include "Walker/Options/InitPolicy.hpp"
#include "Walker/Options/CoeffPolicy.hpp"
#include "Walker/InputDeck/InputDeck.hpp"

namespace tk { namespace ctr { struct Term; } }

namespace walker {

extern ctr::InputDeck g_inputdeck_defaults;
extern ctr::InputDeck g_inputdeck;

//! WalkerPrint : tk::RNGPrint
class WalkerPrint : public tk::RNGPrint {

  public:
    //! Constructor
    //! \param[in] screen Screen output filename
    //! \param[in,out] str Verbose stream
    //! \param[in] mode Open mode for screen output file, see
    //!   http://en.cppreference.com/w/cpp/io/ios_base/openmode
    //! \param[in,out] qstr Quiet stream
    //! \see tk::RNGPrint::RNGPrint and tk::Print::Print
    explicit WalkerPrint( const std::string& screen,
                          std::ostream& str = std::clog,
                          std::ios_base::openmode mode = std::ios_base::out,
                          std::ostream& qstr = std::cout ) :
      RNGPrint( screen, str, mode, qstr ) {}

    //! Print section only if differs from its default
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
    template< typename... tags >
    void Item( const std::string& name ) const {
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
