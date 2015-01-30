//******************************************************************************
/*!
  \file      src/Main/QuinoaPrint.h
  \author    J. Bakosi
  \date      Wed 28 Jan 2015 12:40:58 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Quinoa-specific pretty printer functionality
  \details   Quinoa-specific pretty printer functionality.
*/
//******************************************************************************
#ifndef QuinoaPrint_h
#define QuinoaPrint_h

#include <algorithm>

#include <boost/algorithm/string/replace.hpp>

#include <RNGPrint.h>
#include <Quinoa/Types.h>
#include <Quinoa/InputDeck/InputDeck.h>

namespace quinoa {

extern ctr::InputDeck g_inputdeck_defaults;
extern ctr::InputDeck g_inputdeck;

//! QuinoaPrint : tk::RNGPrint
class QuinoaPrint : public tk::RNGPrint {

  public:
    //! Constructor
    //! \param[inout] str Verbose stream
    //! \param[inout] qstr Quiet stream
    //! \see tk::RNGPrint::RNGPrint and tk::Print::Print
    //! \author J. Bakosi
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
    template< typename Option, typename... tags >
    void Item() const {
      if (g_inputdeck.get< tags... >() !=
            g_inputdeck_defaults.get< tags... >() ) {
        Option opt;
        m_stream << m_item_name_value_fmt
                    % m_item_indent % opt.group()
                    % opt.name( g_inputdeck.get< tags... >() );
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

} // quinoa::

#endif // QuinoaPrint_h
