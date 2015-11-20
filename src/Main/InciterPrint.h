//******************************************************************************
/*!
  \file      src/Main/InciterPrint.h
  \author    J. Bakosi
  \date      Tue 17 Nov 2015 01:06:42 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Inciter-specific pretty printer functionality
  \details   Inciter-specific pretty printer functionality.
*/
//******************************************************************************
#ifndef InciterPrint_h
#define InciterPrint_h

#include <iostream>
#include <string>

#include <boost/format.hpp>
#include <boost/optional.hpp>

#include "Print.h"
#include "RNGPrint.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck_defaults;
extern ctr::InputDeck g_inputdeck;

//! InciterPrint : tk::RNGPrint
class InciterPrint : public tk::RNGPrint {

  public:
    //! Constructor
    //! \param[inout] str Verbose stream
    //! \param[inout] qstr Quiet stream
    //! \see tk::RNGPrint::RNGPrint and tk::Print::Print
    //! \author J. Bakosi
    explicit InciterPrint( std::ostream& str = std::clog,
                           std::ostream& qstr = std::cout ) :
      RNGPrint( str, qstr ) {}

//     //! Print control option: 'group : option' only if differs from its default
//     template< typename Option, typename... tags >
//     void Section() const {
//       if (g_inputdeck.get< tags... >() !=
//             g_inputdeck_defaults.get< tags... >() ) {
//         Option opt;
//         auto& group = opt.group();
//         auto& value = opt.name( g_inputdeck.get< tags... >() );
//         m_stream << m_section_title_value_fmt % m_section_indent
//                                               % m_section_bullet
//                                               % group
//                                               % value;
//         m_stream << m_section_underline_fmt
//                     % m_section_indent
//                     % std::string( m_section_indent.size() + 3 +
//                                    group.size() + value.size(), '-' );
//       }
//     }
// 
//     //! Print item: 'name : value' only if differs from its default
//     //! \param[in] name Name of item
//     //! \author J. Bakosi
//     template< typename... tags >
//     void Item( const std::string& name ) const {
//       if (g_inputdeck.get< tags... >() !=
//             g_inputdeck_defaults.get< tags... >() )
//         m_stream << m_item_name_value_fmt
//                     % m_item_indent % name % g_inputdeck.get<tags...>();
//     }

    //! Print control option: 'group : option'
    //! \author J. Bakosi
    template< typename Option, typename... tags >
    void Item() const {
      Option opt;
      m_stream << m_item_name_value_fmt
                  % m_item_indent % opt.group()
                  % opt.name( g_inputdeck.get< tags... >() );
    }

    //! Print time integration header
    void inthead( const std::string& title, const std::string& name,
                  const std::string& legend, const std::string& head ) const;
};

} // inciter::

#endif // InciterPrint_h
