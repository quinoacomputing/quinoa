//******************************************************************************
/*!
  \file      src/Main/InciterPrint.h
  \author    J. Bakosi
  \date      Fri 05 Feb 2016 07:59:50 AM MST
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
#include "ContainerUtil.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck_defaults;
extern ctr::InputDeck g_inputdeck;

//! InciterPrint : tk::Print
class InciterPrint : public tk::Print {

  public:
    //! Constructor
    //! \param[inout] str Verbose stream
    //! \param[inout] qstr Quiet stream
    //! \see tk::RNGPrint::RNGPrint and tk::Print::Print
    //! \author J. Bakosi
    explicit InciterPrint( std::ostream& str = std::clog,
                           std::ostream& qstr = std::cout ) :
      Print( str, qstr ) {}

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

    // Helper class for compact output of PDE policies
    //! \author J. Bakosi
    class Policies {
      public:
        // Default constructor
        explicit Policies() {}
        // Initializer constructor
        explicit Policies( const std::string& p ) :
          prob(p) {}
        // Operator += for adding up two Policies structs
        Policies& operator+= ( const Policies& p ) {
          prob += p.prob;
          return *this;
        }
        // Output unique policies to output stream
        friend std::ostream& operator<< ( std::ostream& os, const Policies& p )
        {
          Policies copy( p );     // copy policies
          copy.unique();          // get rid of duplicate policies
          os << "p:" << copy.prob;
          return os;
        }

      private:
        // Make all policies unique
        void unique() {
          tk::unique( prob );
        }

        std::string prob;
    };

    //! Print equation list with policies
    //! \param[in] title Section title
    //! \param[in] factory Factory to get equation data from
    //! \param[in] ntypes Unique equation types
    //! \author J. Bakosi
    template< class Factory >
    void eqlist( const std::string& title,
                 const Factory& factory,
                 std::size_t ntypes ) const {
      if (!factory.empty()) {
        section( title );
        item( "Unique equation types", ntypes );
        item( "With all policy combinations", factory.size() );
        raw( '\n' );
        raw( m_item_indent + "Legend: equation name : supported policies\n" );
        raw( '\n' );
        raw( m_item_indent + "Policy codes:\n" +
             m_item_indent + " * p: problem configuration:\n" +
             m_item_indent + "   " +
               kw::user_defined::info::name() + " - user-defined\n" +
             m_item_indent + "   " +
               kw::shear_diff::info::name() + " - shear diffusion\n" +
             m_item_indent + "   " +
               kw::slot_cyl::info::name() + " - slotted cylinder\n\n" );
        // extract eqname and supported policies
        const auto p = ctr::Problem();
        std::map< std::string, Policies > eqs;      // eqname : policies
        for (const auto& f : factory)
          eqs[ PDEName( f.first ) ] +=
            Policies( p.name( f.first.template get< tag::problem >() ) );
        // output eqname and supported policies
        for (const auto& e : eqs)
          m_stream << m_item_name_value_fmt % m_item_indent % e.first % e.second;
      }
    }

    //! Print configuration of a stack of partial differential equations
    void pdes( const std::string& title,
      const std::vector< std::vector< std::pair< std::string, std::string > > >&
        info ) const;

    //! Print time integration header
    void inthead( const std::string& title, const std::string& name,
                  const std::string& legend, const std::string& head ) const;

  private:
    //! Return partial differential equation name
    //! \param[in] key Equation key
    //! \return Partial differential equation name based on key
    template< class Key >
    std::string PDEName ( const Key& key ) const
    { return ctr::PDE().name( key.template get< tag::pde >() ); }
};

} // inciter::

#endif // InciterPrint_h
