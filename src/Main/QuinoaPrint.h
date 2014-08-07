//******************************************************************************
/*!
  \file      src/Main/QuinoaPrint.h
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 08:48:46 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Quinoa's printer
  \details   Quinoa's printer
*/
//******************************************************************************
#ifndef QuinoaPrint_h
#define QuinoaPrint_h

#include <RNGPrint.h>
#include <Quinoa/Options/SDE.h>
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

    //! Print configuration of a model
    template< typename Option, typename... tags >
    void Model( const quinoa::Model& model ) const {
      if (g_inputdeck.get< tags... >() !=
            g_inputdeck_defaults.get< tags... >() ) {
        Option opt;
        subsection( opt.group() + ": " +
                    opt.name( g_inputdeck.get< tags... >() ) );
        printModel( model );
      }
    }

    //! Print configuration of a model in vector
    template< typename Option, typename... tags >
    void Model( const quinoa::Model& model, std::size_t i ) const {
      if (g_inputdeck.get< tags... >() !=
            g_inputdeck_defaults.get< tags... >() ) {
        Option opt;
        subsection( opt.group() + ": " +
                    opt.name( g_inputdeck.get< tags... >()[i] ) );
        printModel( model );
      }
    }

    //! Print list: name: option names with policies...
    template< class Factory >
    void optionlist( const std::string& title, const Factory& factory ) const {
      if (!factory.empty()) section( title );
      for (const auto& f : factory) {
        std::vector< std::string > entries = SDEPolicyNames( f.first );
        std::stringstream ss;
        for (const auto& e : entries) ss << e << " ";
        m_stream << m_item_name_value_fmt
                    % m_item_indent % SDEName( f.first ) % ss.str();
      }
    }

    //! Echo requested statistics if differs from default.
    //! Fields of vector<vector< struct{field, name, plot} >> must exist.
    //! See src/Control/Quinoa/InputDeck/Types.h for the definition of operator
    //! <<= for outputing requested Term and vector<Term>.
    void RequestedStats(const std::string& msg) const;
 
    //! Echo estimated statistics if differs from default.
    //! Fields of vector<vector< struct{field, name, plot} >> must exist.
    //! See src/Control/Quinoa/InputDeck/Types.h for the definition of operator
    //! << for outputing estimated Term and vector<Term>.
    void EstimatedStats(const std::string& msg) const;

  private:
    //! Print configuration of a model
    void printModel( const quinoa::Model& model ) const;

    //! Return SDE name
    std::string SDEName( const ctr::SDEKey& key ) const
    { return ctr::SDE().name( key.get<tag::sde>() ); }

    //! Return SDE policies names
    std::vector< std::string > SDEPolicyNames( const ctr::SDEKey& key ) const;
};

} // quinoa::

#endif // QuinoaPrint_h
