//******************************************************************************
/*!
  \file      src/Base/QuinoaPrint.h
  \author    J. Bakosi
  \date      Thu 13 Feb 2014 09:27:43 PM CET
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's printer
  \details   Quinoa's printer
*/
//******************************************************************************
#ifndef QuinoaPrint_h
#define QuinoaPrint_h

#include <Print.h>
#include <Quinoa/InputDeck/InputDeck.h>
#include <Quinoa/Options/SDE.h>

namespace quinoa {

//! QuinoaPrint : Print
class QuinoaPrint : public tk::Print {

  public:
    //! Bring vanilla overloads from base into scope in case local overloads fail
    //using Print::section;
    //using Print::item;

    //! Constructor
    explicit QuinoaPrint(const std::unique_ptr< ctr::InputDeck >& control) :
      m_ctr(*control) {}

    //! Destructor
    ~QuinoaPrint() override = default;

    //! Print control option: 'group : option' only if differs from its default
    template<typename OptionType, typename... tags>
    void Section() const {
      if (m_ctr.get<tags...>() != ctr::InputDeckDefaults.get<tags...>()) {
        tk::Option< OptionType > opt;
        auto& group = opt.group();
        auto& value = opt.name(m_ctr.get< tags... >());
        m_stream << m_section_title_value_fmt % m_section_indent
                                              % m_section_bullet
                                              % group
                                              % value;
        m_stream << m_section_underline_fmt
                    % m_section_indent
                    % std::string( m_section_indent_size + 3 +
                                   group.size() + value.size(), '-' );
      }
    }

    //! Print item: 'name : value' only if differs from its default
    template<typename... tags>
    void Item(const std::string& name) const {
      if (m_ctr.get<tags...>() != ctr::InputDeckDefaults.get<tags...>())
        m_stream << m_item_name_value_fmt % m_item_indent
                                          % name
                                          % m_ctr.get<tags...>();
    }

    //! Print control option: 'group : option' only if differs from its default
    template<typename OptionType, typename... tags>
    void Item() const {
      if (m_ctr.get<tags...>() != ctr::InputDeckDefaults.get<tags...>()) {
        tk::Option<OptionType> opt;
        m_stream << m_item_name_value_fmt % m_item_indent
                                          % opt.group()
                                          % opt.name(m_ctr.get<tags...>());
      }
    }

    //! Print configuration of a model
    template< typename OptionType, typename... tags >
    void Model( const quinoa::Model* const model ) const {
      if (m_ctr.get<tags...>() != ctr::InputDeckDefaults.get<tags...>()) {
        tk::Option< OptionType > opt;
        // Echo option "group : value" as subsection title
        subsection( opt.group() + ": " + opt.name( m_ctr.get<tags...>() ) );
        // Echo equation type and RNG if model is stochastic
        if (model->stochastic()) {
          m_stream << m_item_name_value_fmt % m_item_indent
                                            % "Equation"
                                            % "stochastic";
          tk::Option< tk::ctr::RNG > rng;
          m_stream << m_item_name_value_fmt % m_item_indent
                                            % rng.group()
                                            % rng.name( model->rng() );
        } else {
          // Echo equation type if model is deterministic
          m_stream << m_item_name_value_fmt % m_item_indent
                                            % "Equation"
                                            % "deterministic";
        }
        // Echo initialization policy
        m_stream << m_item_name_value_fmt % m_item_indent
                                          % "Init policy"
                                          % model->initPolicy();
        // Echo coefficients policy
        m_stream << m_item_name_value_fmt % m_item_indent
                                          % "Coefficients policy"
                                          % model->coeffPolicy();
        // Echo number of components
        m_stream << m_item_name_value_fmt % m_item_indent
                                          % "Number of components"
                                          % model->ncomp();
      }
    }

    //! Print list: name: option names with policies...
    template< class Factory >
    void optionlist( const std::string& title, const Factory& factory ) const {
      if (!factory.empty()) section( title );
      for (const auto& f : factory) {
        std::vector< std::string > entries = SDEPolicyNames( f.first );
        std::stringstream ss;
        for (const auto& e : entries) {
          ss << e << " ";
        }
        m_stream << m_item_name_value_fmt % m_item_indent
                                          % SDEName( f.first )
                                          % ss.str();
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
    //! Don't permit copy constructor
    QuinoaPrint(const QuinoaPrint&) = delete;
    //! Don't permit copy assigment
    QuinoaPrint& operator=(const QuinoaPrint&) = delete;
    //! Don't permit move constructor
    QuinoaPrint(QuinoaPrint&&) = delete;
    //! Don't permit move assigment
    QuinoaPrint& operator=(QuinoaPrint&&) = delete;

    //! Return SDE name
    std::string SDEName( const ctr::SDEKey& key ) const {
      return ctr::SDE().name( key.get<tag::sde>() );
    }

    //! Return SDE policies names
    std::vector< std::string > SDEPolicyNames( const ctr::SDEKey& key ) const;

    const ctr::InputDeck& m_ctr;         //!< Parsed control
};

} // quinoa::

#endif // QuinoaPrint_h
