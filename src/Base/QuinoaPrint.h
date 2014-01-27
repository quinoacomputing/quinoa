//******************************************************************************
/*!
  \file      src/Base/QuinoaPrint.h
  \author    J. Bakosi
  \date      Mon 27 Jan 2014 09:25:53 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's printer
  \details   Quinoa's printer
*/
//******************************************************************************
#ifndef QuinoaPrint_h
#define QuinoaPrint_h

#include <Print.h>
#include <Quinoa/InputDeck/InputDeck.h>

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
        tk::Option<OptionType> opt;
        auto& group = opt.group();
        auto& value = opt.name(m_ctr.get<tags...>());
        std::cout << m_section_title_value_fmt % m_section_indent
                                               % m_section_bullet
                                               % group
                                               % value;
        std::cout << m_section_underline_fmt
                     % m_section_indent
                     % std::string(m_section_indent_size + 3 +
                                   group.size() + value.size(), '-');
      }
    }

    //! Print item: 'name : value' only if differs from its default
    template<typename... tags>
    void Item(const std::string& name) const {
      if (m_ctr.get<tags...>() != ctr::InputDeckDefaults.get<tags...>())
        std::cout << m_item_name_value_fmt % m_item_indent
                                           % name
                                           % m_ctr.get<tags...>();
    }

    //! Print control option: 'group : option' only if differs from its default
    template<typename OptionType, typename... tags>
    void Item() const {
      if (m_ctr.get<tags...>() != ctr::InputDeckDefaults.get<tags...>()) {
        tk::Option<OptionType> opt;
        std::cout << m_item_name_value_fmt % m_item_indent
                                           % opt.group()
                                           % opt.name(m_ctr.get<tags...>());
      }
    }

    //! Echo requested statistics if differs from default.
    //! Fields of vector<vector< struct{field, name, plot} >> must exist.
    //! See src/Control/Quinoa/InputDeck/Types.h for the definition of operator
    //! <<= for outputing requested Term and vector<Term>.
    void RequestedStats(const std::string& msg) const {
      if (m_ctr.get<tag::stat>() != ctr::InputDeckDefaults.get<tag::stat>()) {
        std::cout << m_item_name_fmt % m_item_indent % msg;
        for (auto& v : m_ctr.get<tag::stat>()) std::cout <<= v;
        std::cout << '\n';
      }
    }

    //! Echo estimated statistics if differs from default.
    //! Fields of vector<vector< struct{field, name, plot} >> must exist.
    //! See src/Control/Quinoa/InputDeck/Types.h for the definition of operator
    //! << for outputing estimated Term and vector<Term>.
    void EstimatedStats(const std::string& msg) const {
      if (m_ctr.get<tag::stat>() != ctr::InputDeckDefaults.get<tag::stat>()) {
        std::cout << m_item_name_fmt % m_item_indent % msg;
        for (auto& v : m_ctr.get<tag::stat>()) std::cout << v;
        std::cout << '\n';
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
          std::cout << m_item_name_value_fmt % m_item_indent
                                             % "Equation"
                                             % "stochastic";
          tk::Option< tk::ctr::RNG > rng;
          std::cout << m_item_name_value_fmt % m_item_indent
                                             % rng.group()
                                             % rng.name( model->rng() );
        } else {
          // Echo equation type if model is deterministic
          std::cout << m_item_name_value_fmt % m_item_indent
                                             % "Equation"
                                             % "deterministic";
        }
        // Echo initialization policy
        std::cout << m_item_name_value_fmt % m_item_indent
                                           % "Init policy"
                                           % model->initPolicy();
        // Echo coefficients policy
        std::cout << m_item_name_value_fmt % m_item_indent
                                           % "Coefficients policy"
                                           % model->coeffPolicy();
        // Echo number of components
        std::cout << m_item_name_value_fmt % m_item_indent
                                           % "Number of components"
                                           % model->ncomp();
      }
    }

  private:
    //! Don't permit copy constructor
    QuinoaPrint(const QuinoaPrint&) = delete;
    //! Don't permit copy assigment
    QuinoaPrint& operator=(const QuinoaPrint&) = delete;
    //! Don't permit move constructor
    QuinoaPrint(QuinoaPrint&&) = delete;
    //! Don't permit move assigment
    QuinoaPrint& operator=(QuinoaPrint&&) = delete;

    const ctr::InputDeck& m_ctr;         //!< Parsed control
};

} // quinoa::

#endif // QuinoaPrint_h
