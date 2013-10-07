//******************************************************************************
/*!
  \file      src/Base/QuinoaPrint.h
  \author    J. Bakosi
  \date      Mon Oct  7 09:08:35 2013
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
    using Print::section;
    using Print::item;

    //! Constructor
    explicit QuinoaPrint(const InputDeck& control) : m_ctr(control) {}

    //! Destructor
    ~QuinoaPrint() noexcept override = default;

    //! Print control option: 'group : option' only if differs from its default
    template<typename OptionType, typename... tags>
    void section() const {
      if (m_ctr.get<tags...>() != QuinoaDefaults.get<tags...>()) {
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
    void item(const std::string& name) const {
      if (m_ctr.get<tags...>() != QuinoaDefaults.get<tags...>())
        std::cout << m_item_name_value_fmt % m_item_indent
                                           % name
                                           % m_ctr.get<tags...>();
    }

    //! Print control option: 'group : option' only if differs from its default
    template<typename OptionType, typename... tags>
    void item() const {
      if (m_ctr.get<tags...>() != QuinoaDefaults.get<tags...>()) {
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
    void requestedStats(const std::string& msg) const {
      if (m_ctr.get<ctr::stat>() != QuinoaDefaults.get<ctr::stat>()) {
        std::cout << m_item_name_fmt % m_item_indent % msg;
        for (auto& v : m_ctr.get<ctr::stat>()) std::cout <<= v;
        std::cout << '\n';
      }
    }

    //! Echo estimated statistics if differs from default.
    //! Fields of vector<vector< struct{field, name, plot} >> must exist.
    //! See src/Control/Quinoa/InputDeck/Types.h for the definition of operator
    //! << for outputing estimated Term and vector<Term>.
    void estimatedStats(const std::string& msg) const {
      if (m_ctr.get<ctr::stat>() != QuinoaDefaults.get<ctr::stat>()) {
        std::cout << m_item_name_fmt % m_item_indent % msg;
        for (auto& v : m_ctr.get<ctr::stat>()) std::cout << v;
        std::cout << '\n';
      }
    }

//     //! Echo vector of Option.names
//     template<class OptionType, typename... tags>
//     void echoVecOptName(const std::string& msg) const {
//       ctr::Option<OptionType> opt;
//       std::cout << "   - " << msg << ": {";
//       for (auto& v : get<tags...>()) {
//         std::cout << " " << opt.name(v);
//       }
//       std::cout << " }" << std::endl;
//     }

  private:
    //! Don't permit copy constructor
    QuinoaPrint(const QuinoaPrint&) = delete;
    //! Don't permit copy assigment
    QuinoaPrint& operator=(const QuinoaPrint&) = delete;
    //! Don't permit move constructor
    QuinoaPrint(QuinoaPrint&&) = delete;
    //! Don't permit move assigment
    QuinoaPrint& operator=(QuinoaPrint&&) = delete;

    const InputDeck& m_ctr;         //!< Parsed control
};

} // quinoa::

#endif // QuinoaPrint_h
