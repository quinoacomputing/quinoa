//******************************************************************************
/*!
  \file      src/Base/QuinoaPrint.h
  \author    J. Bakosi
  \date      Fri Sep 27 09:08:12 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's printer
  \details   Quinoa's printer
*/
//******************************************************************************
#ifndef QuinoaPrint_h
#define QuinoaPrint_h

#include <Print.h>
#include <QuinoaControl.h>

namespace quinoa {

//! QuinoaPrint : Print
class QuinoaPrint : public Print {

  public:
    //! Bring also Print::item() into scope so if local overload fails
    //! Print::item() is available
    using Print::item;

    //! Constructor
    explicit QuinoaPrint(const QuinoaControl& control,
                         const QuinoaControl& defctr) :
      m_ctr(control), m_def(defctr) {}

    //! Destructor
    ~QuinoaPrint() noexcept override {}

    //! Print item: name : control's value only if differs from its default
    template<typename... tags>
    void item(const std::string& name) const {
      if (m_ctr.get<tags...>() != m_def.get<tags...>())
        std::cout << m_item_name_value_fmt % m_item_indent
                                           % name
                                           % m_ctr.get<tags...>();
    }

    //! Print control option: group : option only if differs from its default
    template<typename OptionType, typename... tags>
    void option() const {
      if (m_ctr.get<tags...>() != m_def.get<tags...>()) {
        ctr::Option<OptionType> opt;
        std::cout << m_item_name_value_fmt % m_item_indent
                                           % opt.group()
                                           % opt.name(m_ctr.get<tags...>());
      }
    }

    //! Echo vector of vector of element names
    //! Fields of vector<vector< struct{field, name, plot} >> must exist
    //! See src/Control/QuinoaControlTypes.h for the definitions of operator <<
    //! for outputing Term and vector<Term>, and operator <<= for outputing
    //! requested (i.e., plotted) Term
    template<typename... tags>
    void vecvecNames(const QuinoaControl& ctr,
                     const std::string& msg,
                     const bool req = false) const {
      std::cout << m_item_name_fmt % m_item_indent % msg;
      if (req) for (auto& v : ctr.get<tags...>()) std::cout <<= v;
      else for (auto& v : ctr.get<tags...>()) std::cout << v;
      std::cout << '\n';
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

    const QuinoaControl& m_ctr;         //!< Parsed control
    const QuinoaControl& m_def;         //!< Default control
};

} // namespace quinoa

#endif // QuinoaPrint_h
