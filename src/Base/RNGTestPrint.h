//******************************************************************************
/*!
  \file      src/Base/RNGTestPrint.h
  \author    J. Bakosi
  \date      Fri 29 Nov 2013 10:03:59 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNGTest's printer
  \details   RNGTest's printer
*/
//******************************************************************************
#ifndef RNGTestPrint_h
#define RNGTestPrint_h

#include <Macro.h>
#include <Print.h>
#include <RNGTest/InputDeck/InputDeck.h>

namespace rngtest {

//! RNGTestPrint : Print
class RNGTestPrint : public tk::Print {

  public:
    //! Bring vanilla overloads from base into scope in case local overloads fail
    using Print::section;
    using Print::item;

    //! Constructor
    explicit RNGTestPrint(const std::unique_ptr< ctr::InputDeck >& control) :
      m_ctr(*control) {}

    //! Destructor
    ~RNGTestPrint() override = default;

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

    //! Print all fields of MKL RNG parameters
    template< typename RNG, typename UniformMethod, typename GaussianMethod,
              typename MapType >
    void Mklparams( const MapType& map ) const {
      tk::Option< RNG > rng;
      tk::Option< UniformMethod > um;
      tk::Option< GaussianMethod > gm;
      for (auto& m : map) {
        subsection( rng.name(m.first) );
        std::cout << m_item_name_value_fmt
                     % m_item_indent
                     % "seed"
                     % m.second.template get<quinoa::ctr::seed>();
        std::cout << m_item_name_value_fmt
                     % m_item_indent
                     % um.group()
                     % um.name(
                         m.second.template get<quinoa::ctr::uniform_method>() );
        std::cout << m_item_name_value_fmt
                     % m_item_indent
                     % gm.group()
                     % gm.name(
                         m.second.template get<quinoa::ctr::gaussian_method>() );
        endsubsection();
      }
    }

    //! Print failed statistical tests' names and p-values
    template< class StatTest, class TestContainer >
    void failed( const std::string& name,
                 const std::vector< typename StatTest::Pvals >& pvals,
                 const TestContainer& tests ) const {
      using Psize = typename StatTest::Pvals::size_type;
      using Pval = typename StatTest::Pvals::value_type;
      using Tsize = typename TestContainer::size_type;
      section( name );
      Tsize ntest = tests.size();
      for (Tsize i=0; i<ntest; ++i) {
        Psize npval = pvals[i].size();
        for (Psize p=0; p<npval; ++p) {
          if ( fabs(pvals[i][p]+1.0) > std::numeric_limits<Pval>::epsilon() ) {
            std::cout << m_item_name_value_fmt
                         % m_item_indent
                         % tests[i]->name(p)
                         % pvals[i][p];
          }
        }
      }
    }

  private:
    //! Don't permit copy constructor
    RNGTestPrint(const RNGTestPrint&) = delete;
    //! Don't permit copy assigment
    RNGTestPrint& operator=(const RNGTestPrint&) = delete;
    //! Don't permit move constructor
    RNGTestPrint(RNGTestPrint&&) = delete;
    //! Don't permit move assigment
    RNGTestPrint& operator=(RNGTestPrint&&) = delete;

    const ctr::InputDeck& m_ctr;         //!< Parsed control
};

} // rngtest::

#endif // RNGTestPrint_h
