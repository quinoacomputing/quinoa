//******************************************************************************
/*!
  \file      src/Base/RNGTestPrint.h
  \author    J. Bakosi
  \date      Thu 26 Dec 2013 02:12:15 PM MST
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
    template< class OptionType, typename... tags >
    void Section() const {
      if (m_ctr.get< tags... >() != ctr::InputDeckDefaults.get< tags... >()) {
        tk::Option< OptionType > opt;
        auto& group = opt.group();
        auto& value = opt.name( m_ctr.get< tags... >() );
        std::cout << m_section_title_value_fmt % m_section_indent
                                               % m_section_bullet
                                               % group
                                               % value;
        std::cout << m_section_underline_fmt
                     % m_section_indent
                     % std::string( m_section_indent_size + 3 +
                                    group.size() + value.size(), '-' );
      }
    }

    //! Print battery option: 'group : option (info)' only if differs from def.
    template< class Psize, class Tsize >
    void battery( const Tsize& ntest, const Psize& nstat ) const {
      if (m_ctr.get< ctr::selected, ctr::battery >() !=
          ctr::InputDeckDefaults.get< ctr::selected, ctr::battery >() ) {
        tk::Option< ctr::Battery > bat;
        auto& group = bat.group();
        auto& value = bat.name( m_ctr.get< ctr::selected, ctr::battery >() );
        std::stringstream ss;
        ss << value << " (" << ntest << " tests, " << nstat << " stats)";
        std::cout << m_section_title_value_fmt % m_section_indent
                                               % m_section_bullet
                                               % group
                                               % ss.str();
        std::cout << m_section_underline_fmt
                     % m_section_indent
                     % std::string( m_section_indent_size + 3 +
                                    group.size() + value.size(), '-' );
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
    template< class MapType >
    void MKLParams( const std::vector< quinoa::ctr::RNGType >& vec,
                    const MapType& map ) const
    {
      quinoa::ctr::RNG rng;
      for (const auto& r : vec) {
        if (rng.lib(r) == quinoa::ctr::RNGLibType::MKL) {
          subsection( rng.name(r) );
          const auto& m = map.find(r);
          if (m == map.end()) {   // no parameter map entry, print defaults
            echoMKLParams( quinoa::ctr::MKLRNGParam() );
          } else {
            echoMKLParams( m->second );
          }
          endsubsection();
        }
      }
    }

    //! Print all fields of RNGSSE parameters
    template< class MapType >
    void RNGSSEParams( const std::vector< quinoa::ctr::RNGType >& vec,
                       const MapType& map ) const
    {
      quinoa::ctr::RNG rng;
      for (const auto& r : vec) {
        if (rng.lib(r) == quinoa::ctr::RNGLibType::RNGSSE) {
          subsection( rng.name(r) );
          const auto& m = map.find(r);
          if (m == map.end()) {   // no parameter map entry, print defaults
            echoRNGSSEParams( quinoa::ctr::RNGSSEParam() );
          } else {
            echoRNGSSEParams( m->second );
          }
          endsubsection();
        }
      }
    }

    //! Print statistical test names
    template< class StatTest, class TestContainer >
    void names( const TestContainer& tests,
                const typename TestContainer::size_type& ntest ) const
    {
      using Psize = typename StatTest::Psize;
      using Tsize = typename TestContainer::size_type;
      for (Tsize i=0; i<ntest; ++i) {
        Psize npval = tests[i]->nstat();
        for (Psize p=0; p<npval; ++p) {
          std::string name( tests[i]->name(p) );
          if (p>0) name += " *";
          std::cout << m_list_item_fmt % m_item_indent % name;
        }
      }
      raw( "\n" );
      raw( m_item_indent + "Note: Tests followed by an asterisk (*) are\n" +
           m_item_indent + "statistics computed from the preceding test.\n" );
    }

    //! Print statistical tests header (with legend)
    void statshead(const std::string& title) const {
      std::cout << m_section_title_fmt % m_section_indent
                                       % m_section_bullet
                                       % title;
      std::cout << m_section_underline_fmt
                   % m_section_indent
                   % std::string(m_section_indent_size + 2 + title.size(),'-');
      raw( m_item_indent +
           "Legend: [done/total/failed] Test, RNG : p-value\n\n" );
    }

    //! Print a single test name, RNG and pass or fail + p-value
    template< class StatTest, class TestContainer >
    void test( const typename StatTest::Psize& ncomplete,
               const typename StatTest::Rsize& nfail,
               const typename StatTest::Psize& npval,
               const typename TestContainer::value_type& test,
               const typename StatTest::Psize& p ) const
    {
      // Construct info-line
      tk::Option< quinoa::ctr::RNG > rng;
      std::stringstream ss;
      ss << "[" << ncomplete << "/" << npval << "/" << nfail << "] "
         << test->name(p) << ", " << rng.name(test->rng());
      std::string pvalstr("pass");
      // Put in p-value if test failed
      if (test->fail(p)) pvalstr = "fail, p-value = " + test->pvalstr(p);
      // Output
      std::cout << m_item_widename_value_fmt
                   % m_item_indent
                   % ss.str()
                   % pvalstr;
    }

    //! Print failed statistical test names, RNGs, and p-values
    template< class StatTest, class TestContainer >
    void failed(
      const std::string& name,
      const typename std::vector< typename StatTest::Pvals >::size_type& total,
      const typename std::vector< typename StatTest::Pvals >::size_type& failed,
      const TestContainer& tests ) const
    {
      using Psize = typename StatTest::Psize;
      using Tsize = typename TestContainer::size_type;
      std::stringstream ss;
      ss << name << " (" << failed << "/" << total << ")";
      section( ss.str() );
      raw( m_item_indent + "The following tests gave p-values outside "
                           "[0.001, 0.999]\n" +
           m_item_indent + "List groupped by RNG, in the order given in the "
                           "input file\n" +
           m_item_indent + "Legend: Test, RNG : p-value\n" +
           m_item_indent + "(eps  means a value < 1.0e-300)\n" +
           m_item_indent + "(eps1 means a value < 1.0e-15)\n\n" );
      Tsize ntest = tests.size();
      std::string oldname;
      for (Tsize i=0; i<ntest; ++i) {
        Psize npval = tests[i]->nstat();
        for (Psize p=0; p<npval; ++p) {
          if ( tests[i]->fail(p) ) {
            tk::Option< quinoa::ctr::RNG > rng;
            std::string newname( rng.name( tests[i]->rng() ) );
            std::string rngname( newname == oldname ? "" : (", " + newname) );
            oldname = newname;
            std::cout << m_item_widename_value_fmt
                         % m_item_indent
                         % (tests[i]->name(p) + rngname)
                         % tests[i]->pvalstr(p);
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

    void echoMKLParams( const quinoa::ctr::MKLRNGParam& p ) const {
      tk::Option< quinoa::ctr::MKLUniformMethod > um;
      tk::Option< quinoa::ctr::MKLGaussianMethod > gm;
      std::cout << m_item_name_value_fmt
                   % m_item_indent
                   % "seed"
                   % p.get<quinoa::ctr::seed>();
      std::cout << m_item_name_value_fmt
                   % m_item_indent
                   % um.group()
                   % um.name( p.get<quinoa::ctr::uniform_method>() );
      std::cout << m_item_name_value_fmt
                   % m_item_indent
                   % gm.group()
                   % gm.name( p.get<quinoa::ctr::gaussian_method>() );
    }

    void echoRNGSSEParams( const quinoa::ctr::RNGSSEParam& p ) const {
      tk::Option< quinoa::ctr::RNGSSESeqLen > seq;
      std::cout << m_item_name_value_fmt
                   % m_item_indent
                   % "seed"
                   % p.get<quinoa::ctr::seed>();
      std::cout << m_item_name_value_fmt
                   % m_item_indent
                   % seq.group()
                   % seq.name( p.get<quinoa::ctr::seqlen>() );
    }

    const ctr::InputDeck& m_ctr;         //!< Parsed control
};

} // rngtest::

#endif // RNGTestPrint_h
