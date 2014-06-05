//******************************************************************************
/*!
  \file      src/Base/RNGTestPrint.h
  \author    J. Bakosi
  \date      Thu 29 May 2014 10:58:49 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNGTest's printer
  \details   RNGTest's printer
*/
//******************************************************************************
#ifndef RNGTestPrint_h
#define RNGTestPrint_h

#include <Types.h>
#include <RNGPrint.h>
#include <RNGTest/InputDeck/InputDeck.h>

namespace rngtest {

//! RNGTestPrint : RNGPrint
class RNGTestPrint : public tk::RNGPrint {

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

    //! Print battery option: 'group : option (info)' only if differs from def.
    void battery( std::size_t ntest, std::size_t nstat ) const {
      if (m_ctr.get< tag::selected, tag::battery >() !=
          ctr::InputDeckDefaults.get< tag::selected, tag::battery >() ) {
        tk::Option< ctr::Battery > bat;
        auto& group = bat.group();
        auto& value = bat.name( m_ctr.get< tag::selected, tag::battery >() );
        std::stringstream ss;
        ss << value << " (" << ntest << " tests, " << nstat << " stats)";
        m_stream << m_section_title_value_fmt % m_section_indent
                                              % m_section_bullet
                                              % group
                                              % ss.str();
        m_stream << m_section_underline_fmt
                    % m_section_indent
                    % std::string( m_section_indent.size() + 3 +
                                   group.size() + value.size(), '-' );
      }
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

    //! Print statistical test names
    template< class StatTest, class TestContainer >
    void names( const TestContainer& tests, std::size_t ntest ) const
    {
      for (std::size_t i=0; i<ntest; ++i) {
        auto npval = tests[i]->nstat();
        for (std::size_t p=0; p<npval; ++p) {
          std::string name( tests[i]->name(p) );
          if (p>0) name += " *";
          m_stream << m_list_item_fmt % m_item_indent % name;
        }
      }
      raw( '\n' );
      raw( m_item_indent + "Note: Tests followed by an asterisk (*) are\n" +
           m_item_indent + "statistics computed from the preceding test.\n" );
    }

    //! Print statistical tests header (with legend)
    void statshead(const std::string& title) const {
      m_stream << m_section_title_fmt % m_section_indent
                                      % m_section_bullet
                                      % title;
      m_stream << m_section_underline_fmt
                  % m_section_indent
                  % std::string(m_section_indent.size() + 2 + title.size(),'-');
      raw( m_item_indent + "Legend: [done/total/failed] Test, RNG : p-value\n" +
           m_item_indent + "(eps  means a value < 1.0e-300)\n" +
           m_item_indent + "(eps1 means a value < 1.0e-15)\n\n" );

    }

    //! Print a single test name, RNG and pass or fail + p-value
    template< class StatTest, class TestContainer >
    void test( std::size_t ncomplete,
               std::size_t nfail,
               std::size_t npval,
               const typename TestContainer::value_type& tst,
               std::size_t p ) const
    {
      // Construct info-line
      tk::Option< tk::ctr::RNG > rng;
      std::stringstream ss;
      ss << "[" << ncomplete << "/" << npval << "/" << nfail << "] "
         << tst->name(p) << ", " << rng.name(tst->rng());
      std::string pvalstr("pass");
      // Put in p-value if test failed
      if (tst->fail(p)) pvalstr = "fail, p-value = " + tst->pvalstr(p);
      // Output
      m_stream << m_item_widename_value_fmt
                  % m_item_indent
                  % ss.str()
                  % pvalstr;
    }

    //! Print failed statistical test names, RNGs, and p-values
    template< class StatTest, class TestContainer >
    void failed( const std::string& name, std::size_t total, std::size_t fail,
                 const TestContainer& tests ) const
    {
      std::stringstream ss;
      ss << name << " (" << fail << "/" << total << ")";
      section( ss.str() );
      raw( m_item_indent + "The following tests gave p-values outside "
                           "[0.001, 0.999]\n" +
           m_item_indent + "List groupped by RNG, in the order given in the "
                           "input file\n" +
           m_item_indent + "Legend: Test, RNG : p-value\n" +
           m_item_indent + "(eps  means a value < 1.0e-300)\n" +
           m_item_indent + "(eps1 means a value < 1.0e-15)\n\n" );
      std::size_t ntest = tests.size();
      std::string oldname;
      for (std::size_t i=0; i<ntest; ++i) {
        auto npval = tests[i]->nstat();
        for (std::size_t p=0; p<npval; ++p) {
          if ( tests[i]->fail(p) ) {
            tk::Option< tk::ctr::RNG > rng;
            std::string newname( rng.name( tests[i]->rng() ) );
            std::string rngname( newname == oldname ? "" : (", " + newname) );
            oldname = newname;
            m_stream << m_item_widename_value_fmt
                        % m_item_indent
                        % (tests[i]->name(p) + rngname)
                        % tests[i]->pvalstr(p);
          }
        }
      }
    }

    //! Print RNGs and their measured runtimes
    //! (taking a copy of rngtimes for sorting)
    void cost( const std::string& name,
               const std::string& costnote,
               std::vector< std::pair< tk::real, std::string > > rngtimes )
    const {
      std::sort( begin(rngtimes), end(rngtimes) );
      section( name );
      raw( m_item_indent + costnote + "\n\n" );
      tk::real fastest = rngtimes[0].first;
      for (const auto& t : rngtimes) {
        std::stringstream ss;
        ss << t.first << "  ("
           << std::setprecision(3) << t.first/fastest << "x)";
        item( t.second, ss.str() );
      }
    }

    //! Print RNGs and their number of failed tests
    //! (taking a copy of rngnfail for sorting)
    void rank( const std::string& name,
               const std::string& ranknote,
               std::vector< std::pair< std::size_t, std::string > > rngnfail )
    const {
      std::sort( begin(rngnfail), end(rngnfail) );
      section( name );
      raw( m_item_indent + ranknote + "\n\n" );
      for (const auto& t : rngnfail) {
        item( t.second, t.first );
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
