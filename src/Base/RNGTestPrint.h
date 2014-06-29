//******************************************************************************
/*!
  \file      src/Base/RNGTestPrint.h
  \author    J. Bakosi
  \date      Sat 28 Jun 2014 10:13:57 PM MDT
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

extern ctr::InputDeck g_inputdeck_defaults;
extern ctr::InputDeck g_inputdeck;

//! RNGTestPrint : RNGPrint
class RNGTestPrint : public tk::RNGPrint {

  public:
    //! Bring vanilla overloads from base into scope in case local overloads fail
    using Print::section;
    using Print::item;

    //! Print control option: 'group : option' only if differs from its default
    template< class OptionType, typename... tags >
    void Section() const {
      if (g_inputdeck.get< tags... >() !=
            g_inputdeck_defaults.get< tags... >()) {
        tk::Option< OptionType > opt;
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

    //! Print control option: 'group : option' only if differs from its default
    template<typename OptionType, typename... tags>
    void Item() const {
      if (g_inputdeck.get<tags...>() != g_inputdeck_defaults.get<tags...>()) {
        tk::Option<OptionType> opt;
        m_stream << m_item_name_value_fmt % m_item_indent
                                          % opt.group()
                                          % opt.name(g_inputdeck.get<tags...>());
      }
    }

    //! Print battery option: 'group : option (info)' only if differs from def.
    void battery( std::size_t ntest, std::size_t nstat ) const {
      if (g_inputdeck.get< tag::selected, tag::battery >() !=
          g_inputdeck_defaults.get< tag::selected, tag::battery >() ) {
        tk::Option< ctr::Battery > b;
        auto& group = b.group();
        auto& value = b.name( g_inputdeck.get< tag::selected, tag::battery >() );
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

    //! Print statistical test name(s)
    void names( const std::vector< std::string >& testnames ) const {
      for (const auto& n : testnames)
        m_stream << m_list_item_fmt % m_item_indent % n;
    }

    //! Print statistical tests header (with legend)
    void statshead( const std::string& title ) const {
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

    //! Print one-liner info for test. Columns:
    //! [done/total/failed]
    //!   - done: number of tests completed so far (note that a test may produce
    //!     more than one statistics and thus p-values)
    //!   - total: total number of tests: number of tests in the suite times the
    //!     number of RNGs tested (note that a test may produce more than one
    //!     statistics and thus p-values)
    //!   - failed: number of failed tests by a given RNG so far
    //! name of the statistical test name
    //! name of RNG
    //! result of test: "pass" or "fail, p-value = ..."
    void test( std::size_t ncomplete,
               std::size_t ntest,
               std::map< std::string, std::size_t >& nfail,
               const std::vector< std::vector< std::string > >& status ) const
    {
      // Assumed fields for status:
      // status[0]: vector of name(s) of the test(s),
      //            length: number of p-values
      // status[1]: vector of p-value strings: "pass" or "fail, p-value = ...",
      //            length: number of p-values
      // status[2]: vector of length 1: RNG name used to run the test

      const auto& numfail = nfail.find( status[2][0] );
      Assert( numfail != nfail.end(), tk::ExceptType::FATAL, "Cannot find RNG" );

      // Lambda to count number of failed tests
      auto nfailed = [ &numfail, &status ]() {
        for (const auto& p : status[1]) if (p.size() > 4) ++numfail->second;
        return numfail->second;
      };

      // Construct and echo info-line for all statistics resulted from test
      for (std::size_t t=0; t<status[0].size(); ++t) {
        std::stringstream ss;
        ss << "[" << ncomplete << "/" << ntest << "/" << nfailed()
           << "] " << status[0][t];
        if (t==0) ss << ", " << status[2][0];
        m_stream << m_item_widename_value_fmt
                    % m_item_indent
                    % ss.str()
                    % status[1][t];
      }
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
};

} // rngtest::

#endif // RNGTestPrint_h
