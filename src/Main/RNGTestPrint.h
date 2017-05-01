// *****************************************************************************
/*!
  \file      src/Main/RNGTestPrint.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     RNGTest-specific pretty printer functionality
  \details   RNGTest-specific pretty printer functionality.
*/
// *****************************************************************************
#ifndef RNGTestPrint_h
#define RNGTestPrint_h

#include "Types.h"
#include "RNGPrint.h"
#include "RNGTest/InputDeck/InputDeck.h"
#include "Flip_map.h"

namespace rngtest {

extern ctr::InputDeck g_inputdeck_defaults;
extern ctr::InputDeck g_inputdeck;

//! RNGTestPrint : tk::RNGPrint
class RNGTestPrint : public tk::RNGPrint {

  public:
    //! Constructor
    //! \param[in,out] str Verbose stream
    //! \param[in,out] qstr Quiet stream
    //! \see tk::RNGPrint::RNGPrint and tk::Print::Print
    //! \author J. Bakosi
    explicit RNGTestPrint( std::ostream& str = std::clog,
                           std::ostream& qstr = std::cout ) :
      RNGPrint( str, qstr ) {}

    //! Bring vanilla overloads from base into scope in case local overloads fail
    using Print::section;
    using Print::item;

    //! Print section only if differs from default
    //! \author J. Bakosi
    template< class Option, typename... tags >
    void Section() const {
      if (g_inputdeck.get< tags... >() !=
            g_inputdeck_defaults.get< tags... >()) {
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

    //! Print control option: 'group : option' only if differs from its default
    //! \author J. Bakosi
    template< class Option, typename... tags>
    void Item() const {
      if (g_inputdeck.get<tags...>() != g_inputdeck_defaults.get<tags...>()) {
        Option opt;
        m_stream << m_item_name_value_fmt % m_item_indent
                                          % opt.group()
                                          % opt.name(g_inputdeck.get<tags...>());
      }
    }

    //! Print battery only if differs from default
    //! \param[in] ntest Number tests in battery
    //! \param[in] nstat Number statistics in battery
    //! \author J. Bakosi
    void battery( std::size_t ntest, std::size_t nstat ) const {
      if (g_inputdeck.get< tag::selected, tag::battery >() !=
          g_inputdeck_defaults.get< tag::selected, tag::battery >() ) {
        ctr::Battery b;
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
                                   group.size() + ss.str().size(), '-' );
      }
    }

    //! Print statistical test name(s)
    //! \param[in] testnames Names of tests
    //! \author J. Bakosi
    void names( const std::vector< std::string >& testnames ) const {
      for (const auto& n : testnames)
        m_stream << m_list_item_fmt % m_item_indent % n;
    }

    //! Print statistical tests header (with legend)
    //! \param[in] t String to use as title
    //! \param[in] npval Number of p-values from tests
    //! \param[in] ntest Number of tests
    //! \author J. Bakosi
    void statshead( const std::string& t, std::size_t npval, std::size_t ntest )
    const {
      std::stringstream ss;
      ss << t << " (" << npval << " stats from " << ntest << " tests)";
      m_stream << m_section_title_fmt % m_section_indent
                                      % m_section_bullet
                                      % ss.str();
      m_stream << m_section_underline_fmt
                  % m_section_indent
                  % std::string( m_section_indent.size() + 2 + ss.str().size(),
                                '-');
      raw( m_item_indent + "Legend: [done/total/failed] Test, RNG : p-value\n" +
           m_item_indent + "(eps  means a value < 1.0e-300)\n" +
           m_item_indent + "(eps1 means a value < 1.0e-15)\n\n" );

    }

    //! \brief Print one-liner info for test
    //! \details Columns: [done/total/failed]
    //!   - done: number of tests completed so far (note that a test may produce
    //!     more than one statistics and thus p-values)
    //!   - total: total number of tests: number of tests in the suite times the
    //!     number of RNGs tested (note that a test may produce more than one
    //!     statistics and thus p-values)
    //!   - failed: number of failed tests by a given RNG so far
    //! name of the statistical test
    //! name of RNG
    //! result of test: "pass" or "fail, p-value = ..."
    //! \param[in] ncomplete Number of completed tests
    //! \param[in] ntest Total number of tests
    //! \param[in] nfail Number of failed tests for RNG
    //! \param[in] status Vector of vector of string with the following assumed
    //!   structure:
    //!   - status[0]: vector of name(s) of the test(s),
    //!                length: number of p-values
    //!   - status[1]: vector of p-value strings: "pass" or "fail, p-value = ...",
    //!                length: number of p-values
    //!   - status[2]: vector of length 1: RNG name used to run the test
    //! \author J. Bakosi
    void test( std::size_t ncomplete,
               std::size_t ntest,
               std::map< std::string, std::size_t >& nfail,
               const std::vector< std::vector< std::string > >& status ) const
    {
      const auto& numfail = nfail.find( status[2][0] );
      Assert( numfail != nfail.end(), "Cannot find RNG" );

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
        (status[1][t] == "pass" ? m_stream : m_qstream) <<
          m_item_widename_value_fmt % m_item_indent % ss.str() % status[1][t];
      }
    }

    //! \brief Print failed statistical test names, RNGs, and p-values
    //! \details Requirements on the template argument, class Failed: must have
    //!    public fields test, rng, and pval.
    //! \param[in] t String to use as title
    //! \param[in] npval Number of p-values from tests
    //! \param[in] nfail Number of failed tests for RNG
    //! \author J. Bakosi
    template< class Failed >
    void failed( const std::string& t,
                 std::size_t npval,
                 const std::vector< Failed >& nfail ) const
    {
      std::stringstream ss;
      ss << t << " (" << nfail.size() << "/" << npval << ")";
      section( ss.str() );
      raw( m_item_indent + "The following tests gave p-values outside "
                           "[0.001, 0.999]\n" +
           m_item_indent + "List groupped by RNG, in the order given in the "
                           "input file\n" +
           m_item_indent + "Legend: Test, RNG : p-value\n" +
           m_item_indent + "(eps  means a value < 1.0e-300)\n" +
           m_item_indent + "(eps1 means a value < 1.0e-15)\n\n" );
      std::string oldname;
      for (const auto& f : nfail) {
        std::string newname( f.rng );
        std::string rngname( newname == oldname ? "" : (", " + newname) );
        oldname = newname;
        m_stream << m_item_widename_value_fmt
                    % m_item_indent
                    % (f.test + rngname)
                    % f.pval;
      }
    }

    //! Print RNGs and their measured run times
    //! \param[in] name Section name
    //! \param[in] costnote A note on how to interpret the costs
    //! \param[in] c Costs for RNGs
    //! \author J. Bakosi
    void cost( const std::string& name,
               const std::string& costnote,
               std::map< std::string, tk::real > c ) const
    {
      Assert( !c.empty(), "Empty map passed to cost()" );
      std::multimap< tk::real, std::string > times = tk::flip_map( c );
      section< tk::QUIET >( name );
      raw< tk::QUIET >( m_item_indent + costnote + "\n\n" );
      tk::real fastest = times.begin()->first;
      for (const auto& t : times) {
        std::stringstream ss;
        ss << t.first << "  (" << std::setprecision(3) << t.first/fastest
           << "x)";
        item< tk::QUIET >( t.second, ss.str() );
      }
    }

    //! Print RNGs and their number of failed tests
    //! \param[in] name Section name
    //! \param[in] ranknote A note on how to interpret ranks
    //! \param[in] f Ranks for RNGs
    //! \author J. Bakosi
    void rank( const std::string& name,
               const std::string& ranknote,
               std::map< std::string, std::size_t > f ) const
    {
      Assert( !f.empty(), "Empty map passed to rank()" );
      std::multimap< std::size_t, std::string > nfail = tk::flip_map( f );
      section< tk::QUIET >( name );
      raw< tk::QUIET >( m_item_indent + ranknote + "\n\n" );
      for (const auto& t : nfail) item< tk::QUIET >( t.second, t.first );
    }
};

} // rngtest::

#endif // RNGTestPrint_h
