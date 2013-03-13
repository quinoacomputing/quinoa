//******************************************************************************
/*!
  \file      src/Parser/Parser.C
  \author    J. Bakosi
  \date      Tue 12 Mar 2013 10:18:48 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Parser base
  \details   Parser base
*/
//******************************************************************************

#include <pegtl.hh>

#include <Grammar.h>
#include <Control.h>
#include <Defaults.h>
#include <Parser.h>
#include <IOException.h>
#include <ParserException.h>

using namespace Quinoa;

Parser::Parser(const string& filename, Control* const control) :
  m_filename(filename), m_control(control)
//******************************************************************************
//  Constructor
//! \param[in]  filename      File to parse
//! \param[in]  control       Control category where parsed data are stored
//! \author  J. Bakosi
//******************************************************************************
{
  // Check if control file exists, throw exception if not
  m_q.open(m_filename, ifstream::in);
  Assert(m_q.good(), IOException,FATAL,IO_FAILED_OPEN,m_filename);
  m_q.close();
  Assert(!m_q.fail(), IOException,FATAL,IO_FAILED_CLOSE,m_filename);
}

void
Parser::parse()
//******************************************************************************
//  Parse
//! \author  J. Bakosi
//******************************************************************************
{
  // Initialize new bundle for parsed data with defaults
  control::Bundle stack(control::DEFAULTS);
  // Initialize new bool bundle for indicating what data is set in bundle
  control::BoolBundle
    boolstack(tuple_size<decltype(control::DEFAULTS)>::value, false);

  //cout << "==== PARSE START ====" << endl;
#ifdef NDEBUG
  pegtl::dummy_parse_file<grammar::read_file>(m_filename, stack, boolstack);
#else  // NDEBUG
  pegtl::basic_parse_file<grammar::read_file>(m_filename, stack, boolstack);
#endif // NDEBUG
  //cout << "==== PARSE END ====" << endl << endl;

  // Filter out repeated statistics
  unique(get<control::STATISTICS>(stack));

  // Store off parsed bundles
  m_control->set(stack);
  m_control->set(boolstack);
}

void
Parser::unique(vector<control::Product>& statistics)
//******************************************************************************
//  Unique: Make requested statistics unique
//! \param[in,out]  statistics  Vector of statistics
//! \author  J. Bakosi
//******************************************************************************
{
  std::sort(statistics.begin(), statistics.end());
  auto it = std::unique(statistics.begin(), statistics.end());
  statistics.resize( std::distance(statistics.begin(), it) );
}

void
Parser::echo()
//******************************************************************************
//  Echo information on stuff parsed
//! \author  J. Bakosi
//******************************************************************************
{
  cout << "Parsed from " << m_filename << ":\n" << setfill('-')
       << setw(13+m_filename.length()) << "-" << endl;

  if (m_control->set<control::TITLE>())
    cout << " * Title: " << m_control->get<control::TITLE>() << endl;

  if (m_control->set<control::PHYSICS>())
    cout << " * Physics: " << m_control->physicsName() << endl;

  if (m_control->set<control::HYDRO>())
    cout << " * Hydrodynamics: " << m_control->hydroName() << endl;

  if (m_control->set<control::MIX>()) {
    cout << " * Material mixing: " << m_control->mixName() << endl;

    if (m_control->set<control::NSTEP>())
      cout << "   - Number of time steps: " << m_control->get<control::NSTEP>()
           << endl;

    if (m_control->set<control::TERM>())
      cout << "   - Terminate time: " << m_control->get<control::TERM>()
           << endl;

    if (m_control->set<control::DT>())
      cout << "   - Time step size: " << m_control->get<control::DT>() << endl;

    if (m_control->set<control::NSCALAR>())
      cout << "   - Number of mixing scalars: "
           << m_control->get<control::NSCALAR>() << endl;

    if (m_control->set<control::NPAR>())
      cout << "   - Number of particles: " << m_control->get<control::NPAR>()
           << endl;

    if (m_control->set<control::TTYI>())
      cout << "   - TTY output interval: " << m_control->get<control::TTYI>()
           << endl;

    if (m_control->set<control::DUMP>())
      cout << "   - Dump output interval = " << m_control->get<control::DUMP>()
           << endl;

    if (m_control->set<control::PLTI>())
      cout << "   - Plot output interval = " << m_control->get<control::PLTI>()
           << endl;

    if (m_control->set<control::PDFI>())
      cout << "   - PDF output interval = " << m_control->get<control::PDFI>()
           << endl;

    if (m_control->set<control::GLOB>())
      cout << "   - Glob output interval = " << m_control->get<control::GLOB>()
           << endl;

    if (m_control->set<control::B>()) {
      cout << "   - Parameter vector b = {";
      for (auto& v : m_control->get<control::B>()) cout << " " << v;
      cout << " }" << endl;
    }
    if (m_control->set<control::S>()) {
      cout << "   - Parameter vector S = {";
      for (auto& v : m_control->get<control::S>()) cout << " " << v;
      cout << " }" << endl;
    }
    if (m_control->set<control::KAPPA>()) {
      cout << "   - Parameter vector kappa = {";
      for (auto& v : m_control->get<control::KAPPA>()) cout << " " << v;
      cout << " }" << endl;
    }
    if (m_control->set<control::C>()) {
      cout << "   - Parameter vector c = {";
      for (auto& v : m_control->get<control::C>()) cout << " " << v;
      cout << " }" << endl;
    }
    if (m_control->set<control::STATISTICS>()) {
      cout << "   - Requested statistics = {";
      for (auto& product : m_control->get<control::STATISTICS>()) {
        if (product[0].plot) {  // only output user-requested stats
          cout << " <";
          for (auto& term : product) cout << term.name;
          cout << ">";
        }
      }
      cout << " }" << endl;
    }
    if (m_control->set<control::STATISTICS>()) {
      cout << "   - Estimated statistics = {";
      for (auto& product : m_control->get<control::STATISTICS>()) {
        cout << " <";
        for (auto& term : product) cout << term.name;
        cout << ">";
      }
      cout << " }" << endl;
    }
  }

  cout << endl;
}
