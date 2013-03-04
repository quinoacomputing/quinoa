//******************************************************************************
/*!
  \file      src/Parser/Parser.C
  \author    J. Bakosi
  \date      Sun 03 Mar 2013 09:55:56 PM MST
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
using namespace control;

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
  Bundle stack(DEFAULTS);
  // Initialize new bool bundle for indicating what data is set in bundle
  BoolBundle boolstack(tuple_size<decltype(DEFAULTS)>::value, false);

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
  m_control->set(JPDF_FILENAME_BASE);
}

void
Parser::unique(vector<Product>& statistics)
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

  if (m_control->set<TITLE>())
    cout << " * Title: " << m_control->get<TITLE>() << endl;

  if (m_control->set<PHYSICS>())
    cout << " * Physics: " << m_control->physicsName() << endl;

  if (m_control->set<HYDRO>())
    cout << " * Hydrodynamics: " << m_control->hydroName() << endl;

  if (m_control->set<MIX>()) {

    cout << " * Material mixing: " << m_control->mixName() << endl;

    if (m_control->set<NSTEP>())
      cout << "   - Number of time steps: " << m_control->get<NSTEP>() << endl;

    if (m_control->set<TERM>())
      cout << "   - Terminate time: " << m_control->get<TERM>() << endl;

    if (m_control->set<DT>())
      cout << "   - Time step size: " << m_control->get<DT>() << endl;

    if (m_control->set<NSCALAR>())
      cout << "   - Number of mixing scalars: " << m_control->get<NSCALAR>()
           << endl;

    if (m_control->set<NPAR>())
      cout << "   - Number of particles: " << m_control->get<NPAR>() << endl;

    if (m_control->set<TTYI>())
      cout << "   - TTY output interval: " << m_control->get<TTYI>() << endl;

    if (m_control->set<DUMP>())
      cout << "   - Dump output interval = " << m_control->get<DUMP>() << endl;

    if (m_control->set<PLTI>())
      cout << "   - Plot output interval = " << m_control->get<PLTI>() << endl;

    if (m_control->set<PDFI>())
      cout << "   - PDF output interval = " << m_control->get<PDFI>() << endl;

    if (m_control->set<B>()) {
      cout << "   - Parameter vector b = {";
      for (auto& v : m_control->get<B>()) cout << " " << v;
      cout << " }" << endl;
    }
    if (m_control->set<S>()) {
      cout << "   - Parameter vector S = {";
      for (auto& v : m_control->get<S>()) cout << " " << v;
      cout << " }" << endl;
    }
    if (m_control->set<KAPPA>()) {
      cout << "   - Parameter vector kappa = {";
      for (auto& v : m_control->get<KAPPA>()) cout << " " << v;
      cout << " }" << endl;
    }
    if (m_control->set<C>()) {
      cout << "   - Parameter vector c = {";
      for (auto& v : m_control->get<C>()) cout << " " << v;
      cout << " }" << endl;
    }
    if (m_control->set<STATISTICS>()) {
      cout << "   - Requested statistics = {";
      for (auto& product : m_control->get<STATISTICS>()) {
        cout << " <";
        for (auto& term : product) {
          //cout << " " << term.field << " " << term.quantity << " " << term.moment;
          cout << term.name;
        }
        cout << ">";
      }
      cout << " }" << endl;
    }
  }

  cout << endl;
}
