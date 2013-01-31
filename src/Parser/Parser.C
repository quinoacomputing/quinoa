//******************************************************************************
/*!
  \file      src/Parser/Parser.C
  \author    J. Bakosi
  \date      Thu 31 Jan 2013 06:57:27 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Parser base
  \details   Parser base
*/
//******************************************************************************

#include <pegtl.hh>

#include <Control.h>
#include <Grammar.h>
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
  cout << "==== PARSE START ====" << endl;
  grammar::stack_type stack;
#ifdef NDEBUG
  pegtl::dummy_parse_file< grammar::read_file >( m_filename, stack );
#else  // NDEBUG
  pegtl::basic_parse_file< grammar::read_file >( m_filename, stack );
#endif // NDEBUG
  cout << "==== PARSE END ====" << endl << endl;

  // Store off stuff parsed
  m_control->setPhysics(get<1>(stack));
  m_control->setHydro(get<2>(stack));

//   for (auto& s : stack) {
//     switch (s.first) {
//       case grammar::TITLE : m_control->setTitle(s.second); break;
//       case grammar::HYDRO :
//         try {
//           //pegtl::basic_parse_string< grammar::match_hydro >( s.second );
//         } catch (exception&) {
//           Throw(ParserException,FATAL,UNKNOWN_HYDRO,s.second);
//         }
//         //m_control->setHydro(s.second);
//         break;
//     }
//   }

  cout << static_cast<int>(m_control->physics()) << endl;
  cout << static_cast<int>(m_control->hydro()) << endl;

//   for (auto& s : stack) {
//     cout << s.first << ": " << s.second << endl;
//   }
//   cout << endl;
}

void
Parser::echo()
//******************************************************************************
//  Echo information on stuff parsed
//! \author  J. Bakosi
//******************************************************************************
{
}
