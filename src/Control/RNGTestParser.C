//******************************************************************************
/*!
  \file      src/Control/RNGTestParser.C
  \author    J. Bakosi
  \date      Thu Aug 29 17:00:05 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator test suite parser
  \details   Random number generator test suite parser
*/
//******************************************************************************

#include <pegtl.hh>

#include <RNGTestParser.h>
#include <RNGTestGrammar.h>
#include <Control.h>

using namespace rngtest;

void
RNGTestParser::parse()
//******************************************************************************
//  Parse random number generator test suite control file
//! \author  J. Bakosi
//******************************************************************************
{
//   // Initialize new bundle for parsed data with defaults
//   control::Bundle stack(control::DEFAULTS);
//   // Initialize new bool bundle for indicating what data is set in bundle
//   control::BoolBundle
//     boolstack(std::tuple_size<decltype(control::DEFAULTS)>::value, false);
// 
//   //std::cout << "==== PARSE START ====" << std::endl;
// #ifdef NDEBUG
//   pegtl::dummy_parse_file<grammar::read_file>(m_filename, stack, boolstack);
// #else  // NDEBUG
//   pegtl::basic_parse_file<grammar::read_file>(m_filename, stack, boolstack);
// #endif // NDEBUG
//   //std::cout << "==== PARSE END ====" << std::endl << std::endl;
// 
//   // Filter out repeated statistics
//   unique(std::get<control::STATISTICS>(stack));
// 
//   // Store off parsed bundles
//   m_control->set(stack);
//   m_control->set(boolstack);
}

void
RNGTestParser::echo() const
//******************************************************************************
//  Echo parsed information from random number generator test suite control
//! \author  J. Bakosi
//******************************************************************************
{
//   std::cout << "Parsed from " << m_filename << ":\n" << std::setfill('-')
//             << std::setw(13+m_filename.length()) << "-" << std::endl;
// 
//   if (m_control->set<control::TITLE>())
//     std::cout << " * Title: " << m_control->get<control::TITLE>() << std::endl;
// 
//   if (m_control->set<control::RNGTEST>()) {
//     std::cout << " * RNG test suite: "
//               << grammar::RNGTest.name(m_control->get<control::RNGTEST>())
//               << std::endl;
//     m_control->echoVecOptName<control::RNGS, select::RNG>("Test RNGs");
// 
//   }
// 
//   std::cout << std::endl;
}
