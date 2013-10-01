//******************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/Parser.C
  \author    J. Bakosi
  \date      Mon 30 Sep 2013 10:44:58 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator test suite input deck parser
  \details   Random number generator test suite input deck parser
*/
//******************************************************************************

#include <RNGTest/InputDeck/Parser.h>
#include <RNGTest/InputDeck/Grammar.h>

using namespace rngtest;

void
InputDeckParser::parse()
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
//   pegtl::dummy_parse_file<grm::read_file>(m_filename, stack, boolstack);
// #else  // NDEBUG
//   pegtl::basic_parse_file<grm::read_file>(m_filename, stack, boolstack);
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
