//******************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/Grammar.h
  \author    J. Bakosi
  \date      Mon 30 Sep 2013 09:57:00 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator test suite grammar definition
  \details   Random number generator test suite input deck grammar definition.
  We use the Parsing Expression Grammar Template Library (PEGTL) to create the
  grammar and the associated parser. Credit goes to Colin Hirsch (pegtl@cohi.at)
  for PEGTL. Word of advice: read from the bottom up.
*/
//******************************************************************************
#ifndef RNGTestInputDeckGrammar_h
#define RNGTestInputDeckGrammar_h

namespace rngtest {
//! Grammar definition: state, actions, grammar definition
namespace grm {

  //using namespace pegtl;

  // State

//   //! RNG test (suite) options
//   static ctr::Option<sel::RNGTest> RNGTest;
//   //! RNG options
//   static ctr::LibOption<sel::RNG> RNG;

//   // specialize convert to RNG
//   template<>
//   sel::RNGType convert(const std::string& str) {
//     return RNG.value(str);
//   }

  // Grammar

//   // common to all RNG test suites
//   struct rngtest_common :
//          sor< process< kw::suite,
//                        store_option<sel::RNGTest, ctr::RNGTEST> >,
//               list<kw::rngs, push<ctr::RNGS, sel::RNGType>, rng>
//             > {};
//

//   // rngtest block
//   struct rngtest :
//          ifmust< parse< kw::rngtest,
//                         store_option<sel::Physics, ctr::PHYSICS> >,
//                  block< rngtest_common >
//                > {};
// 

} // grm::
} // rngtest::

#endif // RNGTestInputDeckGrammar_h
