//******************************************************************************
/*!
  \file      src/Control/RNGTestGrammar.h
  \author    J. Bakosi
  \date      Thu Aug 29 14:50:27 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator test suite grammar definition
  \details   Grammar definition. We use the Parsing Expression Grammar Template
             Library (PEGTL) to create the grammar and the associated parser.
             Credit goes to Colin Hirsch (pegtl@cohi.at) for PEGTL. Word of
             advice: read from the bottom up.
*/
//******************************************************************************
#ifndef RNGTestGrammar_h
#define RNGTestGrammar_h

namespace rngtest {

namespace grammar {

  //using namespace pegtl;

  // State

//   //! RNG test (suite) options
//   static control::Option<select::RNGTest> RNGTest;
//   //! RNG options
//   static control::LibOption<select::RNG> RNG;

//   // specialize convert to RNG
//   template<>
//   select::RNGType convert(const std::string& str) {
//     return RNG.value(str);
//   }

  // Grammar

//   // common to all RNG test suites
//   struct rngtest_common :
//          sor< process< keyword::suite,
//                        store_option<select::RNGTest, control::RNGTEST> >,
//               list<keyword::rngs, push<control::RNGS, select::RNGType>, rng>
//             > {};
//

//   // rngtest block
//   struct rngtest :
//          ifmust< parse< keyword::rngtest,
//                         store_option<select::Physics, control::PHYSICS> >,
//                  block< rngtest_common >
//                > {};
// 

} // namespace grammar

} // namespace Quinoa

#endif // RNGTestGrammar_h
