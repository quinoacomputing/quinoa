//******************************************************************************
/*!
  \file      src/Control/RNGTestGrammar.h
  \author    J. Bakosi
  \date      Thu Sep 19 09:35:05 2013
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
namespace grm {

  //using namespace pegtl;

  // State

//   //! RNG test (suite) options
//   static ctr::Option<select::RNGTest> RNGTest;
//   //! RNG options
//   static ctr::LibOption<select::RNG> RNG;

//   // specialize convert to RNG
//   template<>
//   select::RNGType convert(const std::string& str) {
//     return RNG.value(str);
//   }

  // Grammar

//   // common to all RNG test suites
//   struct rngtest_common :
//          sor< process< kw::suite,
//                        store_option<select::RNGTest, ctr::RNGTEST> >,
//               list<kw::rngs, push<ctr::RNGS, select::RNGType>, rng>
//             > {};
//

//   // rngtest block
//   struct rngtest :
//          ifmust< parse< kw::rngtest,
//                         store_option<select::Physics, ctr::PHYSICS> >,
//                  block< rngtest_common >
//                > {};
// 

} // grm::
} // rngtest::

#endif // RNGTestGrammar_h
