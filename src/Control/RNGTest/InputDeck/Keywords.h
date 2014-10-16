//******************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/Keywords.h
  \author    J. Bakosi
  \date      Sun 10 Nov 2013 06:38:27 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \details   All keywords recognized by Quinoa's random number generator (RNG)
  test suite input deck parser. The keywords are defined by specializing struct
  'keyword', defined in Control/Keyword.h. Introducing a new keyword requires a
  more human readable (bust still short) name as well as a short, few-line,
  help-like description.
*/
//******************************************************************************
#ifndef RNGTestInputDeckKeywords_h
#define RNGTestInputDeckKeywords_h

//! Signal to compiler that we are building a list of keywords. This is used by
//! the inline includes, such as *Keywords.h, below (if any) to make sure they
//! get included in the correct namespace and not polluting the global one.
#define Keywords

#include <Keyword.h>
#include <SharedKeywords.h>

namespace rngtest {
//! List of keywords the parser understands
namespace kw {

using namespace pegtl::ascii;
using tk::kw::keyword;

// Keyword 'smallcrush'
struct smallcrush_info {
  static const char* name() { return "SmallCrush"; }
  static const char* help() { return
    "This keyword is used to introduce the description of the random number "
    "generator test suite, i.e., battery, 'SmallCrush'. SmallCrush is a "
    "battery of relatively small number, O(10), of tests, defined in TestU01, "
    "a library for the empirical testing of random number generators. For more "
    "info, see http://www.iro.umontreal.ca/~simardr/testu01/tu01.html.";
  }
};
using smallcrush = keyword< smallcrush_info, s,m,a,l,l,c,r,u,s,h >;

// Keyword 'crush'
struct crush_info {
  static const char* name() { return "Crush"; }
  static const char* help() { return
    "This keyword is used to introduce the description of the random number "
    "generator test suite, i.e., battery, 'Crush'. Crush is a suite of "
    "stringent statistical tests, O(100), defined in TestU01, a library for "
    "the empirical testing of random number generators. For more info, see "
    "http://www.iro.umontreal.ca/~simardr/testu01/tu01.html.";
  }
};
using crush = keyword< crush_info, c,r,u,s,h >;

// Keyword 'bigcrush'
struct bigcrush_info {
  static const char* name() { return "BigCrush"; }
  static const char* help() { return
    "This keyword is used to introduce the description of the random number "
    "generator test suite, i.e., battery, 'BigCrush'. BigCrush is a "
    "suite of very stringent statistical tests, O(100), defined in TestU01, a "
    "library for the empirical testing of random number generators. For more "
    "info, see http://www.iro.umontreal.ca/~simardr/testu01/tu01.html.";
  }
};
using bigcrush = keyword< bigcrush_info, b,i,g,c,r,u,s,h >;

} // kw::
} // rngtest::

#undef Keywords

#endif // RNGTestInputDeckKeywords_h
