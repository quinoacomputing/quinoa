// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Base/TestPrint.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for Base/Print.h
  \details   Unit tests for Base/Print.h
*/
// *****************************************************************************
#ifndef test_Print_h
#define test_Print_h

#include "NoWarning/tut.h"

#include "Print.h"
#include "RNGTest/Options/Battery.h"
#include "UnitTest/CmdLine/CmdLine.h"

namespace tut {

//! All tests in group inherited from this base
struct Print_common {
  Print_common() : verb(), quiet(), prd(), prv( verb ), prq( verb, quiet ) {}
  std::stringstream verb;
  std::stringstream quiet;
  tk::Print prd;                //!< for testing default-constructed object
  tk::Print prv;                //!< for testing verbose stream
  tk::Print prq;                //!< for testing quiet stream
};

//! Test group shortcuts
using Print_group = test_group< Print_common, MAX_TESTS_IN_GROUP >;
using Print_object = Print_group::object;

//! Define test group
static Print_group Print( "Base/Print" );

//! Test definitions for group

//! Test operator<< to default-constructed Print object's verbose stream
//! \author J. Bakosi
template<> template<>
void Print_object::test< 1 >() {
  set_test_name( "operator<< to default-ctor'd verbose stream" );

  // This test most likely never fails as a test (unless it throws an
  // exception), but if there is output, that's not the intended behaviour (as
  // the default should be no verbose output), so it complains on the screen if
  // that's not the case.
  prd << "Testing operator<< to default-constructed tk::Print object, which by "
         "default, should be a /dev/null-like sink, i.e., this should not be "
         "printed.";
}

//! Test operator<< to verbose stream of Print object
//! \author J. Bakosi
template<> template<>
void Print_object::test< 2 >() {
  set_test_name( "operator<< to verbose stream" );

  prv << "blah";
  ensure_equals( "output to verbose stream", verb.str(), "blah" );
}

//! Test operator% to quiet stream of Print object
//! \author J. Bakosi
template<> template<>
void Print_object::test< 3 >() {
  set_test_name( "operator% to quiet stream" );

  prq % "blah";
  ensure_equals( "output to quiet stream", quiet.str(), "blah" );
}

//! Test saving pointer to verbose stream
//! \author J. Bakosi
template<> template<>
void Print_object::test< 4 >() {
  set_test_name( "save/reset pointer to verbose stream" );

  prv << "blah";        // write to verbose stream
  auto p = prv.save();  // save verbose stream pointer
  tk::Print pr;         // instantiate new printer
  // Reset verbose stream in new printer to saved one, has "blah"
  pr.reset( p );
  // Write to new verbose stream, now should have "blah blah"
  pr << " blah";
  ensure_equals( "initialize new object with previously saved verbose stream "
                 "of another object", verb.str(), "blah blah" );
}

//! Test saving pointer to quiet stream
//! \author J. Bakosi
template<> template<>
void Print_object::test< 5 >() {
  set_test_name( "save/reset pointer to quiet stream" );

  prq % "blah";                      // write to quiet stream
  auto p = prq.save< tk::QUIET >();  // save quiet stream pointer
  // Instantiate new printer using two new stringstreams. Note, if we use the
  // default constructor arguments here, the quiet stream will be std::cout,
  // which then will be overwritten by the stream saved out of prq above,
  // containing "blah" so far. Overwriting std::cout results in hijacking
  // std::cout and thus we would see no screen output after the reset call
  // below.
  std::stringstream v, q;
  tk::Print pr( v, q );
  // Reset quiet stream in new printer to saved one, has "blah"
  pr.reset< tk::QUIET >( p );
  // Write to new quiet stream, now should have "blah blah"
  pr % " blah";
  ensure_equals( "initialize new object with previously saved quiet stream "
                 "of another object", quiet.str(), "blah blah" );
}

//! Test that tk::Print::part() runs without throwing an exception
//! \author J. Bakosi
template<> template<>
void Print_object::test< 6 >() {
  set_test_name( "part() does not throw" );
  prv.part( "title" );
}

//! Test that tk::Print::section(title) runs without throwing an exception
//! \author J. Bakosi
template<> template<>
void Print_object::test< 7 >() {
  set_test_name( "section(title) does not throw" );
  prv.section( "title" );
}

//! Test that tk::Print::section(name,value) runs without throwing an exception
//! \author J. Bakosi
template<> template<>
void Print_object::test< 8 >() {
  set_test_name( "section(name,value) does not throw" );
  prv.section( "name", "value" );
}

//! Test that tk::Print::subsection(title) does not throw an exception
//! \author J. Bakosi
template<> template<>
void Print_object::test< 9 >() {
  set_test_name( "subsection(title) does not throw" );
  prv.subsection( "title" );
}

//! Test that tk::Print::title(title) runs without throwing an exception
//! \author J. Bakosi
template<> template<>
void Print_object::test< 10 >() {
  set_test_name( "title(title) does not throw" );
  prv.title( "value" );
  prv.title( "value\n\n asdgfads g f\n\t\t afg   " );
}

//! Test that tk::Print::item(name) does not throw an exception
//! \author J. Bakosi
template<> template<>
void Print_object::test< 11 >() {
  set_test_name( "item(name) does not throw" );
  prv.item( "name" );
}

//! Test that tk::Print::item(name,value) does not throw an exception
//! \author J. Bakosi
template<> template<>
void Print_object::test< 12 >() {
  set_test_name( "item(name,value) does not throw" );
  prv.item( "name", "value" );
  prv.item( "name", 1 );
  prv.item( "name", 3.14 );
}

//! Test that tk::Print::item(name,watch) does not throw an exception
//! \author J. Bakosi
template<> template<>
void Print_object::test< 13 >() {
  set_test_name( "item(name,watch) does not throw" );
  prv.item( "name", tk::Timer::Watch() );
}

//! Test that tk::Print::list(name,entries) does not throw an exception
//! \author J. Bakosi
template<> template<>
void Print_object::test< 14 >() {
  set_test_name( "list(name,entries) does not throw" );
  std::vector< std::string > vec{ "blah1", "blah2", "blah3" };
  prv.list( "name", vec );
  std::list< std::string > list{ "blah1", "blah2", "blah3" };
  prv.list( "name", list );
  std::set< std::string > set{ "blah1", "blah2", "blah3" };
  prv.list( "name", set );
}

//! Test that tk::Print::list(name,factory) does not throw an exception
//! \author J. Bakosi
template<> template<>
void Print_object::test< 15 >() {
  set_test_name( "list(name,factory) does not throw" );

  // Create a fake factory. Only the key is used in list(name,factory), so the
  // value is just an int here.
  using rngtest::ctr::BatteryType;
  std::map< BatteryType, int > factory;
  factory.emplace( BatteryType::SMALLCRUSH, 1 );
  factory.emplace( BatteryType::CRUSH, 2 );
  factory.emplace( BatteryType::BIGCRUSH, 3 );
  prv.list< rngtest::ctr::Battery >( "factory list", factory );
}

//! Test that tk::Print::time(title,clocks) does not throw an exception
//! \author J. Bakosi
template<> template<>
void Print_object::test< 16 >() {
  set_test_name( "time(title,clocks) does not throw" );

  std::vector< std::pair< std::string, tk::Timer::Watch > > timesw;
  timesw.emplace_back( "first timer", tk::Timer::Watch() );
  timesw.emplace_back( "second timer", tk::Timer::Watch() );
  prv.time( "timings as Watch", timesw );

  std::vector< std::pair< std::string, tk::real > > timesd;
  timesd.emplace_back( "first timer", tk::Timer().dsec() );
  timesd.emplace_back( "second timer", tk::Timer().dsec() );
  prv.time( "timings as dsec", timesd );
}

//! Test that tk::Print::note() does not throw an exception
//! \author J. Bakosi
template<> template<>
void Print_object::test< 17 >() {
  set_test_name( "note() does not throw" );
  prv.note( "some note" );
}

//! Test that tk::Print::help() does not throw an exception
//! \author J. Bakosi
template<> template<>
void Print_object::test< 18 >() {
  set_test_name( "help() does not throw" );
  prv.help( "executable", unittest::ctr::CmdLine().get< tag::cmdinfo >(),
            "Command-line Parameters:", "-" );
}

//! Test that tk::Print::helpkw() does not throw an exception
//! \author J. Bakosi
template<> template<>
void Print_object::test< 19 >() {
  set_test_name( "helpkw() does not throw" );

  // Default-construct command line object (will have info on keywords)
  const unittest::ctr::CmdLine cmdline;
  // Get info on all command-line keywords
  const auto& cmdinfo = cmdline.get< tag::cmdinfo >();
  // Find info for keyword 'help'
  auto it = cmdinfo.find( "help" );
  if (it != cmdinfo.end())
    prv.helpkw( "executable", tk::ctr::HelpKw{ it->first, it->second, true } );
  else
    fail( "Couldn't test tk::Print::helpkw(), because couldn't find keyword "
          "'help' in unittest::ctr::CmdLine().get< tag::cmdinfo >()" );
}

//! Test that tk::Print::endpart() does not throw an exception
//! \author J. Bakosi
template<> template<>
void Print_object::test< 20 >() {
  set_test_name( "endpart() does not throw" );
  prv.endpart();
}

//! Test that tk::Print::endsubsection() does not throw an exception
//! \author J. Bakosi
template<> template<>
void Print_object::test< 21 >() {
  set_test_name( "endsubsection() does not throw" );
  prv.endsubsection();
}

//! Test that tk::Print::raw() does not throw an exception
//! \author J. Bakosi
template<> template<>
void Print_object::test< 22 >() {
  set_test_name( "raw() does not throw" );
  prv.raw( "blah" );
}

//! Test tk::Print::stream to access verbose stream
//! \author J. Bakosi
template<> template<>
void Print_object::test< 23 >() {
  set_test_name( "verbose stream access" );

  std::ostream_iterator< int > i( prv.stream< tk::VERBOSE >(), "; " );
  std::fill_n( i, 5, 12 );      // assign the value 12 five times to iterator
  ensure_equals( "repeated values assigned using std::fill_n() to verbose "
                 "stream of tk::Print via std::ostream_iterator created with "
                 "the implicit conversion operator to ostream&",
                 verb.str(), "12; 12; 12; 12; 12; " );
}

//! Test tk::Print::stream to access quiet stream
//! \author J. Bakosi
template<> template<>
void Print_object::test< 24 >() {
  set_test_name( "quiet stream access" );

  std::ostream_iterator< int > i( prq.stream< tk::QUIET >(), "; " );
  std::fill_n( i, 5, 12 );      // assign the value 12 five times to iterator
  ensure_equals( "repeated values assigned using std::fill_n() to quiet "
                 "stream of tk::Print via std::ostream_iterator created with "
                 "the implicit conversion operator to ostream&",
                 quiet.str(), "12; 12; 12; 12; 12; " );
}

//! Test that tk::Print::headerInciter() does not throw an exception
//! \author J. Bakosi
template<> template<>
void Print_object::test< 25 >() {
  set_test_name( "headerInciter() does not throw" );
  prv.headerInciter();
}

//! Test that tk::Print::headerRNGTest() does not throw an exception
//! \author J. Bakosi
template<> template<>
void Print_object::test< 26 >() {
  set_test_name( "headerRNGTest() does not throw" );
  prv.headerRNGTest();
}

//! Test that tk::Print::headerUnitTest() does not throw an exception
//! \author J. Bakosi
template<> template<>
void Print_object::test< 27 >() {
  set_test_name( "headerUnitTest() does not throw" );
  prv.headerUnitTest();
}

//! Test that tk::Print::headerMeshConv() does not throw an exception
//! \author J. Bakosi
template<> template<>
void Print_object::test< 28 >() {
  set_test_name( "headerMeshConv() does not throw" );
  prv.headerMeshConv();
}

//! Test that tk::Print::headerWalker() does not throw an exception
//! \author J. Bakosi
template<> template<>
void Print_object::test< 29 >() {
  set_test_name( "headerWalker() does not throw" );
  prv.headerWalker();
}

} // tut::

#endif // test_Print_h
