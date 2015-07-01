#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/qt/TestRunner.h>
#include <cppunit/Outputter.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestFailure.h>
#include <cppunit/Test.h>

#include <qapplication.h>

int main( int argc, char* argv[] )
{
  QApplication app( argc, argv );
  
  CppUnit::QtUi::TestRunner runner;
  CppUnit::Test* test;
  
  test = CppUnit::TestFactoryRegistry::getRegistry("Unit").makeTest();
  runner.addTest( test );
  test = CppUnit::TestFactoryRegistry::getRegistry("Regression").makeTest();
  runner.addTest( test );
  
  runner.run();
  return 0;
}
