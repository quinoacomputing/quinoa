/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/Outputter.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestFailure.h>
#include <cppunit/Test.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <vector>
#include <iostream>
using namespace std;

#include "Mesquite_MsqFPE.hpp"

class CPPUNIT_API SummaryOutput : public CppUnit::Outputter 
{
  public:
    SummaryOutput( FILE* file, CppUnit::TestResultCollector* result ) 
      : file_(file), results_(result) {}
    void write();
  private:
    FILE* file_;
    CppUnit::TestResultCollector* results_;
};


int main(int argc, char **argv)
{
  CppUnit::Test* test;
  vector<CppUnit::Test*> test_list;
  CppUnit::TextUi::TestRunner runner;
  int firsttest = 1;
  bool list = false;
  Mesquite::MsqFPE trap_fpe(true);

    // Check for command line arguments
  if (argc > 2 && !strcmp(argv[1],"-s"))
  {
    FILE* file = fopen(argv[2],"w");
    if (!file) 
    {
      perror( argv[2] );
      exit(1);
    }
    runner.setOutputter( new SummaryOutput( file, &runner.result() ) );
    firsttest += 2;
  }
  else if (argc > 1 && !strcmp(argv[1],"-l"))
  {
    ++firsttest;
    list = true;
  }
  
    // If the user requested a specific test...
  if (argc > firsttest)
  {
    while (argc > firsttest)
    {
      argc--;
      CppUnit::TestFactoryRegistry &registry =
        CppUnit::TestFactoryRegistry::getRegistry(argv[argc]);
      test = registry.makeTest();
      if (!test->countTestCases()) {
        std::cerr << argv[argc] << ": does not match any test or group" << std::endl;
        return 1;
      }
      test_list.push_back( test );
    }
    
  }
    // Otherwise do Unit and Regression suites
  else
  {
     test = CppUnit::TestFactoryRegistry::getRegistry("Unit").makeTest();
     test_list.push_back( test );
     test = CppUnit::TestFactoryRegistry::getRegistry("Regression").makeTest();
     test_list.push_back( test );
  }
  
    // If user just wants list of tests
  if (list) {
    for (vector<CppUnit::Test*>::iterator i = test_list.begin();
         i != test_list.end(); ++i) {
      CppUnit::TestSuite* suite = dynamic_cast<CppUnit::TestSuite*>(*i);
      if (!suite) {
        cout << (*i)->getName() << endl;
        continue;
      }
      const vector<CppUnit::Test*>& list = suite->getTests();
      for (vector<CppUnit::Test*>::const_iterator j = list.begin();
         j != list.end(); ++j) 
        cout << (*j)->getName() << endl;
    }
  }
    // Otherwise run the tests
  else {
    for (vector<CppUnit::Test*>::iterator i = test_list.begin();
         i != test_list.end(); ++i) 
      runner.addTest( *i );
    return !runner.run();
  }
  
    // Return 0 if there were no errors
  return 0;
}

void SummaryOutput::write()
{
  CppUnit::TestResultCollector::TestFailures fails = results_->failures();
  CppUnit::TestResultCollector::Tests tests = results_->tests();
  
  CppUnit::TestResultCollector::TestFailures::const_iterator f_iter = fails.begin();
  CppUnit::TestResultCollector::Tests::const_iterator t_iter;
  
  fprintf(file_,"****Tests Run:\n");
  for (t_iter = tests.begin(); t_iter != tests.end(); ++t_iter)
    fprintf(file_, "%s\n", (*t_iter)->getName().c_str());

  fprintf(file_,"****Tests Failed:\n");
  for (f_iter = fails.begin(); f_iter != fails.end(); ++f_iter)
    fprintf(file_, "%s\n", (*f_iter)->failedTestName().c_str());
 
}

