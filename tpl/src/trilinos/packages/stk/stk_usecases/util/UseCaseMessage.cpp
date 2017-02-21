// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <iostream>
#include <fstream>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>

using namespace use_case;

int
use_case_message(
  UseCaseEnvironment &  use_case_environment)
{
  int proc_rank = stk::parallel_machine_rank(use_case_environment.m_comm);
//  int proc_size = stk::parallel_machine_size(use_case_environment.m_comm);
  
// Direct calls to report warning and doomeds
  stk::report_symmetric_warning("Direct call warning message");
  stk::report_symmetric_doomed("Direct call doomed message");

// Single statement warning and doomeds reports
  stk::RuntimeWarningSymmetric() << "Single statement warning message";
  stk::RuntimeDoomedSymmetric() << "Single statement doomed message";

// Multiple statement warning and doomeds reports
  bool error = true;
  if (error) {
    stk::RuntimeWarningSymmetric x;
      
    x << "Assembled warning message";
    for (int i = 0; i < 10; ++i) 
      x << " " << i;
  }
  
  if (error) {
    stk::RuntimeDoomedSymmetric x;

    x << "Assembled doomed message";
    for (int i = 0; i < 10; ++i) 
      x << " " << i;
  }

// Message code for limited display    
  for (int i = 0; i < 100; ++i) {
    static stk::MessageCode code(2);
      
    stk::RuntimeWarningSymmetric(code) << "Limited display warning, count = " << i;
  }

// Message code for limited display    
  for (int i = 0; i < 100; ++i) {
    static stk::MessageCode code(2);
      
    stk::RuntimeDoomedSymmetric(code) << "Limited display doomed, count = " << i;
  }
    
  if (proc_rank == 0) {
    sierra::out() << "There were " << stk::get_warning_count() << " warnings" << std::endl;
    sierra::out() << "There were " << stk::get_doomed_count() << " fatal errors" << std::endl << std::endl;
  }
        
// Doomed exception thrown when max exceeded (disabled.  Is this the behavior we want?)
  stk::set_max_doomed_count(2);
  stk::reset_doomed_count();
  try {
    for (int i = 0; i < 100; ++i) {
      stk::RuntimeDoomedSymmetric() << "Throw exception when doomed count is 2, count = " << i;
    }
  }
  catch (std::exception &x) {
    sierra::out() << "Caught exception: " << x.what() << std::endl;
  }

// Message aggregation with limited display
  stk::set_max_doomed_count(1000);
  stk::set_max_warning_count(1000);
  stk::reset_doomed_count();
  stk::reset_warning_count();
    
  for (int i = 0; i < 100; ++i)
    if (error) {
      static stk::MessageCode code(2);
      stk::RuntimeWarningSymmetric x(code);
      
      x << "Every processor may have something to contribute: ";

      std::ostringstream s;
      if (proc_rank%2)
        s << proc_rank;
      
      stk::aggregate_messages(use_case_environment.m_comm, s, ", ");

      x << s.str() << " contributed";
    }

// Print summary    
  if (proc_rank == 0) {
    sierra::out() << "There were " << stk::get_warning_count() << " warnings" << std::endl;
    sierra::out() << "There were " << stk::get_doomed_count() << " fatal errors" << std::endl << std::endl;
  }
        
// Deferred message
  stk::set_max_doomed_count(1000);
  stk::set_max_warning_count(1000);
  stk::reset_doomed_count();
  stk::reset_warning_count();
    
  if (proc_rank%2) {
    static stk::MessageCode code(2);
    stk::RuntimeWarningDeferred x(code);
        
    x << "Deferred warnings from processor ";
    x.aggregate << proc_rank;
      
    static stk::MessageCode code2(2);
    stk::RuntimeDoomedDeferred x2(code2);

    x2 << "Deferred dooms from processor ";
    x2.aggregate << proc_rank;
  }
  stk::report_deferred_messages(use_case_environment.m_comm);

  if (proc_rank == 0) {
    sierra::out() << "There were " << stk::get_warning_count() << " warnings" << std::endl;
    sierra::out() << "There were " << stk::get_doomed_count() << " fatal errors" << std::endl << std::endl;
  }
        
// Deferred warning message with varying key
  for (int i = 0; i < 100; ++i) {
    if (proc_rank%2 == i%3) {
      for (int j = 0; j < 7; ++j) {
          
        static stk::MessageCode code(2);
        stk::RuntimeWarningDeferred x(code);
        
        x << "Deferred warnings with varying key " << i << " from processors ";
        x.aggregate << proc_rank;
      }
    }
  }
  stk::report_deferred_messages(use_case_environment.m_comm);

// Deferred warning message with limited display
  stk::set_max_doomed_count(2);
  stk::set_max_warning_count(2);
  stk::reset_doomed_count();
  stk::reset_warning_count();
  for (int i = 0; i < 100; ++i) {
    if (proc_rank%2 == i%3) {
      static stk::MessageCode code(2);
      stk::RuntimeWarningDeferred x(code);
        
      x << "Deferred warnings with limited display from processors ";
      x.aggregate << proc_rank;
    }
    stk::report_deferred_messages(use_case_environment.m_comm);
  }
  
  sierra::out() << "There were " << stk::get_warning_count() << " warnings" << std::endl;
  sierra::out() << "There were " << stk::get_doomed_count() << " fatal errors" << std::endl;
  sierra::out() << "Done" << std::endl;
  
  return 0;
}
