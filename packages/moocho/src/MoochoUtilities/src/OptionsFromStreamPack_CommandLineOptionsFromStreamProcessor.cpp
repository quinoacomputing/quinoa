// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "OptionsFromStreamPack_CommandLineOptionsFromStreamProcessor.hpp"

// Define this if you want to debug the parser
//#define PRINT_COMMAND_LINE_OPTIONS_FROM_STREAM_PROCESSOR_TRACE

namespace OptionsFromStreamPack {

CommandLineOptionsFromStreamProcessor::CommandLineOptionsFromStreamProcessor(
  const std::string  &options_file_name_opt_name
  ,const std::string &options_file_name_opt_doc
  ,const std::string &options_file_name
  ,const std::string &extra_options_str_opt_name
  ,const std::string &extra_options_str_opt_doc
  ,const std::string &extra_options_str
  )
  :options_file_name_opt_name_(options_file_name_opt_name)
  ,options_file_name_opt_doc_(options_file_name_opt_doc)
  ,options_file_name_(options_file_name)
  ,extra_options_str_opt_name_(extra_options_str_opt_name)
  ,extra_options_str_opt_doc_(extra_options_str_opt_doc)
  ,extra_options_str_(extra_options_str)
{}

void CommandLineOptionsFromStreamProcessor::set_options(
  Teuchos::RCP<OptionsFromStream> const& options
  )
{
  options_ = options;
}

Teuchos::RCP<OptionsFromStream>
CommandLineOptionsFromStreamProcessor::get_options() const
{
  return options_;
}

void CommandLineOptionsFromStreamProcessor::setup_commandline_processor(
  Teuchos::CommandLineProcessor *clp
  )
{
  clp->setOption(options_file_name_opt_name().c_str(),&options_file_name_,options_file_name_opt_doc().c_str());
  clp->setOption(extra_options_str_opt_name().c_str(),&extra_options_str_,extra_options_str_opt_doc().c_str());
}

void CommandLineOptionsFromStreamProcessor::process_options()
{
  // Process the file options first
  if(options_file_name_.length()) {
    std::ifstream options_in(options_file_name_.c_str());
    if(options_in) {
      if(!options_.get())
        options_ = Teuchos::rcp(new OptionsFromStream());
      options_->read_options(options_in);
    }
  }
  // Process the extra commandline options
  const int len = extra_options_str_.length();
  if(len) {
    if(!options_.get())
      options_ = Teuchos::rcp(new OptionsFromStream());
    const char colon = ':';
    const std::string::size_type npos = std::string::npos;
    std::ostringstream ooptsstream;
    ooptsstream << "\nbegin_options\n\n";
    std::string::size_type start_i = 0, last_i = 0;
    while(true) {
      last_i = extra_options_str_.find(colon,start_i);
      std::string optgroup = extra_options_str_.substr(start_i,last_i-start_i);
#ifdef PRINT_COMMAND_LINE_OPTIONS_FROM_STREAM_PROCESSOR_TRACE
      std::cout << "\nstart_i = " << start_i;
      std::cout << "\nlast_i = " << last_i;
      std::cout << "\noptgroup (before replacement) = \""<<optgroup<<"\"\n";
#endif
      std::replace( optgroup.begin(), optgroup.end(), ',',';' ); // See above!
#ifdef PRINT_COMMAND_LINE_OPTIONS_FROM_STREAM_PROCESSOR_TRACE
      std::cout << "\noptgroup (after replacement) = \""<<optgroup<<"\"\n";
#endif
      ooptsstream << "options_group " << optgroup << "\n";
      if(last_i == npos) break;
      start_i = last_i + 1;
    }
    ooptsstream << "\nend_options\n";
    const std::string options_str = ooptsstream.str();
#ifdef PRINT_COMMAND_LINE_OPTIONS_FROM_STREAM_PROCESSOR_TRACE
    std::cout << "options_str:\n" << options_str;
#endif
    std::istringstream ioptsstream(options_str);
    options_->read_options(ioptsstream);
  }
}

Teuchos::RCP<OptionsFromStream>
CommandLineOptionsFromStreamProcessor::process_and_get_options()
{
  process_options();
  return get_options();
}

}	// end namespace OptionsFromStreamPack
