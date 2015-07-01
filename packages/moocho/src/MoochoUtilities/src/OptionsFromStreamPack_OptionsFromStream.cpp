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

#include <string>
#include <ostream>
#include <istream>

#include "OptionsFromStreamPack_OptionsFromStream.hpp"
#include "InputStreamHelperPack_EatInputComment.hpp"
#include "Teuchos_Assert.hpp"

// Define this if you want to debug the parser
//#define PRINT_OPTIONS_FROM_STREAM_TRACE

namespace {

// Remove white space from beginning and end.
inline
void clip_ws( std::string* str ) {
  // Remove the first non ' ' characters
  {
    std::string::iterator itr;
    for( itr = str->begin(); itr != str->end(); ++itr )
      if( *itr != ' ' ) break;
    str->erase( str->begin(), itr ); 
  }
  // Remove the last non ' ' characters
  {
    std::string::iterator itr;
    for( itr = str->end() - 1; itr > str->begin() - 1; --itr )
      if( *itr != ' ' ) break;
    str->erase( itr + 1, str->end() ); 
  }
}

}	// end namespace

namespace OptionsFromStreamPack {

namespace OptionsFromStreamUtilityPack {

// class OptionsGroup

std::string OptionsGroup::option_does_not_exist_;

}	// end namespace OptionsFromStreamUtilityPack

// class OptionsFromStream

OptionsFromStream::options_group_t
OptionsFromStream::options_group( const std::string& options_group_name )
{
  options_group_map_t::iterator itr = options_group_map_.find( options_group_name );
  if( itr != options_group_map_.end() ) {
    (*itr).second.second.set(true); // flag that we have accessed this options group
    return options_group_t(&(*itr).second.first);
  }
  return options_group_t(NULL);
}

const OptionsFromStream::options_group_t
OptionsFromStream::options_group( const std::string& options_group_name ) const
{
  options_group_map_t::const_iterator itr = options_group_map_.find( options_group_name );
  if( itr != options_group_map_.end() ) {
    const_cast<false_bool_t&>((*itr).second.second).set(true); // flag that we have accessed this options group
    return options_group_t(const_cast<option_to_value_map_t*>(&(*itr).second.first));
  }
  return options_group_t(NULL);
}

void OptionsFromStream::reset_unaccessed_options_groups()
{
  for( iterator og_itr = begin() ; og_itr != end(); ++og_itr ) {
    (*og_itr).second.second.set(false);	// rest to not accessed yet.
  }
}

void OptionsFromStream::print_unaccessed_options_groups( std::ostream& out ) const
{
  const_iterator og_itr = begin();
  for( ; og_itr != end(); ++og_itr ) {
    if( (*og_itr).second.second == false ) // if options group was not accessed
      out << "options_group " << options_group_name( og_itr ) << " {}\n";
  }
}

void OptionsFromStream::read_options( std::istream& in )
{
  using std::getline;
  
  using InputStreamHelperPack::eat_comment_lines;

  std::string curr_word;

  if(!in)
    return; // No options to read!

#ifdef PRINT_OPTIONS_FROM_STREAM_TRACE
  std::cout << "\n*** Entering OptionsFromStream::read_options(...)!\n\n";
#endif

  // Eat words until you get to begin_options.
  while( curr_word != "begin_options" && !in.eof() )
    in >> curr_word;
#ifdef PRINT_OPTIONS_FROM_STREAM_TRACE
  std::cout << "Found begin_options, start parsing options!\n";
#endif
  // Loop through each options group.
  while( !in.eof() ) {
     eat_comment_lines(in,'*');
    // read options_group
    in >> curr_word;
#ifdef PRINT_OPTIONS_FROM_STREAM_TRACE
    std::cout << "\ncurr_word = \""<<curr_word<<"\"\n";
#endif
    if( curr_word == "}" ) {
#ifdef PRINT_OPTIONS_FROM_STREAM_TRACE
      std::cout << "Found \'}\', Moving on to the next options group or the end!\n";
#endif
      break;
    }
    if( curr_word == "end_options" ) {
#ifdef PRINT_OPTIONS_FROM_STREAM_TRACE
      std::cout << "Found \'end_options\', stoping parsing options!\n";
#endif
      break;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      curr_word != "options_group", InputStreamError
      ,"OptionsFromStream::read_options(...) : "
      "Error, curr_word = \'" << curr_word << " != \'options_group\'" );
    // read the name of the options group up to {
    std::getline(in,curr_word,'{');
    clip_ws( &curr_word );
    const std::string optgroup_name = curr_word;
#ifdef PRINT_OPTIONS_FROM_STREAM_TRACE
    std::cout << "\noptgroup_name = \"" << optgroup_name << "\"\n";
#endif
    // Access the options and values map for this options group.
    option_to_value_map_t& optval = options_group_map_[optgroup_name].first;
    // Grap all of the options for this options group
     eat_comment_lines(in,'*');
    std::string optgroup_options;
    getline(in,optgroup_options,'}');
#ifdef PRINT_OPTIONS_FROM_STREAM_TRACE
    std::cout << "optgroup_options = \"" << optgroup_options << "\"\n";
#endif
    std::istringstream optgroup_options_in(optgroup_options);
    // Loop through and add the options.
    while(true) {
      eat_comment_lines(optgroup_options_in,'*');
      // Get an option and its value
      std::string option_and_value;
      getline(optgroup_options_in,option_and_value,';');
      // Note: above If ';' is missing it will take the rest of the string to
      // the end of '}' for the end of the options group.  These means that if
      // there is no comments after the last option=value pair then the last
      // semicolon is optional!  This turns out to work nicely for the
      // CommandLineOptionsFromStreamProcessor class so this is good behavior!
      clip_ws(&option_and_value);
#ifdef PRINT_OPTIONS_FROM_STREAM_TRACE
      std::cout << "  option_and_value = \"" << option_and_value << "\"\n";
#endif
      if(!option_and_value.length())
        break;
      // Process the option and value
      const std::string::size_type equal_idx = option_and_value.find('=',0);
      TEUCHOS_TEST_FOR_EXCEPTION(
        equal_idx==std::string::npos, std::logic_error,
        "Error, for the option group \"" << optgroup_name << "\""
        << " the option value string \"" << option_and_value << "\" is missing the \"=\" separator!"
        );
      std::string option = option_and_value.substr(0,equal_idx);
      std::string value = option_and_value.substr(equal_idx+1);
      clip_ws(&option);
      clip_ws(&value);
#ifdef PRINT_OPTIONS_FROM_STREAM_TRACE
      std::cout << "    option = \"" << option << "\"\n";
      std::cout << "    value  = \"" << value << "\"\n";
#endif
      optval[option] = value;
    }
  }
#ifdef PRINT_OPTIONS_FROM_STREAM_TRACE
  std::cout << "\n*** Leaving OptionsFromStream::read_options(...)!\n\n";
#endif
}

void OptionsFromStream::print_options( std::ostream& out ) const {
  out << "\nbegin_options\n";
  const_iterator og_itr = begin();
  for( ; og_itr != end(); ++og_itr ) {
    const options_group_t optgrp = OptionsFromStreamPack::options_group( og_itr );
    options_group_t::const_iterator itr = optgrp.begin();
    if(itr == optgrp.end()) continue;
    out << "\noptions_group " << options_group_name( og_itr ) << " {\n";
    for( ; itr != optgrp.end(); ++itr ) {
      const std::string
        &name  = option_name(itr),
        &value = option_value(itr);
      out	<< "    " << name << " = " << value << ";\n";
    }
    out << "}\n";
  }
  out << "\nend_options\n\n";
}

}	// end namespace OptionsFromStreamPack
