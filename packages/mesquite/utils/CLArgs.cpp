/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2007 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

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
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2008) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file CLArgs.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "CLArgs.hpp"
#include "CLArgFlag.hpp"
#include "ManPage.hpp"
#include <iostream>
#include <sstream>
#include <limits>
#include <stdlib.h>

const char HELP_FLAG = 'h';
const char MAN_FLAG = 'M';

class CLArgImpl
{
  private:

    std::vector<CLArgFlag*> mFlags;
    std::vector< std::string > reqArgNames, optArgNames;
    std::vector< std::string > reqArg, optArg;
    std::string progName, shortDesc, longDesc;
    
  public:
  
    CLArgImpl( const char* progname, const char* brief, const char* desc )
      : progName(progname), shortDesc(brief), longDesc(desc)
      {}
  
    void add_req_arg( const char* name )
      { reqArgNames.push_back( name ); }
  
    void add_opt_arg( const char* name )
      { optArgNames.push_back( name ); }
  
    bool add_parsed_arg( const char* arg );
  
    bool add_flag( CLArgFlag* flag );
    
    CLArgFlag* find_flag( char flag );
    const CLArgFlag* find_flag( char flag ) const
      { return const_cast<CLArgImpl*>(this)->find_flag(flag); }
    
    bool have_required_args() const
      { return reqArgNames.size() == reqArg.size(); }
      
    void get_args( std::vector< std::string >& result ) 
      { result = reqArg;
        result.resize( reqArg.size() + optArg.size() );
        std::copy( optArg.begin(), optArg.end(), result.begin()+reqArg.size() );
      }
      
    void print_help( std::ostream& stream );
    void print_brief_help( std::ostream& stream );
    void print_man( std::ostream& stream );
    
    void print_arg_names( std::ostream& stream );
};

bool CLArgImpl::add_flag( CLArgFlag* arg )
{
  if (find_flag(arg->flag()))
    return false;
  
  mFlags.push_back( arg );
  return true;
}

bool CLArgImpl::add_parsed_arg( const char* arg )
{
  if (reqArg.size() < reqArgNames.size())
    reqArg.push_back( arg );
  else if (optArg.size() < optArgNames.size())
    optArg.push_back( arg );
  else
    return false;
  
  return true;
}

CLArgFlag* CLArgImpl::find_flag( char flag )
{
  for (unsigned i = 0; i < mFlags.size(); ++i)
    if (mFlags[i]->flag() == flag)
      return mFlags[i];
  return 0;
}

void CLArgImpl::print_help( std::ostream& stream )
{
  stream << progName << " : " << shortDesc << std::endl;
  stream  << std::endl << longDesc << std::endl << std::endl;
  print_brief_help( stream );
  stream << '-' << HELP_FLAG << " : Print this help text." << std::endl;
  stream << '-' << MAN_FLAG << " : Print man page text to standard output stream." << std::endl;
  stream << "--" << " : Treat all subsequent arguments as non-flag arguments." << std::endl;
  for (unsigned i = 0; i < mFlags.size(); ++i) {
    stream << '-' << mFlags[i]->flag() << " : " << mFlags[i]->desc();
    std::string extra = mFlags[i]->callback()->desc_append();
    if (!extra.empty()) 
      stream << " " << extra;
    std::string defval = mFlags[i]->callback()->default_str();
    if (!defval.empty()) 
      stream << " (default: " << defval << ")";
    stream << std::endl;
  }
}

void CLArgImpl::print_brief_help( std::ostream& stream )
{
  stream << progName;
  for (unsigned i = 0; i < mFlags.size(); ++i) {
    std::string str = mFlags[i]->callback()->brief();
    if (!str.empty()) {
      stream << "[-" << mFlags[i]->flag() << ' ' << str << "]";
    }
    else {
      str = mFlags[i]->brief();
      if (!str.empty()) 
        stream << " [" << str << "]";
    }
  }
  print_arg_names( stream );
  stream << std::endl;
  stream << progName << " -" << HELP_FLAG << std::endl;
  stream << progName << " -" <<  MAN_FLAG << std::endl;
}

void CLArgImpl::print_man( std::ostream& stream )
{
  ManPage::begin_manpage( stream, progName, 1 );
  
  ManPage::begin_section( stream, "NAME" );
  ManPage::begin_paragraph( stream );
  stream << progName << " - " << shortDesc << std::endl << std::endl;
  
  ManPage::begin_section( stream, "SYNOPSIS" );
  ManPage::begin_hanging_paragraph( stream );
  ManPage::bold( stream, progName );
  for (unsigned i = 0; i < mFlags.size(); ++i) {
    std::string s = mFlags[i]->callback()->manstr();
    if (!s.empty()) {
      stream << '[';
      ManPage::begin_bold( stream );
      stream << '-' << mFlags[i]->flag();
      ManPage::end_bold(stream);
      stream << s << ']';
    }
    else {
      s = mFlags[i]->manstr();
      if (!s.empty()) 
        stream << " [" << s << "]";
    }
  }
  print_arg_names( stream );
  stream << std::endl;
  ManPage::begin_hanging_paragraph( stream );
  ManPage::bold( stream, progName + " -h" );
  ManPage::begin_hanging_paragraph( stream );
  ManPage::bold( stream, progName + " -M" );
  
  ManPage::begin_section( stream, "DESCRIPTION" );
  ManPage::write_text( stream, false, longDesc );
  
  ManPage::begin_section( stream, "OPTIONS" );
  for (unsigned i = 0; i < mFlags.size(); ++i) {
    std::string s = mFlags[i]->callback()->manstr();
    if (!s.empty()) {
      char tmp[] = { '-', mFlags[i]->flag(), ' ', '\0' };
      s = std::string(tmp) + s;
    }
    else {
      s = mFlags[i]->manstr();
      if (s.empty())
        continue;
    }
    
    ManPage::begin_hanging_paragraph( stream );
    stream << s;
    ManPage::begin_indent( stream );
    ManPage::begin_paragraph( stream );
    stream << mFlags[i]->desc();
    s = mFlags[i]->callback()->desc_append();
    if (!s.empty())
      stream << " " << s;
    std::string defval = mFlags[i]->callback()->default_str();
    if (!defval.empty()) 
      stream << " (default: " << defval << ")";
    ManPage::end_indent( stream );
  }
  
}

void CLArgImpl::print_arg_names( std::ostream& stream )
{
  unsigned i;
  for (i = 0; i < reqArgNames.size(); ++i)
    stream << " <" << reqArgNames[i] << ">";
  for (i = 0; i < optArgNames.size(); ++i)
    stream << " [" << optArgNames[i] << "]";
}

CLArgs::CLArgs( const char* progname, const char* brief, const char* desc )
  { impl = new CLArgImpl( progname, brief, desc ); }

CLArgs::~CLArgs()
  { delete impl; }

bool CLArgs::is_flag_available( char fl ) const
{
  return (fl != HELP_FLAG) && (fl != MAN_FLAG) && !(impl->find_flag(fl));
}
    
bool CLArgs::str_flag( char fl, 
                       const char* name, 
                       const char* desc, 
                       CLArgs::StringArgI* callback )
{
  if (!is_flag_available(fl))
    return false;
  return impl->add_flag( new CLArgString( fl, name, desc, callback ) );
}
                      
bool CLArgs::int_flag( char fl, 
                       const char* name, 
                       const char* desc, 
                       CLArgs::IntArgI* callback )
{
  if (!is_flag_available(fl))
    return false;
  return impl->add_flag( new CLArgInt( fl, name, desc, callback ) );
}

bool CLArgs::long_flag( char fl, 
                        const char* name, 
                        const char* desc, 
                        CLArgs::LongArgI* callback )
{
  if (!is_flag_available(fl))
    return false;
  return impl->add_flag( new CLArgLong( fl, name, desc, callback ) );
}

bool CLArgs::double_flag( char fl, 
                          const char* name, 
                          const char* desc, 
                          CLArgs::DoubleArgI* callback )
{
  if (!is_flag_available(fl))
    return false;
  return impl->add_flag( new CLArgDouble( fl, name, desc, callback ) );
}

bool CLArgs::toggle_flag( char on_flag, 
                          char off_flag, 
                          const char* desc, 
                          CLArgs::ToggleArgI* callback )
{
  if (!(is_flag_available(on_flag) && is_flag_available(off_flag)))
    return false;

  CLArgToggle* t1 = new CLArgToggle( on_flag, desc, true, callback );
  impl->add_flag( t1 );
  impl->add_flag( new CLArgToggle( off_flag, desc, t1 ) );
  return true;
}


bool CLArgs::toggle_flag( char fl, const char* desc, CLArgs::ToggleArgI* callback )
{
  if (!is_flag_available(fl))
    return false;
  return impl->add_flag( new CLArgToggle( fl, desc, true, callback ) );
}
  
bool CLArgs::id_list_flag( char fl, const char* desc, CLArgs::IntListArgI* callback )
{
  if (!is_flag_available(fl))
    return false;
  return impl->add_flag( new CLArgIDList( fl, desc, callback ) );
}
  

bool CLArgs::int_list_flag( char fl, 
                            const char* desc, 
                            CLArgs::IntListArgI* callback )
{
  if (!is_flag_available(fl))
    return false;
  return impl->add_flag( new CLArgIntList( fl, desc, callback ) );
}

bool CLArgs::double_list_flag( char fl, 
                               const char* desc, 
                               CLArgs::DoubleListArgI* callback )
{
  if (!is_flag_available(fl))
    return false;
  return impl->add_flag( new CLArgDoubleList( fl, desc, callback ) );
}

bool CLArgs::limit_list_flag( char fl,
                              int num_values,
                              const char* const* value_names )
{
  CLArgFlag* f = impl->find_flag( fl );
  return f ? f->add_set( num_values, value_names ) : false;
}

void CLArgs::add_required_arg( const char* name )
{
  impl->add_req_arg( name );
}

void CLArgs::add_optional_arg( const char* name )
{
  impl->add_opt_arg( name );
}

bool CLArgs::parse_options( int argc, 
                            char** argv,
                            std::vector< std::string >& args_out,
                            std::ostream& error_stream )
{
  std::vector<CLArgFlag*> pending;
  bool no_more_flags = false;
  for (int i = 1; i < argc; ++i) {
    if (!pending.empty()) {
      CLArgFlag* flag = pending.front();
      pending.erase( pending.begin() );
      if (!flag->parse( argv[i] )) {
        error_stream << argv[0] << ": invalid value for flag: -" << flag->flag() << " \"" << argv[i] << '"' << std::endl;
        return false;
      }
    }
    else if (!no_more_flags && argv[i][0] == '-' && argv[i][1] !='\0') {
      for (int j = 1; argv[i][j]; ++j) {
        if (argv[i][j] == HELP_FLAG) {
          print_help( std::cout );
          exit( 0 );
        }
        else if (argv[i][j] == MAN_FLAG) {
          print_man_page( std::cout );
          exit( 0 );
        }
      
        CLArgFlag* flag = impl->find_flag( argv[i][j] );
        if (!flag) {
          error_stream << argv[0] << ": invalid flag: -" << argv[i][j] << std::endl;
          return false;
        }
        else if (!flag->is_toggle()) {
          pending.push_back( flag );
        }
        else if (!flag->parse( NULL )) {
          error_stream << argv[0] << ": conflicting flag: -" << argv[i][j] << std::endl;
          return false;
        }
      }
    }
    else if (!impl->add_parsed_arg( argv[i] )) {
      error_stream << argv[0] << ": unexpected argument: \"" <<argv[i] << '"' << std::endl;
      return false;
    }
  }
  
  impl->get_args( args_out );
  
  if (!pending.empty()) {
    error_stream << argv[0] << ": expected argument following flag: -" << pending.front()->flag() << std::endl;
    return false;
  }
  if (!impl->have_required_args()) {
    error_stream << argv[0] << ": insufficient arguments" << std::endl;
    return false;
  }
  
  return true;
}

void CLArgs::print_help( std::ostream& stream ) const
{
  impl->print_help( stream );
}

void CLArgs::print_man_page( std::ostream& stream ) const
{
  impl->print_man( stream );
}

void CLArgs::print_usage( std::ostream& stream ) const
{
  impl->print_brief_help( stream );
}



void CLArgs::KeyWordArg::initialize( const char* keyword_list[], int list_length )
{
  mKeyWords.resize( list_length );
  std::copy( keyword_list, keyword_list + list_length, mKeyWords.begin() );
}

bool CLArgs::KeyWordArg::value( const std::string& val )
{
  std::vector< std::string >::const_iterator i;
  for (i = mKeyWords.begin(); i != mKeyWords.end(); ++i) 
    if (compare_no_case( i->c_str(), val.c_str() )) {
      return value(*i);
    }
  
  return false;
}

std::string CLArgs::KeyWordArg::brief() const
{
  std::ostringstream ss;
  std::vector< std::string >::const_iterator i = mKeyWords.begin();
  if (i == mKeyWords.end())
    return std::string();
  
  ss << '{' << *i;
  for (++i; i != mKeyWords.end(); ++i)
    ss << '|' << *i;
  ss << '}';
  return ss.str();
}

std::string CLArgs::KeyWordArg::manstr() const
{
  if (mKeyWords.empty())
    return std::string();

  std::ostringstream ss;
  ManPage::bold( ss, mKeyWords[0].c_str() );
  for (unsigned i = 1; i < mKeyWords.size(); ++i) {
    ss << "|";
    ManPage::bold( ss, mKeyWords[i].c_str() );
  }
  return ss.str();
}

bool CLArgs::KeyWordArg::compare_no_case( const char* s1, const char* s2 )
{
  for (; *s1; ++s1, ++s2)
    if (toupper(*s1) != toupper(*s2))
      return false;
  return !*s2;
}

CLArgs::IntRange::IntRange( const int* min, const int* max )
  : mMin( min ? *min : std::numeric_limits<int>::min() ),
    mMax( max ? *max : std::numeric_limits<int>::max() )
  {}

bool CLArgs::IntRange::is_valid( int val ) const
  { return val >= mMin && val <= mMax; }

std::string CLArgs::IntRange::desc_append()  const
{
  std::ostringstream ss;
  ss << "[" << mMin << "," << mMax << "]";
  return ss.str();
}

bool CLArgs::IntRangeArg::value( const int& val )
{
  if (!mRange.is_valid(val))
    return false;
  return IntArg::value(val);
}

bool CLArgs::IntListRangeArg::value( const std::vector<int>&  val )
{
  for (std::vector<int>::const_iterator i = val.begin(); i != val.end(); ++i)
    if (!mRange.is_valid(*i))
      return false;
  return IntListArg::value(val);
}

CLArgs::DoubleRange::DoubleRange( const double* min, const double* max, bool inclusive )
  : haveMin( min != NULL ),
    haveMax( max != NULL ),
    mInclusive( inclusive )
{
  if (haveMin)
    mMin = *min;
  if (haveMax)
    mMax = *max;
}

bool CLArgs::DoubleRange::is_valid( double val ) const
{
  if (mInclusive) 
    return (!haveMin || val >= mMin)
        && (!haveMax || val <= mMax);
  else
    return (!haveMin || val >  mMin)
        && (!haveMax || val <  mMax);
}

std::string CLArgs::DoubleRange::desc_append()  const
{
  std::ostringstream ss;
  if (mInclusive && haveMin)
    ss << '[';
  else 
    ss << '(';
  if (haveMin)
    ss << mMin;
  else
    ss << "-inf";
  ss << ",";
  if (haveMax)
    ss << mMax;
  else
    ss << "inf";
  if (mInclusive && haveMax)
    ss << ']';
  else
    ss << ')';
  return ss.str();
}

bool CLArgs::DoubleRangeArg::value( const double& val )
{
  if (!mRange.is_valid(val))
    return false;
  return DoubleArg::value(val);
}

bool CLArgs::DoubleListRangeArg::value( const std::vector<double>&  val )
{
  for (std::vector<double>::const_iterator i = val.begin(); i != val.end(); ++i)
    if (!mRange.is_valid(*i))
      return false;
  return DoubleListArg::value(val);
}


  
