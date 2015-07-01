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


/** \file CLArgFlag.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "CLArgFlag.hpp"
#include "ManPage.hpp"
#include <ctype.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>

template <typename T> static std::string stringify( T value ) {
  std::ostringstream ss;
  ss << value;
  return ss.str();
}

/************************* CLArgFlag ********************************/

CLArgFlag::~CLArgFlag() {}

std::string CLArgFlag::make_man_string( int count, const char* names[] ) const
{
  std::ostringstream ss;
  ManPage::begin_bold(ss) << "-" << flag();
  ManPage::end_bold(ss);
  if (count) {
    ss << "<" << names[0];
    for (int i = 1; i < count; ++i) {
      ss << "|";
    }
    ss << ">";
  }
  return ss.str();
}

std::string CLArgFlag::make_literal_man_string( int count, const char* args[] ) const
{
  std::ostringstream ss;
  ManPage::begin_bold(ss) << "-" << flag();
  ManPage::end_bold(ss);
  if (count) {
    ss << args[0];
    for (int i = 1; i < count; ++i)
      ManPage::bold( ss, args[i] );
  }
  return ss.str();
}  

bool CLArgFlag::is_toggle() const
  { return false; }

bool CLArgFlag::add_set( int, const char* const* )
  { return false; }

/************************* CLArgToggle ********************************/

bool CLArgToggle::parse( const char* ) const
{
  if (!mCallback->value( mValue ))
    return false;
  mCallback->set_seen();
  return true;
}

std::string CLArgToggle::brief() const
{
  if (!mOpposite) {
    const char result[] = { '-', flag(), '\0' };
    return result;
  }
  else if (mOpposite->flag() < flag()) {
    return std::string();
  }
  else {
    const char result[] = { '-', flag(), '|', '-', mOpposite->flag(), '\0' };
    return result;
  }
}

std::string CLArgToggle::manstr() const
{
  if (!mOpposite) {
    return make_man_string();
  }
  else if (mOpposite->flag() > flag()) {
    std::ostringstream result;
    ManPage::begin_bold(result);
    result << "-" << flag();
    ManPage::end_bold(result);
    result << "|";
    ManPage::begin_bold(result);
    result << "-" << mOpposite->flag();
    ManPage::end_bold(result);
    return result.str();
  }
  else {
    return std::string();
  }
}
    

/************************* CLArgString ***********************************/

std::string CLArgString::brief() const
{
  std::ostringstream ss;
  ss << '-' << flag() << "  <" << mName << '>';
  return ss.str();
}

std::string CLArgString::manstr() const
  { return make_man_string( mName.c_str() );  }

bool CLArgString::parse( const char* option ) const
{
  if (!mCallback->value( option ))
    return false;
  mCallback->set_seen();
  return true;
}

/************************* CLArgLong ***********************************/

std::string CLArgLong::brief() const
{
  const char str[] = { '-', flag(), ' ', '<', '\0' };
  std::string result( str );
  result += mName;
  result += ">";
  return result;
}
std::string CLArgLong::manstr() const
{
  return make_man_string( &mName[0] );
}

bool CLArgLong::parse( const char* option ) const
{
  long l;
  char* end_ptr;
  l = strtol( option, &end_ptr, 0 );
  if (*end_ptr || !mCallback->value(l))
    return false;
  mCallback->set_seen();
  return true;
}

/************************* CLArgInt ***********************************/

std::string CLArgInt::brief() const
{
  const char str[] = { '-', flag(), ' ', '<', '\0' };
  std::string result( str );
  result += mName;
  result += ">";
  return result;
}
std::string CLArgInt::manstr() const
{
  return make_man_string( &mName[0] );
}

bool CLArgInt::parse( const char* option ) const
{
  long l;
  char* end_ptr;
  l = strtol( option, &end_ptr, 0 );
  int i = (int)l;
  if (*end_ptr || (long)i != l || !mCallback->value(i))
    return false;
  mCallback->set_seen();
  return true;
}



/************************* CLArgDouble ********************************/

std::string CLArgDouble::brief() const
{
  const char str[] = { '-', flag(), ' ', '<', '\0' };
  std::string result( str );
  result += mName;
  result += ">";
  return result;
}
std::string CLArgDouble::manstr() const
{
  return make_man_string( &mName[0] );
}

bool CLArgDouble::parse( const char* option ) const
{
  char* endptr;
  double val = strtod( option, &endptr );
  if (*endptr || !mCallback->value(val) )
    return false;
  mCallback->set_seen();
  return true;
}

/************************* CLArgIDList ********************************/

std::string CLArgIDList::brief() const
{
  const char str[] = { '-', flag(), '\0' };
  std::string result( str );
  result += " <id_list>";
  return result;
}
std::string CLArgIDList::manstr() const
{
  return make_man_string( "id_list" );
}

bool CLArgIDList::parse( const char* str ) const
{
  std::vector<int> list;
  long v;
  int i, j;
  char* endptr;
  for (;;) {
    while (isspace(*str))
      ++str;
    if (!*str)
      break;
      
    v = strtol( str, &endptr, 10 );
    if (endptr == str || v <= 0)
      return false;
    
    i = (int)v;
    if ((long)i != v)
      return false;
    
    list.push_back(i);
    str = endptr;
    
    while (isspace(*str))
      ++str;
    
    if (*str == '-') {
      do { ++str; } while (isspace(*str));
      v = strtol( str, &endptr, 10 );
      if (endptr == str || v < i)
        return false;
      j = (int)v;
      if ((long)j != v)
        return false;
      for (++i; i < j; ++j)
        list.push_back(i);
      
      str = endptr;
      while (isspace(*str))
        ++str;
    } 
      
    if (!*str)
      break;
    if (*str != ',')
      return false;
    ++str;
  }
  
  if (list.empty())
    return false;
  
  if (!mCallback->value( list ))
    return false;
  mCallback->set_seen();
  return true;
}

/************************* CLArgListData ********************************/


bool CLArgListData::add_set( int size, const char* const* names )
{
  if (size < 0 || (!mSets.empty() && acceptable_length(size)))
    return false;
  
  mSets.resize( mSets.size() + 1 );
  for (int i = 0; i < size; ++i) 
    mSets.back().push_back( names[i] );
  return true;
}

bool CLArgListData::acceptable_length( unsigned len ) const
{
  for (unsigned i = 0; i < mSets.size(); ++i)
    if (len == mSets[i].size())
      return true;
  return mSets.empty();
}

std::string CLArgListData::set_string( int i ) const
{
  if (mSets[i].empty())
    return std::string();
    
  std::ostringstream result;
  result << mSets[i][0];
  for (unsigned j = 1; j < mSets[i].size(); ++j)
    result << "," << mSets[i][j];
  return result.str();
}

std::string CLArgListData::brief() const
{
  if (mSets.empty())
    return std::string();
  
  std::ostringstream result;
  for (unsigned i = 0; i < mSets.size(); ++i) {
    if (i)
      result << '|';
    result << '{' << set_string(i) << "}";
  }
  return result.str();
}

std::string CLArgListData::manstr( char type_char, const CLArgFlag& f ) const
{
  if (mSets.empty()) {
    const char argstr[] = { type_char, '1', ',', type_char, '2', ',', '.', '.', '.', '\0' };
    return f.make_man_string( argstr );
  }
  else {
    const char flagstr[] = { '-', f.flag(), '\0' };
    std::ostringstream ss;
    for (unsigned i = 0; i < mSets.size(); ++i) {
      if (i)
        ss << "|";
      ManPage::begin_bold( ss );
      ss << flagstr;
      ManPage::end_bold( ss );
      ss << "<" << set_string(i) << ">";
    }
    return ss.str();
  }
}


/************************* CLArgIntList ********************************/

std::string CLArgIntList::brief() const
{
  const char str[] = { '-', flag(), ' ', '\0' };
  std::string result( str );
  
  if (listData.accept_any_length()) 
    result += "<n1,n2,...>";
  else
    result += listData.brief();
  
  return result;
}

std::string CLArgIntList::manstr() const
  { return listData.manstr('n', *this); }

bool CLArgIntList::parse( const char* str ) const
{
  long i;
  char* endptr;
  std::vector<int> list;
  for (;;) {
    while (isspace(*str))
      ++str;
    
    if (!*str)
      break;
    
    i = strtol( str, &endptr, 0 );
    list.push_back(i);
    
    while (isspace(*str))
      ++str;

    if (!*str)
      break;
    
    if (*str != ',')
      return false;
    
    ++str;
  }
  
  if (!listData.acceptable_length( list.size() ) || !mCallback->value( list ))
    return false;
  
  mCallback->set_seen();
  return true;
}
    

/************************* CLArgDoubleList ********************************/

std::string CLArgDoubleList::brief() const
{
  const char str[] = { '-', flag(), ' ', '\0' };
  std::string result( str );
  
  if (listData.accept_any_length()) 
    result += "<r1,r2,...>";
  else
    result += listData.brief();
  
  return result;
}

std::string CLArgDoubleList::manstr() const
  { return listData.manstr('r', *this); }

bool CLArgDoubleList::parse( const char* str ) const
{
  std::vector<double> list;
  char* endptr;
  for (;;) {
    while (isspace(*str))
      ++str;
    if (!*str)
      break;

    double d = strtod( str, &endptr );
    list.push_back(d);
    str = endptr;

    while (isspace(*str))
      ++str;
    if (!*str)
      break;

    if (*str != ',')
      return false;
    ++str;
  }
  
  if (!listData.acceptable_length( list.size() ) || !mCallback->value( list ))
    return false;
  
  mCallback->set_seen();
  return true;
}

