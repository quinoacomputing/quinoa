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


/** \file ManPage.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "ManPage.hpp"


std::ostream& ManPage::write_text( std::ostream& str, bool hanging_indent, const std::string& text )
{
  std::string::size_type i = 0, j;
  if (hanging_indent)
    begin_hanging_paragraph( str );
  else
    begin_paragraph( str );
  const std::string::size_type n = text.size();
  while (i != n) {
    if (text[i] == '\n') {
      begin_paragraph( str );
      ++i;
      continue;
    }
    j = text.find( "\n", i );
    if (j == std::string::npos)
      j = n;
    if (text[i] == '.')
      str << '\\';
    str << text.substr( i, j - i );
    i = j;
  }
  return str;
}
