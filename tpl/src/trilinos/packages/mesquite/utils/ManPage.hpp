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


/** \file ManPage.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_MAN_PAGE_HPP
#define MSQ_MAN_PAGE_HPP

#include "Mesquite.hpp"
#include <iostream>
#include <string>

class ManPage
{
public:
  static std::ostream& begin_bold( std::ostream& str )
    { return str << std::endl << ".B" << std::endl; }
  static std::ostream& end_bold( std::ostream& str )
    { return str << std::endl; }
  static std::ostream& bold( std::ostream& str, const std::string& s )
    { return end_bold( begin_bold(str) << s ); }

  static std::ostream& begin_italic( std::ostream& str )
    { return str << std::endl << ".I" << std::endl; }
  static std::ostream& end_italic( std::ostream& str )
    { return str << std::endl; }
  static std::ostream& italic( std::ostream& str, const std::string& s )
    { return end_italic( begin_italic(str) << s ); }
    
  static std::ostream& begin_section( std::ostream& str, const std::string& name )
    { return str << std::endl << ".SH " << name << std::endl; }
    
  static std::ostream& begin_subsection( std::ostream& str, const std::string& name )
    { return str << std::endl << ".SS " << name << std::endl; }
    
  static std::ostream& begin_paragraph( std::ostream& str )
    { return str << std::endl << ".P " << std::endl; }
    
  static std::ostream& begin_hanging_paragraph( std::ostream& str )
    { return str << std::endl << ".HP " << std::endl; }
    
  static std::ostream& begin_indent( std::ostream& str )
    { return str << std::endl << ".RS " << std::endl; }
  static std::ostream& end_indent( std::ostream& str )
    { return str << std::endl << ".RE " << std::endl; }
    
  static std::ostream& begin_manpage( std::ostream& str, const std::string& name, int section )
    { return str << std::endl << ".TH " << name << " " << section << std::endl; }
    
  static std::ostream& write_text( std::ostream& str, bool hanging_indent, const std::string& text );
};

#endif
