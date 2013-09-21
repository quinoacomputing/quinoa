//******************************************************************************
/*!
  \file      src/Control/BaseKeywords.h
  \author    J. Bakosi
  \date      Sat 21 Sep 2013 07:30:24 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Basic keywords recognized by all parsers
  \details   Basic keywords recognized by all parsers
*/
//******************************************************************************
#ifndef Keywords
#error "BaseKeywords.h should only be included within a *Keywords.h"
#endif

#ifndef BaseKeywords_h
#define BaseKeywords_h

//! Keyword 'title'
const char title_name[] = "Analysis title";
const char title_help[] =
  "The analysis title may be specified in the input file using the 'title' "
  "keyword. The 'title' keyword must be followed by a doubly-quoted string "
  "specifying the analysis title. Example: title \"Example problem\".";
using title = keyword<title_name, title_help, t,i,t,l,e >;

// Keyword 'end'
const char end_name[] = "End of block";
const char end_help[] =
  "The end of a block is given by the 'end' keyword. Example:\n"
  "\tstatistics\n"
  "\t  ...\n"
  "\tend";
using end = keyword<end_name, end_help, e,n,d >;

#endif // BaseKeywords_h
