//******************************************************************************
/*!
  \file      src/Control/BaseKeywords.h
  \author    J. Bakosi
  \date      Thu Aug 29 14:43:47 2013
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

// Title
using title = pegtl::string< t,i,t,l,e >;

// End of block
using end = pegtl::string< e,n,d >;

#endif // BaseKeywords_h
