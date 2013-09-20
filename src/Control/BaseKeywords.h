//******************************************************************************
/*!
  \file      src/Control/BaseKeywords.h
  \author    J. Bakosi
  \date      Fri Sep 20 12:21:06 2013
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
using title = keyword< t,i,t,l,e >;

// End of block
using end = keyword< e,n,d >;

#endif // BaseKeywords_h
