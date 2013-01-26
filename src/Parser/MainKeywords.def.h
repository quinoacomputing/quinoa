//******************************************************************************
/*!
  \file      src/Parser/MainKeywords.def.h
  \author    J. Bakosi
  \date      Fri 25 Jan 2013 07:39:48 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Main keywords
  \details   Keywords recognized outside of any block
*/
//******************************************************************************

namespace keyword {

  using title = pegtl::string<t,i,t,l,e>;
  using end   = pegtl::string<e,n,d>;

  struct main :
         sor< title, end > {};

} // namespace keyword
