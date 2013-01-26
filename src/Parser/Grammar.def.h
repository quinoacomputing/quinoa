//******************************************************************************
/*!
  \file      src/Parser/Grammar.def.h
  \author    J. Bakosi
  \date      Sat 26 Jan 2013 09:48:04 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Grammar definition
  \details   Grammar definition
*/
//******************************************************************************
#ifndef Grammar_def_h
#define Grammar_def_h

namespace grammar {

  using namespace pegtl;
  using namespace pegtl::ascii;

  // Keywords

  #include <Keywords.def.h>

  // State

  // Actions

  struct do_comment : action_base< do_comment > {
    static void apply(const std::string& m) {
      std::cout << "COMMENT: \"" << m << "\"" << endl;
    }
  };

  struct do_input : action_base< do_input > {
    static void apply(const std::string& m) {
      std::cout << "INPUT  : \"" << m << "\"" << endl;
    }
  };

  struct parse_title : action_base< parse_title > {
    static void apply(const std::string& title) {
      //m_control->setTitle( title );
      std::cout << "TITLE  : \"" << title << "\"" << endl;
    }
  };

  // Utilities

  template< class what, class erased >
  struct trim :
         seq< what, until< at<erased> > > {};

  template< class what >
  struct read :
         pad< trim<what, space>, blank, space > {};

  // Grammar

  struct quoted :
         trim< not_one<'"'>, one<'"'> > {};

  struct read_title :
         ifmust< one<'"'>, ifapply< quoted, parse_title >, one<'"'>, space > {};

  struct process_title :
         ifmust< read<keyword::title>, read_title > {};

  struct keyword_main :
         sor< process_title,
              read< keyword::end > > {};

  struct keyword_physics :
         sor< read< keyword::HomogeneousDirichlet >,
              read< keyword::HomogeneousGeneralizedDirichlet >,
              read< keyword::SPINSFlow > > {};

  struct keyword_any :
         sor< keyword_main,
              keyword_physics > {};

  struct input :
         pad< ifapply< trim<alnum, space>, do_input >, blank, space > {};

  struct comment :
         //pad< ifapply< trim<one<'#'>,eol>, do_comment >, blank, eol > {};
         pad< trim<one<'#'>,eol>, blank, eol> {};

  struct ignore :
         sor< comment, eol > {};

  struct read_file :
         until< eof, sor<keyword_any, input, ignore> > {};

} // namespace grammar

#endif // Grammar_def_h
