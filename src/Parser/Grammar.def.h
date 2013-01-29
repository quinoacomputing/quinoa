//******************************************************************************
/*!
  \file      src/Parser/Grammar.def.h
  \author    J. Bakosi
  \date      Mon 28 Jan 2013 10:14:19 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Grammar definition
  \details   Grammar definition
*/
//******************************************************************************
#ifndef Grammar_def_h
#define Grammar_def_h

#include <unordered_map>

#include <ParserException.h>

namespace grammar {

  using namespace std;
  using namespace pegtl;
  using namespace pegtl::ascii;
  using namespace Quinoa;

  // Keywords

  #include <Keywords.def.h>

  // State

  enum key_type { TITLE=0,
                  HYDRO,
                  MIX
  };

  using value_type = std::string;
  using stack_type = std::unordered_map< key_type, value_type, hash<int> >;

  // Actions

//   struct do_comment : action_base< do_comment > {
//     static void apply(const std::string& m, stack_type& stack) {
//       std::cout << "COMMENT: \"" << m << "\"" << endl;
//     }
//   };
//
//   struct unknown : action_base< unknown > {
//     static void apply(const std::string& m, stack_type& stack) {
//       Throw(ParserException, FATAL, UNKNOWN_KEYWORD);
//       //std::cout << "UNKNOWN: \"" << m << "\"" << endl;
//     }
//   };

  // insert value to stack at 'key'
  template< key_type key >
  struct insert : action_base< insert<key> > {
    static void apply(const std::string& value, stack_type& stack) {
      stack[key] = value;
    }
  };

  // Utilities

  // read 'token' until 'erased' trimming 'erased'
  template< class token, class erased >
  struct trim :
         seq< token, until< at<erased> > > {};

  // read 'token' padded by blank at left and space at right
  template< class token >
  struct read :
         pad< trim<token, space>, blank, space > {};

  // parse input padded by blank at left and space at right and apply 'action'
  template< class action >
  struct parse :
         pad< ifapply< trim<alnum, space>, action >, blank, space > {};

  // parse block of 'tokens' until 'end' keyword
  template< typename ... tokens >
  struct block :
         until< read<keyword::end>, sor<tokens ...> > {};

  // read 'keyword' and call its 'insert' action
  template< class keyword, class insert >
  struct process :
         ifmust< read<keyword>, parse<insert> > {};

  // Grammar

  // title: within double quotes
  struct quoted :
         trim< not_one<'"'>, one<'"'> > {};

  struct parse_title :
         ifmust< one<'"'>, ifapply<quoted, insert<TITLE>>, one<'"'>, space > {};

  struct process_title :
         ifmust< read<keyword::title>, parse_title > {};

  // spinsflow block
  struct spinsflow :
         ifmust< read<keyword::spinsflow>,
                 block< process<keyword::hydro, insert<HYDRO>>,
                        process<keyword::mix, insert<MIX>> > > {};

  // physics keywords
  struct physics :
         sor< read< keyword::homdir >,
              read< keyword::homgendir >,
              spinsflow > {};

  // keywords
  struct keywords :
         sor< process_title,
              physics > {};

  // comment: start with '#' until eol
  struct comment :
         //pad< ifapply< trim<one<'#'>,eol>, do_comment >, blank, eol > {};
         pad< trim<one<'#'>,eol>, blank, eol> {};

  // ignore comments and empty lines
  struct ignore :
         sor< comment, until<eol, space> > {};

  // parser entry point: parse keywords and ignores until eof
  struct read_file :
         until< eof, sor<keywords, ignore> > {};

} // namespace grammar

#endif // Grammar_def_h
