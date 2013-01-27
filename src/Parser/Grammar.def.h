//******************************************************************************
/*!
  \file      src/Parser/Grammar.def.h
  \author    J. Bakosi
  \date      Sat 26 Jan 2013 11:14:44 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Grammar definition
  \details   Grammar definition
*/
//******************************************************************************
#ifndef Grammar_def_h
#define Grammar_def_h

#include <unordered_map>

namespace grammar {

  using namespace std;
  using namespace pegtl;
  using namespace pegtl::ascii;

  // Keywords

  #include <Keywords.def.h>

  // State

  using key_type = std::string;
  using value_type = std::string;
  using stack_type = std::unordered_map< key_type, value_type >;

  // Actions

  struct do_comment : action_base< do_comment > {
    static void apply(const std::string& m, stack_type& stack) {
      std::cout << "COMMENT: \"" << m << "\"" << endl;
    }
  };

  struct do_input : action_base< do_input > {
    static void apply(const std::string& m, stack_type& stack) {
      std::cout << "INPUT  : \"" << m << "\"" << endl;
    }
  };

  struct insert_title : action_base< insert_title > {
    static void apply(const std::string& value, stack_type& stack) {
      stack["title"] = value;
    }
  };

  struct insert_hydro : action_base< insert_hydro > {
    static void apply(const std::string& value, stack_type& stack) {
      stack["hydro"] = value;
    }
  };

  struct insert_mix : action_base< insert_mix > {
    static void apply(const std::string& value, stack_type& stack) {
      stack["mix"] = value;
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

  // Grammar

  // title
  struct quoted :
         trim< not_one<'"'>, one<'"'> > {};

  struct parse_title :
         ifmust< one<'"'>, ifapply< quoted, insert_title >, one<'"'>, space > {};

  struct process_title :
         ifmust< read<keyword::title>, parse_title > {};

  // spinsflow
  struct process_hydro :
         ifmust< read<keyword::hydro>, parse<insert_hydro> > {};

  struct process_mix :
         ifmust< read<keyword::mix>, parse<insert_mix> > {};

  struct read_spinsflow :
         until< read<keyword::end>, sor< process_hydro, process_mix > > {};

  struct process_spinsflow :
         ifmust< read<keyword::spinsflow>, read_spinsflow > {};

  // physics
  struct physics :
         sor< read< keyword::homdir >,
              read< keyword::homgendir >,
              process_spinsflow > {};

  struct keywords :
         sor< process_title,
              physics > {};

  struct comment :
         //pad< ifapply< trim<one<'#'>,eol>, do_comment >, blank, eol > {};
         pad< trim<one<'#'>,eol>, blank, eol> {};

  struct ignore :
         sor< comment, eol > {};

  struct read_file :
         until< eof, sor<keywords, parse<do_input>, ignore> > {};

} // namespace grammar

#endif // Grammar_def_h
