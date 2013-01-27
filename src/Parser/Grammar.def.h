//******************************************************************************
/*!
  \file      src/Parser/Grammar.def.h
  \author    J. Bakosi
  \date      Sat 26 Jan 2013 08:45:39 PM MST
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

  // Grammar

  struct quoted :
         trim< not_one<'"'>, one<'"'> > {};

  struct read_title :
         ifmust< one<'"'>, ifapply< quoted, insert_title >, one<'"'>, space > {};
  
  struct read_hydro :
         pad< ifapply< trim<alnum, space>, insert_hydro >, blank, space > {};

  struct read_mix :
         pad< ifapply< trim<alnum, space>, insert_mix >, blank, space > {};

  struct process_hydro :
         ifmust< read<keyword::hydro>, read_hydro > {};

  struct process_mix :
         ifmust< read<keyword::mix>, read_mix > {};

  struct read_spinsflow :
         until< read<keyword::end>, sor< process_hydro, process_mix > > {};

  struct process_title :
         ifmust< read<keyword::title>, read_title > {};

  struct process_spinsflow :
         ifmust< read<keyword::spinsflow>, read_spinsflow > {};

  struct keyword_main :
         process_title {};

  struct keyword_physics :
         sor< read< keyword::homdir >,
              read< keyword::homgendir >,
              process_spinsflow > {};

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
