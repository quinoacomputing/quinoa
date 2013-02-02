//******************************************************************************
/*!
  \file      src/Parser/Grammar.h
  \author    J. Bakosi
  \date      Sat 02 Feb 2013 06:55:09 AM MST
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

  using namespace pegtl;
  using namespace pegtl::ascii;
  using namespace Quinoa;

  // Keywords

  #include <Keywords.h>

  // State

  using stack_type = tuple< std::string,     //!< 0: title
                            PhysicsType,     //!< 1: physics
                            HydroType,       //!< 2: hydro
                            MixType,         //!< 3: mix
                            int,             //!< 4: nstep
                            real,            //!< 5: term
                            real,            //!< 6: dt
                            int,             //!< 7: nscalar
                            int,             //!< 8: npar
                            int >;           //!< 9: echo
  // Actions

//   struct do_comment : action_base< do_comment > {
//     static void apply(const std::string& m, stack_type& stack) {
//       cout << "COMMENT: \"" << m << "\"" << endl;
//     }
//   };
//
//   struct unknown : action_base< unknown > {
//     static void apply(const std::string& m, stack_type& stack) {
//       Throw(ParserException, FATAL, UNKNOWN_KEYWORD);
//       //cout << "UNKNOWN: \"" << m << "\"" << endl;
//     }
//   };

  // convert string to 'type'
  template< class type >
  static type convert(const std::string& str) {
    stringstream ss(str);
    type num;
    ss >> num;
    return num;
  }

  // store value in state at 'position'
  template< class type, size_t position >
  struct store : action_base< store<type, position> > {
    static void apply(const std::string& value, stack_type& stack) {
      get<position>(stack) = convert<type>(value);
    }
  };

  // store selected title
  struct store_title : action_base< store_title > {
    static void apply(const std::string& value, stack_type& stack) {
      get<0>(stack) = value;
    }
  };

  // store selected physics
  struct store_physics : action_base< store_physics > {
    static void apply(const std::string& value, stack_type& stack) {
      get<1>(stack) = associate::Physics[value];
    }
  };

  // store selected hydrodynamics model
  struct store_hydro : action_base< store_hydro > {
    static void apply(const std::string& value, stack_type& stack) {
      get<2>(stack) = associate::Hydro[value];
    }
  };

  // store selected material mix model
  struct store_mix : action_base< store_mix > {
    static void apply(const std::string& value, stack_type& stack) {
      get<3>(stack) = associate::Mix[value];
    }
  };

  // Grammar

  // read 'token' until 'erased' trimming 'erased'
  template< class token, class erased >
  struct trim :
         seq< token, until< at<erased> > > {};

  // read 'token' padded by blank at left and space at right
  template< class token >
  struct read :
         pad< trim<token, space>, blank, space > {};

  // match all accepted as hydro
  struct hydro : sor< keyword::slm,
                      keyword::glm > {};

  // match all accepted as mix
  struct mix : sor< keyword::iem,
                    keyword::iecm,
                    keyword::dir,
                    keyword::gendir > {};

  // parse input padded by blank at left and space at right and if it matches
  // 'keywords', apply 'action'
  template< class action, class keywords >
  struct parse :
         pad< ifapply< trim<keywords, space>, action >, blank, space > {};

  // comment: start with '#' until eol
  struct comment :
         //pad< ifapply< trim<one<'#'>,eol>, do_comment >, blank, eol > {};
         pad< trim<one<'#'>,eol>, blank, eol> {};

  // plow through block of 'tokens' until 'end' keyword
  template< typename ... tokens >
  struct block :
         until< read<keyword::end>, sor<comment, tokens ...> > {};

  // process 'keyword' and call its 'insert' action if matches 'keywords'
  template< class keyword, class insert, class keywords = alnum >
  struct process :
         ifmust< read<keyword>, parse<insert, keywords> > {};

  // title: within double quotes
  struct quoted :
         trim< not_one<'"'>, one<'"'> > {};

  struct parse_title :
         ifmust< one<'"'>,
                 ifapply< quoted, store_title >,
                 one<'"'>, space > {};

  struct title :
         ifmust< read<keyword::title>, parse_title > {};

  // spinsflow block
  struct spinsflow :
         ifmust< parse<store_physics, keyword::spinsflow>,
                 block< process<keyword::hydro, store_hydro, hydro>,
                        process<keyword::mix, store_mix, mix> > > {};

  // homdir block
  struct homdir :
         ifmust< parse<store_physics, keyword::homdir>,
                 block< process<keyword::nstep, store<int,4>>,
                        process<keyword::term, store<real,5>>,
                        process<keyword::dt, store<real,6>>,
                        process<keyword::nscalar, store<int,7>>,
                        process<keyword::npar, store<int,8>>,
                        process<keyword::echo, store<int,9>> > > {};

  // physics block
  struct physics :
         sor< homdir,
              read< keyword::homgendir >,
              spinsflow > {};

  // keywords
  struct keywords :
         sor< title,
              physics > {};

  // ignore comments and empty lines
  struct ignore :
         sor< comment, until<eol, space> > {};

  // Entry point: parse keywords and ignores until eof
  struct read_file :
         until< eof, sor<keywords, ignore> > {};

} // namespace grammar

#endif // Grammar_def_h
