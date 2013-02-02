//******************************************************************************
/*!
  \file      src/Parser/Grammar.h
  \author    J. Bakosi
  \date      Sat 02 Feb 2013 12:43:24 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Grammar definition
  \details   Grammar definition
*/
//******************************************************************************
#ifndef Grammar_h
#define Grammar_h

#include <Type.h>
#include <FwdAssociate.h>
#include <ParserException.h>

namespace Quinoa {

namespace grammar {

  using namespace pegtl;
  using namespace pegtl::ascii;

  // Keywords

  #include <Keywords.h>

  // State

  using stack_type = control::Bundle;
  stack_type dummy_stack;       // dummy instant for decltype in store()

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
  template< control::BundlePosition at >
  struct store : action_base< store<at> > {
    static void apply(const std::string& value, stack_type& stack) {
      using type = typename std::tuple_element<at, decltype(dummy_stack)>::type;
      get<at>(stack) = convert<type>(value);
    }
  };

  // store selected title
  struct store_title : action_base< store_title > {
    static void apply(const std::string& value, stack_type& stack) {
      get<control::TITLE>(stack) = value;
    }
  };

  // store selected physics
  struct store_physics : action_base< store_physics > {
    static void apply(const std::string& value, stack_type& stack) {
      get<control::PHYSICS>(stack) = associate::PhysicsEnum[value];
    }
  };

  // store selected hydrodynamics model
  struct store_hydro : action_base< store_hydro > {
    static void apply(const std::string& value, stack_type& stack) {
      get<control::HYDRO>(stack) = associate::HydroEnum[value];
    }
  };

  // store selected material mix model
  struct store_mix : action_base< store_mix > {
    static void apply(const std::string& value, stack_type& stack) {
      get<control::MIX>(stack) = associate::MixEnum[value];
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
                 block< process<keyword::nstep, store<control::NSTEP>>,
                        process<keyword::term, store<control::TERM>>,
                        process<keyword::dt, store<control::DT>>,
                        process<keyword::nscalar, store<control::NSCALAR>>,
                        process<keyword::npar, store<control::NPAR>>,
                        process<keyword::echo, store<control::ECHO>> > > {};

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

} // namespace Quinoa

#endif // Grammar_h
