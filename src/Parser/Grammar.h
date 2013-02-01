//******************************************************************************
/*!
  \file      src/Parser/Grammar.h
  \author    J. Bakosi
  \date      Fri 01 Feb 2013 06:06:55 AM MST
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

  struct stack_type {
    std::string title;
    PhysicsType physics;
    HydroType hydro;
    MixType mix;
    int nstep;
    real term;
    real dt;
    int nscalar;
    int npar;
    int echo;
  };
 
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

  // store selected title
  struct store_title : action_base< store_title > {
    static void apply(const std::string& value, stack_type& stack) {
      stack.title = value;
    }
  };

  // store selected physics
  struct store_physics : action_base< store_physics > {
    static void apply(const std::string& value, stack_type& stack) {
      stack.physics = associate::Physics[value];
    }
  };

  // store selected hydrodynamics model
  struct store_hydro : action_base< store_hydro > {
    static void apply(const std::string& value, stack_type& stack) {
      stack.hydro = associate::Hydro[value];
    }
  };

  // store selected material mix model
  struct store_mix : action_base< store_mix > {
    static void apply(const std::string& value, stack_type& stack) {
      stack.mix = associate::Mix[value];
    }
  };

  // store selected number of time steps
  struct store_nstep : action_base< store_nstep > {
    static void apply(const std::string& value, stack_type& stack) {
      stack.nstep = 0;
    }
  };

  // store selected terminate time
  struct store_term : action_base< store_term > {
    static void apply(const std::string& value, stack_type& stack) {
      stack.term = 0.0;
    }
  };

  // store selected time step size
  struct store_dt : action_base< store_dt > {
    static void apply(const std::string& value, stack_type& stack) {
      stack.dt = 0.0;
    }
  };

  // store selected number of mixing scalars
  struct store_nscalar : action_base< store_nscalar > {
    static void apply(const std::string& value, stack_type& stack) {
      stack.nscalar = 0;
    }
  };

  // store selected number of particles
  struct store_npar : action_base< store_npar > {
    static void apply(const std::string& value, stack_type& stack) {
      stack.npar = 0.0;
    }
  };

  // store selected one-liner info frequency
  struct store_echo : action_base< store_echo > {
    static void apply(const std::string& value, stack_type& stack) {
      stack.echo = 0;
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
         ifmust< one<'"'>, ifapply<quoted, store_title>, one<'"'>, space > {};

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
                 block< process<keyword::nstep, store_nstep>,
                        process<keyword::term, store_term>,
                        process<keyword::nscalar, store_nscalar>,
                        process<keyword::npar, store_npar>,
                        process<keyword::echo, store_echo> > > {};

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
