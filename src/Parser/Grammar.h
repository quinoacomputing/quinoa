//******************************************************************************
/*!
  \file      src/Parser/Grammar.h
  \author    J. Bakosi
  \date      Thu May 30 08:27:42 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Grammar definition
  \details   Grammar definition. We use the Parsing Expression Grammar Template
             Library (PEGTL) to create the grammar and the associated parser.
             Credit goes to Colin Hirsch (pegtl@cohi.at) for PEGTL. See
             src/ThirdParty/PEGTL/documentation.wiki for more details. Word of
             advice: read from the bottom up.
*/
//******************************************************************************
#ifndef Grammar_h
#define Grammar_h

#include <algorithm>

#include <Macro.h>
#include <ControlTypes.h>
#include <Option.h>
#include <PhysicsOptions.h>
#include <PositionOptions.h>
#include <MassOptions.h>
#include <HydroOptions.h>
#include <MixOptions.h>

namespace Quinoa {

namespace grammar {

  using namespace pegtl;
  using namespace pegtl::ascii;

  // Keywords

  #include <Keywords.h>

  // State

  //! Bundle is where everything is stored during parsing
  using Stack = control::Bundle;
  //! BoolBundle indicates stored fields (true or false)
  using BoolStack = control::BoolBundle;
  //! Out-of-struct storage of field ID for push_term
  static int field = 0;
  //! Physics options
  static control::Option<select::Physics> Physics;
  //! Position options
  static control::Option<select::Position> Position;
  //! Mass options
  static control::Option<select::Mass> Mass;
  //! Hydro options
  static control::Option<select::Hydro> Hydro;
  //! Mix options
  static control::Option<select::Mix> Mix;

  // Actions

  // convert string to 'type'
  template< class type >
  static type convert(const std::string& str) {
    std::stringstream ss(str);
    type num;
    ss >> num;
    return num;
  }

  // convert 'type' to string
  template< class type >
  static std::string convert(const type& val) {
    std::stringstream ss;
    ss << val;
    return ss.str();
  }

  // convert & store value in state 'at' position
  template< control::BundlePosition at >
  struct cstore : action_base< cstore<at> > {
    static void apply(const std::string& value,
                      Stack& stack,
                      BoolStack& boolstack) {
      // Convert to correct type and store element at position 'at'
      std::get<at>(stack) =
        convert< typename std::tuple_element<at,Stack>::type >(value);
      boolstack[at] = true;
    }
  };

  // store value in state 'at' position
  template< control::BundlePosition at >
  struct store : action_base< store<at> > {
    static void apply(const std::string& value,
                      Stack& stack,
                      BoolStack& boolstack) {
      std::get<at>(stack) = value;
      boolstack[at] = true;
    }
  };

  // push value in state 'at' position
  template< control::BundlePosition at >
  struct push : action_base< push<at> > {
    static void apply(const std::string& value,
                      Stack& stack,
                      BoolStack& boolstack) {
      // Convert to correct type and push element to vector at position 'at'
      std::get<at>(stack).push_back(convert<real>(value));
      boolstack[at] = true;
    }
  };

  // start new product in vector of statistics
  struct start_product : action_base< start_product > {
    static void apply(const std::string& value,
                      Stack& stack,
                      BoolStack& boolstack) {
      std::get<control::STATISTICS>(stack).push_back(std::vector<control::Term>());
      boolstack[control::STATISTICS] = true;
      IGNORE(value);        // suppress compiler warning on unused variable
    }
  };

  // add matched value as Term into vector of Product in vector of statistics
  template< control::Quantity quantity, control::Moment moment, char name='\0' >
  struct push_term : action_base< push_term<quantity, moment, name> > {
    static void apply(const std::string& value,
                      Stack& stack,
                      BoolStack& boolstack) {

      // if name is given, push name, otherwise push first char of value
      char n(name ? name : value[0]);
      // if name is given, it is triggered not user-requested
      bool plot(name ? false : true);
      // Use stats for shorthand of reference in bundle
      std::vector<control::Product>& stats = std::get<control::STATISTICS>(stack);
      // Push term into current product
      stats.back().push_back(control::Term(field, quantity, moment, n, plot));

      // If central moment, trigger mean
      if (moment == control::Moment::CENTRAL) {
        control::Term term(field,
                           quantity,
                           control::Moment::ORDINARY,
                           toupper(n),
                           false);
        stats.insert(stats.end()-1, control::Product(1,term));
      }

      field = 0;            // reset default field
      IGNORE(boolstack);    // suppress compiler warning on unused variable
    }
  };

  // save field ID so push_term can pick it up
  struct save_field : action_base< save_field > {
    static void apply(const std::string& value,
                      Stack& stack,
                      BoolStack& boolstack) {
      field = convert<int>(value) - 1;  // numbering of field IDs start from 0
      IGNORE(stack);        // suppress compiler warning on unused variable
      IGNORE(boolstack);    // suppress compiler warning on unused variable
    }
  };

  // store selected physics
  struct store_physics : action_base< store_physics > {
    static void apply(const std::string& value,
                      Stack& stack,
                      BoolStack& boolstack) {
      // Issue warning if overwrite
      if (boolstack[control::PHYSICS]) {
        std::cout << ">>> PARSER WARNING: Multiple physics defined in input file"
                  << std::endl << ">>> Overwriting \""
                  << Physics.name(std::get<control::PHYSICS>(stack))
                  << "\" with \""
                  << Physics.name(Physics.option(value)) << "\""
                  << std::endl;
      }
      std::get<control::PHYSICS>(stack) = Physics.option(value);
      boolstack[control::PHYSICS] = true;
    }
  };

  // store selected position model
  struct store_position : action_base< store_position > {
    static void apply(const std::string& value,
                      Stack& stack,
                      BoolStack& boolstack) {
      // Issue warning if overwrite
      if (boolstack[control::POSITION]) {
        std::cout << ">>> PARSER WARNING: Multiple position models defined in input"
                " file" << std::endl << ">>> Overwriting \""
             << Position.name(std::get<control::POSITION>(stack))
             << "\" with \""
             << Position.name(Position.option(value))
             << "\"" << std::endl;
      }
      std::get<control::POSITION>(stack) = Position.option(value);
      boolstack[control::POSITION] = true;
    }
  };

  // store selected mass model
  struct store_mass : action_base< store_mass > {
    static void apply(const std::string& value,
                      Stack& stack,
                      BoolStack& boolstack) {
      // Issue warning if overwrite
      if (boolstack[control::MASS]) {
        std::cout << ">>> PARSER WARNING: Multiple mass models defined in input file"
             << std::endl << ">>> Overwriting \""
             << Mass.name(std::get<control::MASS>(stack))
             << "\" with \""
             << Mass.name(Mass.option(value)) << "\""
             << std::endl;
      }
      std::get<control::MASS>(stack) = Mass.option(value);
      boolstack[control::MASS] = true;
    }
  };

  // store selected hydrodynamics model
  struct store_hydro : action_base< store_hydro > {
    static void apply(const std::string& value,
                      Stack& stack,
                      BoolStack& boolstack) {
      // Issue warning if overwrite
      if (boolstack[control::HYDRO]) {
        std::cout << ">>> PARSER WARNING: Multiple hydro models defined in input "
                "file" << std::endl << ">>> Overwriting \""
             << Hydro.name(std::get<control::HYDRO>(stack))
             << "\" with \""
             << Hydro.name(Hydro.option(value)) << "\""
             << std::endl;
      }
      std::get<control::HYDRO>(stack) = Hydro.option(value);
      boolstack[control::HYDRO] = true;
    }
  };

  // store selected material mix model
  struct store_mix : action_base< store_mix > {
    static void apply(const std::string& value,
                      Stack& stack,
                      BoolStack& boolstack) {
      // Issue warning if overwrite
      if (boolstack[control::MIX]) {
        std::cout << ">>> PARSER WARNING: Multiple mix models defined in input file"
             << std::endl << ">>> Overwriting \""
             << Mix.name(std::get<control::MIX>(stack))
             << "\" with \""
             << Mix.name(Mix.option(value)) << "\"" << std::endl;
      }
      std::get<control::MIX>(stack) = Mix.option(value);
      boolstack[control::MIX] = true;
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

  // match all accepted as position models
  struct position : sor< keyword::invpos,
                         keyword::vispos > {};

  // match all accepted as hydro models
  struct hydro : sor< keyword::slm,
                      keyword::glm > {};

  // match all accepted as mix models
  struct mix : sor< keyword::iem,
                    keyword::iecm,
                    keyword::dir,
                    keyword::gendir > {};

  // parse input padded by blank at left and space at right and if it matches
  // 'keywords', apply 'actions'
  template< class keywords, typename ... actions >
  struct parse :
         pad< ifapply< trim<keywords, space>, actions ... >, blank, space > {};

  // comment: start with '#' until eol
  struct comment :
         pad< trim<one<'#'>,eol>, blank, eol> {};

  // plow through block of 'tokens' until 'end' keyword
  template< typename ... tokens >
  struct block :
         until< read<keyword::end>, sor<comment, tokens ...> > {};

  // plow through list of values between keywords 'key' and "end", calling
  // 'insert' for each if matches and allowing comments between values
  template< class key, class insert, class value = digit >
  struct list :
         ifmust< read<key>,
                 until< read<keyword::end>, sor<comment, parse<value,insert>> >
               > {};

  // process 'keyword' and call its 'insert' action if matches 'keywords'
  template< class keyword, class insert, class keywords = alnum >
  struct process :
         ifmust< read<keyword>, parse<keywords,insert> > {};

  // moment: 'keyword' optionally followed by a digit, pushed to vector of terms
  template< class keyword, control::Quantity q, control::Moment m >
  struct moment :
         sor < ifapply<seq<keyword, ifapply<digit,save_field>>, push_term<q,m>>,
               ifapply<keyword, push_term<q,m>>
             > {};

  // terms recognized within an expectation and their mapping
  struct terms :
         sor< moment<keyword::transported_scalar,
                     control::Quantity::SCALAR,
                     control::Moment::ORDINARY>,
              moment<keyword::transported_scalar_fluctuation,
                     control::Quantity::SCALAR,
                     control::Moment::CENTRAL>,
              moment<keyword::velocity_x,
                     control::Quantity::VELOCITY_X,
                     control::Moment::ORDINARY>
            > {};

  // plow through terms in expectation until character 'rbound'
  template< char rbound >
  struct expectation :
         until< one<rbound>, terms > {};

  // plow through expectations between characters 'lbound' and 'rbound'
  template< char lbound, char rbound >
  struct parse_expectations :
         read< ifmust< one<lbound>, apply<start_product>, expectation<rbound> >
             > {};

  // scan title between characters 'lbound' and 'rbound'
  template< char lbound, char rbound >
  struct parse_title :
         ifmust< one<lbound>,
                 ifapply< trim<not_one<lbound>, one<rbound>>,
                          store<control::TITLE> >,
                 one<rbound>
               > {};

  // title
  struct title :
         ifmust< read<keyword::title>, parse_title<'"','"'> > {};

  // dir block
  struct dir :
         ifmust< parse<keyword::dir, store_mix>,
                 block< process<keyword::nscalar, cstore<control::NSCALAR>>,
                        list<keyword::dir_B, push<control::B>>,
                        list<keyword::dir_S, push<control::S>>,
                        list<keyword::dir_kappa, push<control::KAPPA>>
                      >
               > {};

  // gendir block
  struct gendir :
         ifmust< parse<keyword::gendir, store_mix>,
                 block< process<keyword::nscalar, cstore<control::NSCALAR>>,
                        list<keyword::dir_B, push<control::B>>,
                        list<keyword::dir_S, push<control::S>>,
                        list<keyword::dir_kappa, push<control::KAPPA>>,
                        list<keyword::gendir_C, push<control::C>>
                      >
               > {};

  // statistics block
  struct statistics :
         ifmust< read< keyword::statistics >,
                 block< parse_expectations<'<','>'> >
               > {};

  // slm block
  struct slm :
         ifmust< parse< keyword::slm, store_hydro,
                        // trigger estimating the diagonal of Reynolds-stress
                        start_product,
                        push_term<control::Quantity::VELOCITY_X,
                                  control::Moment::CENTRAL, 'u'>,
                        push_term<control::Quantity::VELOCITY_X,
                                  control::Moment::CENTRAL, 'u'>,
                        start_product,
                        push_term<control::Quantity::VELOCITY_Y,
                                  control::Moment::CENTRAL, 'v'>,
                        push_term<control::Quantity::VELOCITY_Y,
                                  control::Moment::CENTRAL, 'v'>,
                        start_product,
                        push_term<control::Quantity::VELOCITY_Z,
                                  control::Moment::CENTRAL, 'w'>,
                        push_term<control::Quantity::VELOCITY_Z,
                                  control::Moment::CENTRAL, 'w'>
                      >,
                 block< process<keyword::SLM_C0, cstore<control::C0>>,
                        process<keyword::nvelocity, cstore<control::NVELOCITY>>
                      >
               > {};

  // beta block
  struct beta :
         ifmust< parse< keyword::beta, store_mass >,
                 block< process<keyword::ndensity, cstore<control::NDENSITY>>,
                        process<keyword::Beta_At, cstore<control::AT>> >
               > {};

  // common to all physics
  struct physics_common :
         sor< process<keyword::nstep, cstore<control::NSTEP>>,
              process<keyword::term, cstore<control::TERM>>,
              process<keyword::dt, cstore<control::DT>>,
              process<keyword::npar, cstore<control::NPAR>>,
              process<keyword::pdfname, store<control::PDFNAME>>,
              process<keyword::globname, store<control::GLOBNAME>>,
              process<keyword::plotname, store<control::PLOTNAME>>,
              process<keyword::glob, cstore<control::GLOB>>,
              process<keyword::pdfi, cstore<control::PDFI>>,
              process<keyword::plti, cstore<control::PLTI>>,
              process<keyword::ttyi, cstore<control::TTYI>>,
              process<keyword::dump, cstore<control::DUMP>>
            > {};

  // hommix block
  struct hommix :
         ifmust< parse<keyword::hommix, store_physics>,
                 block< physics_common,
                        dir, gendir, statistics > > {};

  // homrt block
  struct homrt :
         ifmust< parse<keyword::homrt, store_physics>,
                 block< physics_common,
                        dir, gendir, slm, beta, statistics > > {};

  // homhydro block
  struct homhydro :
         ifmust< parse<keyword::homhydro, store_physics>,
                 block< physics_common,
                        slm, statistics > > {};

  // spinsflow block
  struct spinsflow :
         ifmust< parse<keyword::spinsflow, store_physics>,
                 block< physics_common,
                        slm, dir, gendir, beta > > {};

  // physics
  struct physics :
         sor< hommix,
              homhydro,
              homrt,
              spinsflow > {};

  // main keywords
  struct keywords :
         sor< title,
              physics > {};

  // ignore: comments and empty lines
  struct ignore :
         sor< comment, until<eol, space> > {};

  // entry point: parse keywords and ignores until eof
  struct read_file :
         until< eof, sor<keywords, ignore> > {};

} // namespace grammar

} // namespace Quinoa

#endif // Grammar_h
