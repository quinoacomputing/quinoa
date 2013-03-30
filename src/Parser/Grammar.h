//******************************************************************************
/*!
  \file      src/Parser/Grammar.h
  \author    J. Bakosi
  \date      Sat 30 Mar 2013 01:28:35 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Grammar definition
  \details   Grammar definition
*/
//******************************************************************************
#ifndef Grammar_h
#define Grammar_h

#include <algorithm>

#include <Macro.h>
#include <ControlTypes.h>
#include <FwdAssociate.h>
#include <ParserException.h>

namespace Quinoa {

namespace grammar {

  using namespace pegtl;
  using namespace pegtl::ascii;

  // Keywords

  #include <Keywords.h>

  // State

  //! Bundles is where everything is stored during parsing
  using stack_type = control::Bundle;
  //! BoolBundle indicates stored fields (true or false)
  using boolstack_type = control::BoolBundle;
  //! Dummy Bundle instant for decltype in cstore()
  static stack_type dummy_stack;
  //! Field ID for statistics
  static int field;

  // Actions

  // convert string to 'type'
  template< class type >
  static type convert(const std::string& str) {
    stringstream ss(str);
    type num;
    ss >> num;
    return num;
  }

  // convert & store value in state 'at' position
  template< control::BundlePosition at >
  struct cstore : action_base< cstore<at> > {
    static void apply(const std::string& value,
                      stack_type& stack,
                      boolstack_type& boolstack) {
      // Figure out element type at position 'at'
      using type = typename std::tuple_element<at, decltype(dummy_stack)>::type;
      // Convert to correct type and store element at position 'at'
      get<at>(stack) = convert<type>(value);
      boolstack[at] = true;
    }
  };

  // store value in state 'at' position
  template< control::BundlePosition at >
  struct store : action_base< store<at> > {
    static void apply(const std::string& value,
                      stack_type& stack,
                      boolstack_type& boolstack) {
      get<at>(stack) = value;
      boolstack[at] = true;
    }
  };

  // push value in state 'at' position
  template< control::BundlePosition at >
  struct push : action_base< push<at> > {
    static void apply(const std::string& value,
                      stack_type& stack,
                      boolstack_type& boolstack) {
      // Convert to correct type and push element to vector at position 'at'
      get<at>(stack).push_back(convert<real>(value));
      boolstack[at] = true;
    }
  };

  // start new product in vector of statistics
  struct start_product : action_base< start_product > {
    static void apply(const std::string& value,
                      stack_type& stack,
                      boolstack_type& boolstack) {
      get<control::STATISTICS>(stack).push_back(vector<control::Term>());
      boolstack[control::STATISTICS] = true;
      IGNORE(value);        // suppress compiler warning on unused variable
    }
  };

  // add Term into vector of Product in vector of statistics
  template< control::Quantity quantity, control::Moment moment >
  struct push_term : action_base< push_term<quantity, moment> > {
    static void apply(const std::string& value,
                      stack_type& stack,
                      boolstack_type& boolstack) {
      control::Term term = {field, quantity, moment, value, true};
      vector<control::Product>& stats = get<control::STATISTICS>(stack);
      // Push term into current product
      stats.back().push_back(term);
      // If central moment, trigger mean
      if (moment == control::CENTRAL) {
        // Convert name to upper-case for human-readable name
        std::string upper(value);
        std::for_each(upper.begin(), upper.end(),
                      [](char& c){ c = static_cast<char>(std::toupper(c)); } );
        // Put in request for mean
        control::Term term = {field, quantity, control::ORDINARY, upper, false};
        control::Product mean(1,term);
        stats.insert(stats.end()-1, mean);
      }
      IGNORE(boolstack);    // suppress compiler warning on unused variable
    }
  };

  // save field ID so push_term can pick it up
  struct save_field : action_base< save_field > {
    static void apply(const std::string& value,
                      stack_type& stack,
                      boolstack_type& boolstack) {
      field = convert<int>(value) - 1;  // numbering of field IDs start from 0
      IGNORE(stack);        // suppress compiler warning on unused variable
      IGNORE(boolstack);    // suppress compiler warning on unused variable
    }
  };

  // store selected physics
  struct store_physics : action_base< store_physics > {
    static void apply(const std::string& value,
                      stack_type& stack,
                      boolstack_type& boolstack) {
      get<control::PHYSICS>(stack) = associate::PhysicsEnum[value];
      boolstack[control::PHYSICS] = true;
    }
  };

  // store selected hydrodynamics model
  struct store_hydro : action_base< store_hydro > {
    static void apply(const std::string& value,
                      stack_type& stack,
                      boolstack_type& boolstack) {
      get<control::HYDRO>(stack) = associate::HydroEnum[value];
      boolstack[control::HYDRO] = true;
    }
  };

  // store selected material mix model
  struct store_mix : action_base< store_mix > {
    static void apply(const std::string& value,
                      stack_type& stack,
                      boolstack_type& boolstack) {
      get<control::MIX>(stack) = associate::MixEnum[value];
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
                 until< read<keyword::end>, sor<comment, parse<insert,value>> >
               > {};

  // process 'keyword' and call its 'insert' action if matches 'keywords'
  template< class keyword, class insert, class keywords = alnum >
  struct process :
         ifmust< read<keyword>, parse<insert, keywords> > {};

  // moment: 'keyword' followed by a digit, pushed to vector of terms
  template< class keyword, control::Quantity q, control::Moment m >
  struct moment :
         ifapply< seq<keyword, ifapply<digit,save_field>>, push_term<q, m> > {};

  // terms recognized within an expectation and their mapping
  struct terms :
         sor< moment<keyword::transported_scalar,
                     control::TRANSPORTED_SCALAR,
                     control::ORDINARY>,
              moment<keyword::transported_scalar_fluctuation,
                     control::TRANSPORTED_SCALAR,
                     control::CENTRAL>
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
         ifmust< parse<store_mix, keyword::dir>,
                 block< process<keyword::nscalar, cstore<control::NSCALAR>>,
                        list<keyword::B, push<control::B>>,
                        list<keyword::S, push<control::S>>,
                        list<keyword::kappa, push<control::KAPPA>>
                      >
               > {};

  // gendir block
  struct gendir :
         ifmust< parse<store_mix, keyword::gendir>,
                 block< process<keyword::nscalar, cstore<control::NSCALAR>>,
                        list<keyword::B, push<control::B>>,
                        list<keyword::S, push<control::S>>,
                        list<keyword::kappa, push<control::KAPPA>>,
                        list<keyword::C, push<control::C>>
                      >
               > {};

  // statistics block
  struct statistics :
         ifmust< read< keyword::statistics >,
                 block< parse_expectations<'<','>'> >
               > {};

  // hommix block
  struct hommix :
         ifmust< parse<store_physics, keyword::hommix>,
                 block< process<keyword::nstep, cstore<control::NSTEP>>,
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
                        process<keyword::dump, cstore<control::DUMP>>,
                        dir, gendir, statistics
                      >
               > {};

  // homhydro block
  struct homhydro :
         ifmust< parse<store_physics, keyword::homhydro>,
                 block< process<keyword::nstep, cstore<control::NSTEP>>,
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
                        process<keyword::dump, cstore<control::DUMP>>,
                        process<keyword::hydro, store_hydro, hydro>,
                        statistics
                      >
               > {};

  // spinsflow block
  struct spinsflow :
         ifmust< parse<store_physics, keyword::spinsflow>,
                 block< process<keyword::hydro, store_hydro, hydro>,
                        process<keyword::mix, store_mix, mix> > > {};

  // physics
  struct physics :
         sor< hommix,
              homhydro,
              spinsflow > {};

  // keywords
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
