// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef OPTIONS_FROM_STREAM_H
#define OPTIONS_FROM_STREAM_H

#include "Moocho_ConfigDefs.hpp"

// /////////////////////////////////////////////////
// Forward declarations.

namespace OptionsFromStreamPack {
  namespace OptionsFromStreamUtilityPack {

    // Simple class for a boolean variable that is false on construction.
    class false_bool_t {
    public:
      false_bool_t() : val_(false) {}
      void set(bool val) { val_ = val; }
      operator bool() const { return val_; }
    private:
      bool val_;
    };

    // Implementation type for the map for looking up a value given
    // a an options name.
    typedef std::map< std::string , std::string >			option_to_value_map_t;

    // Implementation type for an options group's options and
    // a boolean variable to determine if this option has been visited yet.
    typedef std::pair< option_to_value_map_t, false_bool_t >	options_group_pair_t;

    // Implementation type of the map for looking up a set of options
    // given the option groups name.
    typedef std::map< std::string, options_group_pair_t >	options_group_map_t;

    // The above declarations sets up the data structure:
    // 
    //       map<string,options_group_pair_t> (options_group_map_t)
    //                  |
    //                  | value_type
    //                 `.'
    //       pair<string,options_group_pair_t> 
    //                  |
    //       -----------------------              
    //      |                       |         
    //      | first                 | second
    //     `.'                     `.'
    //    string         pair<option_to_value_map_t,false_bool_t> (options_group_pair_t)
    //      |                      |
    // "solver_options"            |
    //                             |
    //        ------------------------------------------------
    //       |                                                |
    //       | first                                          | second
    //      `.'                                              `.'
    //    map<string,string> (option_to_value_map_t)       false_bool_t
    //       |                                                |
    //       | value_type                                   "true"
    //      `.'
    // pair<string,string>
    //       |
    //     --------------
    //    | first        | second
    //   `.'            `.'
    //  string         string
    //    |              |    
    //  "tol"          "1e-6"
    //  

    class OptionsGroup;
  }
  class OptionsFromStream;
}

// //////////////////////////////////////////////////
// Class and function declarations

namespace OptionsFromStreamPack {

/* \defgroup OptionsFromStreamPack_grp Untilities for maintianing a simple options database.
 * \ingroup Misc_grp
 * 
 */
// @ {

namespace OptionsFromStreamUtilityPack {

/** \brief Class used to encapsulate options belonging an options group.
  *
  * This class offers two ways to access the options and their values in
  * an options group.  These two methods differ in their convience,
  * the cost of performing them and the level of validation possible.
  * In the following discussion, n_total
  * is the total number of options being checked for and n_set is the
  * number of options set.
  *
  * 1) Looking up the values directly.
  *
  * Here the software client just looks up the value directly.  If the option with the given
  * name does not exist then the function <tt>option_exists( ... )</tt> will return <tt>false</tt>.
  *
  * In the following example, the lookup of the option named "tol" is attempted
  * on the <tt>OptionsGroup</tt> object <tt>optgrp</tt>.
  \verbatim

  const std::string& val = optgrp.option_value( "tol" );
  if( OptionsGroup::option_exists( val ) )
      std::cout << "\ntol = " << val;

  \endverbatim
  * The total cost of this way of looking up option values is:
  *
  * O(n_total*log(n_set)).
  *
  * The disadvantage of using this method is that if the user misspelled the
  * name of the option then the option will not be found or set.  In the
  * above example, if the user misspelled the option "tol" as "tal" then
  * <tt>OptionsGroup::option_exists( val )</tt> would return false and the option
  * value would not be accessed.
  *
  * 2) Setting options by iterating through options present.
  *
  * In this method the software client must keep a map to translate from an options
  * name to an id representing the option.  Then a switch can be used to take
  * action.
  *
  * In the following example (see the example from \Ref{OptionsFromStream}
  * , three possible options are looked for.
  * Here the client has defined an enum for the options and initialized
  * a \Ref{StringToIntMap}.
  *
  \verbatim

  const char optgrp_name[] = "MySolverOptions";
  const int num_opt = 3;
  enum EOptions {
    TOL
    ,MAX_ITER
    ,PROB_TYPE
  };
  const char* SOptions[num_opt] = {
    "tol"
    ,"max_iter"
    ,"prob_type"
  };
  StringToIntMap	opt_map( optgrp_name, num_opt, SOptions );
  OptionsGroup::const_iterator
    itr     = optgrp.begin(),
    itr_end = optgrp.end();
  for( ; itr != itr_end; ++itr ) {
    switch( (EOptions)opt_map(option_name(itr)) ) {
      case TOL:
        std::cout << "\noption tol = " << option_value(itr);
        break;
      case MAX_ITER:
        std::cout << "\noption max_iter = " << option_value(itr);
        break;
      case PROB_TYPE:
        std::cout << "\noption prob_type = " << option_value(itr);
        break;
      default:
        std::cout << "\nThe option " << option_name(itr) << " is not valid";
        exit(-1);	            
    }
  }

  \endverbatim
  *
  * The above code terminates the program in case the option with
  * the name returned from <tt>option_name(itr)</tt> is not expected.  This option
  * could be ignored (as would have been in method 1) or some other action
  * could be taken.  It is a good idea to use the above method to validate
  * that all of the options that the user sets for an options_group are
  * known by the client software so that misspellings are caught.
  * For example, if the user misspelled "tol" as "tal" then the
  * "opt_map" object above would not find this option name an the
  * <tt>default:</tt> case above would be executed.
  *
  * The total cost of this way off looking up option values is:
  *
  * O(n_set*log(n_total))
  *
  * As it can be clearly seen, method 2 will be faster since n_total >= n_set.
  * Method 1 takes less overhead for the client (no map or enum needed)
  * and is easier to maintain.  The real advantage of
  * Method 2 though is that if a user misspells an option then Method 1
  * will not catch this while Method 2 will.  For this reason alone Method
  * 2 is to be preferred.
  */
class OptionsGroup {
public:

  // friends
  friend class OptionsFromStream;

  /** @name Public Types */
  //@{

  /** \brief . */
  typedef option_to_value_map_t::iterator			iterator;
  /** \brief . */
  typedef option_to_value_map_t::const_iterator	const_iterator;

  //@}

  /** @name Constructors.
    *
    * Default constructor and assignment operator are not publicly allowed.
    * Default copy constructor is allowed.
    */
  //@{

  //@}

  /** @name Lookup an an value of an option given the option's name.
    *
    * If the option does not exist then <tt>option_exists( option_value( option_name ) )</tt>
    * will return <tt>false</tt> otherwise the option exists and can be read.
    */
  //@{

  /** \brief . */
  std::string& option_value( const std::string& option_name );
  /** \brief . */
  const std::string& option_value( const std::string& option_name ) const;

  //@}

  /** \brief . */
  static bool option_exists( const std::string& option_value );

  /// Returns true if this options groups exists.
  bool options_group_exists() const;

  /** @name Iterator access */
  //@{

  /** \brief . */
  int				num_options() const;
  /** \brief . */
  iterator		begin();
  /** \brief . */
  iterator		end();
  /** \brief . */
  const_iterator	begin() const;
  /** \brief . */
  const_iterator	end() const;

  //@}

private:
  
  option_to_value_map_t*	option_to_value_map_;	// 0 used for no options group.
  static std::string		option_does_not_exist_;

  // Not defined and not to be called
  OptionsGroup();
  OptionsGroup& operator=(const OptionsGroup&);

public:
  // Is given zero by OptionsFromStream to signify the option group does't exist.
  // It is put here to keep it away from the eyes of the general user.
  OptionsGroup( option_to_value_map_t* option_to_value_map )
    : option_to_value_map_(option_to_value_map)
  {}

};	// end class OptionsGroup

inline
/** \brief . */
const std::string& option_name( OptionsGroup::const_iterator& itr ) {
  return (*itr).first;
}

inline
/** \brief . */
const std::string& option_value( OptionsGroup::const_iterator& itr ) {
  return (*itr).second;
}


// @ }

}	// end namespace OptionsFromStreamUtilityPack 


/** \brief Extracts options from a text stream and then allows
  * convenient access to them.
  *
  * The basic idea is that options are read in from a stream (which can be
  * a file, C++ string etc.) and then parsed and stored in a format so
  * that options can be efficiently looked up by client software.
  *
  * The syntax for the file (or any C++ istream) is as follows:
  \verbatim

  begin_options

  *** These are my solver options
  options_group MySolverOptions {
      tol       = 1e-5; *** Convergence tolerance
      max_iter  = 100;  *** Maximum number of iterations
      prob_type = LINEAR;
      *prob_type = NON_LINEAR;  *** Comment this line out
  }

  *** Options for another solver
  options_group YourSolverOptions {
      tol       = 1e-4;
      *** These options determine the type of problem solved
  *    type_prob = LP;
      type_prob = QP;
  }

  *** Reset the tolerance
  options_group MySolverOptions {
      tol = 1e-8; *** Reset to a tighter tolerance
  }

  end_options

  \endverbatim

  * The text stream will be read up to the <tt>end_options</tt> line.
  * Options groups will be read starting with the <tt>begin_options</tt> line.
  * Options groups can not be nested.  The names for the option groups
  * or the option names themselves can not contain any white space.  The
  * text for the option values however can contain white space.  The <tt>=</tt> must
  * separate each option from its value.  The value for an option begins
  * with the first non-whitespace character after the <tt>=</tt> and ends with the
  * last non-whitespace character before the <tt>;</tt>.  For the option and value
  * pair "tol = + 1e-5 ;" the option value would be
  * stored as "+ 1e-5".  Comment lines starting with
  * <tt>*</tt> can be placed anywhere in the stream and will be ignored.
  * Also comments starting with <tt>*</tt> after the <tt>;</tt> for an option and value
  * pair can occur on the same line as shown above.
  *
  * <b>Warning!</b> Do not use the char '}' in any comment within an options
  * group!  This will break the parser.
  *
  * The options groups are also reentrant which means that they may be included
  * more than once as shown above.  Therefore options may be set and reset in the
  * same or another declaration of the options groups.
  * In the above example, the second declaration of the options group declaration
  * for <tt>OptionsGroup1</tt> resets the value of <tt>OptionsGroup1::option1</tt> from
  * <tt>value1</tt> to <tt>another_value</tt>.  This feature provides much more flexibility
  * in letting options be changed easily without having to edit the text stream
  * before it is read in.
  * 
  * Now, what if the user misspells an options group name?  One strategy is
  * to require that the options group exist and then for the client to throw
  * exceptions if the options group does not exist.  However, where will be
  * occasions where and options group may be "optional" and you don't want the user
  * to have to specify every options group.  Therefore we don't want to make all
  * of the options groups mandatory.  However, what if a user thinks he/she is
  * setting options in an "optional" options group but misspells the name of the
  * options group?  In this case the options group would never be read by the
  * client software and the user may be perplexed as to why nothing changed.
  * To help deal with this problem, after all of the option groups have been
  * accessed, the client can call the function <tt>print_unaccessed_options_groups(...)</tt>
  * to print the list of options groups that have not been accessed.  This way
  * the user can see if they may have not spelled and "optional" options group
  * correctly.
  *
  * Careful use of this simple class for specifying and setting options has the
  * following advantages/features:
  * <ul>
  * <li> Options are read from a stream, not a file, and can therefore represent
  *     any kind of source (C++ strings, files, console etc.).
  * <li> Options are inclosed in namespaces (called options_group) and therefore
  *     can handle options from many different sources.
  * <li> An OptionsFromStream object can be passed around in a software application
  *     and interested sub-objects can search for 'their' options group and can pull
  *     out options that apply to them.
  * <li> Reasonable mechanisms are available to validate that users have spelled the
  *     names of options groups and their options names or the allowable value of
  *     the option values properly.
  * <li> The implementation is very light weight and uses very little code.
  *     Most of the major functionality comes from the ISO standard C++ library
  *     and therefore should be portable to any up to date C++ compiler.
  *</ul>
  */
class OptionsFromStream {
public:

  /** @name Public Types */
  //@{

  /// const iterator through options group access options
  typedef OptionsFromStreamUtilityPack::options_group_map_t::iterator			iterator;
  /// non-const iterator through options group access options
  typedef OptionsFromStreamUtilityPack::options_group_map_t::const_iterator	const_iterator;

  /// \Ref{OptionsGroup} typedef
  typedef	OptionsFromStreamUtilityPack::OptionsGroup							options_group_t;

  /// Thrown if there is an input error
  class InputStreamError : public std::logic_error
  {public: InputStreamError(const std::string& what_arg) : std::logic_error(what_arg) {}};

  //@}

  /** @name Constructors / Initializes.
    *
    * The default constructor, copy constructor and assignment operator functions
    * are allowed.
    */
  //@{

  /// Construct with no options set.
  OptionsFromStream();

  /** \brief Construct initialized from a text stream.
    *
    * This is equivalent to calling the default constructor and then calling
    * <tt>read_options(in)</tt>.
    */
  explicit OptionsFromStream( std::istream& in );

  /// Clear all the options
  void clear_options();

  /** \brief Add / modify options read in from a text stream.
    *
    * The format of the text stream is described in the introduction.
    *
    * The options read in from <tt>in</tt> will either be added anew or will
    * overwrite options already present.
    *
    * If the format of the stream is not correct then a
    * <tt>InputStreamError</tt> exception will be thrown.
    *
    * <b>Warning!</b> Do not use the char '}' in any comment within an options
    * group!  This will break the parser.
    */
  void read_options( std::istream& in );

  //@}

    /** \brief Print the options to an output stream.
    *
      * This is useful for debugging and also to record exactly what options have been set.
    */
  void print_options( std::ostream& out ) const;

  /** @name Get an options group access object given its name.
    *
    * If the option group does not exist then\\
    * <tt>options_group_exists( this->options_group( options_group_name ) ) == false</tt>\\
    * where <tt>options_group_name</tt> is the string name of the option group. 
    */
  //@{

  /** \brief . */
  options_group_t options_group( const std::string& options_group_name );
  /** \brief . */
  const options_group_t options_group( const std::string& options_group_name ) const;

  //@}

  /** \brief . */
  static bool options_group_exists( const options_group_t& options_group );

  /** @name Determine what options groups where not accessed.
    *
    * The only the options groups accessed through the this->options_group(...)
    * functions are maked as accessed.  When the options groups are
    * accessed through the iterator access, it is assumed that the client
    * will not need this other information.  Note that all of the
    * flags are false by default.
    */
  //@{

  /// Reset the flags to false for if the options groups was accessed.
  void reset_unaccessed_options_groups();

  /// Print a list of options groups never accessed (accessed flag is falsed).
  void print_unaccessed_options_groups( std::ostream& out ) const;

  //@}

  /** @name Iterator access to options groups. */
  //@{

  /** \brief . */
  int				num_options_groups() const;
  /** \brief . */
  iterator		begin();
  /** \brief . */
  iterator		end();
  /** \brief . */
  const_iterator	begin() const;
  /** \brief . */
  const_iterator	end() const;

  //@}

private:
  typedef OptionsFromStreamUtilityPack::false_bool_t			false_bool_t;
  typedef OptionsFromStreamUtilityPack::option_to_value_map_t	option_to_value_map_t;
  typedef OptionsFromStreamUtilityPack::options_group_map_t	options_group_map_t;
  options_group_map_t											options_group_map_;

};	// end class OptionsFromStream

/** \brief . */
inline
const std::string&
options_group_name( OptionsFromStream::const_iterator& itr )
{
  return (*itr).first;
}

/** \brief . */
inline
const OptionsFromStream::options_group_t
options_group( OptionsFromStream::const_iterator& itr )
{
  return OptionsFromStream::options_group_t(
    const_cast<OptionsFromStreamUtilityPack::option_to_value_map_t*>(&(*itr).second.first) );
}


// @ }

/* @name Return the name or value of an option from an OptionsGroup iterator. */
// @ {

//
using OptionsFromStreamUtilityPack::option_name;
//
using OptionsFromStreamUtilityPack::option_value;

// @ }

//	end namespace OptionsFromStreamPack 
// @ }

}	// end namespace OptionsFromStreamPack

// //////////////////////////////////////////////////////////////////////////////////////////
// Inline definitions.

namespace OptionsFromStreamPack {

namespace OptionsFromStreamUtilityPack {

// class OptionsGroup.

inline
std::string& OptionsGroup::option_value( const std::string& option_name ) {
  option_to_value_map_t::iterator itr = option_to_value_map_->find( option_name );
  return ( itr != option_to_value_map_->end() ? (*itr).second : option_does_not_exist_ );
}

inline
const std::string& OptionsGroup::option_value( const std::string& option_name ) const {
  option_to_value_map_t::const_iterator itr = option_to_value_map_->find( option_name );
  return ( itr != option_to_value_map_->end() ? (*itr).second : option_does_not_exist_ );
}

inline
bool OptionsGroup::option_exists( const std::string& option_value ) {
  return &option_value != &option_does_not_exist_;
}

inline
bool OptionsGroup::options_group_exists() const {
  return option_to_value_map_ != 0;
}

inline
int	OptionsGroup::num_options() const {
  return option_to_value_map_->size();
}

inline
OptionsGroup::iterator		OptionsGroup::begin() {
  return option_to_value_map_->begin();
}

inline
OptionsGroup::iterator		OptionsGroup::end() {
  return option_to_value_map_->end();
}

inline
OptionsGroup::const_iterator	OptionsGroup::begin() const {
  return option_to_value_map_->begin();
}

inline
OptionsGroup::const_iterator	OptionsGroup::end() const {
  return option_to_value_map_->end();
}

}	// end namespace OptionsFromStreamPack

// class OptionsFromStream

inline
OptionsFromStream::OptionsFromStream()
{}

inline
OptionsFromStream::OptionsFromStream( std::istream& in ) {
  read_options(in);
}

inline
void OptionsFromStream::clear_options() {
  options_group_map_.clear();
}

inline
bool OptionsFromStream::options_group_exists( const options_group_t& options_group )
{
  return options_group.options_group_exists();
}

inline
int	OptionsFromStream::num_options_groups() const {
  return options_group_map_.size();
}

inline
OptionsFromStream::iterator OptionsFromStream::begin() {
  return options_group_map_.begin();
}

inline
OptionsFromStream::iterator OptionsFromStream::end() {
  return options_group_map_.end();
}

inline
OptionsFromStream::const_iterator OptionsFromStream::begin() const {
  return options_group_map_.begin();
}

inline
OptionsFromStream::const_iterator OptionsFromStream::end() const {
  return options_group_map_.end();
}

}	// end namespace OptionsFromStreamPack

#endif	// OPTIONS_FROM_STREAM_H
