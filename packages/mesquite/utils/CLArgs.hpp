/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2008 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    (2008) kraftche@cae.wisc.edu    

  ***************************************************************** */


/** \file CLArgs.hpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#ifndef MSQ_CLARGS_HPP
#define MSQ_CLARGS_HPP


class CLArgImpl;

#include "Mesquite.hpp"
#include <iosfwd>
#include <vector>
#include <string>
#include <sstream>
#include <iterator>

/**\brief Parse command-line arguments
 *
 * This class provides a mechanism for parsing and to some extend validating
 * command-line arguments.  Use of this class can be divided into three steps:
 * 
 * 1) Call *_flag and arg methods to define acceptable command line arguments
 *
 * 2) Call parse_options to parse the command line arguments according to
 *    acceptable flags defined in step 1.
 *
 * 3) Check the values in registerd callback class instances.
 *
 * The '-h' flag is reserved for reqesting a description of the allowable
 * arguments (help).  If it is encountered inside parse_options, the help
 * will be printed and the program will be terminated w/out returning from
 * parse_options.
 *
 * The '-M' flag is similar to '-h', except that the help is written
 * in UNIX 'man page' format.
 */
class CLArgs
{
  public:
  
      /**\brief Base class for callback interface (type-independent functions) */
    class ArgIBase {
      private: bool wasSeen; //!< Keep track of whether or not this flag was encountered
      public: ArgIBase() : wasSeen(false) {} //!< constructor
              virtual ~ArgIBase() {}         //!< virtual destructor for proper cleanup
                /**\brief Get short description string for usage output or empty string for default*/
              virtual std::string brief() const { return std::string(); }
                /**\brief Get short description string for UNIX man page output or empty string for default */
              virtual std::string manstr() const { return std::string(); }
                /**\brief Get optional additional info to print with flag description */
              virtual std::string desc_append() const { return std::string(); }
                /**\brief Get optional string containing default value for option if not specified by user */
              virtual std::string default_str() const { return std::string(); }
                /**\brief Mark this flag as having been specified by the user */
              void set_seen() { wasSeen = true; }
                /**\brief Test if the user specified this flag */
              bool seen() const { return wasSeen; }
    };
  
      /**\brief Interface for type-specific callback classes */
    template <typename T> class ArgTemplateI : public ArgIBase {
      public: virtual bool value( const T& val ) = 0;
    };
      /**\brief Trivial implementation for type-specific classes */
    template <typename T> class ArgTemplate :public ArgTemplateI<T> {
      private: T mValue;         //!< The default or user-specified value for an option.
               bool haveDefault; //!< True if app. provided default value.
      public:  virtual ~ArgTemplate() {}
               virtual bool value( const T& val ) //!< Set value
               { 
                  mValue = val; 
                  ArgTemplateI<T>::set_seen(); 
                  return true; 
               }
               const T& value() const { return mValue; } //!< get value
                /**\brief Initialize with default value */
               ArgTemplate( const T& initial_value ) : mValue(initial_value), haveDefault(true) {}
                /**\brief Initialize without default value */
               ArgTemplate() : mValue(T()), haveDefault(false) {}
                /**\brief Get string representation of default value, or empty string of no default value */
               virtual std::string default_str() const {
                 std::ostringstream ss;
                 if (haveDefault)
                   ss << mValue;
                 return ss.str();
               }    
               
    };

      /**\brief Trivial implementation for type-specific classes */
    template <typename T> class ArgListTemplate : public ArgTemplateI< std::vector<T> > {
      private: std::vector<T> mValue;         //!< The default or user-specified value for an option.
               bool haveDefault; //!< True if app. provided default value.
      public:  virtual ~ArgListTemplate() {}
               virtual bool value( const std::vector<T>& val ) //!< Set value
               { 
                  mValue = val; 
                  ArgTemplateI< std::vector<T> >::set_seen(); 
                  return true; 
               }
               const std::vector<T>& value() const { return mValue; } //!< get value
                /**\brief Initialize with default value */
               ArgListTemplate( const std::vector<T>& initial_value ) : mValue(initial_value), haveDefault(true) {}
                /**\brief Initialize without default value */
               ArgListTemplate() : haveDefault(false) {}
                /**\brief Get string representation of default value, or empty string of no default value */
               virtual std::string default_str() const {
                 std::ostringstream ss;
                 std::copy( mValue.begin(), mValue.end(), 
                   std::ostream_iterator<T>( ss, ", " ) );
                 return ss.str();
               }    
               
    };
    
      /**\brief Callback API for a string argument */
    typedef ArgTemplateI< std::string >     StringArgI;
      /**\brief Callback API for an integer argument */
    typedef ArgTemplateI< int >                 IntArgI;
      /**\brief Callback API for a long integer argument */
    typedef ArgTemplateI< long >                LongArgI;
      /**\brief Callback API for a double-precision floating-point argument */
    typedef ArgTemplateI< double >              DoubleArgI;
      /**\brief Callback API for a Boolean or toggle argument */
    typedef ArgTemplateI< bool >                ToggleArgI;
      /**\brief Callback API for an integer list argument */
    typedef ArgTemplateI< std::vector<int> >    IntListArgI;
      /**\brief Callback API for a double list argument */
    typedef ArgTemplateI< std::vector<double> > DoubleListArgI;

      /**\brief Trivial callback implementation for a string argument */
    typedef ArgTemplate< std::string>      StringArg;
      /**\brief Trivial callback implementation for an integer argument */
    typedef ArgTemplate< int >                 IntArg;
      /**\brief Trivial callback implementation for a long integer argument */
    typedef ArgTemplate< long >                LongArg;
      /**\brief Trivial callback implementation for a ouble-precision floating-point argument */
    typedef ArgTemplate< double >              DoubleArg;
      /**\brief Trivial callback implementation for a Boolean or toggle argument */
    typedef ArgTemplate< bool >                ToggleArg;
      /**\brief Trivial callback implementation for an integer list argument */
    typedef ArgListTemplate< int >    IntListArg;
      /**\brief Trivial callback implementation for a double list argument */
    typedef ArgListTemplate< double > DoubleListArg;

      /**\brief String arugment that is limited to a list of acceptable keywords
       *
       * A specialized string arugment implementation that limits the
       * acceptable string argument to one of a list of keywords.
       * A case-insensitive comparison is done with the allowed keywords.
       * The "value" has the case of the keyword rather than the case
       * used in the literal value specified in the command line argument.
       */
    class KeyWordArg : public StringArg
    {
      private:
        std::vector< std::string > mKeyWords;
        void initialize( const char* keyword_list[], int list_length );
      public:
        KeyWordArg( const char* keyword_list[], int list_length )
          { initialize( keyword_list, list_length ); }
        KeyWordArg( const char* default_val, const char* keyword_list[], int list_length )
          : StringArg( default_val )
          { initialize( keyword_list, list_length ); }
        virtual bool value( const std::string& val );
        virtual std::string brief() const;
        virtual std::string manstr() const;
        static bool compare_no_case( const char* s1, const char* s2 );
    };
    
    class IntRange
    {
      private:
        int mMin, mMax;
      public:
        IntRange( const int* min, const int* max );
        bool is_valid( int val ) const;
        std::string desc_append() const;
    };
    
      /**\brief Integer argument constrained to a range of valid values. */
    class IntRangeArg : public IntArg
    {
      private:
        IntRange mRange;
      public:
        IntRangeArg( const int* min = 0, const int* max = 0 )
          : mRange( min, max ) {}
        IntRangeArg( int default_val, const int* min, const int* max )
          : IntArg(default_val), mRange( min, max ) {}
        bool value( const int& val );
        const int& value() const { return IntArg::value(); }
        std::string desc_append() const
          { return mRange.desc_append(); }
    };
    
      /**\brief Integer list argument constrained to a range of valid values. */
    class IntListRangeArg : public IntListArg 
    {
      private:
        IntRange mRange;
      public:
        IntListRangeArg( const int* min = 0, const int* max = 0 )
          : mRange( min, max ) {}
        bool value( const std::vector<int>& val );
        const std::vector<int>& value() const { return IntListArg::value(); }
        std::string desc_append() const
          { return mRange.desc_append(); }
    };
    
    class DoubleRange
    {
      private:
        bool haveMin, haveMax, mInclusive;
        double mMin, mMax;
      public:
        DoubleRange( const double* min, const double* max, bool inclusive );
        bool is_valid( double value ) const;
        std::string desc_append() const;
    };
    
      /**\brief Double argument constrained to a range of valid values. */
    class DoubleRangeArg : public DoubleArg
    {
      private:
        DoubleRange mRange;
      public:
        DoubleRangeArg( const double* min = 0, 
                        const double* max = 0, 
                        bool inclusive = true ) 
                      : mRange( min, max, inclusive )
                        {}
        DoubleRangeArg( double default_val,
                        const double* min = 0, 
                        const double* max = 0, 
                        bool inclusive = true ) 
                      : DoubleArg(default_val), 
                        mRange( min, max, inclusive )
                        {}
        bool value( const double& val );
        const double& value() const { return DoubleArg::value(); }
        std::string desc_append() const
          { return mRange.desc_append(); }
    };
    
      /**\brief Double list argument constrained to a range of valid values. */
    class DoubleListRangeArg : public DoubleListArg
    {
      private:
        DoubleRange mRange;
      public:
        DoubleListRangeArg( const double* min = 0, 
                            const double* max = 0, 
                            bool inclusive = true ) 
                          : mRange( min, max, inclusive )
                            {}
        bool value( const std::vector<double>& val );
        const std::vector<double>& value() const { return DoubleListArg::value(); }
        std::string desc_append() const
          { return mRange.desc_append(); }
    };
        


  public:
  
    /**\brief Define basic program 
     *
     *\param progname  The program name
     *\param brief_desc A brief description of the purpose of the program.
     *\param desc Program description for documentation.
     */
    CLArgs( const char* progname, const char* brief_desc, const char* desc );
    
    ~CLArgs();
    
      /**\brief Check if flag is undefined */
    bool is_flag_available( char fl ) const;
    
    /**\brief Register a flag that requires a string argument.
     *
     * Define a flag of the form "-f <somestring>".
     *\param fl       The character for the flag. 
     *\param name     The name of the flag argument
     *\param desc     A description of the purpose of the flag and argument.
     *\parma callback Object instance to which to pass the parsed argument value.
     *\return         false if flag is already in use, true otherwise.
     */
    bool str_flag( char fl, 
                   const char* name, 
                   const char* desc, 
                   StringArgI* callback );

    /**\brief Register a flag that requires an integer argument.
     *
     * Define a flag of the form "-f <int>".
     *\param fl       The character for the flag. 
     *\param name     The name of the flag argument
     *\param desc     A description of the purpose of the flag and argument.
     *\parma callback Object instance to which to pass the parsed argument value.
     *\return         false if flag is already in use, true otherwise.
     */
    bool int_flag( char fl, 
                   const char* name, 
                   const char* desc, 
                   IntArgI* callback );

    /**\brief Register a flag that requires an integer argument.
     *
     * Define a flag of the form "-f <int>".
     *\param fl       The character for the flag. 
     *\param name     The name of the flag argument
     *\param desc     A description of the purpose of the flag and argument.
     *\parma callback Object instance to which to pass the parsed argument value.
     *\return         false if flag is already in use, true otherwise.
     */
    bool long_flag( char fl, 
                    const char* name, 
                    const char* desc, 
                    LongArgI* callback );

    /**\brief Register a flag that requires a real umber argument.
     *
     * Define a flag of the form "-f <double>".
     *\param fl        The character for the flag. 
     *\param name      The name of the flag argument
     *\param desc      A description of the purpose of the flag and argument.
     *\parma callback  Object instance to which to pass the parsed argument value.
     *\param inclusive If true, accept 'min' or 'max': [min,max].  If
     *                 false, reject 'min' and 'max' values: (min,max).
     *\return          false if flag is already in use, true otherwise.
     */
    bool double_flag( char fl, 
                      const char* name, 
                      const char* desc, 
                      DoubleArgI* callback );
    
    /**\brief Register a pair of flags that accept no arguments and have
     *        opposing affects.
     *
     * Regstier a flag of the form [-f|-F], where one implies a true
     * state and the other a false state (i.e. enable or disable some
     * functionality.)
     *\param on_flag  Flag corresponding to true or 'on' state.
     *\param off_flag Flag corresponding to false or 'off' state.
     *\param desc     A description of the purpose of the flag and argument.
     *\parma callback Object instance to which to pass the parsed argument value.
     *\return         false if flag is already in use, true otherwise.
     */
    bool toggle_flag( char on_flag, 
                      char off_flag, 
                      const char* desc, 
                      ToggleArgI* callback );
    
    /**\brief Register a flag with no value.
     *
     * Define a flag such that the state of the option is considered
     * to be false unless flag is specified.  If the flag is specified,
     * the option is considered to be true.
     *\param fl       The character for the flag. 
     *\param desc     A description of the purpose of the flag and argument.
     *\parma callback Object instance to which to pass the parsed argument value.
     *\return         false if flag is already in use, true otherwise.
     */
    bool toggle_flag( char fl, const char* desc, ToggleArgI* callback );
    
    
    /**\brief Register a flag that accepts a list of positive integer arguments.
     *
     * Define a flag that accepts a list of ranges of positive integer values
     * separated by commas.  A zero value is rejected.  Ranges can be either
     * a single value or pair of values separated by a dash.  For example:
     * "-f 1,4-10,2,20-25".
     * 
     * Use 'get_int_list' to query values of flag.
     *
     *\param fl       The character for the flag. 
     *\param desc     A description of the purpose of the flag and argument.
     *\parma callback Object instance to which to pass the parsed argument value.
     *\return         false if flag is already in use, true otherwise.
     */
    bool id_list_flag( char fl, const char* desc, IntListArgI* callback );

    /**\brief Register a flag that requires a comma-separated 
     *        list of integer values.
     *
     * Define a flag that has an integer list for its arguments.  The 
     * integers must be specified as a comma-separated list.
     *
     * Use 'limit_list_flag' to limit the number of values accepted
     * in the argument.  If 'limit_list_flag' is never called, any
     * number of argument values will be accepted.
     *
     *\param fl       The character for the flag. 
     *\param desc     A description of the purpose of the flag and argument.
     *\parma callback Object instance to which to pass the parsed argument value.
     *\return         false if flag is already in use, true otherwise.
     */
    bool int_list_flag( char fl, const char* desc, IntListArgI* callback );

    /**\brief Register a flag that requires a comma-separated 
     *        list of double values.
     *
     * Define a flag that has an double list for its arguments.  The 
     * values must be specified as a comma-separated list.
     *
     * Use 'limit_list_flag' to limit the number of values accepted
     * in the argument.  If 'limit_list_flag' is never called, any
     * number of argument values will be accepted.
     *
     *\param fl       The character for the flag. 
     *\param desc     A description of the purpose of the flag and argument.
     *\parma callback Object instance to which to pass the parsed argument value.
     *\return         false if flag is already in use, true otherwise.
     */
    bool double_list_flag( char fl, 
                           const char* desc, 
                           DoubleListArgI* callback );
    
    /**\brief Set a limit on the number of values accepted for a list-type
     *        flag.  May be called multiple times for a flag.
     *
     * Add a limit on the number of values accepted for a list-type flag.
     * This function may be called multiple times if the flag should 
     * accept a finite set of different arugment value counts.
     *
     *\param fl  The flag
     *\param num_values  The number of values to accept
     *\param value_names An array of 'num_value' strings specifying the
     *                   name of each value.
     *\return false if flag is not defined as a list-type flag or this
     *        method has already been called with the same flag AND
     *        number of values.  True otherwise.
     */
    bool limit_list_flag( char fl,
                          int num_values,
                          const char* const* value_names );

    /**\brief Specify that an argument without a flag is expected.
     *
     * Arguments are parsed in the order they are added, with the
     * exception that all optional args are parsed after all required
     * args.
     * \param name 'name' of argument to display in help (e.g. "output_file");
     */
    void add_required_arg( const char* name );
    void add_optional_arg( const char* name );
  
    /**\brief Parse argument list.
     * 
     *\param argc The argument list length passed to the main() routine.
     *            The first value is assumed to be the executable name.
     *            This this value must be at least 1.
     *\param argv The argument list as passed to main().
     *\param args_out  The list of non-flag arguments encountered, as
     *            defined by the 'args' method.  If the 'args' method
     *            has not been called, no non-flag arguments are accepted
     *            and this list will be empty.
     *\param error_stream  stream to which to write error messages.
     *\return true if all arguments were accepted.  false otherwise.
     */
    bool parse_options( int argc, 
                        char* argv[],
                        std::vector< std::string >& args_out,
                        std::ostream& error_stream );
    
    
    
    
    /**\brief Write help
     * 
     * Write help text to passed stream.
     */
    void print_help( std::ostream& stream ) const;
    
    /**\brief Write UNIX man page
     * 
     * Write man page to passed stream.
     */
    void print_man_page( std::ostream& stream ) const;
    
    /**\brief prinint usage (brief help)
     */
    void print_usage( std::ostream& stream ) const;
    
  private:
  
    CLArgImpl* impl;
};

template <typename T> std::ostream& operator<<( std::ostream& str, const std::vector<T>& list )
{ 
  typename std::vector<T>::const_iterator i = list.begin();
  if (i != list.end()) {
    str << *i;
    for (++i; i != list.end(); ++i)
      str << ',' << *i;
  }
  return str;
}

#endif
