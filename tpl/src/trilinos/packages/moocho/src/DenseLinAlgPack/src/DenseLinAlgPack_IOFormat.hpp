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

#ifndef LINALGPACK_IO_FORMAT_H
#define LINALGPACK_IO_FORMAT_H

#include <ios>

#include "DenseLinAlgPack_IOBasic.hpp"

namespace DenseLinAlgPack {
namespace LinAlgPackIO {

/** \brief . */
/* * This class is used to encapsulate a set of bit flags.
  */
class bit_flags {
public:

  /** \brief . */
  typedef LinAlgPackIO::fmtflags fmtflags;	

  /// Initialize the flags to 0x0000
  bit_flags() : flags_((fmtflags)(0x0000)) {}

  /// Get the flags
  fmtflags flags() const {
    return flags_;
  }

  /// Set the flags
  fmtflags flags(fmtflags f) {
    fmtflags f_tmp = flags_;
    flags_ = f;
    return f_tmp;
  }

  /// Set a flag
  fmtflags setf(fmtflags f) {
    return flags( (fmtflags)(flags() | f) );
  }

  /// Set a flag under a mask
  fmtflags setf(fmtflags f, fmtflags mask) {
    return flags( (fmtflags)(flags() | (f&mask)) );
  }

  /// Unset a flag(s)
  void unsetf(fmtflags mask) {
    flags( (fmtflags)(flags() & ~mask) );
  }

private:
  fmtflags flags_;

};	// end class flags

/** \brief . */
/* * Memento class that saves the state of a standard stream.
  *
  * This is a variation of the "Memento" Pattern in Gama et. al.
  *
  * Allow default copy constructor and assignment operator.
  * The default constuctor is private so the only way to 
  * create an object of this class if from another object
  * or from save_format().
  */
class ios_format_memento {
public:
  /// Save a streams format
  static ios_format_memento	save_format(const std::ios& s);

  /// Reset a steams format
  void set_format(std::ios& s) const;

private:
  LinAlgPackIO::fmtflags	flags_;
  int				prec_;
  int				wdt_;
  int				fill_;

  /** \brief . */
  ios_format_memento() : flags_((fmtflags)(0)), prec_(6), wdt_(0)
    , fill_(' ') {}

};	// end class ios_format_memento

template<class T> class bound_format;
template<class T> class const_bound_format;

/** \brief . */
/* * This class is used to allow some flexiblility in inputing
  * and outputing of DVector, DVectorSlice, DMatrix and DMatrixSlice objects.
  * It maintains format information for a #std::ios_base# object and adds
  * additional format information.
  *
  * The idea for this class is based on the user-defined manipulator given in
  * Stroustrup 1997, p. 635.  The main difference is that the non-member funciton
  * bind(format,object) is used to create a bound format object and not that
  * member function operator().  To use operator() would require a template
  * member function and they are not supposed to be allowed in MS VC++ 5.0.
  *
  * Warning: The default constructor sets the formating flags to false but this
  * may not be the default value for the actual ios stream.  It is recommended
  * that you start with the formating settings of a valid ios stream and then
  * change the flags that you want.  You can do this on construction with
  *
  * #format f(cin);#\\
  *
  * or after construction using
  *
  * #format f;#\\
  * #f.copy_format(cin);#\\
  *
  * then set the formating using the manipulator like functions.
  *
  * The default copy constructor and assignment operator functions are allowed.
  */
class format {
public:

  /** \brief . */
  typedef LinAlgPackIO::fmtflags fmtflags;

  /// Sets format to defaults
  format() : prec_(6), wdt_(0), fill_(' ') {}

  /// Copy the formats from a ios stream (This is the suggested method.
  format(const std::ios& s) {
    copy_format(s);
  }

  /* * @name Access format flags objects
    */
  // @{	
  
  /** \brief . */
  bit_flags& ios_base_flags()				{	return ios_base_flags_;	}
  /** \brief . */
  const bit_flags& ios_base_flags() const	{	return ios_base_flags_;	}

  /** \brief . */
  bit_flags& extra_flags()				{	return extra_flags_;	}
  /** \brief . */
  const bit_flags& extra_flags() const	{	return extra_flags_;	}

  // @}

  /* * @name Manipulator like functions for standard floating-point number formating.
    *
    * These member functions are ment to mirror the standard input/output
    * manipulator functions.  Their names are the same and you can set
    * a list of options in one statement like:
    *
    * #format_obj.showpos().setprecision(6).fixed();#\\
    *
    * The primary difference is that #width# presists over all uses in output
    * operations while the standard #width# manipulator only applys to the next
    * operation.  See \Ref{bound_format} and \Ref{operator>>} for details
    * on how objects of this class are used to specify formating for input/output
    * operations.
    */
  // @{

  /** \brief . */
  format& showpoint()		{	ios_base_flags().setf(std::ios_base::showpoint); return *this;	}
  /** \brief . */
  format& noshowpoint()	{	ios_base_flags().unsetf(std::ios_base::showpoint); return *this;	}
  /** \brief . */
  format& showpos()		{	ios_base_flags().setf(std::ios_base::showpos); return *this;	}
  /** \brief . */
  format& noshowpos()		{	ios_base_flags().unsetf(std::ios_base::showpos); return *this;	}
  /** \brief . */
  format& skipws()		{	ios_base_flags().setf(std::ios_base::skipws); return *this;	}
  /** \brief . */
  format& noskipws()		{	ios_base_flags().unsetf(std::ios_base::skipws); return *this;	}
  /** \brief . */
  format& uppercase()		{	ios_base_flags().setf(std::ios_base::uppercase); return *this;	}
  /** \brief . */
  format& nouppercase()	{	ios_base_flags().unsetf(std::ios_base::uppercase); return *this;	}
  /** \brief . */
  format& internal()		{
    ios_base_flags().setf(std::ios_base::internal, std::ios_base::adjustfield);
    return *this;
  }
  /** \brief . */
  format& left()			{
    ios_base_flags().setf(std::ios_base::left, std::ios_base::adjustfield);
    return *this;
  }
  /** \brief . */
  format& right()			{
    ios_base_flags().setf(std::ios_base::right, std::ios_base::adjustfield);
    return *this;
  }
  /** \brief . */
  format& general()		{
    ios_base_flags().setf((fmtflags)0, std::ios_base::floatfield);
    return *this;
  }
  /** \brief . */
  format& fixed()			{
    ios_base_flags().setf(std::ios_base::fixed, std::ios_base::floatfield);
    return *this;
  }
  /** \brief . */
  format& scientific()	{
    ios_base_flags().setf(std::ios_base::scientific, std::ios_base::floatfield);
    return *this;
  }
  /** \brief . */
  format& setfill(int c)	{	fill_ = c; return *this;	}
  /** \brief . */
  format& setprecision(int p)	{	prec_ = p; return *this;	}
  /** \brief . */
  format& setw(int w)		{	wdt_ = w; return *this;	}

  // @}

  /* * @name Manipulator like functions for extra I/O formatin.
    *
    * These member functions are ment to mirror the type of standard input/output
    * manipulator functions accept they are used to set extra flags for 
    * #DenseLinAlgPack::ignore_dim_bit# and #DenseLinAlgPack::no_insert_newlines_bit#.
    * This allows them to be set in the same statement that sets a standard
    * option.  For example you can write:
    *
    * #format_obj.showpos().setprecision(6).fixed().ignore_dim().no_insert_newlines();#\\
    *
    * Like the member functions that simulate the standard manipulators these
    * options presist after input/output operations.
    */
  // @{

  /** \brief . */
  format& ignore_dim()
  {
    extra_flags().setf((fmtflags)(ignore_dim_bit));
    return *this;
  }
  /** \brief . */
  format& no_ignore_dim()
  {
    extra_flags().unsetf((fmtflags)(ignore_dim_bit));
    return *this;
  }
  /** \brief . */
  format& insert_newlines()
  {
    extra_flags().unsetf((fmtflags)(no_insert_newlines_bit));
    return *this;
  }
  /** \brief . */
  format& no_insert_newlines()
  {
    extra_flags().setf((fmtflags)(no_insert_newlines_bit));
    return *this;
  }

  // @}

  /* * @name Other access functions
    */
  // @{

  /** \brief . */
  int precision() const			{	return prec_;	}
  /** \brief . */
  int precision(int p)			{	int tmp = prec_; prec_ = p; return tmp;	}
  /** \brief . */
  int width() const				{	return wdt_;	}
  /** \brief . */
  int width(int w)				{	int tmp = wdt_; wdt_ = w; return tmp;	}
  /** \brief . */
  int fill() const				{	return fill_;	}
  /** \brief . */
  int fill(int c)					{	int tmp = fill_; fill_ = c; return tmp;	}

  // @}

  /* * @name Utility functions for saving and replacing a streams format state
    */
  // @{

  /// Copy a streams format to this object
  void copy_format(const std::ios& s);

  /// Set a streams format to the one in this object
  void set_format(std::ios& s) const;
  
  // @}


// ToDo: Enable these once member templates are supported.  This is much more
// attractive syntasticaly than using bind() and cbind() 
//
//	///
//	template<class T>
//	bound_format<T> operator()(T& obj);
//
//	///
//	template<class T>
//	const_bound_format<T> operator()(const T& obj) const;

private:
  bit_flags ios_base_flags_;
  bit_flags extra_flags_;
  int prec_;
  int wdt_;
  int fill_;

};	// end class LinAlgPackIOFormat


// ////////////////////////////////////////////////////////////////////////////////
// Format binding classes and functions

template<class T> class const_bound_format;

/** \brief . */
/* * Special I/O manipulator class for non-constant objects.
  *
  * This class as a special manipulator and is composed of a composition
  * of a format object and a templated object that is used in an input
  * or output operation.
  *
  * An object of this class should only be created by using 
  * the bind(format,obj) function.  It object can not default constructed
  * or assigned.
  *
  * Objects of this class can be implicitly converted to const_bound_format
  * for use in output operations.
  */
template<class T>
class bound_format {
public:

  /** \brief . */
  bound_format(const format& f, T& obj) : f_(f), obj_(obj) {}

  /** \brief . */
  const format& f() const	{ return f_; }
  /** \brief . */
  T&	obj()				{ return obj_; }
  /** \brief . */
  const T& obj() const	{ return obj_; }
  
private:
  const format&	f_;
  T&				obj_;

  // not defined and not to be called
  bound_format();
  bound_format& operator=(const bound_format&);

};	// end class bound_format

/** \brief . */
/* * Special I/O manipulator class for constant objects.
  *
  * This class as a special manipulator and is composed of a composition
  * of a format object and a const templated object that is used only in
  * an output operation.
  *
  * An object of this class can only be created by using 
  * the bind(f,obj) operation.  It can not be default constructed or assigned.
  * It is only to be used in input or output operations.
  */
template<class T>
class const_bound_format {
public:

  /** \brief . */
  const_bound_format(const format& f, const T& obj) : f_(f), obj_(obj) {}

  /// Allow implicit conversion from a bound_format to a const bound_format
  const_bound_format(const bound_format<T>& bf) : f_(bf.f()), obj_(bf.obj()) {}

  /** \brief . */
  const format& f() const	{ return f_; }
  /** \brief . */
  const T& obj() const	{ return obj_; }

private:
  const format&	f_;
  const T&		obj_;

  // not defined and not to be called
  const_bound_format();
  const_bound_format& operator=(const const_bound_format&);

};	// end class const_bound_format

/** \brief . */
/* * Returns a bound_format<T> object using this format object.
  *
  * Using templated input/output operator functions this function
  * is used to bind a format to a non-const target object for a
  * input/output stream operation.  For example, to bind the format 
  * for a single output operation using the streams set formating 
  * and addition formating use something like:
  *
  * #DVector v(1.0,5);#\\
  * #format f; f.copy_format(cout);#\\
  * #f.width(10).fixed().ignore_dim();#\\
  * #cout << bind(f,v);#\\1
  */
template<class T>
inline bound_format<T> bind(const format& f, T& obj) {
  return bound_format<T>(f,obj);
}

/** \brief . */
/* * Returns a const_bound_format<T> using this format object.
  *
  * This function works the same as the previous \Ref{bind}#()# function
  * accept it is for a const target object.  It could not use
  * the name #bind# due to ambiguity problems.
  */
template<class T>
inline const_bound_format<T> cbind(const format& f, const T& obj) {
  return const_bound_format<T>(f,obj);
}

// ///////////////////////////////////////
// Inline functions for members of format

// Use once member templates are supported.  See class format
//
//template<class T>
//inline bound_format<T> format::operator()(T& obj) {
//	return bind(*this,obj);
//}

//template<class T>
//inline const_bound_format<T> format::operator()(const T& obj) const {
//	return cbind(*this,obj);
//}

}	// end namespace LinAlgPackIO
}	// end namespace DenseLinAlgPack

#endif // LINALGPACK_IO_FORMAT_H
