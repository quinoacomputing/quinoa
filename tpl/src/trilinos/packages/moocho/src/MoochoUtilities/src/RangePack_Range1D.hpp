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

#ifndef RANGE1D_H
#define RANGE1D_H

#include "RTOpPack_Types.hpp"
#include "Teuchos_Range1D.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace RangePack {

/** \brief. One-based subregion index range class.
 * 
 * The class <tt>%Range1D</tt> abstracts a 1-D, 1-based, range of indexes.
 * It is used to index into vectors and matrices and return subregions of them
 * respectively.
 *
 * Constructing using \c Range1D() yields a range that represents the entire dimension
 * of an object <tt>[1, max_ubound]</tt> (an entire vector, all the rows in a matrix,
 * or all the columns in a matrix etc.).
 *
 * Constructing using <tt>\ref Range1D::Range1D "Range1D(INVALID)"</tt> yields an invalid
 * range <tt>[1,0]</tt> with <tt>size() == 0</tt>.  In fact the condition
 * <tt>size() == 0</tt> is the determining flag that a range is not valid.
 * Once constructed with <tt>Range1D(INVALID)</tt>, a <tt>%Range1D</tt> object can
 * pass through many other operations that may change <tt>%lbound()</tt> and <tt>%ubound()</tt>
 * but will never change <tt>size() == 0</tt>.
 *
 * Constructing using <tt>\ref Range1D::Range1D "Range1D(lbound,ubound)"</tt> yields a finite dimensional range.
 * The validity of constructed range will only be checked if \c TEUCHOS_DEBUG is defined.
 *
 * There are many \ref Range1D_funcs_grp "non-member functions" that can be used with <tt>%Range1D</tt> objects.
 *
 * The default copy constructor and assignment operator functions are allowed since they have
 * the correct semantics.
 */
class Range1D {
public:
  /** \brief . */
  typedef RTOpPack::Ordinal  Index;
  /** \brief . */
  enum EInvalidRange { INVALID };
  /// Range1D(INVALID)
  static const Range1D Invalid;
  /** \brief Constructs a range representing the entire range.
   *
   * Postconditions: <ul>
   *	<li> <tt>this->full_range() == true</tt>
   * <li> <tt>this->size() ==</tt> a very large number
   * <li> <tt>this->lbound() == 1</tt>
   * <li> <tt>this->ubound() ==</tt> a very large number
   *	</ul>
   */
  /** \brief . */
  inline Range1D();
  /** \brief Constructs an invalid (zero) range.
   *
   * Postconditions: <ul>
   *	<li> <tt>this->full_range() == false</tt>
   * <li> <tt>this->size() == 0</tt>
   * <li> <tt>this->lbound() == 1</tt>
   * <li> <tt>this->ubound() == 0</tt>
   *	</ul>
   */
  inline Range1D( EInvalidRange );
  /** \brief Constructs a range that represents the range <tt>[lbound, ubound]</tt>.
   *
   * Preconditions: <ul>
   *	<li> <tt>lbound >= 1</tt> (throw \c range_error)
   *	<li> <tt>lbound <= ubound</tt> (throw \c range_error)
   *	</ul>
   *
   * Postconditions: <ul>
   *	<li> <tt>this->full_range() == false</tt>
   * <li> <tt>this->size() == ubound - lbound + 1</tt>
   * <li> <tt>this->lbound() == lbound</tt>
   * <li> <tt>this->ubound() == ubound</tt>
   *	</ul>
   */
  inline Range1D(Index lbound, Index ubound);
  /// Returns \c true if the range represents the entire region (constructed from \c Range1D())
  inline bool full_range() const;
  /// Return lower bound of the range
  inline Index lbound() const;
  /// Return upper bound of the range
  inline Index ubound() const;
  /// Return the size of the range (<tt>ubound() - lbound() + 1</tt>)
  inline Index size() const;
  /// Return true if the index is in range
  inline bool in_range(Index i) const;
  /// Increment the range by a constant
  inline Range1D& operator+=( Index incr );
  /// Deincrement the range by a constant
  inline Range1D& operator-=( Index incr );

private:
  Index lbound_;
  Index ubound_;
  
  // assert that the range is valid
  inline void assert_valid_range(Index lbound, Index ubound) const;
  
}; // end class Range1D
  
/** \brief rng1 == rng2.
 *
 * @return Returns <tt>rng1.lbound() == rng2.ubound() && rng1.ubound() == rng2.ubound()</tt>.
 *
 * \relates Range1D
 */
inline bool operator==(const Range1D& rng1, const Range1D& rng2 )
{
  return rng1.lbound() == rng2.lbound() && rng1.ubound() == rng2.ubound();
}

/** \brief rng_lhs = rng_rhs + i.
 *
 * Increments the upper and lower bounds by a constant.
 *
 * Postcondition: <ul>
 *	<li> <tt>rng_lhs.lbound() == rng_rhs.lbound() + i</tt>
 *	<li> <tt>rng_lhs.ubound() == rng_rhs.ubound() + i</tt>
 *	</ul>
 *
 * \relates Range1D
 */
inline Range1D operator+(const Range1D &rng_rhs, Range1D::Index i)
{
    return Range1D(i+rng_rhs.lbound(), i+rng_rhs.ubound());
}

/** \brief rng_lhs = i + rng_rhs.
 *
 * Increments the upper and lower bounds by a constant.
 *
 * Postcondition: <ul>
 *	<li> <tt>rng_lhs.lbound() == i + rng_rhs.lbound()</tt>
 *	<li> <tt>rng_lhs.ubound() == i + rng_rhs.ubound()</tt>
 *	</ul>
 *
 * \relates Range1D
 */
inline Range1D operator+(Range1D::Index i, const Range1D &rng_rhs)
{
    return Range1D(i+rng_rhs.lbound(), i+rng_rhs.ubound());
}

/** \brief rng_lhs = rng_rhs - i.
 *
 * Deincrements the upper and lower bounds by a constant.
 *
 * Postcondition: <ul>
 *	<li> <tt>rng_lhs.lbound() == rng_rhs.lbound() - 1</tt>
 *	<li> <tt>rng_lhs.ubound() == rng_rhs.ubound() - 1</tt>
 *	</ul>
 *
 * \relates Range1D
 */
inline Range1D operator-(const Range1D &rng_rhs, Range1D::Index i)
{
    return Range1D(rng_rhs.lbound()-i, rng_rhs.ubound()-i);
}

/** \brief Return a bounded index range from a potentially unbounded index range.
 * 
 * Return a index range of lbound to ubound if rng.full_range() == true
 * , otherwise just return a copy of rng.
 *
 * Postconditions: <ul>
 *	<li> [<tt>rng.full_range() == true</tt>] <tt>return.lbound() == lbound</tt>
 *	<li> [<tt>rng.full_range() == true</tt>] <tt>return.ubound() == ubound</tt>
 *	<li> [<tt>rng.full_range() == false</tt>] <tt>return.lbound() == rng.lbound()</tt>
 *	<li> [<tt>rng.full_range() == false</tt>] <tt>return.ubound() == rng.ubound()</tt>
 *	</ul>
 *
 * \relates Range1D
 */
inline Range1D full_range(const Range1D &rng, Range1D::Index lbound, Range1D::Index ubound)
{	return rng.full_range() ? Range1D(lbound,ubound) : rng; }

/** \brief Convert from a 0-based Teuchos::Range1D object to a 1-based RangePack::Range1D object.
 *
 * \relates Range1D
 */
inline Range1D convert( const Teuchos::Range1D &rng )
{
  Range1D rngOut;
  if (rng.full_range()) {
    rngOut = Range1D();
  }
  else if (rng.size() == -1) {
    rngOut = Range1D::Invalid;
  }
  else if (rng.size() == 0) {
    rngOut = Range1D::Invalid;
  }
  else {
    rngOut = Range1D(rng.lbound()+1, rng.ubound()+1);
  }
  return rngOut;
}

/** \brief Convert from a 1-based RangePack::Range1D object to a 0-based Teuchos::Range1D object.
 *
 * \relates Range1D
 */
inline Teuchos::Range1D convert( const Range1D &rng )
{
  Teuchos::Range1D rngOut;
  if (rng.full_range()) {
    rngOut = Teuchos::Range1D();
  }
  else if (rng.size() == 0) {
    rngOut = Teuchos::Range1D::Invalid;
  }
  else {
    rngOut = Teuchos::Range1D(rng.lbound()-1, rng.ubound()-1);
  }
  return rngOut;
}


//@}

// //////////////////////////////////////////////////////////
// Inline members

inline
Range1D::Range1D()
  : lbound_(1), ubound_(std::numeric_limits<Index>::max()-1)
{}

inline
Range1D::Range1D( EInvalidRange )
  : lbound_(1), ubound_(0)
{}


inline
Range1D::Range1D(Index lbound, Index ubound)
  : lbound_(lbound), ubound_(ubound)
{
  assert_valid_range(lbound,ubound);
}

inline
bool Range1D::full_range() const {
  return ubound_ == (std::numeric_limits<Index>::max()-1);
}

inline
Range1D::Index Range1D::lbound() const {
  return lbound_;
}

inline
Range1D::Index Range1D::ubound() const {
  return ubound_;
}

inline
Range1D::Index Range1D::size() const {
  return 1 + ubound_ - lbound_;
}

inline
bool Range1D::in_range(Index i) const {
  return lbound_ <= i && i <= ubound_;
}

inline
Range1D& Range1D::operator+=( Index incr ) {
  assert_valid_range( lbound_ + incr, ubound_ + incr );
  lbound_ += incr;
  ubound_ += incr;
  return *this;
}

inline
Range1D& Range1D::operator-=( Index incr ) {
  assert_valid_range( lbound_ - incr, ubound_ - incr );
  lbound_ -= incr;
  ubound_ -= incr;
  return *this;
}

// See Range1D.cpp
inline
void Range1D::assert_valid_range(Index lbound, Index ubound) const {
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    lbound < 1, std::range_error
    ,"Range1D::assert_valid_range(): Error, lbound ="<<lbound<<" must be greater than 0." );
  TEUCHOS_TEST_FOR_EXCEPTION(
    lbound > ubound, std::range_error
    ,"Range1D::assert_valid_range(): Error, lbound = "<<lbound<<" > ubound = "<<ubound );
#endif
}

} // end namespace RangePack

#endif // end RANGE1D_H
