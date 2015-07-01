/* @HEADER@ */
// ************************************************************************
// 
//                 Playa: Programmable Linear Algebra
//                 Copyright 2012 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */

#ifndef PLAYA_VECTORDECL_HPP
#define PLAYA_VECTORDECL_HPP

#include "PlayaDefs.hpp"
#include "PlayaHandle.hpp"
#include "PlayaVectorBaseDecl.hpp"
#include "PlayaVectorSpaceDecl.hpp"
#include "PlayaVectorFunctorsDecl.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace Playa
{

template <class Scalar, int> class LCN;
  
/** 
 * User-level vector class. 
 *
 * <h2> Vector creation </h2>
 *
 * Ordinarily, you will never construct a Vector directly
 * from a derived type.  Rather, the createMember() method of
 * VectorSpace is used to build a vector of the appropriate
 * type, for example,
 * \code
 * VectorType<double> vecType = new EpetraVectorType();
 * int dimension = 100;
 * VectorSpace<double> space = vecType.createSpace(dimension);
 * Vector<double> x = space.createMember();
 * Vector<double> y = space.createMember();
 * \endcode
 * This hides from you all the ugly
 * details of creating a particular concrete type.
 *
 * You will frequently create an empty vector to be filled in later,
 * for example,
 * \code
 * Vector<double> y;
 * \endcode
 * Note that this vector isn't just empty, it's null. Not only does
 * it have no values assigned, it does not have a concrete type. An
 * call a method on a null vector will result in an error. What you
 * <it>can</it> do with a null vector is
 * <ul>
 * <li> assign another vector to it
 * \code
 * Vector<double> x = space.createVector();
 * Vector<Scalar> y;
 * y = x.copy();
 * \endcode
 * <li> assign the result of a vector operation to it
 * \code
 * Vector<Scalar> z = a*x + b*y;
 * \endcode
 *
 * <h2> Vector creation </h2>
 */
template <class Scalar>
class Vector : public Playa::Handle<VectorBase<Scalar> >
{
public:
  /** \name Constructors, Destructors, and Assignment Operators */
  //@{
  HANDLE_CTORS(Vector<Scalar>, VectorBase<Scalar>);

  /** \brief Construct a vector from an N-term LC */
  template <int N>
  Vector<Scalar>(const LCN<Scalar, N>& x);


  /**  \brief Assign a one-term LC to this vector */
  Vector<Scalar>& operator=(const LCN<Scalar, 1>& x);

  /**  \brief Assign a two-term LC to this vector */
  Vector<Scalar>& operator=(const LCN<Scalar, 2>& x);

  /**  \brief Assign a three-term LC to this vector */
  Vector<Scalar>& operator=(const LCN<Scalar, 3>& x);

  /**  \brief Assign an N-term LC to this vector */
  template <int N>
  Vector<Scalar>& operator=(const LCN<Scalar, N>& x);

  //@}

  /** \name Operators with assignment */
  //@{
  /** Add a vector into this vector */
  Vector<Scalar>& operator+=(const Vector<Scalar>& other);

  /** Add a one-term LC into this vector */
  Vector<Scalar>& operator+=(const LCN<Scalar, 1>& x);

  /** Add a two-term LC into this vector */
  Vector<Scalar>& operator+=(const LCN<Scalar, 2>& x);

  /** Add an N-term LC into this vector */
  template <int N>
  Vector<Scalar>& operator+=(const LCN<Scalar, N>& x);

  /** Subtract a vector from this vector */
  Vector<Scalar>& operator-=(const Vector<Scalar>& other);

  /** Subtract a one-term LC from this vector */
  Vector<Scalar>& operator-=(const LCN<Scalar, 1>& x);

  /** Subtract a two-term LC from this vector */
  Vector<Scalar>& operator-=(const LCN<Scalar, 2>& x);

  /** Subtract an N-term LC from this vector */
  template <int N>
  Vector<Scalar>& operator-=(const LCN<Scalar, N>& x);

  /** Scale by a constant */
  Vector<Scalar>& operator*=(const Scalar& a);

  /** Divide by a constant */
  Vector<Scalar>& operator/=(const Scalar& a);
  //@}

  /** \name Structural information */
  //@{
  /**  \brief My space */
  VectorSpace<Scalar> space() const 
    {return this->ptr()->space();}

  /**  \brief My communicator */
  const MPIComm& comm() const 
    {return this->space().comm();}

  /**  \brief My dimension  */
  int dim() const
    {
      return this->ptr()->space()->dim();
    }
  //@}
      

  /** \name Block operations */
  //@{
  /**  \brief get number of blocks */
  int numBlocks() const ;

  /**  \brief set the i-th block  */
  void setBlock(int i, const Vector<Scalar>& v);
      
  /**  \brief get the i-th block */
  const Vector<Scalar>& getBlock(int i) const;
      
  /**  \brief get the i-th block */
  Vector<Scalar> getNonConstBlock(int i) ;
      
  /**  \brief get the i-th block */
  const Vector<Scalar>& getBlock(const BlockIterator<Scalar>& b) const;
      
  /**  \brief get the i-th block */
  Vector<Scalar> getNonConstBlock(const BlockIterator<Scalar>& b);
  //@}



  /** \name Sequential data accessors */
  //@{
  /**  \brief Get the next chunk of values for read-only access */
  ConstDataChunk<Scalar> nextConstChunk() const ;
    
  /**  \brief Get the next chunk of values for read-write access */
  NonConstDataChunk<Scalar> nextChunk() ;

  /**  \brief Tell whether there are more chunks remaining to be accessed */
  bool hasMoreChunks() const ;

  /**  \brief Reset the data stream back to a state where all chunks are
   * considered unread. */
  void rewind() const ;
  //@}

  /** \name Random access to local elements */
  //@{
  /**  \brief const bracket operator for read-only random access to
   * local elements as specified by
   * a flat index that runs from 0 to space().numLocalElements(). 
   * If the vector does not consist of a single contiguous data chunk,
   * this might be slow (worst case would be O(N), if every element
   * is stored in its own chunk of length 1). 
   */
  const Scalar& operator[](int localIndex) const ;
    
  /**  \brief non-const bracket operator for read-write random access to
   * local elements as specified by
   * a flat index that runs from 0 to space().numLocalElements(). 
   * If the vector does not consist of a single contiguous data chunk,
   * this might be slow (worst case would be O(N), if every element
   * is stored in its own chunk of length 1). 
   */
  Scalar& operator[](int localIndex);

  /** \brief  parentheses operator for read-only random access to
   * local elements as specified by 
   * a block iterator and a flat index indicating the
   * element's location within that block. 
   */
  const Scalar& operator()(const BlockIterator<Scalar>& b,
    int localIndexWithinBlock) const ;
    
/** \brief  parentheses operator for read-write random access to
   * local elements as specified by 
   * a block iterator and a flat index indicating the
   * element's location within that block. 
   */
  Scalar& operator()(const BlockIterator<Scalar>& b,
    int localIndexWithinBlock) ;
  //@}

  /** \name Diagnostic output */
  //@{
  /**  \brief Return a short string description of the vector */
  std::string description() const ;

  /** \brief  Print the vector in some detail */
  void print(std::ostream& os) const ;    
  //@}


  /** \name Math operations */
  //@{
  /** \brief  Multiply this vector by a constant scalar factor 
   * \code
   * this = alpha * this;
   * \endcode
   */
  Vector<Scalar>& scale(const Scalar& alpha);

  /** 
   *  \brief Add a scaled vector to this vector times a constant:
   * \f$ this=\alpha x + \gamma \,this. \f$
   */
  Vector<Scalar>& update(const Scalar& alpha, const Vector<Scalar>& x, 
    const Scalar& gamma=1.0);
  /** 
   * \brief  Add two scaled vectors to this vector times a constant:
   * \f$ this=\alpha x + \beta y + \gamma \, this. \f$
   */
  Vector<Scalar>& update(const Scalar& alpha, const Vector<Scalar>& x, 
    const Scalar& beta, const Vector<Scalar>& y, 
    const Scalar& gamma);
  /** 
   * \brief  Add three scaled vectors to this vector times a constant:
   * \f$ this=\alpha x + \beta y + \gamma x + \delta \, this. \f$
   */
  Vector<Scalar>& update(const Scalar& alpha, const Vector<Scalar>& x, 
    const Scalar& beta, const Vector<Scalar>& y, 
    const Scalar& gamma, const Vector<Scalar>& z, 
    const Scalar& delta);

  /** 
   *  \brief Copy the values of another vector into this vector
   * \code
   * this = x
   * \endcode
   */
  Vector<Scalar>& acceptCopyOf(const Vector<Scalar>& x);

  /** 
   * \brief  Create a new vector that is a copy of this vector 
   */
  Vector<Scalar> copy() const ;

  /** 
   * \brief  In-place element-by-element product (Matlab dot-star operator)
   */
  Vector<Scalar>& selfDotStar(const Vector<Scalar>& other) ;

  /** 
   * \brief  In-place element-by-element division (Matlab dot-slash operator)
   */
  Vector<Scalar>& selfDotSlash(const Vector<Scalar>& other) ;

  /** 
   * \brief  Element-by-element product (Matlab dot-star operator)
   */
  Vector<Scalar> dotStar(const Vector<Scalar>& other) const ;

  /** 
   *  \brief Element-by-element division (Matlab dot-slash operator)
   */
  Vector<Scalar> dotSlash(const Vector<Scalar>& other) const ;

  /** 
   * \brief  Return element-by-element reciprocal as a new vector
   */
  Vector<Scalar> reciprocal() const ;


  /** 
   * \brief  Return element-by-element absolute value as a new vector
   */
  Vector<Scalar> abs() const ;

  /** 
   * \brief  Overwrite self with element-by-element reciprocal
   */
  Vector<Scalar>& selfReciprocal() ;

  /** 
   *  \brief Overwrite self with element-by-element absolute value 
   */
  Vector<Scalar>& selfAbs() ;

  /** 
   *  \brief Overwrite self with random values
   */
  Vector<Scalar>& randomize() ;

  /** 
   *  \brief Set all elements to a constant value
   */
  void setToConstant(const Scalar& alpha) ;

      
  /** 
   *  \brief Take dot product with another vector
   */
  Scalar dot(const Vector<Scalar>& other) const ;

  /** 
   * \brief  Overloaded operator for dot product 
   */
  Scalar operator*(const Vector<Scalar>& other) const ;

  /**
   *  \brief Compute the 1-norm of this vector
   */
  Scalar norm1() const ;

  /**
   *  \brief Compute the 2-norm of this vector
   */
  Scalar norm2() const ;

  /**
   *  \brief Compute the weighted 2-norm of this vector
   */
  Scalar norm2(const Vector<Scalar>& weights) const ;    


  /**
   *  \brief Compute the infinity-norm of this vector
   */
  Scalar normInf() const ;

  /**
   *  \brief Set all elements to zero 
   */
  void zero();


  /**  \brief Return the max element */
  Scalar max() const;

  /**  \brief Return the min element */
  Scalar min()const;


  //@}




  /**  \brief Get a stopwtach for timing vector operations */
  static RCP<Time>& opTimer()
    {
      static RCP<Time> rtn 
        = TimeMonitor::getNewTimer("Low-level vector operations");
      return rtn;
    }

  /** \name Functor application */
  //@{

  /**  \brief Apply a unary functor, overwriting this vector with the results */
  template <class UF> 
  Vector<Scalar>& applyUnaryFunctor(const UF& functor);

  /**  \brief Apply a unary functor to another vector, writing the results
      into this vector. The other vector is unchanged. */
  template <class UF> 
  Vector<Scalar>& acceptUnaryFunctor(const UF& functor,
    const Vector<Scalar>& other);

  /**  \brief Apply a binary functor to this vector and another vector, 
      writing the results
      into this vector. The other vector is unchanged. */
  template <class VF> 
  Vector<Scalar>& applyBinaryFunctor(const VF& functor,
    const Vector<Scalar>& other);

  /** \brief Apply a ternary functor to this vector and two other vectors, 
      writing the results
      into this vector. The other vectors are unchanged. */
  template <class VF> 
  Vector<Scalar>& applyTernaryFunctor(const VF& functor,
    const Vector<Scalar>& x, const Vector<Scalar>& y);

  /** \brief Apply a unary reduction functor */
  template <class RF> 
  typename PlayaFunctors::VectorFunctorTraits<Scalar, RF>::ReturnType 
  applyUnaryReductionFunctor(
    const RF& functor)
    const ;

  /** \brief Apply a binary reduction functor */
  template <class RF> 
  typename PlayaFunctors::VectorFunctorTraits<Scalar, RF>::ReturnType 
  applyBinaryReductionFunctor(const RF& functor, const Vector<Scalar>& other)
    const ;
    

    
  //@}
    
protected:

  /** get a subblock as specified by a deque of indices */
  const Vector<Scalar>& getBlock(const std::deque<int>& iter) const ;

  /** get a non-const subblock as specified by a deque of indices */
  Vector<Scalar> getNonConstBlock(const std::deque<int>& iter) ;
  

private:

};


template <class Scalar> class LoadableVector;
/** \relates Vector \brief Dynamic cast a Vector's underlying pointer to 
a LoadableVector. */
template <class Scalar>
LoadableVector<Scalar>* loadable(Vector<Scalar> vec) ;


/** \relates Vector \brief Return a pointer to the beginning of the vector's
single data chunk. If the underlying VectorBase is not a SingleChunkVector
an exception will be thrown. */
template <class Scalar>
Scalar* dataPtr(Vector<Scalar> vec) ;

/** \relates Vector \brief Return a pointer to the beginning of the vector's
single data chunk. If the underlying VectorBase is not a SingleChunkVector
an exception will be thrown. */
template <class Scalar>
const Scalar* dataPtr(const Vector<Scalar>& vec) ;


}

template <class Scalar> inline
std::ostream& operator<<(std::ostream& os, const Playa::Vector<Scalar>& x) 
{
  x.print(os);
  return os;
}




#endif
