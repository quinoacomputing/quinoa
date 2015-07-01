/*
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
*/

#ifndef RTOPPACK_OLD_TYPES_HPP
#define RTOPPACK_OLD_TYPES_HPP

#include "RTOpPack_Types.hpp"
#include "RTOpPack_SparseSubVectorT.hpp"
#include "RTOp.h"

namespace RTOpPack {

/** \brief Class for a non-mutable sub-vector.
 */
template<class Scalar>
class SubVectorT1B {
public:
  /** \brief . */
  SubVectorT1B() : globalOffset_(0), subDim_(0), values_(Teuchos::null), stride_(0) {}
  /** \brief . */
  SubVectorT1B(Teuchos_Ordinal globalOffset, Teuchos_Ordinal subDim, const Scalar *values, ptrdiff_t stride)
    :globalOffset_(globalOffset), subDim_(subDim), values_(values,0,subDim*stride,false), stride_(stride) 
    {}
  /** \brief . */
  SubVectorT1B( const SubVectorT1B<Scalar>& sv )
    :globalOffset_(sv.globalOffset()), subDim_(sv.subDim()), values_(sv.arcp_values()), stride_(sv.stride()) 
    {}
  /** \brief . */
  SubVectorT1B( const ConstSubVectorView<Scalar>& sv )
    :globalOffset_(sv.globalOffset()), subDim_(sv.subDim()), values_(sv.values()), stride_(sv.stride()) 
    {}
  /** \brief . */
  operator ConstSubVectorView<Scalar>()
    { return ConstSubVectorView<Scalar>(globalOffset(),subDim(),arcp_values(),stride()); }
  /** \brief . */
  void initialize(Teuchos_Ordinal globalOffset, Teuchos_Ordinal subDim, const Scalar *values, ptrdiff_t stride)
    { globalOffset_=globalOffset; subDim_=subDim; values_=Teuchos::arcp(values,0,subDim*stride,false); stride_=stride; }
  /** \brief . */
  void set_uninitialized()
    { globalOffset_ = 0; subDim_=0; values_=Teuchos::null; stride_ = 0; }
  /** \brief . */
  void setGlobalOffset(Teuchos_Ordinal globalOffset) { globalOffset_ = globalOffset; } 
  /** \brief . */
  Teuchos_Ordinal globalOffset() const { return globalOffset_; }
  /** \brief . */
  Teuchos_Ordinal subDim() const { return subDim_; }
  /** \brief . */
  const Scalar* values() const { return values_.get(); }
  /** \brief . */
  const Teuchos::ArrayRCP<const Scalar> arcp_values() const { return values_; }
  /** \brief . */
  ptrdiff_t stride() const { return stride_; }
  /// Zero-based indexing (Preconditions: <tt>values()!=NULL && (0 <= i < subDim())</tt>)
  const Scalar& operator[](Teuchos_Ordinal i) const
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        !( 0 <= i && i < subDim_ ), std::logic_error
        ,"Error, index i="<<i<<" does not fall in the range [0,"<<(subDim_-1)<<"]!"
        );
#endif
      return values_[ stride_*i ];
    }
  /// One-based indexing (Preconditions: <tt>values()!=NULL && (1 <= i <= subDim())</tt>)
  const Scalar& operator()(Teuchos_Ordinal i) const { return (*this)[i-1]; }
private:
  Teuchos_Ordinal globalOffset_;
  Teuchos_Ordinal subDim_;
  Teuchos::ArrayRCP<const Scalar> values_;
  ptrdiff_t stride_;
};

/** \brief Class for a mutable sub-vector.
 */
template<class Scalar>
class MutableSubVectorT1B : public SubVectorT1B<Scalar> {
public:
  /** \brief . */
  MutableSubVectorT1B() {}
  /** \brief . */
  MutableSubVectorT1B(Teuchos_Ordinal globalOffset, Teuchos_Ordinal subDim, Scalar *values, ptrdiff_t stride)
    :SubVectorT1B<Scalar>(globalOffset, subDim, values, stride)
    {}
  /** \brief . */
  MutableSubVectorT1B( const MutableSubVectorT1B<Scalar> & sv)
    :SubVectorT1B<Scalar>(sv)
    {}
  /** \brief . */
  MutableSubVectorT1B( const SubVectorView<Scalar>& sv )
    :SubVectorT1B<Scalar>(ConstSubVectorView<Scalar>(sv))
    {}
  /** \brief . */
  operator SubVectorView<Scalar>()
    { return SubVectorView<Scalar>(this->globalOffset(),this->subDim(),this->arcp_values(),this->stride()); }
  /** \brief . */
  void initialize(Teuchos_Ordinal globalOffset, Teuchos_Ordinal subDim, Scalar *values, ptrdiff_t stride)
    { SubVectorT1B<Scalar>::initialize(globalOffset, subDim, values, stride); }
  /** \brief . */
  void set_uninitialized()
    { SubVectorT1B<Scalar>::set_uninitialized(); }
  /** \brief . */
  Scalar* values() const { return const_cast<Scalar*>(SubVectorT1B<Scalar>::values()); }
  /** \brief . */
  const Teuchos::ArrayRCP<Scalar> arcp_values() const { return Teuchos::arcp_const_cast<Scalar>(SubVectorT1B<Scalar>::arcp_values()); }
  /// Zero-based indexing (Preconditions: <tt>values()!=NULL && (0 <= i < subDim())</tt>)
  Scalar& operator[](Teuchos_Ordinal i) const { return const_cast<Scalar&>(SubVectorT1B<Scalar>::operator[](i)); } // Is range changed in subclass!
  /// One-based indexing (Preconditions: <tt>values()!=NULL && (1 <= i <= subDim())</tt>)
  Scalar& operator()(Teuchos_Ordinal i) const { return (*this)[i-1]; }
};

template<class Scalar>
void assign_entries( const MutableSubVectorT1B<Scalar> *msv, const SubVectorT1B<Scalar> &sv )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(msv==NULL);
  TEUCHOS_TEST_FOR_EXCEPT(msv->subDim() != sv.subDim());
#endif
  for( int i = 1; i <= sv.subDim(); ++i ) {
    (*msv)(i) = sv(i);
  }
}

/** \brief Class for a non-mutable sub-multi-vector (submatrix).
 */
template<class Scalar>
class SubMultiVectorT1B {
public:
  /** \brief . */
  SubMultiVectorT1B()
    :globalOffset_(0), subDim_(0), colOffset_(0), numSubCols_(0)
    ,values_(NULL), leadingDim_(0)
    {}
  /** \brief . */
  SubMultiVectorT1B(
    Teuchos_Ordinal globalOffset, Teuchos_Ordinal subDim
    ,Teuchos_Ordinal colOffset, Teuchos_Ordinal numSubCols
    ,const Scalar *values, Teuchos_Ordinal leadingDim
    )
    :globalOffset_(globalOffset), subDim_(subDim)
    ,colOffset_(colOffset), numSubCols_(numSubCols)
    ,values_(values), leadingDim_(leadingDim)
    {}
  /** \brief . */
  SubMultiVectorT1B( const SubMultiVectorT1B<Scalar>& smv )
    :globalOffset_(smv.globalOffset()), subDim_(smv.subDim())
    ,colOffset_(smv.colOffset()), numSubCols_(smv.numSubCols())
    ,values_(smv.values()), leadingDim_(smv.leadingDim())
    {}
/*
  SubMultiVectorT1B( const ConstSubMultiVectorView<Scalar>& smv )
  :globalOffset_(smv.globalOffset()), subDim_(smv.subDim())
  ,colOffset_(smv.colOffset()), numSubCols_(smv.numSubCols())
  ,values_(smv.values()), leadingDim_(smv.leadingDim())
  {}
  operator ConstSubMultiVectorView<Scalar>()
  { return ConstSubMultiVectorView<Scalar>(globalOffset(),subDim(),colOffset(),numSubCols(),values(),leadingDim()); }
*/
  /** \brief . */
  void initialize(
    Teuchos_Ordinal globalOffset, Teuchos_Ordinal subDim
    ,Teuchos_Ordinal colOffset, Teuchos_Ordinal numSubCols
    ,const Scalar *values, Teuchos_Ordinal leadingDim
    )
    { globalOffset_=globalOffset; subDim_=subDim; colOffset_=colOffset; numSubCols_=numSubCols;
      values_=values; leadingDim_=leadingDim; }
  /** \brief . */
  void set_uninitialized()
    { globalOffset_ = 0; subDim_=0; colOffset_=0, numSubCols_=0; values_=NULL; leadingDim_=0; }
  /** \brief . */
  void setGlobalOffset(Teuchos_Ordinal globalOffset) { globalOffset_ = globalOffset; } 
  /** \brief . */
  Teuchos_Ordinal   globalOffset()   const { return globalOffset_; }
  /** \brief . */
  Teuchos_Ordinal   subDim()         const { return subDim_; }
  /** \brief . */
  Teuchos_Ordinal   colOffset()      const { return colOffset_; }
  /** \brief . */
  Teuchos_Ordinal   numSubCols()     const { return numSubCols_; }
  /** \brief . */
  const Scalar*   values()         const { return values_; }
  /** \brief . */
  Teuchos_Ordinal   leadingDim()     const { return leadingDim_;  }
  /// One-based indexing (Preconditions: <tt>values()!=NULL && (1<=i<=subDim()) && (1<=j<= numSubCols()</tt>)
  const Scalar& operator()(Teuchos_Ordinal i, Teuchos_Ordinal j) const
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        !( 1 <= i && i < subDim_ ), std::logic_error
        ,"Error, index i="<<i<<" does not fall in the range [1,"<<(subDim_-1)<<"]!"
        );
      TEUCHOS_TEST_FOR_EXCEPTION(
        !( 1 <= j && j <= numSubCols_ ), std::logic_error
        ,"Error, index j="<<j<<" does not fall in the range [1,"<<(numSubCols_-1)<<"]!"
        );
#endif
      return values_[ (i-1) + leadingDim_*(j-1) ];
    }
  /// Return a <tt>SubVectorT1B</tt> view of the jth sub-column (Preconditions: <tt>values()!=NULL (1<=j<=numSubCols()</tt>)
  SubVectorT1B<Scalar> col( const Teuchos_Ordinal j ) const
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        !( 1 <= j && j <= numSubCols_ ), std::logic_error
        ,"Error, index j="<<j<<" does not fall in the range [1,"<<(numSubCols_-1)<<"]!"
        );
#endif
      return SubVectorT1B<Scalar>(globalOffset(),subDim(),values()+(j-1)*leadingDim(),1);
    }
private:
  Teuchos_Ordinal     globalOffset_;
  Teuchos_Ordinal     subDim_;
  Teuchos_Ordinal     colOffset_;
  Teuchos_Ordinal     numSubCols_;
  const Scalar        *values_;
  Teuchos_Ordinal     leadingDim_;
};

/** \brief Class for a mutable sub-vector.
 */
template<class Scalar>
class MutableSubMultiVectorT1B : public SubMultiVectorT1B<Scalar> {
public:
  /** \brief . */
  MutableSubMultiVectorT1B() {}
  /** \brief . */
  MutableSubMultiVectorT1B(
    Teuchos_Ordinal globalOffset, Teuchos_Ordinal subDim
    ,Teuchos_Ordinal colOffset, Teuchos_Ordinal numSubCols
    ,const Scalar *values, Teuchos_Ordinal leadingDim
    )
    :SubMultiVectorT1B<Scalar>(globalOffset,subDim,colOffset,numSubCols,values,leadingDim)
    {}
  /** \brief . */
  MutableSubMultiVectorT1B( const MutableSubMultiVectorT1B<Scalar> & smv)
    :SubMultiVectorT1B<Scalar>(smv)
    {}
/*
  MutableSubMultiVectorT1B( const SubMultiVectorView<Scalar>& smv )
  :SubMultiVectorT1B<Scalar>(
  MutableSubMultiVectorT1B<Scalar>(
  smv.globalOffset(),smv.subDim(),smv.colOffset(),smv.numSubCols()
  ,smv.values(),smv.leadingDim()
  )
  )
  {}
  operator SubMultiVectorView<Scalar>()
  { return SubMultiVectorView<Scalar>(this->globalOffset(),this->subDim(),this->colOffset(),this->numSubCols(),this->values(),this->leadingDim()); }
*/
  /** \brief . */
  void initialize(
    Teuchos_Ordinal globalOffset, Teuchos_Ordinal subDim
    ,Teuchos_Ordinal colOffset, Teuchos_Ordinal numSubCols
    ,const Scalar *values, Teuchos_Ordinal leadingDim
    )
    { SubMultiVectorT1B<Scalar>::initialize(globalOffset,subDim,colOffset,numSubCols,values,leadingDim); }
  /** \brief . */
  void set_uninitialized()
    { SubMultiVectorT1B<Scalar>::set_uninitialized(); }
  /** \brief . */
  Scalar* values() const { return const_cast<Scalar*>(SubMultiVectorT1B<Scalar>::values());  }
  /// One-based indexing (Preconditions: <tt>values()!=NULL && (1<=i<= subDim()) && (1<=j<=numSubCols()</tt>)
  Scalar& operator()(Teuchos_Ordinal i, Teuchos_Ordinal j) const
    { return const_cast<Scalar&>(SubMultiVectorT1B<Scalar>::operator()(i,j)); } // Is range checked in subclass
  /// Return a <tt>MutableSubVectorT1B</tt> view of the jth sub-column (Preconditions: <tt>values()!=NULL && (1<=j<=numSubCols()</tt>)
  MutableSubVectorT1B<Scalar> col( const Teuchos_Ordinal j ) const
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(
        !( 1 <= j && j <= this->numSubCols() ), std::logic_error
        ,"Error, index j="<<j<<" does not fall in the range [1,"<<(this->numSubCols())<<"]!"
        );
#endif
      return MutableSubVectorT1B<Scalar>(this->globalOffset(),this->subDim(),values()+(j-1)*this->leadingDim(),1);
    }
};

template<class Scalar>
void assign_entries( const MutableSubMultiVectorT1B<Scalar> *msmv, const SubMultiVectorT1B<Scalar> &smv )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(msmv==NULL);
  TEUCHOS_TEST_FOR_EXCEPT(msmv->subDim() != smv.subDim());
  TEUCHOS_TEST_FOR_EXCEPT(msmv->numSubCols() != smv.numSubCols());
#endif
  for( Teuchos_Ordinal j = 1; j <= smv.numSubCols(); ++j ) {
    for( Teuchos_Ordinal i = 1; i < smv.subDim(); ++i ) {
      (*msmv)(i,j) = smv(i,j);
    }
  }
}

//
// Typedefs
//

/** \brief . */
typedef SubVectorT1B<RTOp_value_type>              SubVector;
/** \brief . */
typedef MutableSubVectorT1B<RTOp_value_type>       MutableSubVector;
/** \brief . */
typedef SparseSubVectorT<RTOp_value_type>          SparseSubVector;
/** \brief . */
typedef SubMultiVectorT1B<RTOp_value_type>         SubMultiVector;
/** \brief . */
typedef MutableSubMultiVectorT1B<RTOp_value_type>  MutableSubMultiVector;
/** \brief . */
typedef RTOpT<RTOp_value_type>                     RTOp;

} // namespace RTOpPack

#endif // RTOPPACK_OLD_TYPES_HPP
