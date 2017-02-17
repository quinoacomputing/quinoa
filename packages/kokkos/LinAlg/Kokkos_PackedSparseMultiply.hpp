//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef KOKKOS_PACKEDSPARSEMULTIPLY_H
#define KOKKOS_PACKEDSPARSEMULTIPLY_H

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CisMatrix.hpp" 
#include "Kokkos_SparseOperation.hpp" 


namespace Kokkos {

//! PackedSparseMultiply: A reference class for computing sparse matrix multiplication operations.

/*! The PackedSparseMultiply provide basic functionality for computing sparse matrix times vector, or
    sparse matrix times multivector operations.  This class is templated on the ordinal (integer) and scalar (floating
    point) types, so it can compute using any reasonable data type.

  <b>Constructing PackedSparseMultiply objects</b>

  Constructing PackedSparseMultiply objects is a multi-step process.  The basic steps are as follows:
  <ol>
  <li> Create PackedSparseMultiply instance:  The constructor takes no arguments.
  <li> Register the structure of a CisMatrix object using initializeStructure(): 
       We provide this method so that derived implementations can
       take advantage of multiple problems that have the same structure.  In this situation, initializeStructure() would
       be called once and then initializeValues() would be called repeatedly, amortizing the cost of setting up the structure.
       This method may be called only once.
  <li> Register the values of a CisMatrix object using initializeValues(): This method is used to pass values to the
       multiply class.  It can be called repeatedly if multiple matrices have the same structure.
  </ol>

  <b> Counting Floating Point Operations </b>

  Each PackedSparseMultiply object keeps track of the number
  of floating point operations performed using the specified object as the \e this argument
  to the function.  The getFlops() function returns this number as a double precision number.  Using this 
  information, in conjunction with the Time class, one can get accurate  performance
  numbers.  The resetFlops() function resets the floating point counter.

*/    

  template<typename OrdinalType, typename ScalarType>
  class PackedSparseMultiply: public virtual SparseOperation<OrdinalType, ScalarType> {
  public:

    //! @name Constructors/Destructor

    //@{

    //! PackedSparseMultiply constuctor with variable number of indices per row.
    PackedSparseMultiply();
  
    //! Copy constructor.
    PackedSparseMultiply(const PackedSparseMultiply& source);
	
    //! PackedSparseMultiply Destructor
    virtual ~PackedSparseMultiply();
    //@}
    //! @name Abstract CisMatrix Interface Initialization Methods

    //@{
 
    //! Initialize structure of matrix
    /*!
      This interface supports matrices that implement the CisMatrix matrix interface.
      \param A (In)  An instance of a class that implements the CisMatrix.  All necessary information
      about the matrix can be obtained via this interface.
      \param willKeepStructure (In) This argument is unused by this implementation of the BaseSparseMultiply
             class since structure and values will be copied.
      \return Integer error code, set to 0 if successful.
    */
    virtual int initializeStructure(const CisMatrix<OrdinalType, ScalarType>& A, bool willKeepStructure = false);
 
    //! Initialize values of matrix
    /*!
      This interface supports matrices that implement the CisMatrix matrix interface.
      \param A (In)  An instance of a class that implements the CisMatrix.  All necessary information
      about the matrix can be obtained via this interface.
      \param willKeepValues (In) This argument is unused by this implementation of the BaseSparseMultiply
             class since structure and values will be copied.
      \param checkStructure (In) If set to true, the structure of A will be checked against the structure of
      the matrix passed in to the initializeStructure() methods.  This parameter is false by default.

      \return Integer error code, set to 0 if successful, returns - 1 if checkStructure is true and structure is changed.
    */
    virtual int initializeValues(const CisMatrix<OrdinalType, ScalarType>& A, bool willKeepValues = false,
				 bool checkStructure = false);
 
    //@}

    //! @name Computational methods

    //@{
	
    //! Returns the result of a Kokkos_PackedSparseMultiply multiplied by multiple vectors in x, results in y.
    /*! 
      \param x (In) A MultiVector to multiply by.
      \param y (Out) A MultiVector containing results.
      \param transA (In) If true, multiply by the transpose of matrix, otherwise just use matrix.
      \param conjA (In) If true, multiply by the conjugate of matrix values, otherwise just use matrix values.
		
      \return Integer error code, set to 0 if successful.
    */
    virtual int apply(const MultiVector<OrdinalType, ScalarType>& x, MultiVector<OrdinalType, ScalarType>& y, 
		      bool transA = false, bool conjA = false) const;
    //@}
	
    //! @name Operator attribute access methods

    //@{

    //! Returns false for this implementation.
    /*! This implementation will not use the user's copy of the matrix structure.
    */
    virtual bool getCanUseStructure() const {return(false);};

    //! Returns false for this implementation.
    /*! This implementation will not use the user's copy of the matrix values.
    */
    virtual bool getCanUseValues() const {return(false);};

    //! Returns a reference to the most recent CisMatrix that was passed into the \e this object.
    virtual const CisMatrix<OrdinalType, ScalarType> & getMatrix() const {
      if (matrixForValues_==0) return(*matrixForStructure_);
      else return(*matrixForValues_);
    };
		
    //@}
  
  protected:

    void copyEntries();
    void deleteStructureAndValues();

    struct EntryStruct {
      OrdinalType index;
      ScalarType value;
    }; 
    typedef struct EntryStruct Entry;
    
    CisMatrix<OrdinalType, ScalarType> * matrixForStructure_;
    CisMatrix<OrdinalType, ScalarType> * matrixForValues_;

    bool isRowOriented_;
    bool haveStructure_;
    bool haveValues_;
    bool hasUnitDiagonal_;
  
    OrdinalType numRows_;
    OrdinalType numCols_;
    OrdinalType numRC_;
    OrdinalType numEntries_;

    OrdinalType * profile_;
    double costOfMatVec_;
    Entry * allEntries_;
  };

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  PackedSparseMultiply<OrdinalType, ScalarType>::PackedSparseMultiply() 
    : matrixForStructure_(0),
      matrixForValues_(0),
      isRowOriented_(true),
      haveStructure_(false),
      haveValues_(false),
      hasUnitDiagonal_(false),
      numRows_(0),
      numCols_(0),
      numRC_(0),
      numEntries_(0),
      profile_(0),
      costOfMatVec_(0.0),
      allEntries_(0) {
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  PackedSparseMultiply<OrdinalType, ScalarType>::PackedSparseMultiply(const PackedSparseMultiply<OrdinalType, ScalarType> &source) 
    : matrixForStructure_(source.matrixForStructure_),
      matrixForValues_(source.matrixForValues_),
      isRowOriented_(source.isRowOriented_),
      haveStructure_(source.haveStructure_),
      haveValues_(source.haveValues_),
      hasUnitDiagonal_(source.hasUnitDiagonal_),
      numRows_(source.numRows_),
      numCols_(source.numCols_),
      numRC_(source.numRC_),
      numEntries_(source.numEntries_),
      profile_(source.profile_),
      costOfMatVec_(source.costOfMatVec_),
      allEntries_(source.allEntries_) {

    copyEntries();
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void PackedSparseMultiply<OrdinalType, ScalarType>::copyEntries() {

    OrdinalType i;

    if (allEntries_!=0) {
      Entry * tmp_entries = new Entry[numEntries_];
      for (i=0; i< numEntries_; i++) {
	tmp_entries[i].index = allEntries_[i].index;
      	tmp_entries[i].value = allEntries_[i].value;
      }
      allEntries_ = tmp_entries;
    }
    return;
  }
  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  void PackedSparseMultiply<OrdinalType, ScalarType>::deleteStructureAndValues() {


    OrdinalType i;

    if (profile_!=0) {
      delete [] profile_;
      profile_ = 0;
    }

    if (allEntries_!=0) {
      delete [] allEntries_;
      allEntries_ = 0;
    }
    return;
  }
  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  PackedSparseMultiply<OrdinalType, ScalarType>::~PackedSparseMultiply(){

    deleteStructureAndValues();

  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int PackedSparseMultiply<OrdinalType, ScalarType>::initializeStructure(const CisMatrix<OrdinalType, ScalarType>& A,
								       bool willKeepStructure) {


    if (haveStructure_) return(-1); // Can only call this one time!

    matrixForStructure_ = const_cast<CisMatrix<OrdinalType, ScalarType> *> (&A);
    OrdinalType i, j;
    isRowOriented_ = A.getIsRowOriented();
    hasUnitDiagonal_ = A.getHasImplicitUnitDiagonal();
    numRows_ = A.getNumRows();
    numCols_ = A.getNumCols();
    numEntries_ = A.getNumEntries();
    numRC_ = numCols_;
    if (isRowOriented_) numRC_ = numRows_;

    profile_ = new OrdinalType[numRC_];

    OrdinalType numRCEntries;
    OrdinalType * indicesRC;

      
    allEntries_ = new Entry[numEntries_]; // Allocate storage for all entries at once
    
    OrdinalType offset = 0;
    for (i=0; i< numRC_; i++) {
      int ierr = A.getIndices(i, numRCEntries, indicesRC);
      if (ierr<0) return(ierr);
      profile_[i] = numRCEntries;
      Entry * curRC = allEntries_+offset;
      for (j=0; j<numRCEntries; j++) curRC[j].index = indicesRC[j];
      offset += numRCEntries;
    }

    costOfMatVec_ = 2.0 * ((double) numEntries_);
    if (hasUnitDiagonal_) costOfMatVec_ += 2.0 * ((double) numRC_);
    haveStructure_ = true;
    return(0);
  }

  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int PackedSparseMultiply<OrdinalType, ScalarType>::initializeValues(const CisMatrix<OrdinalType, ScalarType>& A, 
							       bool willKeepValues, bool checkStructure) {

    if (!haveStructure_) return(-1); // Must have structure first!

    matrixForValues_ = const_cast<CisMatrix<OrdinalType, ScalarType> *> (&A);
    OrdinalType i, j;

    ScalarType * valuesRC;

    OrdinalType offset = 0;
    for (i=0; i<numRC_; i++) {
      int ierr = A.getValues(i, valuesRC);
      if (ierr<0) return(ierr);
      Entry * curRC = allEntries_+offset;
      OrdinalType numRCEntries = profile_[i];
      for (j=0; j<numRCEntries; j++) curRC[j].value = valuesRC[j];
      offset += numRCEntries;
    }
    haveValues_ = true;
    return(0);
  }


  //==============================================================================
  template<typename OrdinalType, typename ScalarType>
  int PackedSparseMultiply<OrdinalType, ScalarType>::apply(const MultiVector<OrdinalType, ScalarType>& x, 
						    MultiVector<OrdinalType, ScalarType> & y,
						    bool transA, bool conjA) const {
    if (!haveValues_) return(-1); // Can't compute without values!
    if (conjA) return(-2); // Unsupported at this time
    if (x.getNumRows()!=numCols_) return(-3); // Number of cols in A not same as number of rows in x
    if (y.getNumRows()!=numRows_) return(-4); // Number of rows in A not same as number of rows in x
    OrdinalType numVectors = x.getNumCols();
    if (numVectors!=y.getNumCols()) return(-5); // Not the same number of vectors in x and y

    OrdinalType i, j, k, curNumEntries;
    Entry * curEntries = allEntries_;

    OrdinalType * profile = profile_;

    ScalarType ** xpp = x.getValues();
    ScalarType ** ypp = y.getValues();

    if ((isRowOriented_ && !transA) ||
	(!isRowOriented_ && transA)) {
      ScalarType sum = 0;
      for(i = 0; i < numRC_; i++) {
	curNumEntries = *profile++;
	for (k=0; k<numVectors; k++) {
	  ScalarType * xp = xpp[k];
	  ScalarType * yp = ypp[k];
	  if (hasUnitDiagonal_)
	    sum = xp[i];
	  else
	    sum = 0.0;
	  for(j = 0; j < curNumEntries; j++)
	    sum += curEntries[j].value * xp[curEntries[j].index];
	  yp[i] = sum;
	}
	curEntries += curNumEntries;
      }
    }
    else {
      
      for (k=0; k<numVectors; k++) {
	ScalarType * yp = ypp[k];
	if (hasUnitDiagonal_) {
	  ScalarType * xp = xpp[k];
	  for(i = 0; i < numRC_; i++)
	    yp[i] = xp[i]; // Initialize y
	}
	else
	  for(i = 0; i < numRC_; i++)
	    yp[i] = 0.0; // Initialize y
      }
      for(i = 0; i < numRC_; i++) {
	curNumEntries = *profile++;
	for (k=0; k<numVectors; k++) {
	  ScalarType * xp = xpp[k];
	  ScalarType * yp = ypp[k];
	  for(j = 0; j < curNumEntries; j++)
	    yp[curEntries[j].index] += curEntries[j].value * xp[i];
	}
	curEntries += curNumEntries;
      }
    }
    updateFlops(this->costOfMatVec_ * ((double) numVectors));
    return(0);
  }

} // namespace Kokkos
#endif /* KOKKOS_PACKEDSPARSEMULTIPLY_H */
