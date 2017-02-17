/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_VectorReducer_hpp_
#define _fei_VectorReducer_hpp_

#include <fei_iosfwd.hpp>
#include <fei_VectorSpace.hpp>
#include <fei_Reducer.hpp>
#include <fei_Vector.hpp>

#undef fei_file
#define fei_file "fei_VectorReducer.hpp"

#include <fei_ErrMacros.hpp>

namespace fei {

  class VectorReducer : public fei::Vector {
  public:

    /** Constructor */
    VectorReducer(fei::SharedPtr<fei::Reducer> reducer,
                  fei::SharedPtr<fei::Vector> target,
	          bool isSolutionVector=false);

    /** Destructor */
    virtual ~VectorReducer();

    /** Query for underlying target vector. */
    fei::SharedPtr<fei::Vector> getTargetVector()
      { return(target_); }

    /** Return a name describing the run-time type
	of this object.
    */
    const char* typeName() const { return(target_->typeName()); }

    /** Update 'this' = b*'this' + a*x
     */
    int update(double a,
	       const fei::Vector* x,
	       double b);

    /** Use data in the underlying non-overlapping decomposition to update
	any shared data in the overlapping decomposition.

	If any data is already held for the shared positions, that data will
	be replaced by the data from the 'owning' processor.
    */
    int scatterToOverlap();

    void setCommSizes() { target_->setCommSizes(); }

    /** Move any shared data from the overlapping decomposition to the
	underlying non-overlapping decomposition.
    */
    int gatherFromOverlap(bool accumulate = true);

    /** Set a specified scalar throughout the vector. */
    int putScalar(double scalar);

    /** Sum values into the vector, adding to any
	that may already exist at the specified indices.
    */
    int sumIn(int numValues, const int* indices, const double* values,
	      int vectorIndex=0);

    /** Copy values into the vector, overwriting any that may already exist
	at the specified indices.
    */
    int copyIn(int numValues, const int* indices, const double* values,
	       int vectorIndex=0);

    /** Obtain the VectorSpace associated with this vector.
     */
    fei::SharedPtr<fei::VectorSpace> getVectorSpace() const
      { return(target_->getVectorSpace()); }

    /** Set the VectorSpace associated with this vector.
     */
    void setVectorSpace(fei::SharedPtr<fei::VectorSpace> vecSpace)
    { target_->setVectorSpace( vecSpace ); }

    /** Sum field data into the vector, adding to any coefficients that may
	already exist at the specified locations.
        If the specified fieldID doesn't exist at one or more of the specified
        IDs, then the corresponding positions in the data array will simply
        not be used.
    */
    int sumInFieldData(int fieldID,
		       int idType,
		       int numIDs,
		       const int* IDs,
		       const double* data,
		       int vectorIndex=0);

    /** Copy field data into the vector, overwriting any coefficients that may
	already exist at the specified locations.
        If the specified fieldID doesn't exist at one or more of the specified
        IDs, then the corresponding positions in the data array will simply
        not be used.
    */
    int copyInFieldData(int fieldID,
			int idType,
			int numIDs,
			const int* IDs,
			const double* data,
			int vectorIndex=0);

    int copyInFieldDataLocalIDs(int fieldID,
			int idType,
			int numIDs,
			const int* localIDs,
			const double* data,
			int vectorIndex=0);

    /** Copy field data out of the vector, into the caller-allocated data
	array.
        If the specified fieldID doesn't exist at one or more of the specified
        IDs, then the corresponding positions in the data array will simply
        not be referenced.
    */
    int copyOutFieldData(int fieldID,
			 int idType,
			 int numIDs,
			 const int* IDs,
			 double* data,
			 int vectorIndex=0);

    int writeToFile(const char* filename,
		    bool matrixMarketFormat=true);

    int writeToStream(FEI_OSTREAM& ostrm,
		      bool matrixMarketFormat=true);

    int copyOut(int numValues,
		const int* indices,
		double* values,
		int vectorIndex=0) const;

  private:
    /** please ignore
     */
    int copyOut_FE(int nodeNumber, int dofOffset, double& value);

    int giveToUnderlyingVector(int numValues,
			       const int* indices,
			       const double* values,
			       bool sumInto=true,
			       int vectorIndex=0);

    int sumIntoFEVector(int blockID,
			int connOffset,
			int numNodes,
			const int* nodeNumbers,
			const int* numIndicesPerNode,
			const double* values);

    fei::SharedPtr<fei::Reducer> reducer_;
    fei::SharedPtr<fei::Vector> target_;
    bool isSolution_;

    int localProc_;
    int numProcs_;
  };//class VectorReducer

} //namespace fei

#endif // _fei_VectorReducer_hpp_

