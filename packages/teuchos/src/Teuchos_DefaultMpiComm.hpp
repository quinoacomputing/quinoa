// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_MPI_COMM_HPP
#define TEUCHOS_MPI_COMM_HPP


#include "Teuchos_Comm.hpp"
#include "Teuchos_CommUtilities.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_OpaqueWrapper.hpp"
#include "Teuchos_MpiReductionOpSetter.hpp"
#include "Teuchos_SerializationTraitsHelpers.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_Assert.hpp"
#include "mpi.h"
#include <iterator>

// This must be defined globally for the whole program!
//#define TEUCHOS_MPI_COMM_DUMP


#ifdef TEUCHOS_MPI_COMM_DUMP
#  include "Teuchos_VerboseObject.hpp"
#endif


namespace Teuchos {

//! Human-readable string version of the given MPI error code.
std::string
mpiErrorCodeToString (const int err);

#ifdef TEUCHOS_MPI_COMM_DUMP
template<typename Ordinal, typename T>
void dumpBuffer(
  const std::string &funcName, const std::string &buffName
  ,const Ordinal bytes, const T buff[]
  )
{
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::OSTab tab(out);
  *out
    << "\n" << funcName << "::" << buffName << ":\n";
  tab.incrTab();
  for( Ordinal i = 0; i < bytes; ++i ) {
    *out << buffName << "[" << i << "] = '" << buff[i] << "'\n";
  }
  *out << "\n";
}
#endif // TEUCHOS_MPI_COMM_DUMP


/// \class MpiCommRequest
/// \brief MPI-specific implementation of \c CommRequest.
///
/// Users would not normally create an instance of this class.  Calls
/// to \c ireceive() and \c isend() return a CommRequest, which for an
/// MPI implementation of \c Comm is an MpiCommRequest.  Users might
/// wish to create an MpiCommRequest directly if they want to
/// encapsulate an MPI_Request returned by an external library or by
/// their own code, and pass it into one of our wrapper functions like
/// \c wait() or \c waitAll().
class MpiCommRequest : public CommRequest {
public:
  //! Constructor (from a raw MPI_Request).
  MpiCommRequest (MPI_Request rawMpiRequest, 
		  const ArrayView<char>::size_type numBytesInMessage) : 
    rawMpiRequest_ (rawMpiRequest), numBytes_ (numBytesInMessage) 
  {}

  /// \brief Return and relinquish ownership of the raw MPI_Request.
  ///
  /// "Relinquish ownership" means that this object sets its raw
  /// MPI_Request to MPI_REQUEST_NULL, but returns the original
  /// MPI_Request.  This effectively gives the caller ownership of the
  /// raw MPI_Request.  This prevents hanging requests.
  MPI_Request releaseRawMpiRequest()
  {
    MPI_Request tmp_rawMpiRequest = rawMpiRequest_;
    rawMpiRequest_ = MPI_REQUEST_NULL;
    return tmp_rawMpiRequest;
  }
  
  //! Whether the raw MPI_Request is MPI_REQUEST_NULL.
  bool isNull() const {
    return rawMpiRequest_ == MPI_REQUEST_NULL;
  }

  /// \brief Number of bytes in the nonblocking send or receive request.
  ///
  /// Remembering this is inexpensive, and is also useful for
  /// debugging (e.g., for detecting whether the send and receive have
  /// matching message lengths).
  ArrayView<char>::size_type numBytes () const {
    return numBytes_;
  }

private:
  //! The raw request (an opaque object).
  MPI_Request rawMpiRequest_;
  //! Number of bytes in the nonblocking send or receive request.
  ArrayView<char>::size_type numBytes_;

  MpiCommRequest(); // Not defined
};

/// \fn mpiCommRequest
/// \brief Nonmember constructor for \c MpiCommRequest.
/// \relates MpiCommRequest
///
/// \param rawMpiRequest [in] The raw MPI_Request opaque object.
/// \param numBytes [in] The number of bytes in the nonblocking
///   send or receive request.
inline RCP<MpiCommRequest>
mpiCommRequest (MPI_Request rawMpiRequest, 
		const ArrayView<char>::size_type numBytes)
{
  return rcp (new MpiCommRequest (rawMpiRequest, numBytes));
}

/// \class MpiCommStatus
/// \brief MPI-specific implementation of \c CommStatus.
///
/// Users would not normally create an instance of this class.  The
/// only time they might wish to do so is to encapsulate an MPI_Status
/// returned by an external library or by their own code, and pass it
/// into one of our functions like \c wait() or \c waitAll().
///
/// \tparam OrdinalType The same template parameter as \c Comm.  Only
///   use \c int here.  We only make this a template class for
///   compatibility with \c Comm.
template<class OrdinalType>
class MpiCommStatus : public CommStatus<OrdinalType> {
public:
  MpiCommStatus (MPI_Status status) : status_ (status) {}

  //! Destructor (declared virtual for memory safety)
  virtual ~MpiCommStatus() {}

  //! The source rank that sent the message.
  OrdinalType getSourceRank () { return status_.MPI_SOURCE; }

  //! The tag of the received message.
  OrdinalType getTag () { return status_.MPI_TAG; }

  //! The error code of the received message.
  OrdinalType getError () { return status_.MPI_ERROR; }

private:
  //! We forbid default construction syntactically.
  MpiCommStatus ();

  //! The raw MPI_Status struct that this class encapsulates.
  MPI_Status status_;
};

/// \fn mpiCommStatus
/// \brief Nonmember constructor for \c MpiCommStatus.
/// \relates MpiCommStatus
template<class OrdinalType>
inline RCP<MpiCommStatus<OrdinalType> >
mpiCommStatus (MPI_Status rawMpiStatus)
{
  return rcp (new MpiCommStatus<OrdinalType> (rawMpiStatus));
}

/** \brief Concrete communicator subclass based on MPI.
 *
 * <b>Assertions:</b><ul>
 * <li><tt>getRawMpiComm().get()!=NULL && *getRawMpiComm()!=MPI_COMM_NULL</tt>
 * <li><tt>getSize() > 0</tt>
 * <li><tt>0 <= getRank() && getRank() < getSize()</tt>
 * </ul>
 *
 * ToDo: Finish documentation!
 */
template<typename Ordinal>
class MpiComm : public Comm<Ordinal> {
public:

  //! @name Constructors 
  //@{

  /** \brief Construct given a wrapped MPI_Comm oqaque object.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>rawMpiComm.get()!=NULL && *rawMpiComm != MPI_COMM_NULL</tt>
   * </ul>
   */
  MpiComm(
    const RCP<const OpaqueWrapper<MPI_Comm> > &rawMpiComm
    );

  /**
   * \brief Construct a communicator with a new context with the same properties
   * as the original.
   *
   * The newly constructed communicator will have a duplicate communication
   * space that has the same properties (e.g. processes, attributes,
   * topologies) as the input communicator.
   *
   * \param other The communicator to copy from.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>
   * other.getRawMpiComm().get() != NULL && *other.getRawMpiComm() != NULL
   * </tt></li>
   * </ul>
   */
  MpiComm(const MpiComm<Ordinal>& other);

  /** \brief Return the embedded wrapped opaque <tt>MPI_Comm</tt> object. */
  RCP<const OpaqueWrapper<MPI_Comm> > getRawMpiComm() const
  {return rawMpiComm_;}

  /// \brief Set the MPI error handler for this communicator.
  ///
  /// MPI lets you set an error handler function specific to each
  /// communicator.  MpiComm wraps this functionality.  Create an
  /// error handler using \c MPI_Errhandler_create(), or use one of
  /// the default error handlers that the MPI standard or your MPI
  /// implementation provides.  You will need to wrap the MPI error
  /// handler in an OpaqueWrapper first.  (See the documentation of
  /// OpaqueWrapper for the rationale behind not using MPI's opaque
  /// objects directly.)
  ///
  /// MpiComm will not attempt to call MPI_Errhandler_free() on the
  /// error handler you provide.  You can always set the RCP's custom
  /// "deallocator" function to free the error handler, if you want it
  /// taken care of automatically.
  ///
  /// \param errHandler [in] The error handler to set.  If null, do
  ///   nothing.
  void setErrorHandler (const RCP<const OpaqueWrapper<MPI_Errhandler> >& errHandler);

  //@}

  //! @name Overridden from Comm 
  //@{

  /** \brief . */
  virtual int getRank() const;
  /** \brief . */
  virtual int getSize() const;
  /** \brief . */
  virtual void barrier() const;
  /** \brief . */
  virtual void broadcast(
    const int rootRank, const Ordinal bytes, char buffer[]
    ) const;
  /** \brief . */
  virtual void gatherAll(
    const Ordinal sendBytes, const char sendBuffer[]
    ,const Ordinal recvBytes, char recvBuffer[]
    ) const;
  /** \brief . */
  virtual void reduceAll(
    const ValueTypeReductionOp<Ordinal,char> &reductOp
    ,const Ordinal bytes, const char sendBuffer[], char globalReducts[]
    ) const;
  /** \brief . */
  virtual void reduceAllAndScatter(
    const ValueTypeReductionOp<Ordinal,char> &reductOp
    ,const Ordinal sendBytes, const char sendBuffer[]
    ,const Ordinal recvCounts[], char myGlobalReducts[]
    ) const;
  /** \brief . */
	virtual void scan(
    const ValueTypeReductionOp<Ordinal,char> &reductOp
    ,const Ordinal bytes, const char sendBuffer[], char scanReducts[]
    ) const;
  /** \brief . */
  virtual void send(
    const Ordinal bytes, const char sendBuffer[], const int destRank
    ) const;
  /** \brief . */
  virtual void ssend(
    const Ordinal bytes, const char sendBuffer[], const int destRank
    ) const;
  /** \brief . */
  virtual int receive(
    const int sourceRank, const Ordinal bytes, char recvBuffer[]
    ) const;
  /** \brief . */
  virtual void readySend(
    const ArrayView<const char> &sendBuffer,
    const int destRank
    ) const;
  /** \brief . */
  virtual RCP<CommRequest> isend(
    const ArrayView<const char> &sendBuffer,
    const int destRank
    ) const;
  /** \brief . */
  virtual RCP<CommRequest> ireceive(
    const ArrayView<char> &Buffer,
    const int sourceRank
    ) const;
  /** \brief . */
  virtual void waitAll(
    const ArrayView<RCP<CommRequest> > &requests
    ) const;
  /** \brief . */
  virtual void 
  waitAll (const ArrayView<RCP<CommRequest> >& requests,
	   const ArrayView<RCP<CommStatus<Ordinal> > >& statuses) const;
  /** \brief . */
  virtual RCP<CommStatus<Ordinal> > 
  wait (const Ptr<RCP<CommRequest> >& request) const;
  /** \brief . */
  virtual RCP< Comm<Ordinal> > duplicate() const;
  /** \brief . */
  virtual RCP< Comm<Ordinal> > split(const int color, const int key) const;
  /** \brief . */
  virtual RCP< Comm<Ordinal> > createSubcommunicator(
    const ArrayView<const int>& ranks) const;
  //@}

  //! @name Overridden from Describable 
  //@{

  /** \brief . */
  std::string description() const;

  //@}

  // These should be private but the PGI compiler requires them be public

  static int const minTag_ = 26000; // These came from Teuchos::MpiComm???
  static int const maxTag_ = 26099; // ""

private:

  // Set internal data members once the rawMpiComm_ data member is valid.
  void setupMembersFromComm();
  static int tagCounter_;

  RCP<const OpaqueWrapper<MPI_Comm> > rawMpiComm_;
  int rank_;
  int size_;
  int tag_;

  //! MPI error handler.  If null, MPI uses the default error handler.
  RCP<const OpaqueWrapper<MPI_Errhandler> > customErrorHandler_;

  void assertRank(const int rank, const std::string &rankName) const;

  // Not defined and not to be called!
  MpiComm();

#ifdef TEUCHOS_MPI_COMM_DUMP
public:
  static bool show_dump;
#endif // TEUCHOS_MPI_COMM_DUMP
	
};


/** \brief Helper function that creates a dynamically allocated
 * <tt>MpiComm</tt> object or returns <tt>Teuchos::null</tt> to correctly
 * represent a null communicator.
 *
 * <b>Postconditions:</b></ul>
 * <li>[<tt>rawMpiComm.get()!=NULL && *rawMpiComm!=MPI_COMM_NULL</tt>]
 *     <tt>return.get()!=NULL</tt>
 * <li>[<tt>rawMpiComm.get()==NULL || *rawMpiComm==MPI_COMM_NULL</tt>]
 *     <tt>return.get()==NULL</tt>
 * </ul>
 *
 * \relates MpiComm
 */
template<typename Ordinal>
RCP<MpiComm<Ordinal> >
createMpiComm(
  const RCP<const OpaqueWrapper<MPI_Comm> > &rawMpiComm
  );


// ////////////////////////
// Implementations


// Static members


template<typename Ordinal>
int MpiComm<Ordinal>::tagCounter_ = MpiComm<Ordinal>::minTag_;


// Constructors


template<typename Ordinal>
MpiComm<Ordinal>::MpiComm(
  const RCP<const OpaqueWrapper<MPI_Comm> > &rawMpiComm
  )
{
  TEUCHOS_TEST_FOR_EXCEPTION(rawMpiComm.get()==NULL, std::invalid_argument, 
    "Teuchos::MpiComm constructor: The input RCP is null.");
  TEUCHOS_TEST_FOR_EXCEPTION(*rawMpiComm == MPI_COMM_NULL, 
    std::invalid_argument, "Teuchos::MpiComm constructor: The given MPI_Comm "
    "is MPI_COMM_NULL.");
  rawMpiComm_ = rawMpiComm;

  // FIXME (mfh 26 Mar 2012) The following is a bit wicked in that it
  // changes the behavior of existing applications that use MpiComm,
  // without warning.  I've chosen to do it because I can't figure out
  // any other way to help me debug MPI_Waitall failures on some (but
  // not all) of the testing platforms.  The problem is that MPI's
  // default error handler is MPI_ERRORS_ARE_FATAL, which immediately
  // aborts on error without returning an error code from the MPI
  // function.  Also, the testing platforms' MPI implementations'
  // diagnostics are not giving me useful information.  Thus, I'm
  // setting the default error handler to MPI_ERRORS_RETURN, so that
  // MPI_Waitall will return an error code.
  //
  // Note that all MpiComm methods check error codes returned by MPI
  // functions, and throw an exception if the code is not MPI_SUCCESS.
  // Thus, this change in behavior will only affect your program in
  // the following case: You call a function f() in the try block of a
  // try-catch, and expect f() to throw an exception (generally
  // std::runtime_error) in a particular case not related to MpiComm,
  // but MpiComm throws the exception instead.  It's probably a bad
  // idea for you to do this, because MpiComm might very well throw
  // exceptions for things like invalid arguments.
  const bool makeMpiErrorsReturn = true;
  if (makeMpiErrorsReturn) {
    RCP<const OpaqueWrapper<MPI_Errhandler> > errHandler = 
      rcp (new OpaqueWrapper<MPI_Errhandler> (MPI_ERRORS_RETURN));
    setErrorHandler (errHandler);
  }

  setupMembersFromComm();
}


template<typename Ordinal>
MpiComm<Ordinal>::MpiComm(const MpiComm<Ordinal>& other)
{
  TEUCHOS_TEST_FOR_EXCEPT(other.getRawMpiComm().get() == NULL);
  TEUCHOS_TEST_FOR_EXCEPT(*other.getRawMpiComm() == MPI_COMM_NULL);
  MPI_Comm newComm;
  const int err = MPI_Comm_dup (*other.getRawMpiComm(), &newComm);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm copy constructor: MPI_Comm_dup failed with error \"" 
    << mpiErrorCodeToString (err) << "\".");
  rawMpiComm_ = opaqueWrapper (newComm);
  setupMembersFromComm();
}


template<typename Ordinal>
void MpiComm<Ordinal>::setupMembersFromComm()
{
  MPI_Comm_size(*rawMpiComm_, &size_);
  MPI_Comm_rank(*rawMpiComm_, &rank_);
  if(tagCounter_ > maxTag_)
    tagCounter_ = minTag_;
  tag_ = tagCounter_++;
}


template<typename Ordinal>
void 
MpiComm<Ordinal>::
setErrorHandler (const RCP<const OpaqueWrapper<MPI_Errhandler> >& errHandler)
{
  if (! is_null (errHandler)) {
    const int err = MPI_Comm_set_errhandler (*getRawMpiComm(), *errHandler);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
      "Teuchos::MpiComm::setErrorHandler: MPI_Comm_set_errhandler() failed with "
      "error \"" << mpiErrorCodeToString (err) << "\".");
  }
  // Wait to set this until the end, in case MPI_Errhandler_set()
  // doesn't succeed.
  customErrorHandler_ = errHandler;
}



// Overridden from Comm

  
template<typename Ordinal>
int MpiComm<Ordinal>::getRank() const
{
  return rank_;
}

  
template<typename Ordinal>
int MpiComm<Ordinal>::getSize() const
{
  return size_;
}

  
template<typename Ordinal>
void MpiComm<Ordinal>::barrier() const
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::barrier()"
    );
  const int err = MPI_Barrier(*rawMpiComm_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::barrier: MPI_Barrier failed with error \"" 
    << mpiErrorCodeToString (err) << "\".");
}

  
template<typename Ordinal>
void MpiComm<Ordinal>::broadcast(
  const int rootRank, const Ordinal bytes, char buffer[]
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::broadcast(...)"
    );
  const int err = MPI_Bcast(buffer,bytes,MPI_CHAR,rootRank,*rawMpiComm_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::broadcast: MPI_Bcast failed with error \"" 
    << mpiErrorCodeToString (err) << "\".");
}

  
template<typename Ordinal>
void MpiComm<Ordinal>::gatherAll(
  const Ordinal sendBytes, const char sendBuffer[],
  const Ordinal recvBytes, char recvBuffer[]
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::gatherAll(...)"
    );
  TEUCHOS_ASSERT_EQUALITY((sendBytes*size_), recvBytes );
  const int err = 
    MPI_Allgather (const_cast<char *>(sendBuffer), sendBytes, MPI_CHAR,
		   recvBuffer, sendBytes, MPI_CHAR, *rawMpiComm_);
  // NOTE: 'sendBytes' is being sent above for the MPI arg recvcount (which is
  // very confusing in the MPI documentation) for MPI_Allgether(...).

  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::gatherAll: MPI_Allgather failed with error \"" 
    << mpiErrorCodeToString (err) << "\".");
}

  
template<typename Ordinal>
void 
MpiComm<Ordinal>::
reduceAll (const ValueTypeReductionOp<Ordinal,char> &reductOp,
	   const Ordinal bytes, 
	   const char sendBuffer[], 
	   char globalReducts[]) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::reduceAll(...)" );

  MpiReductionOpSetter op(mpiReductionOp(rcp(&reductOp,false)));
  MPI_Datatype char_block;

  // TODO (mfh 26 Mar 2012) Check returned error codes of the MPI
  // custom datatype functions.
  MPI_Type_contiguous(bytes, MPI_CHAR, &char_block);
  MPI_Type_commit(&char_block);

  const int err = 
    MPI_Allreduce (const_cast<char*>(sendBuffer), globalReducts, 1, char_block,
		   op.mpi_op(), *rawMpiComm_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::reduceAll (custom op): MPI_Allreduce failed with error \"" 
    << mpiErrorCodeToString (err) << "\".");

  // TODO (mfh 26 Mar 2012) Check returned error codes of the MPI
  // custom datatype functions.
  MPI_Type_free(&char_block);
}


template<typename Ordinal>
void MpiComm<Ordinal>::reduceAllAndScatter(
  const ValueTypeReductionOp<Ordinal,char> &reductOp
  ,const Ordinal sendBytes, const char sendBuffer[]
  ,const Ordinal recvCounts[], char myGlobalReducts[]
  ) const
{

  (void)sendBytes; // Ignore if not in debug mode

  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::reduceAllAndScatter(...)"
    );

#ifdef TEUCHOS_DEBUG
  Ordinal sumRecvBytes = 0;
  for( Ordinal i = 0; i < size_; ++i ) {
    sumRecvBytes += recvCounts[i];
  }
  TEUCHOS_TEST_FOR_EXCEPT(!(sumRecvBytes==sendBytes));
#endif // TEUCHOS_DEBUG

#ifdef TEUCHOS_MPI_COMM_DUMP
  if(show_dump) {
    dumpBuffer<Ordinal,char>(
      "Teuchos::MpiComm<Ordinal>::reduceAllAndScatter(...)",
      "sendBuffer", sendBytes, sendBuffer );
    dumpBuffer<Ordinal,Ordinal>(
      "Teuchos::MpiComm<Ordinal>::reduceAllAndScatter(...)",
      "recvCounts", as<Ordinal>(size_), recvCounts );
    dumpBuffer<Ordinal,char>(
      "Teuchos::MpiComm<Ordinal>::reduceAllAndScatter(...)",
      "myGlobalReducts", as<char>(recvCounts[rank_]), myGlobalReducts );
  }
#endif // TEUCHOS_MPI_COMM_DUMP

  // Create a new recvCount[] if Ordinal!=int
  WorkspaceStore* wss = get_default_workspace_store().get();
  const bool Ordinal_is_int = typeid(int)==typeid(Ordinal);
  Workspace<int> ws_int_recvCounts(wss,Ordinal_is_int?0:size_);
  const int *int_recvCounts = 0;
  if(Ordinal_is_int) {
    int_recvCounts = reinterpret_cast<const int*>(recvCounts);
    // Note: We must do an reinterpet cast since this must
    // compile even if it is not executed.  I could implement
    // code that would not need to do this using template
    // conditionals but I don't want to bother.
  }
  else {
    std::copy(recvCounts, recvCounts+size_, &ws_int_recvCounts[0]);
    int_recvCounts = &ws_int_recvCounts[0];
  }

  // Perform the operation
  MpiReductionOpSetter op(mpiReductionOp(rcp(&reductOp, false)));

  const int err = MPI_Reduce_scatter(
    const_cast<char*>(sendBuffer), myGlobalReducts,
    const_cast<int*>(int_recvCounts),
    MPI_CHAR,
    op.mpi_op(),
    *rawMpiComm_
    );
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::reduceAllAndScatter: MPI_Reduce_scatter failed with "
    "error \"" << mpiErrorCodeToString (err) << "\".");
}


template<typename Ordinal>
void MpiComm<Ordinal>::scan(
  const ValueTypeReductionOp<Ordinal,char> &reductOp
  ,const Ordinal bytes, const char sendBuffer[], char scanReducts[]
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::scan(...)" );

  MpiReductionOpSetter op(mpiReductionOp(rcp(&reductOp,false)));
  const int err = 
    MPI_Scan (const_cast<char*>(sendBuffer), scanReducts, bytes, MPI_CHAR, 
	      op.mpi_op(), *rawMpiComm_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::scan: MPI_Scan() failed with error \"" 
    << mpiErrorCodeToString (err) << "\".");
}


template<typename Ordinal>
void 
MpiComm<Ordinal>::send (const Ordinal bytes, 
			const char sendBuffer[], 
			const int destRank) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::send(...)" );

#ifdef TEUCHOS_MPI_COMM_DUMP
  if(show_dump) {
    dumpBuffer<Ordinal,char>(
      "Teuchos::MpiComm<Ordinal>::send(...)"
      ,"sendBuffer", bytes, sendBuffer
      );
  }
#endif // TEUCHOS_MPI_COMM_DUMP

  const int err = MPI_Send (const_cast<char*>(sendBuffer), bytes, MPI_CHAR,
			    destRank, tag_, *rawMpiComm_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::send: MPI_Send() failed with error \"" 
    << mpiErrorCodeToString (err) << "\".");
}


template<typename Ordinal>
void 
MpiComm<Ordinal>::ssend (const Ordinal bytes, 
			 const char sendBuffer[], 
			 const int destRank) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::ssend(...)" );

#ifdef TEUCHOS_MPI_COMM_DUMP
  if(show_dump) {
    dumpBuffer<Ordinal,char>(
      "Teuchos::MpiComm<Ordinal>::send(...)"
      ,"sendBuffer", bytes, sendBuffer
      );
  }
#endif // TEUCHOS_MPI_COMM_DUMP

  const int err = MPI_Ssend (const_cast<char*>(sendBuffer), bytes, MPI_CHAR,
			     destRank, tag_, *rawMpiComm_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::send: MPI_Ssend() failed with error \"" 
    << mpiErrorCodeToString (err) << "\".");
}


template<typename Ordinal>
void MpiComm<Ordinal>::readySend(
  const ArrayView<const char> &sendBuffer,
  const int destRank
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::readySend" );

#ifdef TEUCHOS_MPI_COMM_DUMP
  if(show_dump) {
    dumpBuffer<Ordinal,char>(
      "Teuchos::MpiComm<Ordinal>::readySend(...)"
      ,"sendBuffer", bytes, sendBuffer
      );
  }
#endif // TEUCHOS_MPI_COMM_DUMP

  const int err = 
    MPI_Rsend (const_cast<char*>(sendBuffer.getRawPtr()), sendBuffer.size(),
	       MPI_CHAR, destRank, tag_, *rawMpiComm_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::readySend: MPI_Rsend() failed with error \"" 
    << mpiErrorCodeToString (err) << "\".");
}


template<typename Ordinal>
int 
MpiComm<Ordinal>::receive (const int sourceRank, 
			   const Ordinal bytes, 
			   char recvBuffer[]) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::receive(...)" );

  // A negative source rank indicates MPI_ANY_SOURCE, namely that we
  // will take an incoming message from any process, as long as the
  // tag matches.
  const int theSrcRank = (sourceRank < 0) ? MPI_ANY_SOURCE : sourceRank;

  MPI_Status status;
  const int err = MPI_Recv (recvBuffer, bytes, MPI_CHAR, theSrcRank, tag_, 
			    *rawMpiComm_, &status);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::receive: MPI_Recv() failed with error \"" 
    << mpiErrorCodeToString (err) << "\".");

#ifdef TEUCHOS_MPI_COMM_DUMP
  if (show_dump) {
    dumpBuffer<Ordinal,char> ("Teuchos::MpiComm<Ordinal>::receive(...)", 
			      "recvBuffer", bytes, recvBuffer);
  }
#endif // TEUCHOS_MPI_COMM_DUMP

  // Returning the source rank is useful in the MPI_ANY_SOURCE case.
  return status.MPI_SOURCE;
}


template<typename Ordinal>
RCP<CommRequest> 
MpiComm<Ordinal>::isend (const ArrayView<const char> &sendBuffer,
			 const int destRank) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::isend(...)" );

  MPI_Request rawMpiRequest = MPI_REQUEST_NULL;
  const int err = 
    MPI_Isend (const_cast<char*>(sendBuffer.getRawPtr()), sendBuffer.size(), 
	       MPI_CHAR, destRank, tag_, *rawMpiComm_, &rawMpiRequest);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::isend: MPI_Isend() failed with error \"" 
    << mpiErrorCodeToString (err) << "\".");

  return mpiCommRequest (rawMpiRequest, sendBuffer.size());
}


template<typename Ordinal>
RCP<CommRequest> 
MpiComm<Ordinal>::ireceive (const ArrayView<char> &recvBuffer,
			    const int sourceRank) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::ireceive(...)" );

  // A negative source rank indicates MPI_ANY_SOURCE, namely that we
  // will take an incoming message from any process, as long as the
  // tag matches.
  const int theSrcRank = (sourceRank < 0) ? MPI_ANY_SOURCE : sourceRank;

  MPI_Request rawMpiRequest = MPI_REQUEST_NULL;
  const int err = 
    MPI_Irecv (const_cast<char*>(recvBuffer.getRawPtr()), recvBuffer.size(), 
	       MPI_CHAR, theSrcRank, tag_, *rawMpiComm_, &rawMpiRequest);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::ireceive: MPI_Irecv() failed with error \"" 
    << mpiErrorCodeToString (err) << "\".");

  return mpiCommRequest (rawMpiRequest, recvBuffer.size());
}


namespace {
  // Called by both MpiComm::waitAll() implementations.
  template<typename Ordinal>
  void 
  waitAllImpl (const ArrayView<RCP<CommRequest> >& requests,
	       const ArrayView<MPI_Status>& rawMpiStatuses)
  {
    typedef ArrayView<RCP<CommRequest> >::size_type size_type;
    const size_type count = requests.size();
    // waitAllImpl() is not meant to be called by users, so it's a bug
    // for the two views to have different lengths.
    TEUCHOS_TEST_FOR_EXCEPTION(rawMpiStatuses.size() != count, 
      std::logic_error, "Teuchos::MpiComm's waitAllImpl: rawMpiStatus.size() = "
      << rawMpiStatuses.size() << " != requests.size() = " << requests.size() 
      << ".  Please report this bug to the Tpetra developers.");
    if (count == 0) {
      return; // No requests on which to wait
    }

    // MpiComm wraps MPI and can't expose any MPI structs or opaque
    // objects.  Thus, we have to unpack requests into a separate array.
    // If that's too slow, then your code should just call into MPI
    // directly.
    //
    // Pull out the raw MPI requests from the wrapped requests.
    // MPI_Waitall should not fail if a request is MPI_REQUEST_NULL, but
    // we keep track just to inform the user.
    bool someNullRequests = false;
    Array<MPI_Request> rawMpiRequests (count, MPI_REQUEST_NULL);
    for (int i = 0; i < count; ++i) {
      RCP<CommRequest> request = requests[i];
      if (! is_null (request)) {
	RCP<MpiCommRequest> mpiRequest = 
	  rcp_dynamic_cast<MpiCommRequest> (request);
	// releaseRawMpiRequest() sets the MpiCommRequest's raw
	// MPI_Request to MPI_REQUEST_NULL.  This makes waitAll() not
	// satisfy the strong exception guarantee.  That's OK because
	// MPI_Waitall() doesn't promise that it satisfies the strong
	// exception guarantee, and we would rather conservatively
	// invalidate the handles than leave dangling requests around
	// and risk users trying to wait on the same request twice.
	rawMpiRequests[i] = mpiRequest->releaseRawMpiRequest();
      }
      else { // Null requests map to MPI_REQUEST_NULL
	rawMpiRequests[i] = MPI_REQUEST_NULL;
	someNullRequests = true;
      }
    }

    // This is the part where we've finally peeled off the wrapper and
    // we can now interact with MPI directly.  
    //
    // One option in the one-argument version of waitAll() is to ignore
    // the statuses completely.  MPI lets you pass in the named constant
    // MPI_STATUSES_IGNORE for the MPI_Status array output argument in
    // MPI_Waitall(), which would tell MPI not to bother with the
    // statuses.  However, we want the statuses because we can use them
    // for detailed error diagnostics in case something goes wrong.
    const int err = MPI_Waitall (count, rawMpiRequests.getRawPtr(), 
				 rawMpiStatuses.getRawPtr());

    // In MPI_Waitall(), an error indicates that one or more requests
    // failed.  In that case, there could be requests that completed
    // (their MPI_Status' error field is MPI_SUCCESS), and other
    // requests that have not completed yet but have not necessarily
    // failed (MPI_PENDING).  We make no attempt here to wait on the
    // pending requests.  It doesn't make sense for us to do so, because
    // in general Teuchos::Comm doesn't attempt to provide robust
    // recovery from failed messages.
    if (err != MPI_SUCCESS) {
      if (err == MPI_ERR_IN_STATUS) {
	//
	// When MPI_Waitall returns MPI_ERR_IN_STATUS (a standard error
	// class), it's telling us to check the error codes in the
	// returned statuses.  In that case, we do so and generate a
	// detailed exception message.
	//
	// Figure out which of the requests failed.
	Array<std::pair<size_type, int> > errorLocationsAndCodes;
	for (size_type k = 0; k < rawMpiStatuses.size(); ++k) {
	  const int curErr = rawMpiStatuses[k].MPI_ERROR;
	  if (curErr != MPI_SUCCESS) {
	    errorLocationsAndCodes.push_back (std::make_pair (k, curErr));
	  }
	}
	const size_type numErrs = errorLocationsAndCodes.size();
	if (numErrs > 0) {
	  // There was at least one error.  Assemble a detailed
	  // exception message reporting which requests failed,
	  // their error codes, and their source 
	  std::ostringstream os;
	  os << "Teuchos::MpiComm::waitAll: MPI_Waitall() failed with error \""
	     << mpiErrorCodeToString (err) << "\".  Of the " << count 
	     << " total request" << (count != 1 ? "s" : "") << ", " << numErrs 
	     << " failed.  Here are the indices of the failed requests, and the "
	    "error codes extracted from their returned MPI_Status objects:" 
	     << std::endl;
	  for (size_type k = 0; k < numErrs; ++k) {
	    const size_type errInd = errorLocationsAndCodes[k].first;
	    os << "Request " << errInd << ": MPI_ERROR = " 
	       << mpiErrorCodeToString (rawMpiStatuses[errInd].MPI_ERROR) 
	       << std::endl;
	  }
	  if (someNullRequests) {
	    os << "  On input to MPI_Waitall, there was at least one MPI_"
	      "Request that was MPI_REQUEST_NULL.  MPI_Waitall should not "
	      "normally fail in that case, but we thought we should let you know "
	      "regardless.";
	  }
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, os.str());
	}
	// If there were no actual errors in the returned statuses,
	// well, then I guess everything is OK.  Just keep going.
      }
      else {
	std::ostringstream os;
	os << "Teuchos::MpiComm::waitAll: MPI_Waitall() failed with error \""
	   << mpiErrorCodeToString (err) << "\".";
	if (someNullRequests) {
	  os << "  On input to MPI_Waitall, there was at least one MPI_Request "
	    "that was MPI_REQUEST_NULL.  MPI_Waitall should not normally fail in "
	    "that case, but we thought we should let you know regardless.";
	}
	TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, os.str());
      }
    }

#ifdef HAVE_TEUCHOS_DEBUG
    if (false) // mfh 02 Apr 2012: This test fails in some cases (e.g., Belos BlockCG), with the MPI_Request reporting 8 bytes and the MPI_Status reporting 0 bytes.  The tests pass otherwise, so I'm disabling this check for now.
    {
      // In debug mode, test whether the requests' message lengths
      // matched the message lengths on completion.
      Array<size_type> nonmatchingIndices;
      Array<std::pair<size_type, size_type> > nonmatchingLengthPairs;
      for (size_type k = 0; k < count; ++k) {
	if (! is_null (requests[k])) {
	  RCP<MpiCommRequest> mpiRequest = 
	    rcp_dynamic_cast<MpiCommRequest> (requests[k]);
	  
	  int statusCount = -1;
	  (void) MPI_Get_count (&rawMpiStatuses[k], MPI_CHAR, &statusCount);
	  if (mpiRequest->numBytes() != as<size_type> (statusCount)) {
	    nonmatchingIndices.push_back (k);
            nonmatchingLengthPairs.push_back (std::make_pair (mpiRequest->numBytes(), Teuchos::as<size_type> (statusCount)));
	  }
	}
      }
      const size_type numNonmatching = nonmatchingIndices.size();
      if (numNonmatching > 0) {
	std::ostringstream os;
	os << "Teuchos::MpiComm::waitAll(): " << numNonmatching << " message "
	  "request" << (numNonmatching != 1 ? "s" : "") << " have a number of "
	  "bytes which does not match the number of bytes in " 
	   << (numNonmatching != 1 ? "their" : "its") << " corresponding status"
	   << (numNonmatching != 1 ? "es" : "") << "." << std::endl;
        os << "Here are the lengths that don't match (from MPI_Request, MPI_Status resp.): " << std::endl;
        for (Array<std::pair<size_type, size_type> >::const_iterator it = nonmatchingLengthPairs.begin(); it != nonmatchingLengthPairs.end(); ++it) {
          os << "(" << it->first << "," << it->second << ") ";
        }
        if (err == MPI_ERR_IN_STATUS) { 
          os << std::endl << "This is that weird case where MPI_Waitall returned MPI_ERR_IN_STATUS, but all of the MPI_Statuses' error codes were MPI_SUCCESS.";
        }
	// This is a bug, so we throw std::logic_error.
	TEUCHOS_TEST_FOR_EXCEPTION(numNonmatching > 0, std::logic_error, os.str());
      }
    }
#endif // HAVE_TEUCHOS_DEBUG

    // Invalidate the input array of requests by setting all entries to
    // null.
    std::fill (requests.begin(), requests.end(), null);
  }
} // namespace (anonymous)



template<typename Ordinal>
void 
MpiComm<Ordinal>::
waitAll (const ArrayView<RCP<CommRequest> >& requests) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::waitAll(requests)" );

  Array<MPI_Status> rawMpiStatuses (requests.size());
  waitAllImpl<Ordinal> (requests, rawMpiStatuses());
}


template<typename Ordinal>
void 
MpiComm<Ordinal>::
waitAll (const ArrayView<RCP<CommRequest> >& requests,
	 const ArrayView<RCP<CommStatus<Ordinal> > >& statuses) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::waitAll(requests, statuses)" );

  typedef ArrayView<RCP<CommRequest> >::size_type size_type;
  const size_type count = requests.size();

  TEUCHOS_TEST_FOR_EXCEPTION(count != statuses.size(), 
    std::invalid_argument, "Teuchos::MpiComm::waitAll: requests.size() = " 
    << count << " != statuses.size() = " << statuses.size() << ".");

  Array<MPI_Status> rawMpiStatuses (count);
  waitAllImpl<Ordinal> (requests, rawMpiStatuses());

  // Repackage the raw MPI_Status structs into the wrappers.
  for (size_type i = 0; i < count; ++i) {
    statuses[i] = mpiCommStatus<Ordinal> (rawMpiStatuses[i]);
  }
}


template<typename Ordinal>
RCP<CommStatus<Ordinal> > 
MpiComm<Ordinal>::wait (const Ptr<RCP<CommRequest> >& request) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::wait(...)" );

  if (is_null(*request)) {
    return null; // Nothing to wait on ...
  }
  const RCP<MpiCommRequest> mpiCommRequest =
    rcp_dynamic_cast<MpiCommRequest>(*request);
  // This function doesn't satisfy the strong exception guarantee,
  // because releaseRawMpiRequest() modifies the MpiCommRequest.
  MPI_Request rawMpiRequest = mpiCommRequest->releaseRawMpiRequest();
  MPI_Status status;
  const int err = MPI_Wait (&rawMpiRequest, &status);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error, 
    "Teuchos::MpiComm::wait: MPI_Wait() failed with error \""
    << mpiErrorCodeToString (err) << "\".");

  *request = null;
  return rcp (new MpiCommStatus<Ordinal> (status));
}


template<typename Ordinal>
RCP< Comm<Ordinal> >
MpiComm<Ordinal>::duplicate() const
{
  return rcp (new MpiComm<Ordinal> (*this));
}


template<typename Ordinal>
RCP< Comm<Ordinal> >
MpiComm<Ordinal>::split(const int color, const int key) const
{
  MPI_Comm newComm;
  const int splitReturn = 
    MPI_Comm_split (*rawMpiComm_, 
		    color < 0 ? MPI_UNDEFINED : color,
		    key,
		    &newComm);
  TEUCHOS_TEST_FOR_EXCEPTION(
    splitReturn != MPI_SUCCESS,
    std::logic_error,
    "Teuchos::MpiComm::split: Failed to create communicator with color " 
    << color << "and key " << key << ".  MPI_Comm_split failed with error \""
    << mpiErrorCodeToString (splitReturn) << "\".");
  if (newComm == MPI_COMM_NULL) {
    return RCP< Comm<Ordinal> >();
  } else {
    return rcp(new MpiComm<Ordinal>(
      rcp_implicit_cast<const OpaqueWrapper<MPI_Comm> >( opaqueWrapper(newComm) )
      ) );
  }
}


template<typename Ordinal>
RCP< Comm<Ordinal> >
MpiComm<Ordinal>::createSubcommunicator(const ArrayView<const int> &ranks) const
{
  int mpiReturn;

  // Get the group that this communicator is in.
  MPI_Group thisGroup;
  mpiReturn = MPI_Comm_group(*rawMpiComm_, &thisGroup);
  TEUCHOS_TEST_FOR_EXCEPTION(mpiReturn != MPI_SUCCESS, std::logic_error,
    "Failed to obtain group.  MPI_Comm_group failed with error \"" 
    << mpiErrorCodeToString (mpiReturn) << "\".");

  // Create a new group with the specified members.
  MPI_Group newGroup;
  mpiReturn = MPI_Group_incl(
    thisGroup, ranks.size(), const_cast<int *>(&ranks[0]), &newGroup);
  TEUCHOS_TEST_FOR_EXCEPTION(mpiReturn != MPI_SUCCESS, std::logic_error,
    "Failed to create subgroup.  MPI_Group_incl failed with error \""
    << mpiErrorCodeToString (mpiReturn) << "\".");

  // Create a new communicator from the new group.
  MPI_Comm newComm;
  mpiReturn = MPI_Comm_create(*rawMpiComm_, newGroup, &newComm);
  TEUCHOS_TEST_FOR_EXCEPTION(mpiReturn != MPI_SUCCESS, std::logic_error,
    "Failed to create subcommunicator.  MPI_Comm_create failed with error \""
    << mpiErrorCodeToString (mpiReturn) << "\".");

  if (newComm == MPI_COMM_NULL) {
    return RCP< Comm<Ordinal> >();
  } else {
    return rcp(new MpiComm<Ordinal>(
      rcp_implicit_cast<const OpaqueWrapper<MPI_Comm> >( opaqueWrapper(newComm) )
      ));
  }
}


// Overridden from Describable


template<typename Ordinal>
std::string MpiComm<Ordinal>::description() const
{
  std::ostringstream oss;
  oss
    << typeName(*this)
    << "{"
    << "size="<<size_
    << ",rank="<<rank_
    << ",rawMpiComm="<<static_cast<MPI_Comm>(*rawMpiComm_)
    <<"}";
  return oss.str();
}


#ifdef TEUCHOS_MPI_COMM_DUMP
template<typename Ordinal>
bool MpiComm<Ordinal>::show_dump = false;
#endif


// private


template<typename Ordinal>
void MpiComm<Ordinal>::assertRank(const int rank, const std::string &rankName) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! ( 0 <= rank && rank < size_ ), std::logic_error
    ,"Error, "<<rankName<<" = " << rank << " is not < 0 or is not"
    " in the range [0,"<<size_-1<<"]!"
    );
}


} // namespace Teuchos


template<typename Ordinal>
Teuchos::RCP<Teuchos::MpiComm<Ordinal> >
Teuchos::createMpiComm(
  const RCP<const OpaqueWrapper<MPI_Comm> > &rawMpiComm
  )
{
  if( rawMpiComm.get()!=NULL && *rawMpiComm != MPI_COMM_NULL )
    return rcp(new MpiComm<Ordinal>(rawMpiComm));
  return Teuchos::null;
}


#endif // TEUCHOS_MPI_COMM_HPP
