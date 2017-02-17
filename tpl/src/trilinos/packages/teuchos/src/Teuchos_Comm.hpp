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

#ifndef TEUCHOS_COMM_HPP
#define TEUCHOS_COMM_HPP

#include "Teuchos_ReductionOp.hpp"
#include "Teuchos_ArrayRCP.hpp"


namespace Teuchos {

/// \class CommRequest
/// \brief Encapsulation of MPI_Request.
///
/// An MPI_Request encapsulates the result of a nonblocking MPI send
/// or receive (which we in turn encapsulate by \c isend() resp. \c
/// ireceive()).  This class in turn wraps MPI_Request.
///
/// This is an opaque object which is meant to be given to \c wait()
/// or \c waitall().
class CommRequest : public Teuchos::Describable {};

/// \class CommStatus
/// \brief Encapsulation of MPI_Status.
///
/// This interface encapsulates the result of a receive (the \c
/// MPI_Status struct).  Its main use is to figure out which process
/// sent you a message, if you received it using \c MPI_ANY_SOURCE.
///
/// \tparam OrdinalType The same template parameter as \c Comm.  Only
///   use \c int here.  We only make this a template class for
///   compatibility with \c Comm.
///
/// \note For now, this class only exposes the rank of the process
///   that sent the message (the "source rank").  Later, we might
///   expose other fields of \c MPI_Status in this interface.  For
///   now, you can attempt a dynamic cast to \c MpiCommStatus to
///   access all three fields (MPI_SOURCE, MPI_TAG, and MPI_ERROR).
template<class OrdinalType>
class CommStatus {
public:
  //! Destructor (declared virtual for memory safety)
  virtual ~CommStatus() {}

  //! The source rank that sent the message.
  virtual OrdinalType getSourceRank () = 0;
};

/** \brief Abstract interface class for a basic communication channel between
 * one or more processes.
 *
 * This interface is templated on the ordinal type but only deals with buffers
 * of untyped data represented as arrays <tt>char</tt> type. All reduction
 * operations that are initiated by the concreate communicator object are
 * performed by user-defined <tt>ReductOpBase</tt> objects.  It is the
 * responsibility of the <tt>ReductOpBase</tt> object to know what the currect
 * data type is, to perform casts or serializations/unserializations to and
 * from <tt>char[]</tt> buffers, and to know how to reduce the objects
 * correctly.  It is strictly up to the client to correctly convert data types
 * to <tt>char[]</tt> arrays but there is a great deal of helper code to make
 * this easy and safe.
 *
 * ToDo: Finish documentation!
 */
template<typename Ordinal>
class Comm : virtual public Describable {
public:
  
  //! @name Query functions 
  //@{

  /** \brief Returns the rank of this process.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>0 <= return && return < this->getSize()</tt>
   * </ul>
   */
  virtual int getRank() const = 0;

  /** \brief Returns the number of processes that make up this communicator.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>return > 0</tt>
   * </ul>
   */
  virtual int getSize() const = 0;
  
  //@}

  //! @name Collective Operations 
  //@{

  /** \brief Pause every process in <tt>*this</tt> communicator until all the
   * processes reach this point.
   */
  virtual void barrier() const = 0;
  
  /** \brief Broadcast values from the root process to the slave processes.
   *
   * \param rootRank [in] The rank of the root process.
   *
   * \param count [in] The number of bytes in <tt>buffer[]</tt>.
   *
   * \param buffer [in/out] Array (length <tt>bytes</tt>) of packed data.
   * Must be set on input on the root processes with rank <tt>root</tt>.  On
   * output, each processs, including the root process contains the data.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>0 <= rootRank && rootRank < this->getSize()</tt>
   * </ul>
   */
  virtual void broadcast(
    const int rootRank, const Ordinal bytes, char buffer[]
    ) const = 0;

  /** \brief Gather values from each process to collect on all processes.
   *
   * \param sendBytes [in] Number of entires in <tt>sendBuffer[]</tt> on
   * input.
   *
   * \param sendBuffer [in] Array (length <tt>sendBytes</tt>) of data being
   * sent from each process.
   *
   * \param recvBytes [in] Number of entires in <tt>recvBuffer[]</tt> which
   * must be equal to <tt>sendBytes*this->getSize()</tt>.  This field is just
   * here for debug checking.
   *
   * \param recvBuffer [out] Array (length <tt>recvBytes</tt>) of all of the
   * entires sent from each processes.  Specifically,
   * <tt>recvBuffer[sendBytes*j+i]</tt>, for <tt>j=0...this->getSize()-1</tt>
   * and <tt>i=0...sendBytes-1</tt>, is the entry <tt>sendBuffer[i]</tt> from
   * process with rank <tt>j</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>recvBytes==sendBytes*this->getSize()</tt>
   * </ul>
   */
  virtual void gatherAll(
    const Ordinal sendBytes, const char sendBuffer[]
    ,const Ordinal recvBytes, char recvBuffer[]
    ) const = 0;

  /** \brief Global reduction.
   *
   * \param reductOp [in] The user-defined reduction operation
   *
   * \param bytes [in] The length of the buffers <tt>sendBuffer[]</tt> and
   * <tt>globalReducts[]</tt>.
   *
   * \param sendBuffer [in] Array (length <tt>bytes</tt>) of the data
   * contributed from each process.
   *
   * \param globalReducts [out] Array (length <tt>bytes</tt>) of the global
   * reduction from each process.
   */
  virtual void reduceAll(
    const ValueTypeReductionOp<Ordinal,char> &reductOp
    ,const Ordinal bytes, const char sendBuffer[], char globalReducts[]
    ) const = 0;

  /** \brief Global reduction combined with a scatter.
   *
   * \param reductOp [in] The user-defined reduction operation that accepts
   * char arrays.
   *
   * \param sendBytes [in] The number of entires in <tt>sendBuffer[]</tt>.
   * This must be the same in each process.
   *
   * \param sendBuffer [in] Array (length <tt>sendBytes</tt>) of the data
   * contributed from each process.
   *
   * \param recvCounts [in] Array (length <tt>this->getSize()</tt>) which
   * gives the number of chars from the global reduction that will be received
   * in each process.
   *
   * \param myGlobalReducts [out] Array (length
   * <tt>blockSize*recvCounts[rank]</tt>) of the global reductions gathered in
   * this process.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>sendBytes == sum(recvCounts[i],i=0...this->getSize()-1)</tt>
   * </ul>
   */
  virtual void reduceAllAndScatter(
    const ValueTypeReductionOp<Ordinal,char> &reductOp
    ,const Ordinal sendBytes, const char sendBuffer[]
    ,const Ordinal recvCounts[], char myGlobalReducts[]
    ) const = 0;

  /** \brief Scan reduction.
   *
   * \param reductOp [in] The user-defined reduction operation
   *
   * \param bytes [in] The length of the buffers <tt>sendBuffer[]</tt> and
   * <tt>scanReducts[]</tt>.
   *
   * \param sendBuffer [in] Array (length <tt>bytes</tt>) of the data
   * contributed from each process.
   *
   * \param scanReducts [out] Array (length <tt>bytes</tt>) of the reduction
   * up to and including this process.
   */
	virtual void scan(
    const ValueTypeReductionOp<Ordinal,char> &reductOp
    ,const Ordinal bytes, const char sendBuffer[], char scanReducts[]
    ) const = 0;

  //! @name Blocking Point-to-Point Operations 
  //@{

  /** \brief Possibly blocking send of data from this process to another process.
   *
   * This routine does not return until you can reuse the send buffer.
   * Whether this routine blocks depends on whether the MPI
   * implementation buffers.
   *
   * \param bytes [in] The number of bytes of data being passed between
   * processes.
   *
   * \param sendBuffer [in] Array (length <tt>bytes</tt>) of data being sent
   * from this process.  This buffer can be immediately destroyed or reused as
   * soon as the function exits (that is why this function is "blocking").
   *
   * \param destRank [in] The rank of the process to receive the data.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>0 <= destRank && destRank < this->getSize()</tt>
   * <li><tt>destRank != this->getRank()</tt>
   * </ul>
   */
  virtual void send(
    const Ordinal bytes, const char sendBuffer[], const int destRank
    ) const = 0;

  /** \brief Always blocking send of data from this process to another process.
   *
   * This routine blocks until the matching receive posts.  After it
   * returns, you are allowed to reuse the send buffer.
   *
   * \param bytes [in] The number of bytes of data being passed between
   * processes.
   *
   * \param sendBuffer [in] Array (length <tt>bytes</tt>) of data being sent
   * from this process.  This buffer can be immediately destroyed or reused as
   * soon as the function exits (that is why this function is "blocking").
   *
   * \param destRank [in] The rank of the process to receive the data.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>0 <= destRank && destRank < this->getSize()</tt>
   * <li><tt>destRank != this->getRank()</tt>
   * </ul>
   */
  virtual void ssend(
    const Ordinal bytes, const char sendBuffer[], const int destRank
    ) const = 0;


  /** \brief Blocking receive of data from this process to another process.
   *
   * \param sourceRank [in] The rank of the process to receive the data from.
   * If <tt>sourceRank < 0</tt> then data will be received from any process.
   *
   * \param bytes [in] The number of bytes of data being passed between
   * processes.
   *
   * \param recvBuffer [out] Array (length <tt>bytes</tt>) of data being
   * received from this process.  This buffer can be immediately used to
   * access the data as soon as the function exits (that is why this function
   * is "blocking").
   *
   * <b>Preconditions:</b><ul>
   * <li>[<tt>sourceRank >= 0] <tt>sourceRank < this->getSize()</tt>
   * <li><tt>sourceRank != this->getRank()</tt>
   * </ul>
   *
   * \return Returns the senders rank.
   */
  virtual int receive(
    const int sourceRank, const Ordinal bytes, char recvBuffer[]
    ) const = 0;


  /** \brief Ready send of data from this process to another process.
   *
   * \param sendBuffer [in] The data to be sent.
   *
   * \param destRank [in] The rank of the process to receive the data.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>0 <= destRank && destRank < this->getSize()</tt>
   * <li><tt>destRank != this->getRank()</tt>
   * </ul>
   */
  virtual void readySend(
    const ArrayView<const char> &sendBuffer,
    const int destRank
    ) const = 0;
  

  //@}


  //! @name Non-blocking Point-to-Point Operations 
  //@{


  /** \brief Non-blocking send. 
   *
   * \param sendBuffer [in] The data buffer to be sent.
   *
   * \param destRank [in] The rank of the process to receive the data.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>0 <= destRank && destRank < this->getSize()</tt>
   * <li><tt>destRank != this->getRank()</tt>
   * </ul>
   */
  virtual RCP<CommRequest> isend(
    const ArrayView<const char> &sendBuffer,
    const int destRank
    ) const = 0;


  /** \brief Non-blocking receive. 
   *
   * \param recvBuffer [out] The location for storing the received data.
   *
   * \param sourceRank [in] The rank of the process to receive the data from.
   * If <tt>sourceRank < 0</tt> then data will be received from any process.
   *
   * <b>Preconditions:</b><ul>
   * <li>[<tt>sourceRank >= 0] <tt>sourceRank < this->getSize()</tt>
   * <li><tt>sourceRank != this->getRank()</tt>
   * </ul>
   *
   * \return Returns the senders rank.
   */
  virtual RCP<CommRequest> ireceive(
    const ArrayView<char> &recvBuffer,
    const int sourceRank
    ) const = 0;


  /** \brief Wait on a set of communication requests.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>requests.size() > 0</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li> <tt>is_null(request[i]))</tt> for <tt>i=0...requests.size()-1</tt>
   * </ul>
   */
  virtual void waitAll(
    const ArrayView<RCP<CommRequest> > &requests
    ) const = 0;

  /// \brief Wait on communication requests, and return their statuses.
  /// 
  /// \pre requests.size() == statuses.size()
  ///
  /// \pre For i in 0, 1, ..., requests.size()-1, requests[i] is
  ///   either null or requests[i] was returned by an ireceive() or
  ///   isend().
  ///
  /// \post For i in 0, 1, ..., requests.size()-1,
  ///   requests[i].is_null() is true.
  ///
  /// \param requests [in/out] On input: the requests on which to
  ///   wait.  On output: all set to null.
  ///
  /// \param statuses [out] The status results of waiting on the
  ///   requests.
  virtual void 
  waitAll (const ArrayView<RCP<CommRequest> >& requests,
	   const ArrayView<RCP<CommStatus<Ordinal> > >& statuses) const = 0;

  /// \brief Wait on a single communication request, and return its status.
  ///
  /// \param request [in/out] On input: request is not null, and
  /// *request is either null (in which case this function does
  /// nothing and returns null) or an RCP of a valid CommRequest
  /// instance representing an outstanding communication request.  On
  /// output: If the communication request completed successfully, we
  /// set *request to null, indicating that the request has completed.
  /// (This helps prevent common bugs like trying to complete the same
  /// request twice.)
  ///
  /// \return If *request is null, this method returns null.
  /// Otherwise this method returns a \c CommStatus instance
  /// representing the result of completing the request.  In the case
  /// of a nonblocking receive request, you can query the \c
  /// CommStatus instance for the process ID of the sending process.
  /// (This is useful for receiving from any process via \c
  /// MPI_ANY_SOURCE.)
  /// 
  /// \pre !is_null(request) (that is, the Ptr is not null).
  /// \post is_null(*request) (that is, the RCP is null).
  ///
  /// This function blocks until the communication operation
  /// associated with the CommRequest object has completed.
  virtual RCP<CommStatus<Ordinal> > 
  wait (const Ptr<RCP<CommRequest> >& request) const = 0;

  //@}

  //! @name Subcommunicator Operations
  //@{

  /**
   * \brief Duplicate this communicator.
   *
   * Make a copy of this communicator with a duplicate communication space.
   * Note that the returned communicator has the same properties as this
   * communicator, but is distinct from the original.
   *
   * \return A new communicator.
   */
  virtual RCP< Comm > duplicate() const = 0;

  /**
   * \brief Split a communicator into subcommunicators based on color
   * and key.
   *
   * Generates a partition of this communicator into multiple disjoint
   * groups.  There will be as many groups as there are globally many
   * distinct values for the color parameter. Within each subset of
   * the partition, the ranks will be ordered according to the key
   * value each rank passed for the key parameter. If multiple ranks
   * pass the same value for the key parameter, then they will be
   * ordered according to their rank value in the original
   * communicator. To return a valid communicator, this function
   * requires a nonnegative value for color. If a negative value
   * is passed for color, a null communicator will be
   * returned.
   *
   * All members of this communicator must call this member if it is
   * called at all.
   *
   * \param color [in] An integer representing the color for the local
   * rank.  
   *
   * \param key [in] A key value to order processes of the same color.
   *
   * \return A partitioned communicator.
   */
  virtual RCP< Comm > split(const int color, const int key) const = 0;

  /**
   * \brief Create a subcommunicator containing the specified ranks of
   * this one.
   *
   * Generates a subcommunicator of this communicator that consists of
   * the specified ranks in the order in which they are listed in the
   * input vector. Ranks that are not specified in the input vector
   * will be given a null communicator.
   *
   * All members of this communicator must call this member if it is
   * called at all.
   *
   * \param ranks A vector containing the ranks to include in the new
   * communicator.
   *
   * \return A new communicator.
   */
  virtual RCP<Comm> createSubcommunicator(
    const ArrayView<const int>& ranks) const = 0;

  //@}
  
}; // class Comm
  
} // namespace Teuchos

#endif // TEUCHOS_COMM_HPP
