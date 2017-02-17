/**
//@HEADER
// ************************************************************************
//
//                   Trios: Trilinos I/O Support
//                 Copyright 2011 Sandia Corporation
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
//Questions? Contact Ron A. Oldfield (raoldfi@sandia.gov)
//
// *************************************************************************
//@HEADER
 */
/*-------------------------------------------------------------------------*/
/**  @file nssi_types_xdr.x
 *
 *   @brief XDR definitions for types used in the NSSI RPC interface.
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov).
 *
 */

#ifdef RPC_HDR
%#include "Trios_xdr.h"
%#include "Trios_nnti_xdr.h"
#endif

#ifdef RPC_XDR
%#include "Trios_xdr.h"
%#include "Trios_nnti_xdr.h"
#endif

/** @addtogroup base_types
 *  @{
 */

/**
 *
 * @brief NSSI Return Codes.
 *
 *  Each NSSI function returns an integer value that corresponds
 *  to an error, a warning, or successfull completion. These values are analogous
 *  to the values from the UNIX include file ``<tt>error.h</tt>''.
 */
enum nssi_return_codes {

    /** @brief The function completed successfully. */
    NSSI_OK = 0,

    /** @brief Unspecified I/O error. */
    NSSI_EIO = 1001,

    /** @brief The size of the message is larger than the supported maximum size. */
    NSSI_EMSGSIZE,

    /** @brief The operation or process has been canceled. */
    NSSI_ECANCELED,

    /** @brief An operation timed out. */
    NSSI_ETIMEDOUT,

    /** @brief The value of the variable or parameter is invalid. */
    NSSI_EINVAL,

    /** @brief  No memory is available. */
    NSSI_ENOMEM,

    /** @brief No such entry. */
    NSSI_ENOENT,

    /** @brief Unsupported operation. */
    NSSI_ENOTSUP,

    /** @brief The item already exists. */
    NSSI_EEXIST,

    /** @brief Unsuccessful RPC operation. */
    NSSI_EBADRPC,

    /** @brief Not initialized. */
    NSSI_ENOTINIT,

    /** @brief Insufficient priveleges to perform operation. */
    NSSI_EPERM,

    /** @brief Error decoding an RPC request. */
    NSSI_EDECODE,

    /** @brief Error encoding an RPC request. */
    NSSI_EENCODE,

    /** @brief An operation would have blocked. */
    NSSI_EWOULDBLOCK,

    /** @brief Operation was interupted, but possibly recoverable. */
    NSSI_EAGAIN

};


/**
 * @brief The <tt>\ref nssi_size</tt> type is used to for ``size'' variables.
 */
typedef uint64_t nssi_size;
typedef uint64_t nssi_ssize;



/**
 * @brief Length of an nssi_name.
 */
const NSSI_HOSTNAME_LEN = NNTI_HOSTNAME_LEN;
const NSSI_URL_LEN      = NNTI_URL_LEN;




/**
 * @brief Enumerator for the different type of transport mechanisms.
 *
 * The <tt>\ref nssi_rpc_transport</tt> enumerator provides integer values
 * to represent the different types of supported transport mechanisms.
 * Initially, the only supported transport mechanism is Portals.
 */
enum nssi_rpc_transport {
    /** @brief No operations permitted. */
    NSSI_RPC_NULL,

    /** @brief Use Portals to transfer rpc requests. */
    NSSI_RPC_PTL,

    /** @brief Use Infiniband to transfer rpc requests. */
    NSSI_RPC_IB,

    /** @brief Use Cray Gemini to transfer rpc requests. */
    NSSI_RPC_GEMINI,

   /** @brief Use Blue Gene/P DCMF Lib to transfer rpc requests. */
    NSSI_RPC_BGPDCMF,

    /** @brief Use Blue Gene/P PAMI Lib to transfer rpc requests. */
    NSSI_RPC_BGQPAMI,

    /** @brief Use MPI to transfer rpc requests. */
    NSSI_RPC_MPI,

    /** @brief Use a local buffer (not a remote operation). */
    NSSI_RPC_LOCAL
};

/**
 * @brief The number of RPC mechanisms supported by NSSI.
 */
const NSSI_RPC_COUNT = 8;


/**
 * @brief Enumerator for the different type of encoding mechanisms.
 *
 * The <tt>\ref nssi_rpc_encode</tt> enumerator provides integer values
 * to represent the different types of supported mechanisms for encoding
 * control messages transferred to/from NSSI servers.
 * Initially, the only supported transport mechanism is XDR.
 */
enum nssi_rpc_encode {
    /** @brief Use XDR to encode/decode rpc requests and results. */
    NSSI_RPC_XDR
};


/**
 * @brief Default timeout for rpc calls (in ms).
 */
const DEFAULT_RPC_TIMEOUT = 10000000;


/**
 * @brief A descriptor for remote services.
 *
 * The <tt>\ref nssi_service</tt> data
 * structure contains information needed by the client
 * to send operation requests to a remote service, including the process ID
 * of the remote host, the encoding mechanism to use for control messages,
 * and the memory location reserved for incoming requests.
 *
 * A client obtains the <tt>\ref nssi_service</tt> structure by calling the
 * <tt>\ref nssi_get_service "nssi_get_service()"</tt> function.
 */
struct nssi_service {

    /** @brief Identifies the RPC mechanism to use for transfering messages. */
    NNTI_transport_id_t transport_id;

    /** @brief Identifies the mechanism to use for encoding messages. */
    nssi_rpc_encode rpc_encode;

    /** @brief The address of the service. */
    NNTI_peer_t svc_host;

    /** @brief The remote memory address reserved for incoming requests. */
    NNTI_buffer_t req_addr;

    /** @brief This service accepts requests up to 'req_size' is length. */
    int32_t req_size;

    /** @brief This service sends results up to 'res_size' is length. */
    int32_t res_size;

    /** @brief The maximum number of requests to process. */
    int32_t max_reqs;

    /** @brief A callback to invoke at intervals defined by progess_callback_timeout. */
    uint64_t progress_callback;

    /** @brief The interval at which progess_callback should be invoked.  Note: This is not a real-time environment, so there is no guarentee of the interval. */
    uint64_t progress_callback_timeout;
};

/**
 * @brief The size of an encoded nssi_request_header buffer.
 */
const NSSI_SHORT_REQUEST_SIZE = NNTI_REQUEST_BUFFER_SIZE;
/**
 * @brief The size of an encoded nssi_result_header buffer.
 */
const NSSI_SHORT_RESULT_SIZE = NNTI_RESULT_BUFFER_SIZE;


/**
 * @brief The request header.
 *
 * The nssi_request_header structure contains details needed by an
 * NSSI server to perform a remote operation.
 *
 * The operation arguments may or may not be sent to the server
 * in the original request. If the field \em fetch_args is
 * \em true, the server will fetch arguments from \em args_addr
 * in a separate operation. Otherwise, the client sends the
 * arguments to the server with the request.
 */
struct nssi_request_header {
    NNTI_transport_header_t transport_header;

    /** @brief ID of the request */
    uint32_t id;

    /** @brief ID of the operation to perform. */
    uint32_t opcode;

    /** @brief A flag that tells the server to fetch args from
     *        <em>\ref args_addr</em>. */
    bool fetch_args;

    /** @brief A flag that tells the server to fetch bulk data from
     *        <em>\ref data_addr</em>. */
    bool fetch_data;

    /** @brief A flag that tells the server if the client expects 
     * a response. */
    bool is_responseless;

    /** @brief The remote memory address reserved for
     *        long arguments. */
    NNTI_buffer_t args_addr;

    /** @brief The remote memory address reserved for bulk
     *        data transfers. */
    NNTI_buffer_t data_addr;

    /** @brief The remote memory address reserved for the short result. */
    NNTI_buffer_t res_addr;
};


/**
 * @brief The result header.
 *
 * The mds_request_header structure contains the details needed by the
 * MDS server to perform a remote operation.
 *
 * The arguments may or may not be sent to the server in the original request.
 * If the field \em fetchargs is \em true, the arguments will be fetched
 * from \em args_portal in a separate operation by the server. Otherwise,
 * the arguments are sent directly to the server.
 *
 * NO POINTERS ALLOWED IN THIS STRUCTURE!
 *
 */
struct nssi_result_header {
    NNTI_transport_header_t transport_header;

    /** @brief ID of the result (should be same as request) */
    uint32_t id;

    /** @brief ID of the operation to perform. */
    uint32_t opcode;

    /** @brief Size of the result. */
    uint32_t result_size;

    /** @brief A flag that tells the client to "get" result from the client. */
    bool fetch_result;

    /** @brief The remote memory address reserved for long results. */
    NNTI_buffer_t result_addr;

    /** @brief The remote memory address reserved for long results ACK. */
    NNTI_buffer_t result_ack_addr;

    /** @brief The return code of the function. */
    uint32_t rc;
};

/** @} */
