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
/**
 *   @file nssi_request.h
 *
 *   @brief Definition of the nssi_request data structure.
 *
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov)
 */

#ifndef _NSSI_REQUEST_H
#define _NSSI_REQUEST_H


struct nssi_request;

#include "Trios_nssi_rpc.h"

#ifdef __cplusplus
extern "C" {
#endif

       /**
         * @ingroup base_types
         *
         * @brief States for a pending RPC request.
         *
         * The <tt>\ref nssi_request_status </tt> enumerator provides
         * integer values that identify the state of a pending NSSI request.
         */
        enum nssi_request_status {
                /** @brief The request is complete with an error. */
                NSSI_REQUEST_ERROR = -1,

                /** @brief The request has no/null status. */
                NSSI_REQUEST_NULL = 0,

                /** @brief The client is sending the request to server. */
                NSSI_SENDING_REQUEST = 1,

                /** @brief The remote server is processing the request. */
                NSSI_PROCESSING_REQUEST,

                /** @brief The client is processing the result. */
                NSSI_PROCESSING_RESULT,

                /** @brief The last wait operation timed out. */
                NSSI_REQUEST_TIMEDOUT,

                /** @brief The request is complete. */
                NSSI_REQUEST_COMPLETE
        };

        typedef enum nssi_request_status nssi_request_status;

        /**
         * @ingroup base_types
         *
         * @brief The request structure.
         *
         * The <tt>\ref nssi_request</tt> structure represents a pending NSSI
         * request.  It contains a unique identifier and pointers to all
         * data structures and buffers needed to send the request to the
         * remote server and to process the result when the server completes.
         *
         * \todo We need to abstract out the implementation-specific portions
         * of this data structure.  These include the fields to encode/decode
         * arguments and the all of the Portals data structures.
         */
        struct nssi_request {
                /** @brief The remote service that is handling this request. */
                const nssi_service *svc;

                int job_type;

                /** @brief An ID for this request. */
                unsigned long id;

                /** @brief The opcode of the remote function. */
                int opcode;

                /** @brief If TRUE, nssi_send_request() will mark this request complete - long args, bulk data and short result are not allowed. */
                int8_t is_responseless;

                /** @brief If TRUE, nssi_send_request() will wait for this request to complete. */
                int8_t is_sync;

                /** @brief If is_sync==TRUE, the timeout for nssi_timedwait(). */
                int8_t sync_timeout;

                /** @brief The error code of request. This value will be \ref NSSI_OK unless the
                 *          request status=\ref NSSI_REQUEST_ERROR . */
                int error_code;

                /** @brief Status of the pending request */
                nssi_request_status status;


                /***********************
                 ** Short Request
                 ***********************/
                /** @brief If TRUE, the app pinned the short request buffer, so do not unpin/free at cleanup. */
                int8_t               app_pinned_short_request;
                /** @brief Handle for the buffer where the short request will be put. */
                NNTI_buffer_t       *short_request_hdl;
                NNTI_buffer_t        short_request;
                /** @brief Work request to track the short request transfer. */
                NNTI_work_request_t  short_request_wr;

                /***********************
                 ** Long Args
                 ***********************/
                /** @brief If TRUE, the app pinned the long args buffer, so do not unpin/free at cleanup. */
                int8_t               app_pinned_long_args;
                int                  use_long_args;
                /** @brief Points to the memory reserved for the args. */
                void                *args;
                /** @brief Points to the XDR function used to encode arguments. */
                xdrproc_t            xdr_encode_args;
                /** @brief Handle for the buffer where the long arguments reside. */
                NNTI_buffer_t       *long_args_hdl;
                NNTI_buffer_t        long_args;
                /** @brief Work request to track the long arguments transfer. */
                NNTI_work_request_t  long_args_wr;

                /***********************
                 ** Bulk Data
                 ***********************/
                /** @brief If TRUE, the app pinned the bulk data buffer, so do not unpin/free at cleanup. */
                int8_t               app_pinned_bulk_data;
                /** @brief Points to the memory reserved for the bulk data transfers (NULL if not used). */
                void                *data;
                /** @brief The max size of the memory used for bulk transfers. */
                uint32_t             data_size;
                /** @brief Points to the XDR function used to encode data. */
                xdrproc_t            xdr_encode_data;
                /** @brief Handle for the buffer where the data reside. */
                NNTI_buffer_t       *bulk_data_hdl;
                NNTI_buffer_t        bulk_data;
                /** @brief Work request to track the bulk data transfer. */
                NNTI_work_request_t  bulk_data_wr;

                /***********************
                 ** Short Result
                 ***********************/
                /** @brief If TRUE, the app pinned the short result buffer, so do not unpin/free at cleanup. */
                int8_t               app_pinned_short_result;
                /** @brief Points to the memory reserved for the result. */
                void               *result;
                /** @brief Points to the XDR function used to decode the result. */
                xdrproc_t           xdr_decode_result;
                /** @brief Handle for the buffer where the short result will be put. */
                NNTI_buffer_t      *short_result_hdl;
                NNTI_buffer_t       short_result;
                /** @brief Work request to track the short result transfer. */
                NNTI_work_request_t short_result_wr;

                /***********************
                 ** Callback
                 ***********************/
                /** @brief A callback function used by the wait() function when a request is complete. */
                int (*callback)(struct nssi_request *);
                /** @brief Optional arguments for the callback function. */
                void *callback_args;

        };

        typedef struct nssi_request nssi_request;

#ifdef __cplusplus
}
#endif

#endif

