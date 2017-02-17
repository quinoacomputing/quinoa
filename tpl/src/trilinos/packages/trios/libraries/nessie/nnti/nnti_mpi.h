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
/**
 * nnti_mpi.h
 *
 *  Created on: Apr 4, 2012
 *      Author: thkorde
 */

#ifndef NNTI_MPI_H_
#define NNTI_MPI_H_

#include "Trios_config.h"

#include "Trios_nnti.h"
#include "nnti_internal.h"


#if defined(HAVE_TRIOS_MPI)
#include <mpi.h>
#endif


#define MPI_SET_MATCH_ANY(p) \
(p)->peer.transport_id = NNTI_TRANSPORT_MPI; \
(p)->peer.NNTI_remote_process_t_u.mpi.rank  = MPI_ANY_SOURCE;


#ifdef __cplusplus
extern "C" {
#endif

NNTI_result_t NNTI_mpi_init (
        const NNTI_transport_id_t  trans_id,
        const char                *my_url,
        NNTI_transport_t          *trans_hdl);

NNTI_result_t NNTI_mpi_get_url (
        const NNTI_transport_t *trans_hdl,
        char                   *url,
        const uint64_t          maxlen);

NNTI_result_t NNTI_mpi_connect (
        const NNTI_transport_t *trans_hdl,
        const char             *url,
        const int               timeout,
        NNTI_peer_t            *peer_hdl);

NNTI_result_t NNTI_mpi_disconnect (
        const NNTI_transport_t *trans_hdl,
        NNTI_peer_t            *peer_hdl);

NNTI_result_t NNTI_mpi_alloc (
        const NNTI_transport_t *trans_hdl,
        const uint64_t          element_size,
        const uint64_t          num_elements,
        const NNTI_buf_ops_t    ops,
        NNTI_buffer_t          *reg_buf);

NNTI_result_t NNTI_mpi_free (
        NNTI_buffer_t    *reg_buf);

NNTI_result_t NNTI_mpi_register_memory (
        const NNTI_transport_t *trans_hdl,
        char                   *buffer,
        const uint64_t          element_size,
        const uint64_t          num_elements,
        const NNTI_buf_ops_t    ops,
        NNTI_buffer_t          *reg_buf);

NNTI_result_t NNTI_mpi_register_segments (
        const NNTI_transport_t *trans_hdl,
        char                  **segments,
        const uint64_t         *segment_lengths,
        const uint64_t          num_segments,
        const NNTI_buf_ops_t    ops,
        NNTI_buffer_t          *reg_buf);

NNTI_result_t NNTI_mpi_unregister_memory (
        NNTI_buffer_t    *reg_buf);

NNTI_result_t NNTI_mpi_send (
        const NNTI_peer_t   *peer_hdl,
        const NNTI_buffer_t *msg_hdl,
        const NNTI_buffer_t *dest_hdl,
        NNTI_work_request_t  *wr);

NNTI_result_t NNTI_mpi_put (
        const NNTI_buffer_t *src_buffer_hdl,
        const uint64_t       src_offset,
        const uint64_t       src_length,
        const NNTI_buffer_t *dest_buffer_hdl,
        const uint64_t       dest_offset,
        NNTI_work_request_t  *wr);

NNTI_result_t NNTI_mpi_get (
        const NNTI_buffer_t *src_buffer_hdl,
        const uint64_t       src_offset,
        const uint64_t       src_length,
        const NNTI_buffer_t *dest_buffer_hdl,
        const uint64_t       dest_offset,
        NNTI_work_request_t  *wr);

NNTI_result_t NNTI_mpi_scatter (
        const NNTI_buffer_t  *src_buffer_hdl,
        const uint64_t        src_length,
        const NNTI_buffer_t **dest_buffer_list,
        const uint64_t        dest_count,
        NNTI_work_request_t  *wr);

NNTI_result_t NNTI_mpi_gather (
        const NNTI_buffer_t **src_buffer_list,
        const uint64_t        src_length,
        const uint64_t        src_count,
        const NNTI_buffer_t  *dest_buffer_hdl,
        NNTI_work_request_t  *wr);

NNTI_result_t NNTI_mpi_atomic_set_callback (
		const NNTI_transport_t *trans_hdl,
		const uint64_t          local_atomic,
		NNTI_callback_fn_t      cbfunc,
		void                   *context);

NNTI_result_t NNTI_mpi_atomic_read (
		const NNTI_transport_t *trans_hdl,
		const uint64_t          local_atomic,
		int64_t                *value);

NNTI_result_t NNTI_mpi_atomic_fop (
		const NNTI_transport_t *trans_hdl,
		const NNTI_peer_t      *peer_hdl,
		const uint64_t          target_atomic,
		const uint64_t          result_atomic,
		const int64_t           operand,
		const NNTI_atomic_op_t  op,
		NNTI_work_request_t    *wr);

NNTI_result_t NNTI_mpi_atomic_cswap (
		const NNTI_transport_t *trans_hdl,
		const NNTI_peer_t      *peer_hdl,
		const uint64_t          target_atomic,
		const uint64_t          result_atomic,
		const int64_t           compare_operand,
		const int64_t           swap_operand,
		NNTI_work_request_t    *wr);

NNTI_result_t NNTI_mpi_create_work_request (
        NNTI_buffer_t        *reg_buf,
        NNTI_work_request_t  *wr);

NNTI_result_t NNTI_mpi_clear_work_request (
        NNTI_work_request_t  *wr);

NNTI_result_t NNTI_mpi_destroy_work_request (
        NNTI_work_request_t  *wr);

NNTI_result_t NNTI_mpi_cancel (
        NNTI_work_request_t *wr);

NNTI_result_t NNTI_mpi_cancelall (
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count);

NNTI_result_t NNTI_mpi_interrupt (
        const NNTI_transport_t *trans_hdl);

NNTI_result_t NNTI_mpi_wait (
        NNTI_work_request_t *wr,
        const int            timeout,
        NNTI_status_t       *status);

NNTI_result_t NNTI_mpi_waitany (
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count,
        const int             timeout,
        uint32_t             *which,
        NNTI_status_t        *status);

NNTI_result_t NNTI_mpi_waitall (
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count,
        const int             timeout,
        NNTI_status_t       **status);

NNTI_result_t NNTI_mpi_fini (
        const NNTI_transport_t *trans_hdl);

#ifdef __cplusplus
}
#endif

#endif /* NNTI_MPI_H_*/
