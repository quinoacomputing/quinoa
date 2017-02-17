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
 * nnti_ptls.c
 *
 *  Created on: Jan 13, 2011
 *      Author: thkorde
 */

#include "Trios_config.h"
#include "Trios_threads.h"
#include "Trios_timer.h"
#include "Trios_signal.h"
#include "Trios_nnti_fprint_types.h"

// MPI is only used to increment the PID
#ifdef HAVE_TRIOS_MPI
#include <mpi.h>
#endif

#include <assert.h>
#include <string.h>

#include <map>
#include <deque>

#include "nnti_ptls.h"
#include "nnti_utils.h"



/* if defined, the RDMA initiator will send an ACK message to the RDMA
 * target when the RDMA op is complete.  the target process must wait
 * on the target buffer in order to get the ACK.  this creates two-sided
 * semantics for RDMA ops.   in this mode, when the wait returns the
 * the RDMA op is complete and status indicates what data was addressed.
 */
#define USE_RDMA_TARGET_ACK
/* if undefined, the ACK message is NOT sent to the RDMA target when
 * the RDMA op is complete.  this creates one-sided semantics for RDMA
 * ops.  in this mode, the target has no idea when the RDMA op is
 * complete and what data was addressed.  NNTI_wait() returns NNTI_EINVAL
 * if passed a target buffer.
 */
#undef USE_RDMA_TARGET_ACK



#define PTL_OP_PUT_INITIATOR  1
#define PTL_OP_GET_INITIATOR  2
#define PTL_OP_PUT_TARGET     3
#define PTL_OP_GET_TARGET     4
#define PTL_OP_SEND           5
#define PTL_OP_NEW_REQUEST    6
#define PTL_OP_RECEIVE        8

/**
 * @brief Track the state of a PtlPut (i'm the initiator).
 */
typedef struct {
    int send_start;
    int send_end;
    int ack;
    int unlink;
} ptl_put_initiator_state_t;
/**
 * @brief Track the state of a PtlGet (i'm the initiator).
 */
typedef struct {
    int send_start;
    int send_end;
    int reply_start;
    int reply_end;
    int unlink;
} ptl_get_initiator_state_t;
/**
 * @brief Track the state of a PtlPut (i'm the target).
 */
typedef struct {
    int put_start;
    int put_end;
} ptl_put_target_state_t;
/**
 * @brief Track the state of a PtlGet (i'm the target).
 */
typedef struct {
    int get_start;
    int get_end;
} ptl_get_target_state_t;

typedef union {
    ptl_put_initiator_state_t put_initiator;
    ptl_get_initiator_state_t get_initiator;
    ptl_put_target_state_t    put_target;
    ptl_get_target_state_t    get_target;
} ptl_op_state_t;

typedef enum {
    REQUEST_BUFFER,
    RECEIVE_BUFFER,
    SEND_BUFFER,
    GET_SRC_BUFFER,
    GET_DST_BUFFER,
    PUT_SRC_BUFFER,
    PUT_DST_BUFFER,
    RDMA_TARGET_BUFFER,
    UNKNOWN_BUFFER
} ptl_buffer_type;

typedef struct portals_work_request {
    NNTI_buffer_t   *reg_buf;
    NNTI_peer_t      peer;
    uint64_t         src_offset;
    uint64_t         dst_offset;
    uint64_t         length;

    ptl_event_t      last_event;

    uint8_t          last_op;
    ptl_op_state_t   op_state;
    uint8_t          is_last_op_complete;

    ptl_handle_me_t  me_h;
    ptl_md_t         md;
    ptl_handle_md_t  md_h;
} portals_work_request;

typedef std::deque<NNTI_work_request_t *>           wr_queue_t;
typedef std::deque<NNTI_work_request_t *>::iterator wr_queue_iter_t;

typedef struct portals_memory_handle {
    ptl_buffer_type  type;

    ptl_handle_eq_t  eq_h;
    ptl_pt_index_t   buffer_id;
    ptl_process_id_t match_id;
    ptl_match_bits_t match_bits;
    ptl_match_bits_t ignore_bits;
    ptl_handle_me_t  me_h;
    ptl_md_t         md;
    ptl_handle_md_t  md_h;

    wr_queue_t wr_queue;
} portals_memory_handle;


#define NUM_REQ_QUEUES 2
typedef struct portals_request_queue_handle {
    NNTI_buffer_t *reg_buf;

    /* incoming queues */
    char *req_queue[NUM_REQ_QUEUES];

    /* each message is no larger than req_size */
    int req_size;

    /* each md can recv reqs_per_queue messages */
    int reqs_per_queue;

    /* keep track of the queue index and count */
    int indices[NUM_REQ_QUEUES];
    int queue_count[NUM_REQ_QUEUES];

    /* for each queue, we need these structs */
    ptl_handle_me_t  me_h[NUM_REQ_QUEUES];
    ptl_md_t         md[NUM_REQ_QUEUES];
    ptl_handle_md_t  md_h[NUM_REQ_QUEUES];
} portals_request_queue_handle;

typedef struct portals_transport_global {

    ptl_handle_ni_t  ni_h;

    ptl_process_id_t me;

    ptl_handle_eq_t  req_eq_h;
    ptl_handle_eq_t  data_eq_h;

    portals_request_queue_handle req_queue;

    bool init_called_mpi_init;

} portals_transport_global;



static nthread_lock_t nnti_ptl_lock;


static void set_req_pid(NNTI_pid *pid);
static const NNTI_work_request_t *decode_event_wr(
        const NNTI_work_request_t *wait_wr,
        const ptl_event_t         *event);
static int process_event(
        const NNTI_work_request_t *reg_buf,
        const ptl_event_t         *event);
//static NNTI_result_t post_recv_work_request(
//        NNTI_buffer_t *reg_buf);
//static NNTI_result_t repost_recv_work_request(
//        NNTI_work_request_t *wr);
static int8_t is_wr_complete(
        portals_work_request *wr);
static int8_t is_any_wr_complete(
        portals_work_request **wr_list,
        const uint32_t         wr_count,
        uint32_t              *which);
static int8_t is_all_wr_complete(
        portals_work_request **wr_list,
        const uint32_t         wr_count);
static int8_t is_wr_complete(
        NNTI_work_request_t *wr);
static int8_t is_any_wr_complete(
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count,
        uint32_t             *which);
static int8_t is_all_wr_complete(
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count);
//static portals_work_request *first_incomplete_wr(
//        portals_memory_handle *ptl_mem_hdl);
//static int8_t is_wr_queue_empty(
//        const NNTI_buffer_t *reg_buf);
//static int8_t is_buf_op_complete(
//        const NNTI_buffer_t *reg_buf);
//static int8_t is_any_buf_op_complete(
//        const NNTI_buffer_t **buf_list,
//        const uint32_t        buf_count,
//        uint32_t             *which);
//static int8_t is_all_buf_ops_complete(
//        const NNTI_buffer_t **buf_list,
//        const uint32_t        buf_count);

static void create_status(
        const NNTI_work_request_t *reg_buf,
        int                        nnti_rc,
        NNTI_status_t             *status);
static void create_peer(
        NNTI_peer_t *peer,
        ptl_nid_t nid,
        ptl_pid_t pid);


#define PTL_MEM_HDL(b) ((portals_memory_handle *)((b)->transport_private))
#define PTL_WORK_REQUEST(wr) ((portals_work_request *)((wr)->transport_private))

static bool ptl_initialized=false;


static portals_transport_global transport_global_data;
static const int MIN_TIMEOUT = 10;  /* in milliseconds */

/**
 * @brief Initialize NNTI to use a specific transport.
 *
 * Enable the use of a particular transport by this process.  <tt>my_url</tt>
 * allows the process to have some control (if possible) over the
 * URL assigned for the transport.  For example, a Portals URL to put
 * might be "ptl://-1,128".  This would tell Portals to use the default
 * network ID, but use PID=128.  If the transport
 * can be initialized without this info (eg. a Portals client), <tt>my_url</tt> can
 * be NULL or empty.
 */
NNTI_result_t NNTI_ptl_init (
        const NNTI_transport_id_t  trans_id,
        const char                *my_url,
        NNTI_transport_t          *trans_hdl)
{
    int rc=0;
    NNTI_result_t nnti_rc=NNTI_OK;

    int max_interfaces;
    ptl_ni_limits_t actual;


    char transport[NNTI_URL_LEN];
    char address[NNTI_URL_LEN];
//    char memdesc[NNTI_URL_LEN];
    char *sep;

//    NNTI_nid nid;
    NNTI_pid pid=-1;

    log_debug(nnti_debug_level, "enter");

    assert(trans_hdl);


    if (!ptl_initialized) {

        nthread_lock_init(&nnti_ptl_lock);

        if (my_url != NULL) {
            if ((nnti_rc=nnti_url_get_transport(my_url, transport, NNTI_URL_LEN)) != NNTI_OK) {
                return(nnti_rc);
            }
            if (0!=strcmp(transport, "ptl")) {
                return(NNTI_EINVAL);
            }

            if ((nnti_rc=nnti_url_get_address(my_url, address, NNTI_URL_LEN)) != NNTI_OK) {
                return(nnti_rc);
            }

            sep=strchr(address, ':');
//            nid=strtol(address, NULL, 0);
            pid=strtol(sep+1, NULL, 0);
        }
        if (pid == -1) {
            set_req_pid(&pid);
        }


        log_debug(nnti_debug_level, "initializing portals");

//        /* The UTCP NAL requires that the application defines where the Portals
//         * API and library should send any output. We'll send the output to
//         * the same file as logger.
//         */
//        utcp_api_out = logger_get_file();
//        utcp_lib_out = logger_get_file();
//
//
//        /* register trace groups (let someone else enable) */
//        trace_counter_gid = trace_get_gid(TRACE_RPC_COUNTER_GNAME);
//        trace_interval_gid = trace_get_gid(TRACE_RPC_INTERVAL_GNAME);


        memset(&transport_global_data, 0, sizeof(portals_transport_global));

        /* initialize the portals library */
        log_debug(nnti_debug_level, "initializing portals library");
        rc = PtlInit(&max_interfaces);
        if (rc) {
            log_fatal(nnti_debug_level,"PtlInit() failed, %s", ptl_err_str[rc]);
            abort();
        }

        /* initialize the portals interface */
        log_debug(nnti_debug_level, "initializing portals network interface - pid=%d", (int)pid);
        rc = PtlNIInit(PTL_IFACE_DEFAULT, pid, NULL, &actual, &transport_global_data.ni_h);
        if ((rc != PTL_OK) && (rc != PTL_IFACE_DUP)) {
            log_fatal(nnti_debug_level, "PtlNIInit() failed, %s", ptl_err_str[rc]);
            abort();
        }

        ptl_process_id_t ptl_id;
        rc = PtlGetId(transport_global_data.ni_h, &ptl_id);
        if (rc != PTL_OK) {
            log_error(nnti_debug_level,"failed %s", ptl_err_str[rc]);
        }

        transport_global_data.me.nid = ptl_id.nid;
        transport_global_data.me.pid = ptl_id.pid;

        /* create an event queue */
        /* TODO: should we share an event queue? */
        nthread_lock(&nnti_ptl_lock);
        rc = PtlEQAlloc(
                transport_global_data.ni_h,
                5000,
                PTL_EQ_HANDLER_NONE,
                &transport_global_data.data_eq_h);
        nthread_unlock(&nnti_ptl_lock);
        if (rc != NNTI_OK) {
            log_error(nnti_debug_level, "failed to allocate eventq");
            abort();
        }
        log_debug(nnti_debug_level, "allocated transport_global_data data eq=%d", transport_global_data.data_eq_h);

        if (logging_info(nnti_debug_level)) {
            fprintf(logger_get_file(), "Portals Initialized: nid=%llu, pid=%llu\n",
                    (unsigned long long)transport_global_data.me.nid,
                    (unsigned long long)transport_global_data.me.pid);
        }

        log_debug(nnti_debug_level, "sizeof(trans_hdl)=%d", trans_hdl);

        create_peer(&trans_hdl->me, transport_global_data.me.nid, transport_global_data.me.pid);

        ptl_initialized = true;
    }


    log_debug(nnti_debug_level, "exit");


    return(nnti_rc);
}


/**
 * @brief Return the URL field of this transport.
 *
 * Return the URL field of this transport.  After initialization, the transport will
 * have a specific location on the network where peers can contact it.  The
 * transport will convert this location to a string that other instances of the
 * transport will recognize.
 *
 * URL format: "transport://address/memory_descriptor"
 *    - transport - (required) identifies how the URL should parsed
 *    - address   - (required) uniquely identifies a location on the network
 *                - ex. "ptl://nid:pid/", "ib://ip_addr:port"
 *    - memory_descriptor - (optional) transport-specific representation of RMA params
 */
NNTI_result_t NNTI_ptl_get_url (
        const NNTI_transport_t *trans_hdl,
        char                   *url,
        const uint64_t          maxlen)
{
    NNTI_result_t nnti_rc=NNTI_OK;

    assert(trans_hdl);
    assert(url);
    assert(maxlen>0);

    strncpy(url, trans_hdl->me.url, maxlen);
    url[maxlen-1]='\0';

    return(nnti_rc);
}


/**
 * @brief Prepare for communication with the peer identified by <tt>url</tt>.
 *
 * Parse <tt>url</tt> in a transport specific way.  Perform any transport specific
 * actions necessary to begin communication with this peer.
 *
// * If the peer is found and responds
// * to a ping, a handle will be allocated and assigned to the pointer.  This
// * handle should be used to move data to/from the peer.
 *
 * Connectionless transport: parse and populate
 * Connected transport: parse, connection and populate
 *
 */
NNTI_result_t NNTI_ptl_connect (
        const NNTI_transport_t *trans_hdl,
        const char             *url,
        const int               timeout,
        NNTI_peer_t            *peer_hdl)
{
    NNTI_result_t nnti_rc=NNTI_OK;

    char transport[NNTI_URL_LEN];
    char address[NNTI_URL_LEN];
//    char memdesc[NNTI_URL_LEN];
    char *sep;

    NNTI_nid nid;
    NNTI_pid pid;

    log_debug(nnti_debug_level, "enter");

    assert(trans_hdl);
    assert(peer_hdl);

    if (url != NULL) {
        if ((nnti_rc=nnti_url_get_transport(url, transport, NNTI_URL_LEN)) != NNTI_OK) {
            return(nnti_rc);
        }
        if (0!=strcmp(transport, "ptl")) {
            /* the peer described by 'url' is not a Portals peer */
            return(NNTI_EINVAL);
        }

        if ((nnti_rc=nnti_url_get_address(url, address, NNTI_URL_LEN)) != NNTI_OK) {
            return(nnti_rc);
        }

        sep=strchr(address, ':');
        nid=strtol(address, NULL, 0);
        pid=strtol(sep+1, NULL, 0);
    } else {
        /*  */
        return(NNTI_EINVAL);
    }

    create_peer(
            peer_hdl,
            nid,
            pid);

    log_debug(nnti_debug_level, "exit");

    return(nnti_rc);
}


/**
 * @brief Terminate communication with this peer.
 *
 * Perform any transport specific actions necessary to end communication with
 * this peer.
 */
NNTI_result_t NNTI_ptl_disconnect (
        const NNTI_transport_t *trans_hdl,
        NNTI_peer_t            *peer_hdl)
{
    NNTI_result_t nnti_rc=NNTI_OK;

    assert(trans_hdl);
    assert(peer_hdl);

    return(nnti_rc);
}


/**
 * @brief Prepare a block of memory for network operations.
 *
 * Wrap a user allocated block of memory in an NNTI_buffer_t.  The transport
 * may take additional actions to prepare the memory for network send/receive.
 * If the memory block doesn't meet the transport's requirements for memory
 * regions, then errors or poor performance may result.
 */
NNTI_result_t NNTI_ptl_alloc (
        const NNTI_transport_t *trans_hdl,
        const uint64_t          element_size,
        const uint64_t          num_elements,
        const NNTI_buf_ops_t    ops,
        NNTI_buffer_t          *reg_buf)
{
    NNTI_result_t nnti_rc=NNTI_OK;

    log_debug(nnti_debug_level, "enter");

    assert(trans_hdl);
    assert(element_size>0);
    assert(num_elements>0);
    assert(ops>0);
    assert(reg_buf);

    char *buf=(char *)malloc(element_size*num_elements);
    assert(buf);

    nnti_rc=NNTI_ptl_register_memory(
            trans_hdl,
            buf,
            element_size,
            num_elements,
            ops,
            reg_buf);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "reg_buf",
                "end of NNTI_ptl_alloc", reg_buf);
    }

    log_debug(nnti_debug_level, "exit");

    return(nnti_rc);
}


/**
 * @brief Cleanup after network operations are complete.
 *
 * Destroy an NNTI_buffer_t that was previously created by NNTI_regsiter_buffer().
 * It is the user's responsibility to release the the memory region.
 */
NNTI_result_t NNTI_ptl_free (
        NNTI_buffer_t    *reg_buf)
{
    NNTI_result_t nnti_rc=NNTI_OK;

    log_debug(nnti_debug_level, "enter");

    assert(reg_buf);

    char *buf=NNTI_BUFFER_C_POINTER(reg_buf);
    assert(buf);

    nnti_rc=NNTI_ptl_unregister_memory(reg_buf);

    free(buf);

    log_debug(nnti_debug_level, "exit");

    return(nnti_rc);
}


/**
 * @brief Prepare a block of memory for network operations.
 *
 * Wrap a user allocated block of memory in an NNTI_buffer_t.  The transport
 * may take additional actions to prepare the memory for network send/receive.
 * If the memory block doesn't meet the transport's requirements for memory
 * regions, then errors or poor performance may result.
 */
NNTI_result_t NNTI_ptl_register_memory (
        const NNTI_transport_t *trans_hdl,
        char                   *buffer,
        const uint64_t          element_size,
        const uint64_t          num_elements,
        const NNTI_buf_ops_t    ops,
        NNTI_buffer_t          *reg_buf)
{
    int rc=0;
    NNTI_result_t nnti_rc=NNTI_OK;
    static uint64_t mbits=1;

    portals_memory_handle *ptl_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter");

    assert(trans_hdl);
    assert(buffer);
    assert(element_size>0);
    assert(num_elements>0);
    assert(ops>0);
    assert(reg_buf);

    ptl_mem_hdl=new portals_memory_handle();
    assert(ptl_mem_hdl);

    reg_buf->transport_id      = trans_hdl->id;
    reg_buf->buffer_owner      = trans_hdl->me;
    reg_buf->ops               = ops;
    reg_buf->payload_size      = element_size;
    reg_buf->payload           = (uint64_t)buffer;
    reg_buf->transport_private = (uint64_t)ptl_mem_hdl;

    log_debug(nnti_debug_level, "rpc_buffer->payload_size=%ld",
            reg_buf->payload_size);

    ptl_mem_hdl->eq_h=PTL_EQ_NONE;
    ptl_mem_hdl->me_h=0;
    ptl_mem_hdl->md_h=0;

    ptl_mem_hdl->match_id.nid = PTL_NID_ANY;
    ptl_mem_hdl->match_id.pid = PTL_PID_ANY;

    reg_buf->buffer_segments.NNTI_remote_addr_array_t_val=(NNTI_remote_addr_t *)calloc(1, sizeof(NNTI_remote_addr_t));
    reg_buf->buffer_segments.NNTI_remote_addr_array_t_len=1;

    reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].transport_id                      = reg_buf->transport_id;
    reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.portals.size = element_size;

    /*
     * Buffer types are divided into four groups here.
     *   Request Buffer - This buffer type requires immediate registration.
     *                    Events are generated, so an EQ is assigned.
     *   Receive Buffer - This buffer type requires immediate registration.
     *                    Events are generated, so an EQ is assigned.
     *   Send/RDMA Initiators - This buffer type allows for lazy registration.
     *                          The buffer will be register when the operation
     *                          is initiated.  Events are generated, so an EQ
     *                          is assigned.
     *   RDMA Targets - This buffer type requires immediate registration.
     *                  Events are NOT generated, so an EQ is NOT assigned.
     */

    if (ops == NNTI_RECV_QUEUE) {
        uint32_t index=0;
        portals_request_queue_handle *q_hdl=&transport_global_data.req_queue;

        ptl_mem_hdl->type=REQUEST_BUFFER;

        ptl_mem_hdl->buffer_id   = NNTI_REQ_PT_INDEX;
        ptl_mem_hdl->match_bits  = 0;
        ptl_mem_hdl->ignore_bits = 0;
        reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.portals.buffer_id  = (NNTI_portals_indices)ptl_mem_hdl->buffer_id;
        reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.portals.match_bits = ptl_mem_hdl->match_bits;

        q_hdl->reg_buf=reg_buf;

        q_hdl->req_size=element_size;
        q_hdl->reqs_per_queue=num_elements/NUM_REQ_QUEUES;

        /* create an event queue */
        nthread_lock(&nnti_ptl_lock);
        rc = PtlEQAlloc(
                transport_global_data.ni_h,
                num_elements*2,
                PTL_EQ_HANDLER_NONE,
                &transport_global_data.req_eq_h);
        nthread_unlock(&nnti_ptl_lock);
        if (rc != PTL_OK) {
            log_error(nnti_debug_level, "PtlEQAlloc() failed");
            nnti_rc=NNTI_ENOMEM;
            goto cleanup;
        }
        ptl_mem_hdl->eq_h=transport_global_data.req_eq_h;
        log_debug(nnti_debug_level, "allocated eq=%d", ptl_mem_hdl->eq_h);

        for (index=0; index<NUM_REQ_QUEUES; index++) {
            /* initialize the indices stored in the MD user pointer */
            q_hdl->indices[index] = index;
            q_hdl->queue_count[index] = 0;

            /* allocate the buffer for the incoming MD */
            q_hdl->req_queue[index] = buffer + (index*q_hdl->reqs_per_queue*q_hdl->req_size);

            /* initialize the buffer */
            memset(q_hdl->req_queue[index], 0, q_hdl->reqs_per_queue*q_hdl->req_size);

            /* initialize the MD */
            memset(&q_hdl->md[index], 0, sizeof(ptl_md_t));
            q_hdl->md[index].start = q_hdl->req_queue[index];
            q_hdl->md[index].length = q_hdl->reqs_per_queue*q_hdl->req_size;
            q_hdl->md[index].threshold = q_hdl->reqs_per_queue;
            q_hdl->md[index].max_size = q_hdl->req_size;
            q_hdl->md[index].options = PTL_MD_OP_PUT | PTL_MD_MAX_SIZE;
            q_hdl->md[index].user_ptr = &q_hdl->indices[index];  /* used to store the index */
            q_hdl->md[index].eq_handle = ptl_mem_hdl->eq_h;

            log_debug(nnti_debug_level, "attaching match entry to index=%d",
                    ptl_mem_hdl->buffer_id);

            /* Attach the match entry to the portal index */
            nthread_lock(&nnti_ptl_lock);
            rc = PtlMEAttach(
                    transport_global_data.ni_h,
                    ptl_mem_hdl->buffer_id,
                    ptl_mem_hdl->match_id,
                    ptl_mem_hdl->match_bits,
                    ptl_mem_hdl->ignore_bits,
                    PTL_RETAIN,
                    PTL_INS_AFTER,
                    &q_hdl->me_h[index]);
            nthread_unlock(&nnti_ptl_lock);
            if (rc != PTL_OK) {
                log_error(nnti_debug_level, "could not attach ME");
                nnti_rc=NNTI_ENOMEM;
                goto cleanup;
            }

            /* Attach the MD to the match entry */
            nthread_lock(&nnti_ptl_lock);
            rc = PtlMDAttach(
                    q_hdl->me_h[index],
                    q_hdl->md[index],
                    PTL_RETAIN,
                    &q_hdl->md_h[index]);
            nthread_unlock(&nnti_ptl_lock);
            if (rc != PTL_OK) {
                log_error(nnti_debug_level, "could not alloc eq: %s",
                        ptl_err_str[rc]);
                nnti_rc=NNTI_ENOMEM;
                goto cleanup;
            }
            log_debug(nnti_debug_level, "attached q_hdl->md_h[%d]: %d", index, q_hdl->md_h[index]);

            reg_buf->payload_size=q_hdl->req_size;
        }

//        post_recv_work_request(reg_buf);

    } else if (ops == NNTI_RECV_DST) {

        ptl_mem_hdl->type=RECEIVE_BUFFER;

        ptl_mem_hdl->buffer_id   = NNTI_RECV_PT_INDEX;
        ptl_mem_hdl->match_bits  = mbits++;
        ptl_mem_hdl->ignore_bits = 0;
        reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.portals.buffer_id  = (NNTI_portals_indices)ptl_mem_hdl->buffer_id;
        reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.portals.match_bits = ptl_mem_hdl->match_bits;

        ptl_mem_hdl->eq_h=transport_global_data.data_eq_h;

        /* create a match entry (unlink with MD) */
        nthread_lock(&nnti_ptl_lock);
        rc = PtlMEAttach(
                transport_global_data.ni_h,
                ptl_mem_hdl->buffer_id,
                ptl_mem_hdl->match_id,
                ptl_mem_hdl->match_bits,
                ptl_mem_hdl->ignore_bits,
                PTL_UNLINK,
                PTL_INS_AFTER,
                &ptl_mem_hdl->me_h);
        nthread_unlock(&nnti_ptl_lock);
        if (rc != PTL_OK) {
            log_error(nnti_debug_level, "failed to attach me");
            nnti_rc=NNTI_ENOMEM;
            goto cleanup;
        }
        log_debug(nnti_debug_level, "allocated me=%d with bufid=%d, match_id(%d,%d), mbits=%d",
                ptl_mem_hdl->me_h, ptl_mem_hdl->buffer_id, ptl_mem_hdl->match_id.nid, ptl_mem_hdl->match_id.pid, ptl_mem_hdl->match_bits);

        /* initialize the md */
        memset(&ptl_mem_hdl->md, 0, sizeof(ptl_md_t));
        ptl_mem_hdl->md.start     = buffer;
        ptl_mem_hdl->md.length    = element_size;
        ptl_mem_hdl->md.threshold = PTL_MD_THRESH_INF;
        ptl_mem_hdl->md.options   = PTL_MD_OP_PUT|PTL_MD_OP_GET|PTL_MD_MANAGE_REMOTE|PTL_MD_TRUNCATE;
        ptl_mem_hdl->md.user_ptr  = reg_buf;
        ptl_mem_hdl->md.eq_handle = ptl_mem_hdl->eq_h;

        /* attach the memory descriptor (manually unlink) */
        nthread_lock(&nnti_ptl_lock);
        rc = PtlMDAttach(
                ptl_mem_hdl->me_h,
                ptl_mem_hdl->md,
                PTL_RETAIN,
                &ptl_mem_hdl->md_h);
        nthread_unlock(&nnti_ptl_lock);
        if (rc != PTL_OK) {
            log_error(nnti_debug_level, "failed to attach md");
            nnti_rc=NNTI_ENOMEM;
            goto cleanup;
        }
        log_debug(nnti_debug_level, "attached ptl_mem_hdl->md_h: %d", ptl_mem_hdl->md_h);

//        post_recv_work_request(reg_buf);

    } else if ((ops == NNTI_SEND_SRC) ||
               (ops == NNTI_GET_DST) ||
               (ops == NNTI_PUT_SRC)) {

        /*
         * These are initiator buffers.
         */

        if (ops == NNTI_SEND_SRC) {
            ptl_mem_hdl->type=SEND_BUFFER;
        } else if (ops == NNTI_GET_DST) {
            ptl_mem_hdl->type=GET_DST_BUFFER;
        } else if (ops == NNTI_PUT_SRC) {
            ptl_mem_hdl->type=PUT_SRC_BUFFER;
        }

        ptl_mem_hdl->buffer_id   = NNTI_DATA_PT_INDEX;
        ptl_mem_hdl->match_bits  = mbits++;
        ptl_mem_hdl->ignore_bits = 0;
        reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.portals.buffer_id  = (NNTI_portals_indices)ptl_mem_hdl->buffer_id;
        reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.portals.match_bits = ptl_mem_hdl->match_bits;

        ptl_mem_hdl->eq_h=transport_global_data.data_eq_h;

    } else if ((ops == NNTI_GET_SRC) ||
               (ops == NNTI_PUT_DST)    ||
               (ops == (NNTI_GET_SRC|NNTI_PUT_DST))) {

        /*
         * These are eventless target buffers.
         */

        if (ops == NNTI_GET_SRC) {
            ptl_mem_hdl->type=GET_SRC_BUFFER;
        } else if (ops == NNTI_PUT_DST) {
            ptl_mem_hdl->type=PUT_DST_BUFFER;
        } else if (ops == (NNTI_GET_SRC|NNTI_PUT_DST)) {
            ptl_mem_hdl->type=RDMA_TARGET_BUFFER;
        }

        ptl_mem_hdl->buffer_id   = NNTI_DATA_PT_INDEX;
        ptl_mem_hdl->match_bits  = mbits++;
        ptl_mem_hdl->ignore_bits = 0;
        reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.portals.buffer_id  = (NNTI_portals_indices)ptl_mem_hdl->buffer_id;
        reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.portals.match_bits = ptl_mem_hdl->match_bits;

#if defined(USE_RDMA_TARGET_ACK)
        ptl_mem_hdl->eq_h = transport_global_data.data_eq_h;
        post_recv_work_request(reg_buf);
#else
        ptl_mem_hdl->eq_h = PTL_EQ_NONE;
#endif

        /* create a match entry (unlink with MD) */
        nthread_lock(&nnti_ptl_lock);
        rc = PtlMEAttach(
                transport_global_data.ni_h,
                ptl_mem_hdl->buffer_id,
                ptl_mem_hdl->match_id,
                ptl_mem_hdl->match_bits,
                ptl_mem_hdl->ignore_bits,
                PTL_UNLINK,
                PTL_INS_AFTER,
                &ptl_mem_hdl->me_h);
        nthread_unlock(&nnti_ptl_lock);
        if (rc != PTL_OK) {
            log_error(nnti_debug_level, "failed to attach me");
            nnti_rc=NNTI_ENOMEM;
            goto cleanup;
        }
        log_debug(nnti_debug_level, "allocated me=%d with bufid=%d, match_id(%d,%d), mbits=%d",
                ptl_mem_hdl->me_h, ptl_mem_hdl->buffer_id, ptl_mem_hdl->match_id.nid, ptl_mem_hdl->match_id.pid, ptl_mem_hdl->match_bits);

        /* initialize the md */
        memset(&ptl_mem_hdl->md, 0, sizeof(ptl_md_t));
        ptl_mem_hdl->md.start     = buffer;
        ptl_mem_hdl->md.length    = element_size;
        ptl_mem_hdl->md.threshold = PTL_MD_THRESH_INF;
        ptl_mem_hdl->md.options   = PTL_MD_OP_PUT|PTL_MD_OP_GET|PTL_MD_MANAGE_REMOTE|PTL_MD_TRUNCATE;
        ptl_mem_hdl->md.user_ptr  = reg_buf;
        ptl_mem_hdl->md.eq_handle = ptl_mem_hdl->eq_h;

        /* attach the memory descriptor (manually unlink) */
        nthread_lock(&nnti_ptl_lock);
        rc = PtlMDAttach(
                ptl_mem_hdl->me_h,
                ptl_mem_hdl->md,
                PTL_RETAIN,
                &ptl_mem_hdl->md_h);
        nthread_unlock(&nnti_ptl_lock);
        if (rc != PTL_OK) {
            log_error(nnti_debug_level, "failed to attach md");
            nnti_rc=NNTI_ENOMEM;
            goto cleanup;
        }
        log_debug(nnti_debug_level, "attached ptl_mem_hdl->md_h: %d", ptl_mem_hdl->md_h);
    }



cleanup:
    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "reg_buf",
                "end of NNTI_ptl_register_memory", reg_buf);
    }

    log_debug(nnti_debug_level, "exit");

    return(nnti_rc);
}


/**
 * @brief Prepare a list of memory segments for network operations.
 *
 * Wrap a list of user allocated memory segments in an NNTI_buffer_t.  The
 * transport may take additional actions to prepare the memory segments for
 * network send/receive.  If the memory segments don't meet the transport's
 * requirements for memory regions, then errors or poor performance may
 * result.
 *
 */
NNTI_result_t NNTI_ptl_register_segments (
        const NNTI_transport_t *trans_hdl,
        char                  **segments,
        const uint64_t         *segment_lengths,
        const uint64_t          num_segments,
        const NNTI_buf_ops_t    ops,
        NNTI_buffer_t          *reg_buf)
{
    return NNTI_OK;
}


/**
 * @brief Cleanup after network operations are complete.
 *
 * Destroy an NNTI_buffer_t that was previously created by NNTI_regsiter_buffer().
 * It is the user's responsibility to release the the memory region.
 */
NNTI_result_t NNTI_ptl_unregister_memory (
        NNTI_buffer_t    *reg_buf)
{
    int rc=0;
    NNTI_result_t nnti_rc=NNTI_OK;
    portals_memory_handle *ptl_mem_hdl=NULL;
    log_level debug_level = nnti_debug_level;

    log_debug(nnti_debug_level, "enter");

    assert(reg_buf);

    ptl_mem_hdl=PTL_MEM_HDL(reg_buf);

    assert(ptl_mem_hdl);

    log_debug(nnti_debug_level, "unregistering reg_buf(%p) buf(%p) md_h(%d) eq_h(%d)",
        reg_buf, reg_buf->payload, ptl_mem_hdl->md_h, ptl_mem_hdl->eq_h);

    if (reg_buf->ops == NNTI_RECV_QUEUE) {
        uint32_t index=0;
        portals_request_queue_handle *q_hdl=&transport_global_data.req_queue;

        for (index=0; index<NUM_REQ_QUEUES; index++) {
            /* unlink the memory descriptor */
            log_debug(debug_level, "unlinking q_hdl->md_h[%d]: %d", index, q_hdl->md_h[index]);
            nthread_lock(&nnti_ptl_lock);
            rc = PtlMDUnlink(q_hdl->md_h[index]);
            nthread_unlock(&nnti_ptl_lock);
            if (rc != PTL_OK) {
                log_warn(debug_level, "unable to unlink memory descriptor for queue %d: %s",
                        index, ptl_err_str[rc]);
                nnti_rc = NNTI_ENOMEM;
                goto cleanup;
            }
        }

        /* free the event queue */
        log_debug(debug_level, "freeing ptl_mem_hdl->eq_h: %d", ptl_mem_hdl->eq_h);
        nthread_lock(&nnti_ptl_lock);
        rc = PtlEQFree(transport_global_data.req_eq_h);
        nthread_unlock(&nnti_ptl_lock);
        transport_global_data.req_eq_h=PTL_EQ_NONE;
        if (rc != PTL_OK) {
            log_fatal(debug_level, "unable to free EQ: %s", ptl_err_str[rc]);
            nnti_rc = NNTI_ENOMEM;
            goto cleanup;
        }
    } else if ((reg_buf->ops == NNTI_RECV_DST) ||
               (reg_buf->ops == NNTI_GET_SRC)  ||
               (reg_buf->ops == NNTI_PUT_DST)  ||
               (reg_buf->ops == (NNTI_GET_SRC|NNTI_PUT_DST))) {

        log_debug(debug_level, "unlinking ptl_mem_hdl->md_h: %d", ptl_mem_hdl->md_h);
        nthread_lock(&nnti_ptl_lock);
        rc = PtlMDUnlink(ptl_mem_hdl->md_h);
        nthread_unlock(&nnti_ptl_lock);
        if (rc != PTL_OK) {
            log_error(debug_level, "failed to unlink MD: %s", ptl_err_str[rc]);
            nnti_rc = NNTI_ENOMEM;
            goto cleanup;
        }

        ptl_mem_hdl->eq_h=PTL_EQ_NONE;
    }


cleanup:

    if (ptl_mem_hdl)
        delete ptl_mem_hdl;
    if (reg_buf->buffer_segments.NNTI_remote_addr_array_t_val)
        free(reg_buf->buffer_segments.NNTI_remote_addr_array_t_val);

    reg_buf->transport_id      = NNTI_TRANSPORT_NULL;
    PORTALS_SET_MATCH_ANY(&reg_buf->buffer_owner);
    reg_buf->ops               = (NNTI_buf_ops_t)0;
    reg_buf->payload_size      = 0;
    reg_buf->payload           = 0;
    reg_buf->transport_private = 0;

    log_debug(debug_level, "Finished unregistering, rc=%d",rc);

    log_debug(nnti_debug_level, "exit");

    return(nnti_rc);
}


/**
 * @brief Send a message to a peer.
 *
 * Send a message (<tt>msg_hdl</tt>) to a peer (<tt>peer_hdl</tt>).  It is expected that the
 * message is small, but the exact maximum size is transport dependent.
 */
NNTI_result_t NNTI_ptl_send (
        const NNTI_peer_t   *peer_hdl,
        const NNTI_buffer_t *msg_hdl,
        const NNTI_buffer_t *dest_hdl,
        NNTI_work_request_t *wr)
{
    int rc=0;
    NNTI_result_t nnti_rc=NNTI_OK;

    portals_memory_handle *ptl_mem_hdl=NULL;
    portals_work_request  *ptl_wr=NULL;
    ptl_process_id_t dest_id;
    ptl_pt_index_t   buffer_id;
    ptl_match_bits_t match_bits;

    log_debug(nnti_debug_level, "enter");

    assert(peer_hdl);
    assert(msg_hdl);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "msg_hdl",
                "NNTI_ptl_send", msg_hdl);
    }
    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "dest_hdl",
                "NNTI_ptl_send", dest_hdl);
    }

    ptl_mem_hdl=PTL_MEM_HDL(msg_hdl);
    assert(ptl_mem_hdl);
    ptl_wr=(portals_work_request *)calloc(1, sizeof(portals_work_request));
    assert(ptl_wr);

    /* create a match entry (unlink with MD) */
    nthread_lock(&nnti_ptl_lock);
    rc = PtlMEAttach(
            transport_global_data.ni_h,
            ptl_mem_hdl->buffer_id,
            ptl_mem_hdl->match_id,
            ptl_mem_hdl->match_bits,
            ptl_mem_hdl->ignore_bits,
            PTL_UNLINK,
            PTL_INS_AFTER,
            &ptl_mem_hdl->me_h);
    nthread_unlock(&nnti_ptl_lock);
    if (rc != PTL_OK) {
        log_error(nnti_debug_level, "failed to attach me");
        nnti_rc=NNTI_ENOMEM;
        goto cleanup;
    }
    log_debug(nnti_debug_level, "allocated me=%d with bufid=%d, match_id(%d,%d), mbits=%d",
            ptl_mem_hdl->me_h, ptl_mem_hdl->buffer_id, ptl_mem_hdl->match_id.nid, ptl_mem_hdl->match_id.pid, ptl_mem_hdl->match_bits);

    /* initialize the md */
    memset(&ptl_mem_hdl->md, 0, sizeof(ptl_md_t));
    ptl_mem_hdl->md.start    =(void*)msg_hdl->payload;
    ptl_mem_hdl->md.length   =msg_hdl->payload_size;
    ptl_mem_hdl->md.threshold=2;
    ptl_mem_hdl->md.options  =PTL_MD_OP_PUT|PTL_MD_OP_GET|PTL_MD_MANAGE_REMOTE|PTL_MD_TRUNCATE;
    ptl_mem_hdl->md.user_ptr =wr;
    ptl_mem_hdl->md.eq_handle=ptl_mem_hdl->eq_h;

    /* attach the memory descriptor (manually unlink) */
    nthread_lock(&nnti_ptl_lock);
    rc = PtlMDAttach(
            ptl_mem_hdl->me_h,
            ptl_mem_hdl->md,
            PTL_UNLINK,
            &ptl_mem_hdl->md_h);
    nthread_unlock(&nnti_ptl_lock);
    if (rc != PTL_OK) {
        log_error(nnti_debug_level, "failed to attach md");
        nnti_rc=NNTI_ENOMEM;
        goto cleanup;
    }


    memset(&ptl_wr->op_state, 0, sizeof(ptl_op_state_t));

    if (dest_hdl == NULL) {
        ptl_wr->peer=*peer_hdl;
        dest_id.nid =peer_hdl->peer.NNTI_remote_process_t_u.portals.nid;
        dest_id.pid =peer_hdl->peer.NNTI_remote_process_t_u.portals.pid;
        buffer_id   =NNTI_REQ_PT_INDEX;
        match_bits  =0;
    } else {
        ptl_wr->peer=dest_hdl->buffer_owner;
        dest_id.nid =dest_hdl->buffer_owner.peer.NNTI_remote_process_t_u.portals.nid;
        dest_id.pid =dest_hdl->buffer_owner.peer.NNTI_remote_process_t_u.portals.pid;
        buffer_id   =dest_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.portals.buffer_id;
        match_bits  =dest_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.portals.match_bits;
    }

    ptl_wr->reg_buf   =(NNTI_buffer_t *)msg_hdl;
    ptl_wr->src_offset=0;
    ptl_wr->dst_offset=0;
    ptl_wr->length    =msg_hdl->payload_size;

    log_debug(nnti_debug_level, "sending to (nid=%d, pid=%d, buffer_id=%d, mbits=%d)", dest_id.nid, dest_id.pid, buffer_id, match_bits);

    rc=PtlPut(
            ptl_mem_hdl->md_h,
            PTL_ACK_REQ,
            dest_id,
            buffer_id,
            0,
            match_bits,
            0,
            0);
    if (rc != PTL_OK) {
        log_error(nnti_debug_level, "failed to send with PUT: %s", ptl_err_str[rc]);
        nnti_rc = NNTI_EBADRPC;
        goto cleanup;
    }

    ptl_wr->last_op=PTL_OP_SEND;

    wr->transport_id     =msg_hdl->transport_id;
    wr->ops              =NNTI_SEND_SRC;
    wr->transport_private=(uint64_t)ptl_wr;

    ptl_mem_hdl->wr_queue.push_back(wr);

cleanup:
    log_debug(nnti_debug_level, "exit");

    return(nnti_rc);
}


/**
 * @brief Transfer data to a peer.
 *
 * Put the contents of <tt>src_buffer_hdl</tt> into <tt>dest_buffer_hdl</tt>.  It is
 * assumed that the destination is at least <tt>src_length</tt> bytes in size.
 *
 */
NNTI_result_t NNTI_ptl_put (
        const NNTI_buffer_t *src_buffer_hdl,
        const uint64_t       src_offset,
        const uint64_t       src_length,
        const NNTI_buffer_t *dest_buffer_hdl,
        const uint64_t       dest_offset,
        NNTI_work_request_t *wr)
{
    int rc=0;
    NNTI_result_t nnti_rc=NNTI_OK;

    portals_memory_handle *ptl_mem_hdl=NULL;
    portals_work_request  *ptl_wr=NULL;
    ptl_process_id_t dest_id;

    log_debug(nnti_debug_level, "enter");

    assert(src_buffer_hdl);
    assert(dest_buffer_hdl);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "src_buffer_hdl",
                "NNTI_ptl_put", src_buffer_hdl);
    }
    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "dest_buffer_hdl",
                "NNTI_ptl_put", dest_buffer_hdl);
    }

    ptl_mem_hdl=PTL_MEM_HDL(src_buffer_hdl);
    assert(ptl_mem_hdl);
    ptl_wr=(portals_work_request *)calloc(1, sizeof(portals_work_request));
    assert(ptl_wr);

    /* create a match entry (unlink with MD) */
    nthread_lock(&nnti_ptl_lock);
    rc = PtlMEAttach(
            transport_global_data.ni_h,
            ptl_mem_hdl->buffer_id,
            ptl_mem_hdl->match_id,
            ptl_mem_hdl->match_bits,
            ptl_mem_hdl->ignore_bits,
            PTL_UNLINK,
            PTL_INS_AFTER,
            &ptl_mem_hdl->me_h);
    nthread_unlock(&nnti_ptl_lock);
    if (rc != PTL_OK) {
        log_error(nnti_debug_level, "failed to attach me");
        nnti_rc=NNTI_ENOMEM;
        goto cleanup;
    }
    log_debug(nnti_debug_level, "allocated me=%d with bufid=%d, match_id(%d,%d), mbits=%d",
            ptl_mem_hdl->me_h, ptl_mem_hdl->buffer_id, ptl_mem_hdl->match_id.nid, ptl_mem_hdl->match_id.pid, ptl_mem_hdl->match_bits);

    /* initialize the md */
    memset(&ptl_mem_hdl->md, 0, sizeof(ptl_md_t));
    ptl_mem_hdl->md.start     = (void*)src_buffer_hdl->payload;
    ptl_mem_hdl->md.length    = src_buffer_hdl->payload_size;
    ptl_mem_hdl->md.threshold = 2;
    ptl_mem_hdl->md.options   = PTL_MD_OP_PUT|PTL_MD_OP_GET|PTL_MD_MANAGE_REMOTE|PTL_MD_TRUNCATE;
    ptl_mem_hdl->md.user_ptr  = wr;
    ptl_mem_hdl->md.eq_handle = ptl_mem_hdl->eq_h;

    /* attach the memory descriptor (manually unlink) */
    nthread_lock(&nnti_ptl_lock);
    rc = PtlMDAttach(
            ptl_mem_hdl->me_h,
            ptl_mem_hdl->md,
            PTL_UNLINK,
            &ptl_mem_hdl->md_h);
    nthread_unlock(&nnti_ptl_lock);
    if (rc != PTL_OK) {
        log_error(nnti_debug_level, "failed to attach md");
        nnti_rc=NNTI_ENOMEM;
        goto cleanup;
    }


    ptl_wr->reg_buf   =(NNTI_buffer_t *)src_buffer_hdl;
    ptl_wr->peer      =dest_buffer_hdl->buffer_owner;
    ptl_wr->src_offset=src_offset;
    ptl_wr->dst_offset=dest_offset;
    ptl_wr->length    =src_length;
    ptl_wr->last_op   =PTL_OP_PUT_INITIATOR;

    memset(&ptl_wr->op_state, 0, sizeof(ptl_op_state_t));

    dest_id.nid=dest_buffer_hdl->buffer_owner.peer.NNTI_remote_process_t_u.portals.nid;
    dest_id.pid=dest_buffer_hdl->buffer_owner.peer.NNTI_remote_process_t_u.portals.pid;

    rc=PtlPutRegion(
            ptl_mem_hdl->md_h,
            src_offset,
            src_length,
            PTL_ACK_REQ,
            dest_id,
            dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.portals.buffer_id,
            0,
            dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.portals.match_bits,
            dest_offset,
            0);
    if (rc != PTL_OK) {
        log_error(nnti_debug_level, "failed to PUT region: %s", ptl_err_str[rc]);
        nnti_rc = NNTI_EBADRPC;
        goto cleanup;
    }

    log_debug(nnti_debug_level, "putting to (%s, eq=%d)", dest_buffer_hdl->buffer_owner.url, ptl_mem_hdl->eq_h);

    wr->transport_id     =src_buffer_hdl->transport_id;
    wr->ops              =NNTI_PUT_SRC;
    wr->transport_private=(uint64_t)ptl_wr;

    ptl_mem_hdl->wr_queue.push_back(wr);

cleanup:
    log_debug(nnti_debug_level, "exit");

    return(nnti_rc);
}


/**
 * @brief Transfer data from a peer.
 *
 * Get the contents of <tt>src_buffer_hdl</tt> into <tt>dest_buffer_hdl</tt>.  It is
 * assumed that the destination is at least <tt>src_length</tt> bytes in size.
 *
 */
NNTI_result_t NNTI_ptl_get (
        const NNTI_buffer_t *src_buffer_hdl,
        const uint64_t       src_offset,
        const uint64_t       src_length,
        const NNTI_buffer_t *dest_buffer_hdl,
        const uint64_t       dest_offset,
        NNTI_work_request_t *wr)
{
    int rc=0;
    NNTI_result_t nnti_rc=NNTI_OK;

    portals_memory_handle *ptl_mem_hdl=NULL;
    portals_work_request  *ptl_wr=NULL;
    ptl_process_id_t src_id;

    log_debug(nnti_debug_level, "enter");

    assert(src_buffer_hdl);
    assert(dest_buffer_hdl);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "src_buffer_hdl",
                "NNTI_ptl_get", src_buffer_hdl);
    }
    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "dest_buffer_hdl",
                "NNTI_ptl_get", dest_buffer_hdl);
    }

    log_debug(nnti_debug_level, "getting from (%s, src_offset=%llu, src_length=%llu, dest_offset=%llu)",
            src_buffer_hdl->buffer_owner.url, src_offset, src_length, dest_offset);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "src_buffer_hdl",
                "NNTI_ptl_get", src_buffer_hdl);
        fprint_NNTI_buffer(logger_get_file(), "dest_buffer_hdl",
                "NNTI_ptl_get", dest_buffer_hdl);
    }

    ptl_mem_hdl=PTL_MEM_HDL(dest_buffer_hdl);
    assert(ptl_mem_hdl);
    ptl_wr=(portals_work_request *)calloc(1, sizeof(portals_work_request));
    assert(ptl_wr);

    /* create a match entry (unlink with MD) */
    nthread_lock(&nnti_ptl_lock);
    rc = PtlMEAttach(
            transport_global_data.ni_h,
            ptl_mem_hdl->buffer_id,
            ptl_mem_hdl->match_id,
            ptl_mem_hdl->match_bits,
            ptl_mem_hdl->ignore_bits,
            PTL_UNLINK,
            PTL_INS_AFTER,
            &ptl_mem_hdl->me_h);
    nthread_unlock(&nnti_ptl_lock);
    if (rc != PTL_OK) {
        log_error(nnti_debug_level, "failed to attach me");
        nnti_rc=NNTI_ENOMEM;
        goto cleanup;
    }
    log_debug(nnti_debug_level, "allocated me=%d with bufid=%d, match_id(%d,%d), mbits=%d",
            ptl_mem_hdl->me_h, ptl_mem_hdl->buffer_id, ptl_mem_hdl->match_id.nid, ptl_mem_hdl->match_id.pid, ptl_mem_hdl->match_bits);

    /* initialize the md */
    memset(&ptl_mem_hdl->md, 0, sizeof(ptl_md_t));
    ptl_mem_hdl->md.start     = (void*)dest_buffer_hdl->payload;
    ptl_mem_hdl->md.length    = dest_buffer_hdl->payload_size;
    ptl_mem_hdl->md.threshold = 1;
    ptl_mem_hdl->md.options   = PTL_MD_OP_PUT|PTL_MD_OP_GET|PTL_MD_MANAGE_REMOTE|PTL_MD_TRUNCATE;
    ptl_mem_hdl->md.user_ptr  = wr;
    ptl_mem_hdl->md.eq_handle = ptl_mem_hdl->eq_h;

    /* attach the memory descriptor (manually unlink) */
    nthread_lock(&nnti_ptl_lock);
    rc = PtlMDAttach(
            ptl_mem_hdl->me_h,
            ptl_mem_hdl->md,
            PTL_UNLINK,
            &ptl_mem_hdl->md_h);
    nthread_unlock(&nnti_ptl_lock);
    if (rc != PTL_OK) {
        log_error(nnti_debug_level, "failed to attach md");
        nnti_rc=NNTI_ENOMEM;
        goto cleanup;
    }


    ptl_wr->reg_buf   =(NNTI_buffer_t *)dest_buffer_hdl;
    ptl_wr->peer      =src_buffer_hdl->buffer_owner;
    ptl_wr->src_offset=src_offset;
    ptl_wr->dst_offset=dest_offset;
    ptl_wr->length    =src_length;
    ptl_wr->last_op   =PTL_OP_GET_INITIATOR;


    memset(&ptl_wr->op_state, 0, sizeof(ptl_op_state_t));

    src_id.nid=src_buffer_hdl->buffer_owner.peer.NNTI_remote_process_t_u.portals.nid;
    src_id.pid=src_buffer_hdl->buffer_owner.peer.NNTI_remote_process_t_u.portals.pid;

    rc=PtlGetRegion(
            ptl_mem_hdl->md_h,
            dest_offset,
            src_length,
            src_id,
            src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.portals.buffer_id,
            0,
            src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.portals.match_bits,
            src_offset);
    if (rc != PTL_OK) {
        log_error(nnti_debug_level, "failed to GET region: %s", ptl_err_str[rc]);
        nnti_rc = NNTI_EBADRPC;
        goto cleanup;
    }

    log_debug(nnti_debug_level, "getting from (%s, eq=%d)", src_buffer_hdl->buffer_owner.url, ptl_mem_hdl->eq_h);

    wr->transport_id     =dest_buffer_hdl->transport_id;
    wr->ops              =NNTI_GET_DST;
    wr->transport_private=(uint64_t)ptl_wr;

    ptl_mem_hdl->wr_queue.push_back(wr);

cleanup:
    log_debug(nnti_debug_level, "exit");

    return(nnti_rc);
}


/**
 * @brief Transfer data to a peer.
 *
 * \param[in] src_buffer_hdl    A buffer containing the data to put.
 * \param[in] src_length        The number of bytes to put.
 * \param[in] dest_buffer_list  A list of buffers to put the data into.
 * \param[in] dest_count        The number of destination buffers.
 * \return A result code (NNTI_OK or an error)
 */
NNTI_result_t NNTI_ptl_scatter (
        const NNTI_buffer_t  *src_buffer_hdl,
        const uint64_t        src_length,
        const NNTI_buffer_t **dest_buffer_list,
        const uint64_t        dest_count,
        NNTI_work_request_t  *wr)
{
    return NNTI_OK;
}


/**
 * @brief Transfer data from a peer.
 *
 * \param[in] src_buffer_list  A list of buffers containing the data to get.
 * \param[in] src_length       The number of bytes to get.
 * \param[in] src_count        The number of source buffers.
 * \param[in] dest_buffer_hdl  A buffer to get the data into.
 * \return A result code (NNTI_OK or an error)
 */
NNTI_result_t NNTI_ptl_gather (
        const NNTI_buffer_t **src_buffer_list,
        const uint64_t        src_length,
        const uint64_t        src_count,
        const NNTI_buffer_t  *dest_buffer_hdl,
        NNTI_work_request_t  *wr)
{
    return NNTI_OK;
}


NNTI_result_t NNTI_ptl_atomic_set_callback (
		const NNTI_transport_t *trans_hdl,
		const uint64_t          local_atomic,
		NNTI_callback_fn_t      cbfunc,
		void                   *context)
{
    return NNTI_ENOTSUP;
}


NNTI_result_t NNTI_ptl_atomic_read (
		const NNTI_transport_t *trans_hdl,
		const uint64_t          local_atomic,
		int64_t                *value)
{
    return NNTI_ENOTSUP;
}


NNTI_result_t NNTI_ptl_atomic_fop (
		const NNTI_transport_t *trans_hdl,
		const NNTI_peer_t      *peer_hdl,
		const uint64_t          target_atomic,
		const uint64_t          result_atomic,
		const int64_t           operand,
		const NNTI_atomic_op_t  op,
		NNTI_work_request_t    *wr)
{
    return NNTI_ENOTSUP;
}


NNTI_result_t NNTI_ptl_atomic_cswap (
		const NNTI_transport_t *trans_hdl,
		const NNTI_peer_t      *peer_hdl,
		const uint64_t          target_atomic,
		const uint64_t          result_atomic,
		const int64_t           compare_operand,
		const int64_t           swap_operand,
		NNTI_work_request_t    *wr)
{
    return NNTI_ENOTSUP;
}


/**
 * @brief Create a receive work request that can be used to wait for buffer
 * operations to complete.
 *
 */
NNTI_result_t NNTI_ptl_create_work_request (
        NNTI_buffer_t        *reg_buf,
        NNTI_work_request_t  *wr)
{
    portals_work_request *ptl_wr=NULL;
    portals_memory_handle *ptl_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    ptl_mem_hdl=PTL_MEM_HDL(reg_buf);
    assert(ptl_mem_hdl);

    ptl_wr=(portals_work_request *)calloc(1, sizeof(portals_work_request));
    assert(ptl_wr);
    ptl_wr->reg_buf = reg_buf;

    memset(&ptl_wr->op_state, 0, sizeof(ptl_op_state_t));

    wr->transport_id     =reg_buf->transport_id;
    wr->reg_buf          =reg_buf;
    wr->ops              =reg_buf->ops;
    wr->transport_private=(uint64_t)ptl_wr;

    ptl_mem_hdl->wr_queue.push_back(wr);

    log_debug(nnti_debug_level, "exit (reg_buf=%p)", reg_buf);

    return(NNTI_OK);
}


/**
 * @brief Disassociates a receive work request from a previous receive
 * and prepares it for reuse.
 *
 */
NNTI_result_t NNTI_ptl_clear_work_request (
        NNTI_work_request_t  *wr)
{
    portals_work_request  *ptl_wr=NULL;
    portals_memory_handle *ptl_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (wr=%p)", wr);

    ptl_wr=PTL_WORK_REQUEST(wr);
    assert(ptl_wr);
    ptl_mem_hdl=PTL_MEM_HDL(ptl_wr->reg_buf);
    assert(ptl_mem_hdl);

    memset(&ptl_wr->op_state, 0, sizeof(ptl_op_state_t));

    ptl_mem_hdl->wr_queue.push_back(wr);

    log_debug(nnti_debug_level, "exit (wr=%p)", wr);

    return(NNTI_OK);
}


/**
 * @brief Disassociates a receive work request from reg_buf.
 *
 */
NNTI_result_t NNTI_ptl_destroy_work_request (
        NNTI_work_request_t  *wr)
{
    portals_work_request  *ptl_wr=NULL;
    portals_memory_handle *ptl_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (wr=%p)", wr);

    ptl_wr=PTL_WORK_REQUEST(wr);
    assert(ptl_wr);
    ptl_mem_hdl=PTL_MEM_HDL(ptl_wr->reg_buf);
    assert(ptl_mem_hdl);

    free(ptl_wr);

    log_debug(nnti_debug_level, "exit (wr=%p)", wr);

    return(NNTI_OK);
}


/**
 * @brief Attempts to cancel an NNTI opertion.
 *
 */
NNTI_result_t NNTI_ptl_cancel (
        NNTI_work_request_t *wr)
{
    return NNTI_OK;
}


/**
 * @brief Attempts to cancel a list of NNTI opertions.
 *
 */
NNTI_result_t NNTI_ptl_cancelall (
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count)
{
    return NNTI_OK;
}


/**
 * @brief Interrupts NNTI_wait*()
 *
 */
NNTI_result_t NNTI_ptl_interrupt (
        const NNTI_transport_t *trans_hdl)
{
    char dummy=0xAA;

    log_debug(nnti_debug_level, "enter");

    log_debug(nnti_debug_level, "exit");

    return NNTI_OK;
}


/**
 * @brief Wait for <tt>remote_op</tt> on <tt>reg_buf</tt> to complete.
 *
 * Wait for <tt>remote_op</tt> on <tt>reg_buf</tt> to complete or timeout
 * waiting.  This is typically used to wait for a result or a bulk data
 * transfer.  The timeout is specified in milliseconds.  A timeout of <tt>-1</tt>
 * means wait forever.  A timeout of <tt>0</tt> means do not wait.
 *
 */
NNTI_result_t NNTI_ptl_wait (
        NNTI_work_request_t  *wr,
        const int             timeout,
        NNTI_status_t        *status)
{
    int rc=PTL_OK;
    NNTI_result_t nnti_rc=NNTI_OK;
    portals_memory_handle *ptl_mem_hdl=NULL;
    portals_work_request  *ptl_wr=NULL;

//    const NNTI_buffer_t  *reg_buf=NULL;

    const NNTI_work_request_t *wait_wr=NULL;

    int elapsed_time=0;
    int timeout_per_call;
    ptl_event_t event;
    int which_eq=0;

    log_level debug_level=nnti_debug_level;

    trios_declare_timer(call_time);
    trios_declare_timer(total_time);

    trios_start_timer(total_time);

    log_debug(debug_level, "enter");

    assert(wr);
    assert(status);

    ptl_wr=PTL_WORK_REQUEST(wr);
    assert(ptl_wr);
    ptl_mem_hdl=PTL_MEM_HDL(ptl_wr->reg_buf);
    assert(ptl_mem_hdl);
//    ptl_wr=first_incomplete_wr(ptl_mem_hdl);
//    assert(ptl_wr);

    if (ptl_mem_hdl->type == REQUEST_BUFFER) {
        memset(&ptl_wr->op_state, 0, sizeof(ptl_op_state_t));
    }

    if (is_wr_complete(ptl_wr) == TRUE) {
        log_debug(debug_level, "work request already complete");
        nnti_rc = NNTI_OK;
    } else {
        log_debug(debug_level, "work request NOT complete");

        timeout_per_call = MIN_TIMEOUT;

        while (1)   {
            if (trios_exit_now()) {
                log_debug(debug_level, "caught abort signal");
                return NNTI_ECANCELED;
            }

            log_debug(debug_level, "waiting on wr(%p) eq_h(%d)", wr, ptl_mem_hdl->eq_h);

            memset(&event, 0, sizeof(ptl_event_t));
            log_debug(debug_level, "lock before poll");
            trios_start_timer(call_time);
            //            nthread_lock(&nnti_ptl_lock);
            rc = PtlEQPoll(&ptl_mem_hdl->eq_h, 1, timeout_per_call, &event, &which_eq);
            //            nthread_unlock(&nnti_ptl_lock);
            trios_stop_timer("NNTI_ptl_wait - PtlEQPoll", call_time);
            log_debug(debug_level, "polling status is %s", ptl_err_str[rc]);

            log_debug(debug_level, "Poll Event= {");
            log_debug(debug_level, "\ttype         = %d", event.type);
            log_debug(debug_level, "\tinitiator    = (%llu, %llu)", (unsigned long long)event.initiator.nid, (unsigned long long)event.initiator.pid);
            log_debug(debug_level, "\tuid          = %d", event.uid);
            log_debug(debug_level, "\tjid          = %d", event.jid);
            log_debug(debug_level, "\tpt_index     = %d", event.pt_index);
            log_debug(debug_level, "\tmatch_bits   = %d", event.match_bits);
            log_debug(debug_level, "\trlength      = %llu", (unsigned long long)event.rlength);
            log_debug(debug_level, "\tmlength      = %llu", (unsigned long long)event.mlength);
            log_debug(debug_level, "\toffset       = %llu", (unsigned long long)event.offset);
            log_debug(debug_level, "\tmd_handle    = %d", event.md_handle);
            log_debug(debug_level, "\tmd.start     = %p", event.md.start);
            log_debug(debug_level, "\tmd.length    = %d", event.md.length);
            log_debug(debug_level, "\tmd.max_size  = %d", event.md.max_size);
            log_debug(debug_level, "\tmd.threshold = %d", event.md.threshold);
            log_debug(debug_level, "\tmd.user_ptr  = %p", event.md.user_ptr);
            log_debug(debug_level, "}");


            /* case 1: success */
            if (rc == PTL_OK) {
                nnti_rc = NNTI_OK;
            }
            /* case 2: success, but some events were dropped */
            else if (rc == PTL_EQ_DROPPED) {
                log_warn(debug_level, "PtlEQPoll dropped some events");
                log_warn(debug_level, "PtlEQPoll succeeded, but at least one event was dropped");
                nnti_rc = NNTI_OK;
            }
            /* case 3: timed out */
            else if (rc == PTL_EQ_EMPTY) {
                elapsed_time += timeout_per_call;

                /* if the caller asked for a legitimate timeout, we need to exit */
                if (((timeout > 0) && (elapsed_time >= timeout))) {
                    log_debug(debug_level, "PtlEQPoll timed out: %s",
                            ptl_err_str[rc]);
                    nnti_rc = NNTI_ETIMEDOUT;
                    break;
                }
                /* continue if the timeout has not expired */
                /* log_debug(debug_level, "timedout... continuing"); */



                continue;
            }
            /* case 4: failure */
            else {
                log_error(debug_level, "PtlEQPoll failed (eq_handle[%d]==%d): %s",
                        which_eq, ptl_mem_hdl->eq_h, ptl_err_str[rc]);
                nnti_rc = NNTI_EIO;
                break;
            }

            wait_wr=decode_event_wr(wr, &event);
            process_event(wait_wr, &event);

            if (is_wr_complete(ptl_wr) == TRUE) {
                break;
            }
        }
    }

    create_status(wr, nnti_rc, status);

    if (logging_debug(debug_level)) {
        fprint_NNTI_status(logger_get_file(), "status",
                "end of NNTI_ptl_wait", status);
    }

    if ((nnti_rc==NNTI_OK) && (ptl_mem_hdl->buffer_id == NNTI_REQ_PT_INDEX)) {
        portals_request_queue_handle *q_hdl=&transport_global_data.req_queue;

        int index = *(int *)ptl_wr->last_event.md.user_ptr;
        /* get the index of the queue */
        q_hdl->queue_count[index]++;

        log_debug(debug_level, "queue_count[%d]: %d", index, q_hdl->queue_count[index]);

        /* if we've processed all we can on this queue, reset */
        if (q_hdl->queue_count[index] >= q_hdl->reqs_per_queue) {

            log_debug(debug_level, "Resetting MD on queue[%d]", index);

            /* Unlink the ME (also unlinks the MD) */
            nthread_lock(&nnti_ptl_lock);
            rc = PtlMEUnlink(q_hdl->me_h[index]);
            nthread_unlock(&nnti_ptl_lock);
            if (rc != PTL_OK) {
                log_error(debug_level, "Could not unlink ME: %s", ptl_err_str[rc]);
                goto cleanup;
            }

            /* Re-attach the match-list entry */
            nthread_lock(&nnti_ptl_lock);
            rc = PtlMEAttach(
                    transport_global_data.ni_h,
                    ptl_mem_hdl->buffer_id,
                    ptl_mem_hdl->match_id,
                    ptl_mem_hdl->match_bits,
                    ptl_mem_hdl->ignore_bits,
                    PTL_RETAIN,
                    PTL_INS_AFTER,
                    &q_hdl->me_h[index]);
            nthread_unlock(&nnti_ptl_lock);
            if (rc != PTL_OK) {
                log_error(debug_level, "Could not reset ME: %s", ptl_err_str[rc]);
                goto cleanup;
            }

            /* Re-attach the MD */
            nthread_lock(&nnti_ptl_lock);
            rc = PtlMDAttach(
                    q_hdl->me_h[index],
                    q_hdl->md[index],
                    PTL_RETAIN,
                    &q_hdl->md_h[index]);
            nthread_unlock(&nnti_ptl_lock);
            if (rc != PTL_OK) {
                log_error(debug_level, "Could not reset MD: %s", ptl_err_str[rc]);
                goto cleanup;
            }

            q_hdl->queue_count[index] = 0;
        }
    }

    if (nnti_rc==NNTI_OK) {
        ptl_mem_hdl=PTL_MEM_HDL(ptl_wr->reg_buf);
        assert(ptl_mem_hdl);
//        ptl_wr=ptl_mem_hdl->wr_queue.front();
//        assert(ptl_wr);
//        ptl_mem_hdl->wr_queue.pop_front();

        switch (ptl_mem_hdl->type) {
            case REQUEST_BUFFER:
            case RECEIVE_BUFFER:
#if defined(USE_RDMA_TARGET_ACK)
            case GET_SRC_BUFFER:
            case PUT_DST_BUFFER:
            case RDMA_TARGET_BUFFER:
#endif
//                repost_recv_work_request(wr);
                break;
            case SEND_BUFFER:
            case GET_DST_BUFFER:
            case PUT_SRC_BUFFER:
                free(ptl_wr);
                break;
            case UNKNOWN_BUFFER:
            default:
                log_error(nnti_debug_level, "unknown buffer type(%llu).", ptl_mem_hdl->type);
                break;
        }
    }

cleanup:
    trios_stop_timer("NNTI_ptl_wait", total_time);
    log_debug(debug_level, "exit");
    return(nnti_rc);
}

/**
 * @brief Wait for <tt>remote_op</tt> on any buffer in <tt>buf_list</tt> to complete.
 *
 * Wait for <tt>remote_op</tt> on any buffer in <tt>buf_list</tt> to complete or timeout
 * waiting.  This is typically used to wait for a result or a bulk data
 * transfer.  The timeout is specified in milliseconds.  A timeout of <tt>-1</tt>
 * means wait forever.  A timeout of <tt>0</tt> means do not wait.
 *
 * Caveats:
 *   1) All buffers in buf_list must be registered with the same transport.
 *   2) You can't wait on the request queue and RDMA buffers in the same call.  Will probably be fixed in the future.
 */
NNTI_result_t NNTI_ptl_waitany (
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count,
        const int             timeout,
        uint32_t             *which,
        NNTI_status_t        *status)
{
    int rc=PTL_OK;
    NNTI_result_t nnti_rc=NNTI_OK;
    portals_memory_handle *ptl_mem_hdl=NULL;
    portals_work_request  *ptl_wr=NULL;

//    const NNTI_buffer_t  *reg_buf=NULL;
    const NNTI_work_request_t *wait_wr=NULL;

    int elapsed_time=0;
    int timeout_per_call;
    ptl_event_t event;
    int which_eq=0;

    log_level debug_level=nnti_debug_level;

    log_debug(debug_level, "enter");

    assert(wr_list);
    assert(wr_count > 0);
//    if (buf_count > 1) {
//        /* if there is more than 1 buffer in the list, none of them can be a REQUEST_BUFFER */
//        for (uint32_t i=0;i<buf_count;i++) {
//            if (buf_list[i] != NULL) {
//                assert(((portals_memory_handle *)buf_list[i]->transport_private)->type != REQUEST_BUFFER);
//            }
//        }
//    }
    assert(status);

    if (wr_count == 1) {
        nnti_rc=NNTI_ptl_wait(wr_list[0], timeout, status);
        *which=0;
        goto cleanup;
    }

    if (is_any_wr_complete(wr_list, wr_count, which) == TRUE) {
        log_debug(debug_level, "work request already complete (which=%u, wr_list[%d]=%p)", *which, *which, wr_list[*which]);
        nnti_rc = NNTI_OK;
    } else {
        log_debug(debug_level, "work request NOT complete (buf_list=%p)", wr_list);

        timeout_per_call = MIN_TIMEOUT;

        while (1)   {
            if (trios_exit_now()) {
                log_debug(debug_level, "caught abort signal");
                return NNTI_ECANCELED;
            }

            log_debug(debug_level, "waiting on eq_h(%d)", transport_global_data.data_eq_h);

            memset(&event, 0, sizeof(ptl_event_t));
            log_debug(debug_level, "lock before poll");
            //        nthread_lock(&nnti_ptl_lock);
            rc = PtlEQPoll(&transport_global_data.data_eq_h, 1, timeout_per_call, &event, &which_eq);
            //        nthread_unlock(&nnti_ptl_lock);
            log_debug(debug_level, "polling status is %s", ptl_err_str[rc]);

            log_debug(debug_level, "Poll Event= {");
            log_debug(debug_level, "\ttype         = %d", event.type);
            log_debug(debug_level, "\tinitiator    = (%llu, %llu)", (unsigned long long)event.initiator.nid, (unsigned long long)event.initiator.pid);
            log_debug(debug_level, "\tuid          = %d", event.uid);
            log_debug(debug_level, "\tjid          = %d", event.jid);
            log_debug(debug_level, "\tpt_index     = %d", event.pt_index);
            log_debug(debug_level, "\tmatch_bits   = %d", event.match_bits);
            log_debug(debug_level, "\trlength      = %llu", (unsigned long long)event.rlength);
            log_debug(debug_level, "\tmlength      = %llu", (unsigned long long)event.mlength);
            log_debug(debug_level, "\toffset       = %llu", (unsigned long long)event.offset);
            log_debug(debug_level, "\tmd_handle    = %d", event.md_handle);
            log_debug(debug_level, "\tmd.start     = %p", event.md.start);
            log_debug(debug_level, "\tmd.length    = %d", event.md.length);
            log_debug(debug_level, "\tmd.max_size  = %d", event.md.max_size);
            log_debug(debug_level, "\tmd.threshold = %d", event.md.threshold);
            log_debug(debug_level, "\tmd.user_ptr  = %p", event.md.user_ptr);


            /* case 1: success */
            if (rc == PTL_OK) {
                nnti_rc = NNTI_OK;
            }
            /* case 2: success, but some events were dropped */
            else if (rc == PTL_EQ_DROPPED) {
                log_warn(debug_level, "PtlEQPoll dropped some events");
                log_warn(debug_level, "PtlEQPoll succeeded, but at least one event was dropped");
                nnti_rc = NNTI_OK;
            }
            /* case 3: timed out */
            else if (rc == PTL_EQ_EMPTY) {
                elapsed_time += timeout_per_call;

                /* if the caller asked for a legitimate timeout, we need to exit */
                if (((timeout > 0) && (elapsed_time >= timeout))) {
                    log_debug(debug_level, "PtlEQPoll timed out: %s",
                            ptl_err_str[rc]);
                    nnti_rc = NNTI_ETIMEDOUT;
                    break;
                }
                /* continue if the timeout has not expired */
                /* log_debug(debug_level, "timedout... continuing"); */



                continue;
            }
            /* case 4: failure */
            else {
                log_error(debug_level, "PtlEQPoll failed (eq_handle[%d]==%d): %s",
                        which_eq, transport_global_data.data_eq_h, ptl_err_str[rc]);
                nnti_rc = NNTI_EIO;
                break;
            }

            wait_wr=decode_event_wr(wr_list[0], &event);
            process_event(wait_wr, &event);

            if (is_any_wr_complete(wr_list, wr_count, which) == TRUE) {
                break;
            }
        }
    }


    create_status(wr_list[*which], nnti_rc, status);

    if (logging_debug(debug_level)) {
        fprint_NNTI_status(logger_get_file(), "status",
                "end of NNTI_ptl_wait", status);
    }

    if (nnti_rc==NNTI_OK) {
        ptl_wr=PTL_WORK_REQUEST(wr_list[*which]);
        assert(ptl_wr);
        ptl_mem_hdl=PTL_MEM_HDL(ptl_wr->reg_buf);
        assert(ptl_mem_hdl);
//        ptl_wr=ptl_mem_hdl->wr_queue.front();
//        assert(ptl_wr);
//        ptl_mem_hdl->wr_queue.pop_front();

        switch (ptl_mem_hdl->type) {
            case REQUEST_BUFFER:
            case RECEIVE_BUFFER:
#if defined(USE_RDMA_TARGET_ACK)
            case GET_SRC_BUFFER:
            case PUT_DST_BUFFER:
            case RDMA_TARGET_BUFFER:
#endif
//                repost_recv_work_request(wr_list[*which]);
                break;
            case SEND_BUFFER:
            case GET_DST_BUFFER:
            case PUT_SRC_BUFFER:
                free(ptl_wr);
                break;
            case UNKNOWN_BUFFER:
            default:
                log_error(nnti_debug_level, "unknown buffer type(%llu).", ptl_mem_hdl->type);
                break;
        }
    }

cleanup:
    log_debug(debug_level, "exit");
    return(nnti_rc);
}

/**
 * @brief Wait for <tt>remote_op</tt> on all buffers in <tt>buf_list</tt> to complete.
 *
 * Wait for <tt>remote_op</tt> on all buffers in <tt>buf_list</tt> to complete or timeout
 * waiting.  This is typically used to wait for a result or a bulk data
 * transfer.  The timeout is specified in milliseconds.  A timeout of <tt>-1</tt>
 * means wait forever.  A timeout of <tt>0</tt> means do not wait.
 *
 * Caveats:
 *   1) All buffers in buf_list must be registered with the same transport.
 *   2) You can't wait on the receive queue and RDMA buffers in the same call.  Will probably be fixed in the future.
 */
NNTI_result_t NNTI_ptl_waitall (
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count,
        const int             timeout,
        NNTI_status_t       **status)
{
    int rc=PTL_OK;
    NNTI_result_t nnti_rc=NNTI_OK;
    portals_memory_handle *ptl_mem_hdl=NULL;
    portals_work_request  *ptl_wr=NULL;

    const NNTI_work_request_t  *wait_wr=NULL;

    int elapsed_time=0;
    int timeout_per_call;
    ptl_event_t event;
    int which_eq=0;

    log_level debug_level=nnti_debug_level;

    log_debug(debug_level, "enter");

    assert(wr_list);
    assert(wr_count > 0);
//    if (buf_count > 1) {
//        /* if there is more than 1 buffer in the list, none of them can be a REQUEST_BUFFER */
//        for (uint32_t i=0;i<buf_count;i++) {
//            if (buf_list[i] != NULL) {
//                assert(((portals_memory_handle *)buf_list[i]->transport_private)->type != REQUEST_BUFFER);
//            }
//        }
//    }
    assert(status);

    if (wr_count == 1) {
        nnti_rc=NNTI_ptl_wait(wr_list[0], timeout, status[0]);
        goto cleanup;
    }

    if (is_all_wr_complete(wr_list, wr_count) == TRUE) {
        log_debug(debug_level, "all buffer ops already complete (buf_list=%p)", wr_list);
        nnti_rc = NNTI_OK;
    } else {
        log_debug(debug_level, "all buffer ops NOT complete (buf_list=%p)", wr_list);

        timeout_per_call = MIN_TIMEOUT;

        while (1)   {
            if (trios_exit_now()) {
                log_debug(debug_level, "caught abort signal");
                return NNTI_ECANCELED;
            }

            log_debug(debug_level, "waiting on eq_h(%d)", transport_global_data.data_eq_h);

            memset(&event, 0, sizeof(ptl_event_t));
            log_debug(debug_level, "lock before poll");
            //        nthread_lock(&nnti_ptl_lock);
            rc = PtlEQPoll(&transport_global_data.data_eq_h, 1, timeout_per_call, &event, &which_eq);
            //        nthread_unlock(&nnti_ptl_lock);
            log_debug(debug_level, "polling status is %s", ptl_err_str[rc]);

            log_debug(debug_level, "Poll Event= {");
            log_debug(debug_level, "\ttype         = %d", event.type);
            log_debug(debug_level, "\tinitiator    = (%llu, %llu)", (unsigned long long)event.initiator.nid, (unsigned long long)event.initiator.pid);
            log_debug(debug_level, "\tuid          = %d", event.uid);
            log_debug(debug_level, "\tjid          = %d", event.jid);
            log_debug(debug_level, "\tpt_index     = %d", event.pt_index);
            log_debug(debug_level, "\tmatch_bits   = %d", event.match_bits);
            log_debug(debug_level, "\trlength      = %llu", (unsigned long long)event.rlength);
            log_debug(debug_level, "\tmlength      = %llu", (unsigned long long)event.mlength);
            log_debug(debug_level, "\toffset       = %llu", (unsigned long long)event.offset);
            log_debug(debug_level, "\tmd_handle    = %d", event.md_handle);
            log_debug(debug_level, "\tmd.start     = %p", event.md.start);
            log_debug(debug_level, "\tmd.length    = %d", event.md.length);
            log_debug(debug_level, "\tmd.max_size  = %d", event.md.max_size);
            log_debug(debug_level, "\tmd.threshold = %d", event.md.threshold);
            log_debug(debug_level, "\tmd.user_ptr  = %p", event.md.user_ptr);


            /* case 1: success */
            if (rc == PTL_OK) {
                nnti_rc = NNTI_OK;
            }
            /* case 2: success, but some events were dropped */
            else if (rc == PTL_EQ_DROPPED) {
                log_warn(debug_level, "PtlEQPoll dropped some events");
                log_warn(debug_level, "PtlEQPoll succeeded, but at least one event was dropped");
                nnti_rc = NNTI_OK;
            }
            /* case 3: timed out */
            else if (rc == PTL_EQ_EMPTY) {
                elapsed_time += timeout_per_call;

                /* if the caller asked for a legitimate timeout, we need to exit */
                if (((timeout > 0) && (elapsed_time >= timeout))) {
                    log_debug(debug_level, "PtlEQPoll timed out: %s",
                            ptl_err_str[rc]);
                    nnti_rc = NNTI_ETIMEDOUT;
                    break;
                }
                /* continue if the timeout has not expired */
                /* log_debug(debug_level, "timedout... continuing"); */



                continue;
            }
            /* case 4: failure */
            else {
                log_error(debug_level, "PtlEQPoll failed (eq_handle[%d]==%d): %s",
                        which_eq, transport_global_data.data_eq_h, ptl_err_str[rc]);
                nnti_rc = NNTI_EIO;
                break;
            }

            wait_wr=decode_event_wr(wr_list[0], &event);
            process_event(wait_wr, &event);

            if (is_all_wr_complete(wr_list, wr_count) == TRUE) {
                break;
            }
        }
    }


    for (uint32_t i=0;i<wr_count;i++) {
        create_status(wr_list[i], nnti_rc, status[i]);

        if (logging_debug(debug_level)) {
            fprint_NNTI_status(logger_get_file(), "status[i]",
                    "end of NNTI_ptl_wait", status[i]);
        }

        ptl_wr=PTL_WORK_REQUEST(wr_list[i]);
        assert(ptl_wr);
        ptl_mem_hdl=PTL_MEM_HDL(ptl_wr->reg_buf);
        assert(ptl_mem_hdl);
//        ptl_wr=ptl_mem_hdl->wr_queue.front();
//        assert(ptl_wr);
//        ptl_mem_hdl->wr_queue.pop_front();

        if (nnti_rc==NNTI_OK) {
            switch (ptl_mem_hdl->type) {
                case REQUEST_BUFFER:
                case RECEIVE_BUFFER:
#if defined(USE_RDMA_TARGET_ACK)
                case GET_SRC_BUFFER:
                case PUT_DST_BUFFER:
                case RDMA_TARGET_BUFFER:
#endif
//                    repost_recv_work_request(wr_list[i]);
                    break;
                case SEND_BUFFER:
                case GET_DST_BUFFER:
                case PUT_SRC_BUFFER:
                    free(ptl_wr);
                    break;
                case UNKNOWN_BUFFER:
                default:
                    log_error(nnti_debug_level, "unknown buffer type(%llu).", ptl_mem_hdl->type);
                    break;
            }
        }
    }

cleanup:
    log_debug(debug_level, "exit");
    return(nnti_rc);
}

/**
 * @brief Disable this transport.
 *
 * Shutdown the transport.  Any outstanding sends, gets and puts will be
 * canceled.  Any new transport requests will fail.
 *
 */
NNTI_result_t NNTI_ptl_fini (
        const NNTI_transport_t *trans_hdl)
{
//    PtlFini();
    nthread_lock_fini(&nnti_ptl_lock);

    ptl_initialized=false;

#if !defined(HAVE_TRIOS_CRAYPORTALS) && defined(HAVE_TRIOS_MPI)
    if (transport_global_data.init_called_mpi_init) {
    	MPI_Finalize();
    }
#endif

    return(NNTI_OK);
}




static void set_req_pid(NNTI_pid *pid)
{
    log_debug(nnti_debug_level, "enter (pid=%d)", *pid);

    transport_global_data.init_called_mpi_init=false;

    *pid=PTL_PID_ANY;

#if !defined(HAVE_TRIOS_CRAYPORTALS) && defined(HAVE_TRIOS_MPI)
    // Schutt's Portals doesn't properly assign a pid when you pass PTL_PID_ANY.  Fix it here.
    int initialized=0;
    int rank;
    MPI_Initialized(&initialized);
    if (!initialized) {
        MPI_Init(NULL, NULL);
        transport_global_data.init_called_mpi_init=true;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    log_debug(nnti_debug_level, "rank=%d", rank);

    *pid = 128 + rank;
#endif

    log_debug(nnti_debug_level, "exit (pid=%d)", *pid);
}


//static portals_work_request *decode_work_request(
//        const ptl_event_t   *event)
//{
//    log_level debug_level = nnti_debug_level;
//
//    const NNTI_buffer_t *event_buf=NULL;
//    portals_memory_handle *ptl_mem_hdl=NULL;
//    NNTI_work_request_t   *wr=NULL;
//    portals_work_request  *ptl_wr=NULL;
//    portals_work_request  *debug_wr=NULL;
//
//    log_debug(debug_level, "enter");
//
//    event_buf=(NNTI_buffer_t *)event->md.user_ptr;
//    assert(event_buf);
//    ptl_mem_hdl=PTL_MEM_HDL(event_buf);
//    assert(ptl_mem_hdl);
//
//    wr_queue_iter_t i;
//    for (i=ptl_mem_hdl->wr_queue.begin(); i != ptl_mem_hdl->wr_queue.end(); i++) {
//        assert(*i);
//        ptl_wr=PTL_WORK_REQUEST(*i);
//        if (is_wr_complete(ptl_wr) == FALSE) {
//            // work request is incomplete, check if it matches this event
//            switch(ptl_mem_hdl->type) {
//                case REQUEST_BUFFER:
//                    if (((ptl_wr)->src_offset == event->offset) &&
//                        ((ptl_wr)->length == event->mlength)) {
//
//                        wr=*i;
//                    } else {
//                        log_debug(debug_level, "work request doesn't match (wr=%p)", ptl_wr);
//                    }
//                    break;
//                case SEND_BUFFER:
//                case PUT_SRC_BUFFER:
//                    if (((ptl_wr)->src_offset == event->offset) &&
//                        ((ptl_wr)->length == event->mlength)) {
//
//                        wr=*i;
//                    } else {
//                        log_debug(debug_level, "work request doesn't match (wr=%p)", ptl_wr);
//                    }
//                    break;
//                case GET_DST_BUFFER:
//                    if (((ptl_wr)->dst_offset == event->offset) &&
//                        ((ptl_wr)->length == event->mlength)) {
//
//                        wr=*i;
//                    } else {
//                        log_debug(debug_level, "work request doesn't match (wr=%p)", ptl_wr);
//                    }
//                    break;
//                case RECEIVE_BUFFER:
//                case GET_SRC_BUFFER:
//                case PUT_DST_BUFFER:
//                case RDMA_TARGET_BUFFER:
//                    wr=*i;
//                    break;
//                default:
//                    log_debug(debug_level, "unknown event type %d (event_buf==%p)", ptl_mem_hdl->type, event_buf);
//                    break;
//            }
//            if (wr) {
//                break;
//            }
//        } else {
//            log_debug(debug_level, "work request is already complete (wr=%p)", ptl_wr);
//        }
//    }
//
//    if (!wr) {
//        for (i=ptl_mem_hdl->wr_queue.begin(); i != ptl_mem_hdl->wr_queue.end(); i++) {
//            debug_wr=PTL_WORK_REQUEST(*i);
//            log_debug(LOG_ALL, "e.offset=%llu, e.rlength=%llu, e.mlength=%llu, wr=%p, wr.length=%llu, wr.src_offset=%llu, wr.dst_offset=%llu, wr.is_complete=%d",
//                    (uint64_t)event->offset, (uint64_t)event->rlength, (uint64_t)event->mlength,
//                    debug_wr,
//                    (uint64_t)debug_wr->length, (uint64_t)debug_wr->src_offset, (uint64_t)debug_wr->dst_offset,
//                    (is_wr_complete(debug_wr)==TRUE));
//        }
//    }
//    assert(wr);
//
//    log_debug(debug_level, "exit (wr==%p)", wr);
//
//    return(PTL_WORK_REQUEST(wr));
//}

static const NNTI_work_request_t *decode_event_wr(
        const NNTI_work_request_t *wait_wr,
        const ptl_event_t         *event)
{
    const NNTI_buffer_t       *event_buf=NULL;
    const NNTI_work_request_t *event_wr=NULL;

    const NNTI_buffer_t       *wait_buf=NULL;
    portals_work_request      *ptl_wr=NULL;

    portals_memory_handle *ptl_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter");

    assert(wait_wr);

    ptl_wr=PTL_WORK_REQUEST(wait_wr);
    wait_buf=ptl_wr->reg_buf;

    switch (event->type) {
        /* these are the events for a PUT initiator */
        case PTL_EVENT_SEND_START:
            log_debug(nnti_debug_level, "got PTL_EVENT_SEND_START - event->user_ptr is a work request");
            event_wr=(NNTI_work_request_t *)event->md.user_ptr;
            break;
        case PTL_EVENT_SEND_END:
            log_debug(nnti_debug_level, "got PTL_EVENT_SEND_END   - event->user_ptr is a work request");
            event_wr=(NNTI_work_request_t *)event->md.user_ptr;
            break;
        case PTL_EVENT_ACK:
            log_debug(nnti_debug_level, "got PTL_EVENT_ACK        - event->user_ptr is a work request");
            event_wr=(NNTI_work_request_t *)event->md.user_ptr;
            break;
        case PTL_EVENT_UNLINK:
            log_debug(nnti_debug_level, "got PTL_EVENT_UNLINK     - event->user_ptr is a work request");
            event_wr=(NNTI_work_request_t *)event->md.user_ptr;
            break;


        /* these are the events for a GET initiator */
        case PTL_EVENT_REPLY_START:
            log_debug(nnti_debug_level,"got PTL_EVENT_REPLY_START - event->user_ptr is a work request");
            event_wr=(NNTI_work_request_t *)event->md.user_ptr;
            break;
        case PTL_EVENT_REPLY_END:
            log_debug(nnti_debug_level,"got PTL_EVENT_REPLY_END   - event->user_ptr is a work request");
            event_wr=(NNTI_work_request_t *)event->md.user_ptr;
            break;


        /* these are the events for a PUT target */
        case PTL_EVENT_PUT_START:
            log_debug(nnti_debug_level, "got PTL_EVENT_PUT_START  - ");
            if (event->pt_index == NNTI_REQ_PT_INDEX) {
                event_wr=wait_wr;
                log_debug(nnti_debug_level, "the wait work request is a REQUEST BUFFER wr, so event.md.user_ptr is the index of the request buffer.");
            } else {
                event_buf=(NNTI_buffer_t *)event->md.user_ptr;
                ptl_mem_hdl=PTL_MEM_HDL(event_buf);
                assert(ptl_mem_hdl);
                event_wr=ptl_mem_hdl->wr_queue.front();
            }
            break;
        case PTL_EVENT_PUT_END:
            log_debug(nnti_debug_level, "got PTL_EVENT_PUT_END    - ");
            if (event->pt_index == NNTI_REQ_PT_INDEX) {
                event_wr=wait_wr;
                log_debug(nnti_debug_level, "the wait work request is a REQUEST BUFFER wr, so event.md.user_ptr is the index of the request buffer.");
            } else {
                event_buf=(NNTI_buffer_t *)event->md.user_ptr;
                ptl_mem_hdl=PTL_MEM_HDL(event_buf);
                assert(ptl_mem_hdl);
                event_wr=ptl_mem_hdl->wr_queue.front();
            }
            break;


        case PTL_EVENT_GET_START:
            log_debug(nnti_debug_level, "got PTL_EVENT_GET_START  - GET targets should NOT generate events");
            break;
        case PTL_EVENT_GET_END:
            log_debug(nnti_debug_level, "got PTL_EVENT_GET_END    - GET targets should NOT generate events");
            break;


        default:
            log_error(nnti_debug_level, "unrecognized event type: %d - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                    event->type, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
//            rc = NNTI_EINVAL;

    }

//    if ((ptl_wr != NULL) &&
//        (PTL_MEM_HDL(ptl_wr->reg_buf)->type == REQUEST_BUFFER)) {
//        /* if the buffer is a request queue, then the wait_wr is the event_wr */
//        event_wr=wait_wr;
//        log_debug(nnti_debug_level, "the wait work request is a REQUEST BUFFER wr, so event.md.user_ptr is the index of the request buffer.");
//    } else if ((event->initiator.nid==transport_global_data.me.nid) && (event->initiator.pid==transport_global_data.me.pid)) {
//        /* if I am the initiator, then the event user_ptr is a work request pointer */
//
//        log_debug(nnti_debug_level, "I am the initiator, so the event user_ptr is a work request pointer");
//
//        event_wr=(NNTI_work_request_t *)event->md.user_ptr;
//        ptl_mem_hdl=PTL_MEM_HDL(PTL_WORK_REQUEST(event_wr)->reg_buf);
//        assert(ptl_mem_hdl);
//
//        if (event_wr == wait_wr) {
//            log_debug(nnti_debug_level, "the wc matches the wait work request (eq=%d, user_ptr=%p, wait_wr=%p)",
//                    ptl_mem_hdl->eq_h, (void *)event->md.user_ptr, wait_wr);
//        } else {
//            log_debug(nnti_debug_level, "the wc does NOT match the wait buffer (eq=%d, user_ptr=%p, wait_wr=%p)",
//                    ptl_mem_hdl->eq_h, (void *)event->md.user_ptr, wait_wr);
//        }
//    } else {
//        /* if I am NOT the initiator, then the event user_ptr is a pointer to the target buffer */
//        event_buf=(NNTI_buffer_t *)event->md.user_ptr;
//        ptl_mem_hdl=PTL_MEM_HDL(event_buf);
//        assert(ptl_mem_hdl);
//
//        if (event_buf == PTL_WORK_REQUEST(wait_wr)->reg_buf) {
//            log_debug(nnti_debug_level, "the wc matches the wait buffer (eq=%d, user_ptr=%p, wait_buf=%p)",
//                    ptl_mem_hdl->eq_h, (void *)event->md.user_ptr, wait_buf);
//        } else {
//            log_debug(nnti_debug_level, "the wc does NOT match the wait buffer (eq=%d, user_ptr=%p, wait_buf=%p)",
//                    ptl_mem_hdl->eq_h, (void *)event->md.user_ptr, wait_buf);
//        }
//    }

    log_debug(nnti_debug_level, "exit (event_wr==%p)", event_wr);

    return(event_wr);
}


static int process_event(
        const NNTI_work_request_t *wr,
        const ptl_event_t         *event)
{
    int rc=NNTI_OK;

    portals_memory_handle *ptl_mem_hdl=NULL;
    portals_work_request  *ptl_wr      =NULL;

    log_level debug_level = nnti_debug_level;

    ptl_wr=PTL_WORK_REQUEST(wr);
    assert(ptl_wr);
    ptl_mem_hdl=PTL_MEM_HDL(ptl_wr->reg_buf);
    assert(ptl_mem_hdl);

//    if (ptl_mem_hdl->type != REQUEST_BUFFER) {
//        ptl_wr=decode_work_request(event);
//    } else {
//        ptl_wr=ptl_mem_hdl->wr_queue.front();
//    }
    if (ptl_mem_hdl->type == REQUEST_BUFFER) {
        NNTI_work_request_t *tmp_wr=NULL;
        tmp_wr=ptl_mem_hdl->wr_queue.front();
        assert(tmp_wr);
        ptl_wr=PTL_WORK_REQUEST(tmp_wr);
    }
    assert(ptl_wr);

    ptl_wr->last_event=*event;

    log_debug(debug_level, "wr=%p; ptl_wr=%p; ptl_wr->last_op=%d", wr, ptl_wr, ptl_wr->last_op);
    switch (ptl_mem_hdl->type) {
        case SEND_BUFFER:
        case PUT_SRC_BUFFER:
            switch (event->type) {
                case PTL_EVENT_SEND_START:
                    log_debug(debug_level, "got PTL_EVENT_SEND_START - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptl_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptl_wr->op_state.put_initiator.send_start = TRUE;
                    break;
                case PTL_EVENT_SEND_END:
                    log_debug(debug_level, "got PTL_EVENT_SEND_END   - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptl_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptl_wr->op_state.put_initiator.send_end = TRUE;
                    break;
                case PTL_EVENT_ACK:
                    log_debug(debug_level, "got PTL_EVENT_ACK        - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptl_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptl_wr->op_state.put_initiator.ack = TRUE;
                    break;
                case PTL_EVENT_UNLINK:
                    log_debug(debug_level, "got PTL_EVENT_UNLINK     - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptl_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptl_wr->op_state.put_initiator.unlink = TRUE;
                    break;
                default:
                    log_error(debug_level, "unrecognized event type: %d - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            event->type, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    rc = NNTI_EINVAL;
                    goto cleanup;
            }
            break;
        case GET_DST_BUFFER:
            switch (event->type) {
                case PTL_EVENT_SEND_START:
                    log_debug(debug_level, "got PTL_EVENT_SEND_START - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptl_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptl_wr->op_state.get_initiator.send_start = TRUE;
                    break;
                case PTL_EVENT_SEND_END:
                    log_debug(debug_level, "got PTL_EVENT_SEND_END   - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptl_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptl_wr->op_state.get_initiator.send_end = TRUE;
                    break;
                case PTL_EVENT_REPLY_START:
                    log_debug(debug_level,"got PTL_EVENT_REPLY_START - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptl_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptl_wr->op_state.get_initiator.reply_start = TRUE;
                    break;
                case PTL_EVENT_REPLY_END:
                    log_debug(debug_level,"got PTL_EVENT_REPLY_END   - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptl_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptl_wr->op_state.get_initiator.reply_end = TRUE;
                    break;
                case PTL_EVENT_UNLINK:
                    log_debug(debug_level, "got PTL_EVENT_UNLINK     - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptl_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptl_wr->op_state.get_initiator.unlink = TRUE;
                    break;
                default:
                    log_error(debug_level, "unrecognized event type: %d - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            event->type, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    rc = NNTI_EINVAL;
                    goto cleanup;
            }
            break;
        case REQUEST_BUFFER:
        case RECEIVE_BUFFER:
            switch (event->type) {
                case PTL_EVENT_PUT_START:
                    log_debug(debug_level, "got PTL_EVENT_PUT_START  - new request - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptl_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);

                    break;
                case PTL_EVENT_PUT_END:
                    log_debug(debug_level, "got PTL_EVENT_PUT_END    - new request - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptl_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptl_wr->op_state.put_target.put_start = TRUE;
                    ptl_wr->op_state.put_target.put_end = TRUE;
                    break;
                default:
                    log_error(debug_level, "unrecognized event type: %d - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            event->type, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    rc = NNTI_EINVAL;
                    goto cleanup;
            }
            break;
        case PUT_DST_BUFFER:
            switch (event->type) {
                case PTL_EVENT_PUT_START:
                    log_debug(debug_level, "got PTL_EVENT_PUT_START  - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptl_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptl_wr->op_state.put_target.put_start = TRUE;
                    break;
                case PTL_EVENT_PUT_END:
                    log_debug(debug_level, "got PTL_EVENT_PUT_END    - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptl_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptl_wr->op_state.put_target.put_end = TRUE;
                    break;
                default:
                    log_error(debug_level, "unrecognized event type: %d - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            event->type, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    rc = NNTI_EINVAL;
                    goto cleanup;
            }
            break;
        case GET_SRC_BUFFER:
            switch (event->type) {
                case PTL_EVENT_GET_START:
                    log_debug(debug_level, "got PTL_EVENT_GET_START  - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptl_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptl_wr->op_state.get_target.get_start = TRUE;
                    break;
                case PTL_EVENT_GET_END:
                    log_debug(debug_level, "got PTL_EVENT_GET_END    - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptl_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptl_wr->op_state.get_target.get_end = TRUE;
                    break;
                default:
                    log_error(debug_level, "unrecognized event type: %d - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            event->type, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    rc = NNTI_EINVAL;
                    goto cleanup;
            }
            break;
        case RDMA_TARGET_BUFFER:
            switch (event->type) {
                case PTL_EVENT_PUT_START:
                    log_debug(debug_level, "got PTL_EVENT_PUT_START  - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptl_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptl_wr->op_state.put_target.put_start = TRUE;
                    break;
                case PTL_EVENT_PUT_END:
                    log_debug(debug_level, "got PTL_EVENT_PUT_END    - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptl_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptl_wr->op_state.put_target.put_end = TRUE;
                    break;
                case PTL_EVENT_GET_START:
                    log_debug(debug_level, "got PTL_EVENT_GET_START  - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptl_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptl_wr->op_state.get_target.get_start = TRUE;
                    break;
                case PTL_EVENT_GET_END:
                    log_debug(debug_level, "got PTL_EVENT_GET_END    - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            ptl_mem_hdl->eq_h, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    ptl_wr->op_state.get_target.get_end = TRUE;
                    break;
                default:
                    log_error(debug_level, "unrecognized event type: %d - event arrived on eq %d - initiator = (%4llu, %4llu, %4d)",
                            event->type, (unsigned long long)event->initiator.nid,(unsigned long long)event->initiator.pid, event->link);
                    rc = NNTI_EINVAL;
                    goto cleanup;
            }
            break;
        case UNKNOWN_BUFFER:
        default:
            break;
    }

    if (event->ni_fail_type != PTL_NI_OK) {
        log_error(debug_level, "failed on put end: ni_fail_type=%d\n",
                event->ni_fail_type);
        rc = event->ni_fail_type;
    }

cleanup:
    return (rc);
}

//static NNTI_result_t post_recv_work_request(
//        NNTI_buffer_t *reg_buf)
//{
//    portals_work_request *ptl_wr=NULL;
//    portals_memory_handle *ptl_mem_hdl=NULL;
//
//    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);
//
//    ptl_mem_hdl=PTL_MEM_HDL(reg_buf);
//    assert(ptl_mem_hdl);
//
//    ptl_wr=(portals_work_request *)calloc(1, sizeof(portals_work_request));
//    assert(ptl_wr);
//    ptl_wr->reg_buf = reg_buf;
//
//    memset(&ptl_wr->op_state, 0, sizeof(ptl_op_state_t));
//
//    ptl_mem_hdl->wr_queue.push_back(ptl_wr);
//
//    log_debug(nnti_debug_level, "exit (reg_buf=%p)", reg_buf);
//
//    return(NNTI_OK);
//}

//static NNTI_result_t repost_recv_work_request(
//        NNTI_work_request_t *wr)
//{
//    portals_work_request  *ptl_wr=NULL;
//    portals_memory_handle *ptl_mem_hdl=NULL;
//
//    log_debug(nnti_debug_level, "enter (wr=%p)", wr);
//
//    ptl_wr=PTL_WORK_REQUEST(wr);
//    assert(ptl_wr);
//    ptl_mem_hdl=PTL_MEM_HDL(ptl_wr->reg_buf);
//    assert(ptl_mem_hdl);
//
//    memset(&ptl_wr->op_state, 0, sizeof(ptl_op_state_t));
//
//    ptl_mem_hdl->wr_queue.push_back(ptl_wr);
//
//    log_debug(nnti_debug_level, "exit (wr=%p)", wr);
//
//    return(NNTI_OK);
//}

static int8_t is_wr_complete(
        portals_work_request *wr)
{
    int rc=FALSE;
    portals_memory_handle *ptl_mem_hdl=NULL;
//    log_level debug_level = nnti_debug_level;

    log_debug(nnti_debug_level, "enter (wr=%p)", wr);

    ptl_mem_hdl=PTL_MEM_HDL(wr->reg_buf);
    assert(ptl_mem_hdl);

    switch (ptl_mem_hdl->type) {
        case SEND_BUFFER:
        case PUT_SRC_BUFFER:
            if ((wr->op_state.put_initiator.send_start==TRUE) &&
                (wr->op_state.put_initiator.send_end==TRUE)   &&
                (wr->op_state.put_initiator.ack==TRUE)        &&
                (wr->op_state.put_initiator.unlink==TRUE))    {
                wr->last_op=PTL_OP_PUT_INITIATOR;
                rc = TRUE;
            }
        case GET_DST_BUFFER:
            /* cray portals */
            if ((wr->op_state.get_initiator.send_start==TRUE)  &&
                (wr->op_state.get_initiator.send_end==TRUE)    &&
                (wr->op_state.get_initiator.reply_start==TRUE) &&
                (wr->op_state.get_initiator.reply_end==TRUE)   &&
                (wr->op_state.get_initiator.unlink==TRUE))  {
                wr->last_op=PTL_OP_GET_INITIATOR;
                rc = TRUE;
                break;
            }
            /* schutt portals */
            if ((wr->op_state.get_initiator.reply_start==TRUE) &&
                (wr->op_state.get_initiator.reply_end==TRUE)   &&
                (wr->op_state.get_initiator.unlink==TRUE))  {
                wr->last_op=PTL_OP_GET_INITIATOR;
                rc = TRUE;
                break;
            }
            break;
        case PUT_DST_BUFFER:
            if ((wr->op_state.put_target.put_start==TRUE) &&
                (wr->op_state.put_target.put_end==TRUE))  {
                wr->last_op=PTL_OP_PUT_TARGET;
                rc = TRUE;
            }
            break;
        case GET_SRC_BUFFER:
            if ((wr->op_state.get_target.get_start==TRUE) &&
                (wr->op_state.get_target.get_end==TRUE))  {
                wr->last_op=PTL_OP_GET_TARGET;
                rc = TRUE;
            }
            break;
        case REQUEST_BUFFER:
            if ((wr->op_state.put_target.put_start==TRUE) &&
                (wr->op_state.put_target.put_end==TRUE)) {
                wr->last_op=PTL_OP_NEW_REQUEST;
                rc = TRUE;
            }
            break;
        case RECEIVE_BUFFER:
            if ((wr->op_state.put_target.put_start==TRUE) &&
                (wr->op_state.put_target.put_end==TRUE)) {
                wr->last_op=PTL_OP_RECEIVE;
                rc = TRUE;
            }
            break;
        case RDMA_TARGET_BUFFER:
            if ((wr->op_state.get_target.get_start==TRUE) &&
                (wr->op_state.get_target.get_end==TRUE))  {
                wr->last_op=PTL_OP_GET_TARGET;
                rc = TRUE;
            }
            if ((wr->op_state.put_target.put_start==TRUE) &&
                (wr->op_state.put_target.put_end==TRUE)) {
                wr->last_op=PTL_OP_PUT_TARGET;
                rc = TRUE;
            }
            break;
        case UNKNOWN_BUFFER:
        default:
            break;
    }

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);

    return(rc);
}

static int8_t is_any_wr_complete(
        portals_work_request **wr_list,
        const uint32_t         wr_count,
        uint32_t              *which)
{
    int8_t rc=FALSE;

    log_debug(nnti_debug_level, "enter");

    for (uint32_t i=0;i<wr_count;i++) {
        if ((wr_list[i] != NULL) &&
            (is_wr_complete(wr_list[i]) == TRUE)) {

            *which=i;
            rc = TRUE;
            break;
        }
    }

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);

    return(rc);
}

static int8_t is_all_wr_complete(
        portals_work_request **wr_list,
        const uint32_t         wr_count)
{
    int8_t rc=TRUE;

    log_debug(nnti_debug_level, "enter");

    for (uint32_t i=0;i<wr_count;i++) {
        if ((wr_list[i] != NULL) &&
            (is_wr_complete(wr_list[i]) == FALSE)) {

            rc = FALSE;
            break;
        }
    }

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);

    return(rc);
}

static int8_t is_wr_complete(
        NNTI_work_request_t *wr)
{
    portals_work_request *ptl_wr=NULL;

    log_debug(nnti_debug_level, "enter (wr=%p)", wr);

    ptl_wr=PTL_WORK_REQUEST(wr);
    assert(ptl_wr);

    return(is_wr_complete(ptl_wr));
}

static int8_t is_any_wr_complete(
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count,
        uint32_t             *which)
{
    int8_t rc=FALSE;

    log_debug(nnti_debug_level, "enter");

    for (uint32_t i=0;i<wr_count;i++) {
        if ((wr_list[i] != NULL) &&
            (is_wr_complete(wr_list[i]) == TRUE)) {

            *which=i;
            rc = TRUE;
            break;
        }
    }

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);

    return(rc);
}

static int8_t is_all_wr_complete(
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count)
{
    int8_t rc=TRUE;

    log_debug(nnti_debug_level, "enter");

    for (uint32_t i=0;i<wr_count;i++) {
        if ((wr_list[i] != NULL) &&
            (is_wr_complete(wr_list[i]) == FALSE)) {

            rc = FALSE;
            break;
        }
    }

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);

    return(rc);
}

//static portals_work_request *first_incomplete_wr(
//        portals_memory_handle *ptl_mem_hdl)
//{
//    NNTI_work_request_t   *wr    =NULL;
//    portals_work_request  *ptl_wr=NULL;
//
//    log_debug(nnti_debug_level, "enter");
//
//    assert(ptl_mem_hdl);
//
//    if (ptl_mem_hdl->wr_queue.empty()) {
//        log_debug(nnti_debug_level, "work request queue is empty");
//    } else {
//        wr_queue_iter_t i;
//        for (i=ptl_mem_hdl->wr_queue.begin(); i != ptl_mem_hdl->wr_queue.end(); i++) {
//            wr=*i;
//            assert(wr);
//            if (is_wr_complete(wr) == FALSE) {
//                break;
//            }
//        }
//    }
//
//    log_debug(nnti_debug_level, "exit (wr=%p)", wr);
//    return(PTL_WORK_REQUEST(wr));
//}

//static int8_t is_wr_queue_empty(
//        const NNTI_buffer_t *reg_buf)
//{
//    int8_t rc=FALSE;
//    portals_memory_handle *ptl_mem_hdl=NULL;
//
//    log_debug(nnti_debug_level, "enter");
//
//    ptl_mem_hdl=PTL_MEM_HDL(reg_buf);
//    assert(ptl_mem_hdl);
//
//    if (ptl_mem_hdl->wr_queue.empty()) {
//        rc=TRUE;
//    }
//
//    log_debug(nnti_debug_level, "exit (rc=%d)", rc);
//    return(rc);
//}
//
//
//static int8_t is_buf_op_complete(
//        const NNTI_buffer_t *reg_buf)
//{
//    int8_t rc=FALSE;
//    portals_memory_handle *ptl_mem_hdl=NULL;
//    NNTI_work_request_t   *wr    =NULL;
//    portals_work_request  *ptl_wr=NULL;
////    log_level debug_level = nnti_debug_level;
//
//    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);
//
//    ptl_mem_hdl=(portals_memory_handle *)reg_buf->transport_private;
//    assert(ptl_mem_hdl);
//
//    if (is_wr_queue_empty(reg_buf) == TRUE) {
//        log_debug(nnti_debug_level, "work request queue is empty - return FALSE");
//        rc=FALSE;
//    } else {
//        wr=ptl_mem_hdl->wr_queue.front();
//        assert(wr);
//
//        rc = is_wr_complete(wr);
//    }
//
//    if (rc==TRUE) {
//        log_debug(nnti_debug_level, "op is complete");
//    }
//    log_debug(nnti_debug_level, "exit (reg_buf=%p)", reg_buf);
//
//    return(rc);
//}
//
//static int8_t is_any_buf_op_complete(
//        const NNTI_buffer_t **buf_list,
//        const uint32_t        buf_count,
//        uint32_t             *which)
//{
//    int8_t rc=FALSE;
//
//    log_debug(nnti_debug_level, "enter");
//
//    for (uint32_t i=0;i<buf_count;i++) {
//        if ((buf_list[i] != NULL) &&
//            (is_wr_queue_empty(buf_list[i]) == FALSE) &&
//            (is_buf_op_complete(buf_list[i]) == TRUE)) {
//
//            *which=i;
//            rc = TRUE;
//            break;
//        }
//    }
//
//    log_debug(nnti_debug_level, "exit (rc=%d)", rc);
//
//    return(rc);
//}
//
//static int8_t is_all_buf_ops_complete(
//        const NNTI_buffer_t **buf_list,
//        const uint32_t        buf_count)
//{
//    int8_t rc=TRUE;
//
//    log_debug(nnti_debug_level, "enter");
//
//    for (uint32_t i=0;i<buf_count;i++) {
//        if ((buf_list[i] != NULL) &&
//            (is_wr_queue_empty(buf_list[i]) == FALSE) &&
//            (is_buf_op_complete(buf_list[i]) == FALSE)) {
//
//            rc = FALSE;
//            break;
//        }
//    }
//
//    log_debug(nnti_debug_level, "exit (rc=%d)", rc);
//
//    return(rc);
//}

static void create_status(
        const NNTI_work_request_t  *wr,
        int                         nnti_rc,
        NNTI_status_t              *status)
{
    portals_memory_handle *ptl_mem_hdl=NULL;
    portals_work_request  *ptl_wr      =NULL;

    status->op     = wr->ops;
    status->result = (NNTI_result_t)nnti_rc;
    if (nnti_rc==NNTI_OK) {
        ptl_wr=(portals_work_request *)wr->transport_private;
        assert(ptl_wr);
        ptl_mem_hdl=(portals_memory_handle *)ptl_wr->reg_buf->transport_private;
        assert(ptl_mem_hdl);

        status->start  = (uint64_t)ptl_wr->last_event.md.start;
        status->offset = ptl_wr->last_event.offset;
        status->length = ptl_wr->last_event.mlength;
        switch (ptl_wr->last_op) {
            case PTL_OP_PUT_INITIATOR:
            case PTL_OP_GET_TARGET:
            case PTL_OP_SEND:
                create_peer(&status->src, transport_global_data.me.nid, transport_global_data.me.pid); // allocates url
                create_peer(&status->dest, ptl_wr->last_event.initiator.nid, ptl_wr->last_event.initiator.pid); // allocates url
                break;
            case PTL_OP_GET_INITIATOR:
            case PTL_OP_PUT_TARGET:
            case PTL_OP_NEW_REQUEST:
            case PTL_OP_RECEIVE:
                create_peer(&status->src, ptl_wr->last_event.initiator.nid, ptl_wr->last_event.initiator.pid); // allocates url
                create_peer(&status->dest, transport_global_data.me.nid, transport_global_data.me.pid); // allocates url
                break;
        }
    }
}

static void create_peer(NNTI_peer_t *peer, ptl_nid_t nid, ptl_pid_t pid)
{
    log_debug(nnti_debug_level, "enter");

    sprintf(peer->url, "ptl://%u:%u/", nid, pid);

    peer->peer.transport_id                        = NNTI_TRANSPORT_PORTALS;
    peer->peer.NNTI_remote_process_t_u.portals.nid = nid;
    peer->peer.NNTI_remote_process_t_u.portals.pid = pid;

    log_debug(nnti_debug_level, "exit");
}
