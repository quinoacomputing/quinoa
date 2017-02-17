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
/**  @file rpc_server.c
 *
 *   @brief Implementation of the \ref rpc_server_api "RPC server API".
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov)
 *   $Revision: 1654 $
 *   $Date: 2007-12-11 22:57:34 -0700 (Tue, 11 Dec 2007) $
 *
 */


#include "Trios_config.h"
#include "Trios_nssi_server.h"

#include "Trios_nssi_types.h"
#include "Trios_nssi_fprint_types.h"
#include "Trios_nnti_fprint_types.h"
#include "Trios_nssi_xdr.h"
#include "Trios_nnti.h"

#include <stdio.h>
#ifdef HAVE_TRIOS_MALLOC_H
#include <malloc.h>
#endif
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <sys/types.h>
#include <signal.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <queue>
#include <map>

#include "Trios_logger.h"
#include "Trios_timer.h"
#include "Trios_signal.h"
#include "Trios_trace.h"
#include "Trios_threads.h"
#include "Trios_nssi_rpc.h"
#include "buffer_queue.h"
#include "Trios_nssi_debug.h"

#include "nssi_opcodes.h"
#include "nssi_trace.h"
#include "nssi_service_args.h"



extern NNTI_transport_t transports[NSSI_RPC_COUNT];
extern nssi_config_t nssi_config;

extern trios_buffer_queue_t send_bq;
extern trios_buffer_queue_t recv_bq;
extern trios_buffer_queue_t rdma_target_bq;
extern trios_buffer_queue_t rdma_get_bq;
extern trios_buffer_queue_t rdma_put_bq;


int rpc_get_service(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const void *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr);

int rpc_kill_service(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const nssi_kill_service_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr);

int rpc_trace_reset(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const nssi_trace_reset_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr);

typedef void (*progress_callback)(bool is_idle);


static bool time_to_die = false;

static nssi_service local_service;

#undef USE_THREADED_SERVERS


static std::map<int, struct nssi_svc_op> supported_ops;
typedef std::map<int, struct nssi_svc_op>::iterator supported_ops_iterator_t;
static nthread_lock_t supported_ops_mutex;


unsigned long max_mem_allowed=0;

static int trace_counter_gid;
static int trace_interval_gid;


#ifdef GNI_PERF
#include <gemini.h>
gemini_state_t gni_state;
#endif

/* Need a struct to encapsulate a IB connection addr/port pair.
 * STL does not like arrays in template defs.
 */
struct caller_reqid {
    NNTI_peer_t caller;
    unsigned long reqid;

    caller_reqid(const NNTI_peer_t *c, const unsigned long r) {
        caller=*c;
        reqid=r;
    }
};
/* Need a comparison operator to pass into the conn_map.
 */
struct caller_reqid_lt
{
    bool operator()(const struct caller_reqid &cr1, const struct caller_reqid &cr2) const
    {
        if (strcmp(cr1.caller.url, cr2.caller.url)<0)   return TRUE;
        if ((strcmp(cr1.caller.url, cr2.caller.url)==0) &&
            (cr1.reqid < cr2.reqid))             return TRUE;

        return FALSE;
    }
};

typedef struct {
    int           opcode;
    unsigned long request_id;
    uint64_t      start_time;
    double        arrival_time;

    NNTI_buffer_t *data_hdl;
    NNTI_buffer_t  shadow_data;
    NNTI_buffer_t *shadow_data_hdl;

    int8_t         is_responseless;
} request_args_t;

static std::map<struct caller_reqid, request_args_t *, caller_reqid_lt> request_args_map;
typedef std::map<struct caller_reqid, request_args_t *, caller_reqid_lt>::iterator request_args_map_iterator_t;
typedef std::pair<struct caller_reqid, request_args_t *> request_args_map_pair_t;

static nthread_lock_t request_args_map_mutex;

typedef struct {
    const NNTI_buffer_t *data_buffer;   // this is the original client data buffer
    const NNTI_buffer_t *shadow_buffer; // this is the shadow buffer on the server
    bool                 data_was_got;
    bool                 data_was_put;
} shadow_buffer_entry;
static std::map<const NNTI_buffer_t *, shadow_buffer_entry *> shadow_buffer_map;
typedef std::map<const NNTI_buffer_t *, shadow_buffer_entry *>::iterator shadow_buffer_map_iterator_t;
static nthread_lock_t shadow_buffer_mutex;

static log_level shadow_debug_level = LOG_UNDEFINED;



//static void print_raw_buf(void *buf, uint32_t size)
//{
//    if (logging_debug(rpc_debug_level)) {
//        FILE* f=logger_get_file();
//        uint64_t print_limit=(size<96) ? size : 96;
//        fprintf(f, "\nbuf (%p)\n", buf);
//        fflush(f);
//        if (buf != NULL) {
//            uint64_t l=0;
//            for (l=0;l<print_limit;l++) {
//                if (l%32 == 0) fprintf(f, "\nbuf (%lu) (offset(%ld)) => ", (uint64_t)buf, l);
//                fprintf(f, "%02hhX", ((char *)buf)[l]);
//            }
//            fprintf(f, "\n");
//        }
//    }
//}



void request_args_add(const NNTI_peer_t *caller, const unsigned long reqid, request_args_t *request_args)
{
    caller_reqid cr(caller, reqid);

    log_debug(rpc_debug_level, "enter - adding caller(%s) reqid(%lu)",
            caller->url, reqid);

    if (nthread_lock(&request_args_map_mutex)) log_warn(rpc_debug_level, "failed to get thread lock");
    request_args_map[cr] = request_args;
    log_debug(rpc_debug_level, "request_args_map.size() == %lu", request_args_map.size());
    nthread_unlock(&request_args_map_mutex);
    log_debug(rpc_debug_level, "end");
}

request_args_t *request_args_get(const NNTI_peer_t *caller, const unsigned long reqid)
{
    caller_reqid cr(caller, reqid);

    log_debug(rpc_debug_level, "enter - looking for caller(%s) reqid(%lu)",
            caller->url, reqid);

    if (nthread_lock(&request_args_map_mutex)) log_warn(rpc_debug_level, "failed to get thread lock");
    request_args_map_iterator_t iter=request_args_map.find(cr);
    request_args_t *request_args=iter->second;
    log_debug(rpc_debug_level, "request_args_map.size() == %lu", request_args_map.size());
    nthread_unlock(&request_args_map_mutex);

    log_debug(rpc_debug_level, "end");

    return(request_args);
}

void request_args_del(const NNTI_peer_t *caller, const unsigned long reqid)
{
    caller_reqid cr(caller, reqid);
    request_args_t *request_args=NULL;

    log_debug(rpc_debug_level, "enter - deleting caller(%s) reqid(%lu)",
            caller->url, reqid);

    if (nthread_lock(&request_args_map_mutex)) log_warn(rpc_debug_level, "failed to get thread lock");
    request_args_map_iterator_t iter=request_args_map.find(cr);
    // it's OK if the caller/reqid is not found.  could have been an error before request args were added.
    if (iter != request_args_map.end()) {
        request_args=iter->second;
        request_args_map.erase(iter);
        log_debug(rpc_debug_level, "request_args_map.size() == %lu", request_args_map.size());
        if (request_args != NULL) {
            free(request_args);
        }
    }
    nthread_unlock(&request_args_map_mutex);

    log_debug(rpc_debug_level, "end");
}


static void print_shadow_buffer_map()
{
    log_level debug_level=shadow_debug_level;

    if (!logging_debug(debug_level)) {
        return;
    }

    shadow_buffer_map_iterator_t i;
    if (nthread_lock(&shadow_buffer_mutex)) log_warn(shadow_debug_level, "failed to get thread lock");
    for (i=shadow_buffer_map.begin(); i != shadow_buffer_map.end(); i++) {
        log_debug(debug_level, "shadow_buffer_map key=%p sbe=%p", i->first, i->second);
    }
    nthread_unlock(&shadow_buffer_mutex);
}
static NNTI_result_t insert_shadow_buffer(shadow_buffer_entry *sbe)
{
    NNTI_result_t  rc=NNTI_OK;

    if (nthread_lock(&shadow_buffer_mutex)) log_warn(shadow_debug_level, "failed to get thread lock");

    log_debug(shadow_debug_level, "adding shadow buffer (sbe=%p)", sbe);

    assert(shadow_buffer_map.find(sbe->shadow_buffer) == shadow_buffer_map.end());
    shadow_buffer_map[sbe->shadow_buffer] = sbe;

    log_debug(shadow_debug_level, "added shadow buffer (sbe=%p)", sbe);

    nthread_unlock(&shadow_buffer_mutex);

    return(rc);

}
static shadow_buffer_entry *get_shadow_buffer(const NNTI_buffer_t *sbuf)
{
    shadow_buffer_entry *sbe=NULL;

    if (nthread_lock(&shadow_buffer_mutex)) log_warn(shadow_debug_level, "failed to get thread lock");

    log_debug(shadow_debug_level, "looking for shadow buffer (sbuf=%p)", sbuf);

    if (shadow_buffer_map.find(sbuf) != shadow_buffer_map.end()) {
        sbe = shadow_buffer_map[sbuf];
    }

    log_debug(shadow_debug_level, "found shadow buffer entry (sbe=%p)", sbe);

    nthread_unlock(&shadow_buffer_mutex);

    if (sbe != NULL) {
        log_debug(shadow_debug_level, "shadow buffer entry found (sbe=%p)", sbe);
        return sbe;
    }

    log_debug(shadow_debug_level, "shadow buffer entry NOT found");

    print_shadow_buffer_map();

    return(NULL);
}
static shadow_buffer_entry *del_shadow_buffer(NNTI_buffer_t *victim)
{
    shadow_buffer_entry *sbe=NULL;

    log_level debug_level = shadow_debug_level;

    if (nthread_lock(&shadow_buffer_mutex)) log_warn(shadow_debug_level, "failed to get thread lock");

    log_debug(debug_level, "deleting shadow buffer (victim=%p)", victim);

    if (shadow_buffer_map.find(victim) != shadow_buffer_map.end()) {
        sbe = shadow_buffer_map[victim];
    }

    if (sbe != NULL) {
        log_debug(shadow_debug_level, "shadow buffer entry found and deleted (sbe=%p ; victim=%p)", sbe, victim);
        shadow_buffer_map.erase(victim);
    } else {
        log_debug(debug_level, "shadow buffer entry NOT found");
    }

    log_debug(debug_level, "deleted (sbe=%p ; victim=%p)", sbe, victim);

    nthread_unlock(&shadow_buffer_mutex);

    print_shadow_buffer_map();

    return(sbe);
}



/* ----------- Implementation of core services ----------- */

/**
  * @brief Return the service description of this service.
  *
  * @param caller @input_type the client's PID
  * @param args @input_type arguments needed to verify the container
  * @param data_addr @input_type address at which the bulk data can be found
  * @param res_addr @input_type address of remote result buffer
  *
  * @returns The service descriptor for this service.
 */
int rpc_get_service(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const void *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NSSI_OK;

    /* copy the service description into the result */
    log_debug(rpc_debug_level, "entered get service");

    /* send the local_service descriptor to the client */
    rc = nssi_send_result(caller, request_id, rc, &local_service, res_addr);
    if (rc != NSSI_OK) {
        log_warn(rpc_debug_level, "Could not send service description to client");
    }

    return rc;
}


/**
  * @brief Schedule this service to be killed.
  *
  * @param caller @input_type the client's PID
  * @param args @input_type arguments needed to verify the container
  * @param data_addr @input_type address at which the bulk data can be found
  * @param res_addr @input_type address of remote result buffer
  *
  * @returns The service descriptor for this service.
 */
int rpc_kill_service(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const nssi_kill_service_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NSSI_OK;

    /* copy the service description into the result */
    log_debug(rpc_debug_level, "killing service");

    /* set the "exit_now" flag. */
    /*nssi_abort();*/

    switch (args->sig) {

    case 0:
        /* Graceful exit. When threads come out of timeout, they will exit */
        time_to_die = true;
        log_debug(rpc_debug_level, "Graceful Abort");
        //trios_abort();  /* sets the exit_now flag */
        break;

    default:
        /* Force exit. When things just must die! */
        log_debug(LOG_ALL, "Forced Abort");
        trios_abort();
        exit(args->sig);
        break;

    }

    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);
    if (rc != NSSI_OK) {
        log_warn(rpc_debug_level, "Unable to send result from kill");
    }

    return rc;
}


/**
  * @brief Reset the tracing library.
  *
  * If the service is currently using the tracing
  * library, this operation will reset the tracing
  * library... forcing a flush of the previous file
  * and creating new file for the trace data. All
  * counts and timers are reset to 0.
  *
  * @param caller @input_type the client's PID
  * @param args @input_type arguments needed to verify the container
  * @param data_addr @input_type address at which the bulk data can be found
  * @param result @output_type no result
  *
  * @returns The service descriptor for this service.
 */
int rpc_trace_reset(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const nssi_trace_reset_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NSSI_OK;

    char *fname;
    const int ftype = args->ftype;
    char *enable;

    if ((!args->fname) || (strcmp(args->fname,"")==0)) {
        fname = NULL;
    }
    else {
        fname = (char *)args->fname;
    }

    if ((!args->enable) || (strcmp(args->enable,"")==0)) {
        enable = NULL;
    }
    else {
        enable = (char *)args->enable;
    }


    /* copy the service description into the result */
    log_debug(rpc_debug_level, "reset tracing(%s, %d)",
            fname, ftype);

    /* set the "exit_now" flag. */
    /*trace_reset(enable_flag, fname, ftype);*/
    //trace_fini();
    //trace_init(fname, ftype);

    trace_reset(fname, ftype, enable);

    /* send result back to client */
    nssi_send_result(caller, request_id, rc, NULL, res_addr);
    return rc;
}


/* ----------- Implementation of the NSSI messaging -------- */

/**
 * @brief  Fetch or extract the operation arguments.
 *
 * If the args were small enough to fit into the short
 * request buffer, we extract the args directly from
 * the buffer.  Otherwise, we GET the arguments from
 * a remote memory descriptor on the client.
 *
 * @param encoded_buf  The encoded short request buffer.
 * @param header  The request header.
 * @param decode_args  The XDR function that decodes the arguments.
 * @param args  Where to place the decoded arguments.
 */
static int fetch_args(
        NNTI_peer_t         *caller,
        nssi_request_header *header,
        xdrproc_t            xdr_decode_args,
        void                *args)
{
    int rc;   /* return code from non-NSSI methods */
    XDR xdrs;

    trios_declare_timer(call_time);

    /* pointer to the decoded buffer for arguments */
    char *buf=NULL;
    NNTI_buffer_t       encoded_args_hdl;
    NNTI_work_request_t encoded_args_wr;
    NNTI_status_t       status;
    nssi_size encoded_args_size = NNTI_BUFFER_SIZE(&header->args_addr);

    /* allocate the decoded buffer */
    rc=NNTI_alloc(
            &transports[caller->peer.transport_id],
            encoded_args_size,
            1,
            NNTI_GET_DST,
            &encoded_args_hdl);
    if (rc != NNTI_OK) {
        log_error(rpc_debug_level, "failed registering long args: %s",
                nnti_err_str(rc));
        goto cleanup;
    }
    buf=NNTI_BUFFER_C_POINTER(&encoded_args_hdl);

    assert(header->fetch_args);

    log_debug(rpc_debug_level,
            "get args from client");

    /* fetch the buffer from the client */
    trios_start_timer(call_time);
    rc=NNTI_get(
            &header->args_addr,
            0,
            encoded_args_size,
            &encoded_args_hdl,
            0,
            &encoded_args_wr);
    trios_stop_timer("NNTI_get - long args", call_time);
    if (rc != NNTI_OK) {
        log_fatal(rpc_debug_level,
                "could not get long args from client");
        goto cleanup;
    }
    trios_start_timer(call_time);
    rc=NNTI_wait(
            &encoded_args_wr,
            -1,
            &status);
    trios_stop_timer("NNTI_wait - long args", call_time);
    if (rc != NNTI_OK) {
        log_error(rpc_debug_level, "failed waiting for long args: %s",
                nnti_err_str(rc));
        goto cleanup;
    }

    /* decode the arguments */
    log_debug(rpc_debug_level,"decoding args, size=%d",
            encoded_args_size);

    /* create an xdr memory stream for the decoded args */
    xdrmem_create(
            &xdrs,
            NNTI_BUFFER_C_POINTER(&encoded_args_hdl),
            NNTI_BUFFER_SIZE(&encoded_args_hdl),
            XDR_DECODE);
    /* decode -- will allocate memory if necessary */
    trios_start_timer(call_time);
    if (!xdr_decode_args(&xdrs, args)) {
        log_fatal(rpc_debug_level,"could not decode args");
        fprint_NNTI_status(logger_get_file(), "status", "FATAL", &status);
//        print_raw_buf(NNTI_BUFFER_C_POINTER(&encoded_args_hdl),
//                NNTI_BUFFER_SIZE(&encoded_args_hdl));
//        fflush(logger_get_file());
        rc = NSSI_EDECODE;
        goto cleanup;
    }
    trios_stop_timer("xdr_decode_args - decode", call_time);

cleanup:
    /* if we had to fetch the args, we need to free the buffer */
    if (buf) {
        int cleanup_rc;
        cleanup_rc=NNTI_free(&encoded_args_hdl);
        if (cleanup_rc != NNTI_OK) {
            log_error(rpc_debug_level, "failed unregistering long args: %s",
                    nnti_err_str(cleanup_rc));
        }
    }

    return rc;
}

/**
 * @brief Send the result back to the client.
 *
 * If the result is small enough to fit inside a small result
 * buffer, the results are sent in one message transfer. If the
 * results are too large, we tell the client to fetch the result
 * (by setting the fetch_result flag of the result header to true)
 * send the result header, then wait for the client to fetch the
 * result.
 *
 * @param dest              @input Where to send the encoded result.
 * @param xdr_encode_result @input function used to encode result.
 * @param return_code       @input The return code of the function.
 * @param result            @input the result of the function.
 */
static int send_result(const NNTI_peer_t   *caller,
                       const unsigned long  request_id,
                       const NNTI_buffer_t *dest_addr,
                       xdrproc_t            xdr_encode_result,
                       const int            return_code,
                       void                *result)
{
    static uint32_t res_counter = 1;

    trios_declare_timer(call_time);

    int rc; /* return code for non-NSSI methods */

    uint32_t hdr_size;
    uint32_t res_size;
    uint32_t res_buf_size;
    uint32_t remaining;
    uint32_t valid_bytes;
    char *buf=NULL;
    NNTI_buffer_t       short_res;
    NNTI_buffer_t      *short_res_hdl=&short_res;
    NNTI_work_request_t short_res_wr;
    NNTI_buffer_t       long_res_hdl;
    NNTI_status_t       wait_status;
    nssi_result_header  header;

    NNTI_buffer_t       long_res_ack_hdl;
    NNTI_work_request_t long_res_ack_wr;
    NNTI_status_t       long_res_ack_status;


    request_args_t *args=request_args_get(caller, request_id);
    int opcode=args->opcode;

    /* xdrs for the header and the result. */
    XDR hdr_xdrs, res_xdrs;

    /* initialize the result header */
    memset(&header, 0, sizeof(nssi_result_header));

    /* --- HANDLE ERROR CASE --- */

    /* If a function returns an error, we return the error code, but
     * no result.  */
    if (return_code != NSSI_OK) {

        /* treat as if there is no result to return */
        xdr_encode_result = (xdrproc_t)&xdr_void;
    }

    /* --- CALCULATE SIZES --- */

    /* Calculate size of the encoded header */
    hdr_size = xdr_sizeof((xdrproc_t)&xdr_nssi_result_header, &header);

    /* Calculate size of the encoded result */
    if (result == NULL) {
        res_size = 0;
    } else {
        res_size = xdr_sizeof(xdr_encode_result, result);
    }

    /* Extract the size of the client-side buffer for the result */
    res_buf_size = NNTI_BUFFER_SIZE(dest_addr);

    /* Calculate space left in the short result buffer */
    remaining = res_buf_size - hdr_size;

    /* allocated an xdr memory stream for the short result buffer */
    if (res_buf_size <= 0) {
        log_error(rpc_debug_level, "********** res_buf_size (%lu) <= 0", res_buf_size);
        fprint_NNTI_buffer(logger_get_file(), "dest_addr", "ERROR res_buf_size<=0 %", dest_addr);
    }
    assert(res_buf_size > 0);

    if (nssi_config.use_buffer_queue) {
        short_res_hdl=trios_buffer_queue_pop(&send_bq);
        assert(short_res_hdl);
    } else {
        rc=NNTI_alloc(
                &transports[caller->peer.transport_id],
                res_buf_size,
                1,
                NNTI_SEND_SRC,
                short_res_hdl);
        if (rc != NNTI_OK) {
            log_error(rpc_debug_level, "failed registering short result: %s",
                    nnti_err_str(rc));
        }
        buf=NNTI_BUFFER_C_POINTER(short_res_hdl);
        memset(buf, 0, res_buf_size);  // address valgrind uninitialized error
    }

    xdrmem_create(
            &hdr_xdrs,
            NNTI_BUFFER_C_POINTER(short_res_hdl),
            NNTI_BUFFER_SIZE(short_res_hdl),
            XDR_ENCODE);

    /* If the result fits in the short result buffer, send it with the header */
    if (res_size < remaining) {

        log_debug(rpc_debug_level, "sending short_result %lu, "
            "result buffer size = %d, header size = %d, available space = %d, result_size = %d",
            request_id, res_buf_size, hdr_size, remaining, res_size);

        /* client needs this information from the header */
        header.fetch_result = FALSE;
        header.id           = request_id;
        header.opcode       = opcode;
        header.result_size  = res_size;
        header.rc           = return_code;

        /* encode the header  */
        log_debug(rpc_debug_level, "encode result header");
        trios_start_timer(call_time);
        if (!xdr_nssi_result_header(&hdr_xdrs, &header)) {
            log_fatal(rpc_debug_level, "failed to encode the result header");
            return NSSI_EENCODE;
        }
        trios_stop_timer("xdr_nssi_result_header - encode", call_time);

        /* encode the result in the header */
        log_debug(rpc_debug_level, "encode result data");
        if (result != NULL) {
            trios_start_timer(call_time);
            if (!xdr_encode_result(&hdr_xdrs, result)) {
                log_fatal(rpc_debug_level, "failed to encode the result");
                return NSSI_EENCODE;
            }
            trios_stop_timer("xdr_encode_result - encode", call_time);
        }
    }

    /* if result does not fit, client has to fetch result */
    else {

        res_counter++;

        log_debug(rpc_debug_level, "sending long result %lu, "
            "available space = %d, result_size = %d", request_id, remaining,
                res_size);

        /* allocate memory for the result
         * structure keeps track of the buffer so it can free
         * the memory later. */
        rc=NNTI_alloc(
                &transports[caller->peer.transport_id],
                res_size,
                1,
                NNTI_GET_SRC,
                &long_res_hdl);
        if (rc != NNTI_OK) {
            log_error(rpc_debug_level, "failed registering long result: %s",
                    nnti_err_str(rc));
        }
        buf=NNTI_BUFFER_C_POINTER(&long_res_hdl);

        header.result_addr=long_res_hdl;

        log_debug(rpc_debug_level, "allocated long_res_buf(%lu) req_id(%lu)", buf, request_id);

        trios_start_timer(call_time);
        rc=NNTI_alloc(
                &transports[caller->peer.transport_id],
                sizeof(int8_t),
                1,
                NNTI_RECV_DST,
                &long_res_ack_hdl);
        trios_stop_timer("NNTI_register_memory - long result ack", call_time);
        if (rc != NNTI_OK) {
            log_error(rpc_debug_level, "failed registering long result ack: %s",
                    nnti_err_str(rc));
        }

        header.result_ack_addr=long_res_ack_hdl;

        /* we want the client to fetch the result */
        /* client needs this information from the header */
        header.fetch_result = TRUE;
        header.id           = request_id;
        header.opcode       = opcode;
        header.result_size  = res_size;
        header.rc           = return_code;

        /* create a xdr memory stream for the encoded args buffer */
        xdrmem_create(
                &res_xdrs,
                NNTI_BUFFER_C_POINTER(&long_res_hdl),
                NNTI_BUFFER_SIZE(&long_res_hdl),
                XDR_ENCODE);

        /* encode the header  */
        log_debug(rpc_debug_level, "encode result %lu header",
                request_id);
        trios_start_timer(call_time);
        if (!xdr_nssi_result_header(&hdr_xdrs, &header)) {
            log_fatal(rpc_debug_level, "failed to encode the result header");
            rc = NSSI_EENCODE;
            goto cleanup;
        }
        trios_stop_timer("xdr_nssi_result_header - encode", call_time);

        /* encode the result  */
        log_debug(rpc_debug_level, "encode result %lu data", request_id);
        trios_start_timer(call_time);
        if (!xdr_encode_result(&res_xdrs, result)) {
            log_fatal(rpc_debug_level, "failed to encode the result");
            rc = NSSI_EENCODE;
            goto cleanup;
        }
        trios_stop_timer("xdr_encode_result - encode", call_time);

//        if (logging_debug(rpc_debug_level)) {
//            u_int64_t print_limit=(long_res_buf->msg_size<90) ? long_res_buf->msg_size : 90;
//            for (int l=0;l<print_limit;l++) {
//                if (l%30 == 0) fprintf(stdout, "\nlong_res_buf (%lu) after encode (req_id(%lu) offset(%d)) => ", (uint64_t)long_res_buf->msg, request_id, l);
//                fprintf(stdout, "%02hhX", ((char *)long_res_buf->msg)[l]);
//            }
//            fprintf(stdout, "\n");
//        }
    }

    if (logging_debug(rpc_debug_level)) {
        fprint_nssi_result_header(logger_get_file(), "header", "nssi_result_header", &header);
    }

    /* send the short result to the client */
    valid_bytes = hdr_size;

    if (!header.fetch_result)
        valid_bytes += res_size;

    if (logging_debug(rpc_debug_level)) {
        log_debug(rpc_debug_level, "send short result %lu "
                "(in xdr bytes:  len=%d bytes: encoded_header=%d bytes, res=%d bytes)",
                request_id, valid_bytes, hdr_size, res_size);

        fprint_NNTI_buffer(logger_get_file(), "nssi_server: 768 -- dest_addr", "%", dest_addr);

        fprint_NNTI_buffer(logger_get_file(), "nssi_server: 768 -- short_res_hdl", "%", short_res_hdl);

        fprint_NNTI_peer(logger_get_file(), "nssi_server: 772 -- caller", "%", caller);
    }


    /* TODO: Handle the timeout case.  This probably means the client died */
    trios_start_timer(call_time);
    rc=NNTI_send(
            caller,
            short_res_hdl,
            dest_addr,
            &short_res_wr);
    trios_stop_timer("NNTI_send - short result", call_time);
    if (rc != NNTI_OK) {
        log_error(rpc_debug_level, "failed sending short result: %s",
                nnti_err_str(rc));
    }
    trios_start_timer(call_time);
    rc=NNTI_wait(
            &short_res_wr,
            -1,
            &wait_status);
    trios_stop_timer("NNTI_wait - short result", call_time);
    if (rc != NNTI_OK) {
        log_error(rpc_debug_level, "failed waiting for short result: %s",
                nnti_err_str(rc));
        goto cleanup;
    }


    /* if the client has to fetch the results, we need to wait for
     * the GET to complete */
    if (header.fetch_result) {
        log_debug(rpc_debug_level, "waiting for client to "
            "ACK request %lu", request_id);

        trios_start_timer(call_time);
        NNTI_create_work_request(
                &long_res_ack_hdl,
                &long_res_ack_wr);
        rc=NNTI_wait(
                &long_res_ack_wr,
                -1,
                &long_res_ack_status);
        NNTI_destroy_work_request(
        		&long_res_ack_wr);
        trios_stop_timer("NNTI_wait - long result ack", call_time);
        if (rc != NNTI_OK) {
            log_error(rpc_debug_level, "failed waiting for client to send long result ack: %s",
                    nnti_err_str(rc));
        }
    }

cleanup:
    if (header.fetch_result) {
        rc=NNTI_free(&long_res_hdl);
        if (rc != NNTI_OK) {
            log_error(rpc_debug_level, "failed unregistering long result: %s",
                    nnti_err_str(rc));
        }

        rc=NNTI_free(&long_res_ack_hdl);
        if (rc != NNTI_OK) {
            log_error(rpc_debug_level, "failed unregistering long result ack: %s",
                    nnti_err_str(rc));
        }
    }

    if (nssi_config.use_buffer_queue) {
        trios_buffer_queue_push(&send_bq, short_res_hdl);
    } else {
        rc=NNTI_free(short_res_hdl);
        if (rc != NNTI_OK) {
            log_error(rpc_debug_level, "failed unregistering short result: %s",
                    nnti_err_str(rc));
        }
    }

    log_debug(rpc_debug_level, "result %lu sent", request_id);

    return rc;
}

/**
 * Lookup an opcode in the list of supported ops.
 */
static int lookup_service_op(
        const int opcode,
        nssi_svc_op *result_op)
{
    int rc = NSSI_OK;

    log_debug(rpc_debug_level, "enter (opcode=%d)", opcode);

    log_debug(rpc_debug_level, "locking ops mutex");
    if (nthread_lock(&supported_ops_mutex)) log_warn(rpc_debug_level, "failed to get thread lock");
    log_debug(rpc_debug_level, "locked ops mutex");
    if (supported_ops.find(opcode) == supported_ops.end()) {
        rc = NSSI_ENOENT;
    } else {
        /* should use the copy constructor */
        *result_op = supported_ops[(int)opcode];
    }
    log_debug(rpc_debug_level, "unlocking ops mutex");
    nthread_unlock(&supported_ops_mutex);
    log_debug(rpc_debug_level, "unlocked ops mutex");

    log_debug(rpc_debug_level, "exit (*result_op=%p)", *result_op);

    return rc;
}


/**
 * @brief Send result back to the client.
 *
 * This function allows the service library to send a result to
 * a remote memory descriptor.  The result will be packaged as
 * an \ref nssi_result, then PUT on the memory descriptor of
 * the remote address.
 *
 * @param opcode  The opcode ID of the server function.
 * @param return_code The return code for the server function.
 * @param result  Pointer to the result data structure (NULL if not used)
 * @param result_addr The remote memory address of the client (where
 * to PUT the result)
 *
 * @returns \ref NSSI_OK if successful.
 */
int nssi_send_result(
        const NNTI_peer_t   *caller,
        const unsigned long  request_id,
        const int            return_code,
        void                *result,
        const NNTI_buffer_t *result_addr)
{
    int rc = NSSI_OK;
    nssi_svc_op op;

    log_debug(rpc_debug_level, "enter");

    request_args_t *args=request_args_get(caller, request_id);

    log_debug(rpc_debug_level, "args=%p", args);

    if (args->is_responseless==FALSE) {
        /* lookup the service description of the opcode */
        rc = lookup_service_op(args->opcode, &op);
        if (rc != NSSI_OK) {
            log_warn(rpc_debug_level, "Invalid opcode=%d", args->opcode);
            return rc;
        }

        rc = send_result(caller, request_id, result_addr, op.encode_res, return_code, result);
        if (rc != NSSI_OK) {
            log_warn(rpc_debug_level, "Unable to send result to client: %s", nssi_err_str(rc));
        }
    }

    return rc;
}


/**
 * @brief Process an rpc service request.
 *
 * Each incoming rpc service request arrives as a chunk of xdr data. This
 * method first decodes the header and arguments from the xdr data, then
 * it calls the appropriate method to process the request. We define the
 * parameters as void pointers because this function may execute as a thread_pool_task.
 *
 *
 * @param args The arguments (rpc_request *).
 */

int nssi_process_rpc_request(nssi_svc_rpc_request *rpc_req)
{
    XDR xdrs;
    int rc;
    nssi_svc_op svc_op;   /* current operation */

    trios_declare_timer(call_time);

    /* space for args and result (these are passed in with the header) */
    void *op_args = NULL;
    void *op_res  = NULL;

    nssi_request_header header;
    int req_count = 0;
    log_level debug_level = rpc_debug_level;

    shadow_buffer_entry *sbe=NULL;
    char          *shadow_data_buf =NULL;
    nssi_size      shadow_data_size=0;

    NNTI_buffer_t *res_addr=NULL;


    request_args_t *req_args=NULL;

    NNTI_peer_t caller      = rpc_req->caller;
    char *req_buf           = rpc_req->req_buf;
    nssi_size short_req_len = rpc_req->short_req_len;
    req_count               = rpc_req->id;

    log_debug(debug_level, "req_buf=%p", req_buf);

    log_debug(debug_level, "Started processing request %d",
            req_count);

    /* memory check - log memory statistics.  if memory in use
     * is greater than maximum memory allowed, then exit.
     */
//    rc = nthread_lock(&meminfo_mutex);
//    log_meminfo(rpc_debug_level);
//    unsigned long main_memory_in_use = main_memory_used();
//    rc = nthread_unlock(&meminfo_mutex);
//
//    log_debug(debug_level,
//            "max memory check (allowed=%lukB, in_use=%lukB, %%free=%f)",
//            max_mem_allowed, main_memory_in_use,
//            100.0-(((float)main_memory_in_use/(float)max_mem_allowed)*100.0));
//    if ((max_mem_allowed > 0) &&
//            (main_memory_in_use > max_mem_allowed)) {
//        rc = NSSI_OK;
//        log_error(rpc_debug_level,
//                "max memory allowed exceeded (allowed=%lukB, in_use=%lukB), exiting",
//                max_mem_allowed, main_memory_in_use);
//        goto cleanup;
//    }

    /* initialize the request header */
    memset(&header, 0, sizeof(nssi_request_header));

    /* create an xdr memory stream from the request buffer */
    xdrmem_create(
            &xdrs,
            req_buf,
            short_req_len,
            XDR_DECODE);

    /* decode the request header */
    log_debug(debug_level, "decoding header for request %d...",
            req_count);
    trios_start_timer(call_time);
    rc = xdr_nssi_request_header(&xdrs, &header);  // this allocates an xdr_string, must use xdr_free
    trios_stop_timer("xdr_nssi_request_header - decode", call_time);
    if (!rc) {
        log_fatal(debug_level, "failed to decode header");
        rc = NSSI_EDECODE;
        goto cleanup;
    }

    log_debug(debug_level, "begin processing request %d with opcode (%u)",
            req_count, header.opcode);

    if (logging_debug(rpc_debug_level)) {
        fprint_nssi_request_header(logger_get_file(), "header", "nssi_request_header", &header);
    }

    /* See if the opcode is in our list of supported ops */
    rc = lookup_service_op(header.opcode, &svc_op);
    if (rc != NSSI_OK) {
        /* if we get here, there is not match */
        log_warn(debug_level, "unrecognized request: opcode=%d",
                header.opcode);
        rc = NSSI_EBADRPC;
        goto cleanup;
    }

    log_debug(LOG_OFF, "header.id=%d", header.id);

    log_info(debug_level, "Found op for opcode=%d", header.opcode);


    op_args = calloc(1, svc_op.sizeof_args);

    op_res = calloc(1, svc_op.sizeof_res);

    /* start interval for decode args */
    trace_start_interval(trace_interval_gid, 0);

    /* initialize args and res */
    log_debug(debug_level, "NNTI_BUFFER_SIZE(&header.res_addr)==%d for request %d",
            NNTI_BUFFER_SIZE(&header.res_addr), req_count);

    /* If the args fit in the header, extract them from the
     * header buffer.  Otherwise, get them from the client.
     */
    if (!header.fetch_args) {
        if (! svc_op.decode_args(&xdrs, op_args)) {
            log_fatal(debug_level,"could not decode args");
            rc = NSSI_EDECODE;
            goto cleanup;
        }
    }
    else {
        /* fetch the operation arguments */
        log_debug(debug_level, "fetching args for request %d",
                req_count);
        // is reentrant??
        trios_start_timer(call_time);
        rc = fetch_args(
                &caller,
                &header,
                svc_op.decode_args,
                op_args);
        trios_stop_timer("fetch_args", call_time);
        if (rc != NSSI_OK) {
            log_fatal(debug_level,
                    "unable to fetch args");
            goto cleanup;
        }
    }

    req_args = (request_args_t *)malloc(sizeof(request_args_t));
    req_args->opcode = header.opcode;
    req_args->request_id = header.id;
    req_args->arrival_time = rpc_req->arrival_time;
    req_args->start_time = trios_get_time_ms();
    req_args->data_hdl=&header.data_addr;
    req_args->is_responseless=header.is_responseless;
    request_args_add(&caller, header.id, req_args);

    shadow_data_size=NNTI_BUFFER_SIZE(&header.data_addr);
    if (header.fetch_data == TRUE) {

        req_args->shadow_data    =header.data_addr;
        req_args->shadow_data_hdl=&req_args->shadow_data;

    } else if ((shadow_data_size > 0) && (header.fetch_data == FALSE)) {

        if ((nssi_config.use_buffer_queue) &&
            (nssi_config.rdma_buffer_queue_buffer_size >= shadow_data_size)) {
            log_debug(rpc_debug_level, "using buffer queue for SHADOW buffer");
            req_args->shadow_data_hdl=trios_buffer_queue_pop(&rdma_target_bq);
            assert(req_args->shadow_data_hdl);
            NNTI_BUFFER_SIZE(req_args->shadow_data_hdl)=shadow_data_size;
        } else {
            log_debug(rpc_debug_level, "allocating buffer for SHADOW buffer");

            req_args->shadow_data_hdl=&req_args->shadow_data;

            trios_start_timer(call_time);
            rc=NNTI_alloc(
                    &transports[rpc_req->svc->transport_id],
                    shadow_data_size,
                    1,
                    (NNTI_buf_ops_t)(NNTI_GET_SRC|NNTI_PUT_DST),
                    req_args->shadow_data_hdl);
            trios_stop_timer("NNTI_register_memory - shadow_data_buf", call_time);
            if (rc != NNTI_OK) {
                log_error(rpc_debug_level, "failed registering shadow_data: %s",
                        nnti_err_str(rc));
                goto cleanup;
            }
        }

        nssi_size data_offset=xdr_getpos(&xdrs);

        log_debug(rpc_debug_level,"extracting data (size=%d, offset=%d) from short request", shadow_data_size, data_offset);
        /* copy small data into the short request */
        memcpy(NNTI_BUFFER_C_POINTER(req_args->shadow_data_hdl), req_buf+data_offset, shadow_data_size);

        sbe=(shadow_buffer_entry *)malloc(sizeof(shadow_buffer_entry));
        sbe->data_buffer  =req_args->data_hdl;
        sbe->shadow_buffer=req_args->shadow_data_hdl;
        sbe->data_was_got =false;
        sbe->data_was_put =true;
        insert_shadow_buffer(sbe);
    }

    if (header.is_responseless == FALSE) {
        res_addr=&header.res_addr;
    }

    /* end the decode args interval */
    trace_end_interval(trace_interval_gid, TRACE_RPC_DECODE,
            0, "decode request");

    /*
     ** Process the request (print warning if method fails), but
     ** don't return error, because some operations are meant to fail
     */
    log_debug(debug_level, "calling the server function"
            " for request %d (id=%lu, opcode=%d, func=%p, obj=%p)",
            req_count, header.id, header.opcode,
            svc_op.func, svc_op.obj);

    // is reentrant??
    trios_start_timer(call_time);
    if (svc_op.func) {
        rc = svc_op.func(header.id, &caller, op_args, req_args->shadow_data_hdl, res_addr);
    } else if (svc_op.obj) {
        rc = svc_op.obj->doRPC(svc_op.opcode, header.id, &caller, op_args, req_args->shadow_data_hdl, res_addr);
    } else {
        rc = NSSI_ENOENT;
    }
    trios_stop_timer("svc_op", call_time);
    if (rc != NSSI_OK) {
        log_info(rpc_debug_level,
                "user op failed: %s",
                nssi_err_str(rc));
    }

    if ((shadow_data_size > 0) && (header.fetch_data == FALSE)) {

        sbe=del_shadow_buffer(req_args->shadow_data_hdl);
        free(sbe);

        if ((nssi_config.use_buffer_queue) &&
            (nssi_config.rdma_buffer_queue_buffer_size >= (uint32_t)shadow_data_size)) {
            trios_buffer_queue_push(&rdma_target_bq, req_args->shadow_data_hdl);
        } else {
            trios_start_timer(call_time);
            rc=NNTI_free(req_args->shadow_data_hdl);
            trios_stop_timer("NNTI_unregister_memory - shadow_data_buf", call_time);
            if (rc != NNTI_OK) {
                log_error(rpc_debug_level, "failed unregistering data: %s",
                        nnti_err_str(rc));
            }
        }
    }

    /* send result back to client */
    log_debug(debug_level, "sending result for request %d "
            "(%lu) back to client", req_count, header.id);

    if (rc != NSSI_OK) {
        log_fatal(debug_level, "unable to send result %lu"
                " for opcode=%d", header.id, header.opcode);
        goto cleanup;
    }

    /* free data structures created for the args and result */
    log_debug(debug_level, "xdr_freeing args for request %d",
            req_count);
    trios_start_timer(call_time);
    xdr_free((xdrproc_t)svc_op.decode_args, (char *)op_args);
    trios_stop_timer("xdr_free - args", call_time);
    log_debug(debug_level, "xdr_freeing result for request %d",
            req_count);


    log_debug(debug_level, "freeing args for request %d",
            req_count);
    free(op_args);
    log_debug(debug_level, "freeing result for request %d",
            req_count);


    /* This can be made the responsibility of the service */
    trios_start_timer(call_time);
    xdr_free((xdrproc_t)svc_op.encode_res, (char *)op_res);
    trios_stop_timer("xdr_free - result", call_time);
    free(op_res);

    log_debug(debug_level, "result freed for request %d",
            req_count);


cleanup:

    // release space allocated by xdr calls:  leak found by valgrind
    xdr_free((xdrproc_t)xdr_nssi_request_header, (char *)&header);

    log_debug(debug_level, "finished processing request %d", req_count);

    request_args_del(&caller, header.id);

    // release the data allocated for the rpc_req buffer (done in server_start)
    delete [] rpc_req->req_buf;

    // delete the actual request
    delete rpc_req;

    return rc;
}


/**
 * @brief An abstract method to get data from a remote memory descriptor.
 *
 * The server stub uses this function to get or put data to a
 * client memory descriptor.
 *
 * @param buf    @input   the buffer for the data.
 * @param len    @input   the maximum length of the buffer.
 * @param src_md @input   the remote memory descriptor.
 */
int nssi_get_data(
        const NNTI_peer_t *caller,
        void *buf,
        const int len,
        const NNTI_buffer_t *data_addr)
{
    int rc = NSSI_OK;
    NNTI_buffer_t       rpc_msg;
    NNTI_buffer_t      *rpc_msg_hdl=NULL;
    NNTI_work_request_t rpc_msg_wr;
    NNTI_status_t       status;
    trios_declare_timer(call_time);

    shadow_buffer_entry *sbe=NULL;

    if (len == 0)
        return rc;

    if (len < 0)
        return NSSI_EINVAL;

    sbe=get_shadow_buffer(data_addr);
    if (sbe != NULL) {
        memcpy(buf, NNTI_BUFFER_C_POINTER(sbe->shadow_buffer), len);
        sbe->data_was_got=true;
        return rc;
    }

    if ((nssi_config.use_buffer_queue) &&
        (nssi_config.rdma_buffer_queue_buffer_size >= (uint32_t)len)) {
        log_debug(rpc_debug_level, "using buffer queue for GET buffer");
        rpc_msg_hdl=trios_buffer_queue_pop(&rdma_get_bq);
        assert(rpc_msg_hdl);
        NNTI_BUFFER_SIZE(rpc_msg_hdl)=len;
    } else {
        log_debug(rpc_debug_level, "using user buffer for GET buffer");
        trios_start_timer(call_time);
        rpc_msg_hdl=&rpc_msg;
        rc=NNTI_register_memory(
                &transports[data_addr->transport_id],
                (char *)buf,
                len,
                1,
                NNTI_GET_DST,
                rpc_msg_hdl);
        trios_stop_timer("register get dest", call_time);
        if (rc != NNTI_OK) {
            log_error(rpc_debug_level, "failed registering data: %s",
                    nnti_err_str(rc));
        }
    }

    trios_start_timer(call_time);
#ifdef GNI_PERF
    gemini_read_counters(MPI_COMM_WORLD, &gni_state);
#endif
    rc=NNTI_get(
            data_addr,
            0,
            len,
            rpc_msg_hdl,
            0,
            &rpc_msg_wr);
#ifdef GNI_PERF
    gemini_read_counters(MPI_COMM_WORLD, &gni_state);
    gemini_print_counters(MPI_COMM_WORLD, &gni_state, "nssi_get_data - NNTI_get");
#endif
    trios_stop_timer("get to get dest", call_time);
    if (rc != NNTI_OK) {
        log_error(rpc_debug_level, "failed getting data: %s",
                nnti_err_str(rc));
    }
    trios_start_timer(call_time);
    rc=NNTI_wait(
            &rpc_msg_wr,
            -1,
            &status);
#ifdef GNI_PERF
    gemini_read_counters(MPI_COMM_WORLD, &gni_state);
    gemini_print_counters(MPI_COMM_WORLD, &gni_state, "nssi_get_data - NNTI_wait");
#endif
    trios_stop_timer("wait for get dest", call_time);
    if (rc != NNTI_OK) {
        log_error(rpc_debug_level, "failed waiting for data: %s",
                nnti_err_str(rc));
    }
    if ((nssi_config.use_buffer_queue) &&
        (nssi_config.rdma_buffer_queue_buffer_size >= (uint32_t)len)) {
        /* copy the RDMA buffer contents into the user buffer */
        trios_start_timer(call_time);
        memcpy(buf, NNTI_BUFFER_C_POINTER(rpc_msg_hdl), len);
        trios_buffer_queue_push(&rdma_get_bq, rpc_msg_hdl);
        trios_stop_timer("memcpy bq to get dest", call_time);
    } else {
        trios_start_timer(call_time);
        rc=NNTI_unregister_memory(rpc_msg_hdl);
        trios_stop_timer("deregister get dest", call_time);
        if (rc != NNTI_OK) {
            log_error(rpc_debug_level, "failed unregistering data: %s",
                    nnti_err_str(rc));
        }
    }

    return rc;
}

/**
 * @brief An abstract method to put data into a remote memory descriptor.
 *
 * The server stub uses this function to put data to a
 * client memory descriptor.
 *
 * @param buf    @input the buffer for the data.
 * @param len    @input   the amount of data to send.
 * @param dest_md   @input the remote memory descriptor.
 */
extern int nssi_put_data(
        const NNTI_peer_t *caller,
        const void *buf,
        const int len,
        const NNTI_buffer_t *data_addr,
        const long timeout)
{
    int rc = NSSI_OK;
    NNTI_buffer_t       rpc_msg;
    NNTI_buffer_t      *rpc_msg_hdl=NULL;
    NNTI_work_request_t rpc_msg_wr;
    NNTI_status_t       status;
    trios_declare_timer(call_time);

    shadow_buffer_entry *sbe=NULL;

    if (len == 0)
        return rc;

    if (len < 0)
        return NSSI_EINVAL;

    /* if there is a shadow data buffer in use, then put the shadow data into the client's data_addr */
    sbe=get_shadow_buffer(data_addr);
    if (sbe != NULL) {
        data_addr = sbe->data_buffer;
    }

    if ((nssi_config.use_buffer_queue) &&
        (nssi_config.rdma_buffer_queue_buffer_size >= (uint32_t)len)) {
        log_debug(rpc_debug_level, "using buffer queue for PUT buffer");
        rpc_msg_hdl=trios_buffer_queue_pop(&rdma_put_bq);
        assert(rpc_msg_hdl);
        NNTI_BUFFER_SIZE(rpc_msg_hdl)=len;
        /* copy the user buffer contents into RDMA buffer */
        trios_start_timer(call_time);
        memcpy(NNTI_BUFFER_C_POINTER(rpc_msg_hdl), buf, len);
        trios_stop_timer("memcpy put src to bq", call_time);
    } else {
        log_debug(rpc_debug_level, "using user buffer for PUT buffer");
        rpc_msg_hdl=&rpc_msg;
        rc=NNTI_register_memory(
                &transports[data_addr->transport_id],
                (char *)buf,
                len,
                1,
                NNTI_PUT_SRC,
                rpc_msg_hdl);
        if (rc != NNTI_OK) {
            log_error(rpc_debug_level, "failed registering data: %s",
                    nnti_err_str(rc));
        }
    }
    trios_start_timer(call_time);
#ifdef GNI_PERF
    gemini_read_counters(MPI_COMM_WORLD, &gni_state);
#endif
    rc=NNTI_put(
            rpc_msg_hdl,
            0,
            len,
            data_addr,
            0,
            &rpc_msg_wr);
#ifdef GNI_PERF
    gemini_read_counters(MPI_COMM_WORLD, &gni_state);
    gemini_print_counters(MPI_COMM_WORLD, &gni_state, "nssi_put_data - NNTI_put");
#endif
    trios_stop_timer("NNTI_put - put to put dest", call_time);
    if (rc != NSSI_OK) {
        log_error(rpc_debug_level, "failed putting data: %s",
                nnti_err_str(rc));
    }
    trios_start_timer(call_time);
    rc=NNTI_wait(
            &rpc_msg_wr,
            -1,
            &status);
#ifdef GNI_PERF
    gemini_read_counters(MPI_COMM_WORLD, &gni_state);
    gemini_print_counters(MPI_COMM_WORLD, &gni_state, "nssi_put_data - NNTI_wait");
#endif
    trios_stop_timer("NNTI_wait - put to put dest", call_time);
    if (rc != NNTI_OK) {
        log_error(rpc_debug_level, "failed waiting for data: %s",
                nnti_err_str(rc));
    }
    if ((nssi_config.use_buffer_queue) &&
        (nssi_config.rdma_buffer_queue_buffer_size >= (uint32_t)len)) {
        trios_buffer_queue_push(&rdma_put_bq, rpc_msg_hdl);
    } else {
        rc=NNTI_unregister_memory(rpc_msg_hdl);
        if (rc != NNTI_OK) {
            log_error(rpc_debug_level, "failed unregistering data: %s",
                    nnti_err_str(rc));
        }
    }

    return rc;
}


/**
 * @brief Register an RPC service.
 *
 * This method creates a named RPC service on the specified
 * registry server.  Along with the name, the server has to
 * specify where (in the form of an \ref nssi_md) the
 * client should "put" requests.
 *
 * @todo For now, the registry is on the same host as the service
 *       (registrcy_id is ignored). At some point, we need to separate
 *       the registry service.
 *
 * @param registry_id @input Process ID of the registry server.
 * @param name        @input Name of the service.
 * @param remote_sd   @input The remote service description to associate with the name.
 * @param req         @output The request handle (used to test for completion).
 */
 /*
int nssi_register_service(
        const nssi_remote_pid registry_id,
        const char *name,
        const nssi_service *svc,
        nssi_request *req)
{
    return NSSI_ERR_NOTSUPP;
}
*/

/**
 * @brief Add an operation to an operation list.
 *
 * This method adds an operation to an operation list used by
 * a registered RPC service. Operations in this list must have
 * the following prototype:
 *
 * \code
 *    int svc_fun(svc_args *args, svc_result *res);
 * \endcode
 *
 * where \b svc_args is a data structure that contains all required
 * arguments for the method, and \b svc_result is a structure
 * that contains the result.
 *
 * @param svc_op      @input The \ref nssi_svc_op operation description.
 * @param op_list     @output The operation list (modified by this function).
 *
 * @returns \ref NSSI_OK if the operation was successfully added.
 */
int nssi_add_svc_op(
        const nssi_svc_op *svc_op,
        nssi_svc_op_list **op_list)
{
    nssi_svc_op_list *new_list = NULL;

    /* create a new entry */
    new_list = (nssi_svc_op_list *) malloc(1*sizeof(nssi_svc_op_list));

    /* initialize the entry */
    new_list->svc_op.opcode = svc_op->opcode;
    new_list->svc_op.func = svc_op->func;
    new_list->svc_op.sizeof_args = svc_op->sizeof_args;
    new_list->svc_op.decode_args = svc_op->decode_args;
    new_list->svc_op.sizeof_res  = svc_op->sizeof_res;
    new_list->svc_op.encode_res  = svc_op->encode_res;

    /* push entry onto the front of the list */
    new_list->next = *op_list;
    *op_list = new_list;

    return NSSI_OK;
}


/**
 * @brief Initialize an RPC server.
 *
 * @ingroup rpc_server_api
 *
 * This method initializes the portals library and allocates space
 * for incoming requests.
 *
 * @param portal_index @input the portals index to use for incoming reqs.
 * @param service      @output the local service descriptor (used by the server).
 */
int nssi_service_init(
        const nssi_rpc_transport rpc_transport,
        const int                short_req_len,
        nssi_service            *service)
{
    int rc = NSSI_OK;

    int reqs_per_queue = 10000;

    nthread_lock_init(&supported_ops_mutex);
    nthread_lock_init(&request_args_map_mutex);
    nthread_lock_init(&shadow_buffer_mutex);

    /* initialize the service descriptors */
    memset(service, 0, sizeof(nssi_service));

    service->transport_id = transports[rpc_transport].id;
    /* use XDR to encode control messages */
    service->rpc_encode = NSSI_RPC_XDR;
    service->svc_host = transports[rpc_transport].me;
    service->req_size = short_req_len;
    service->res_size = short_req_len;

    if (logging_debug(rpc_debug_level)) {
        fprint_NNTI_peer(logger_get_file(), "transports[rpc_transport].me",
                "nssi_service_init", &transports[rpc_transport].me);
        fprint_NNTI_peer(logger_get_file(), "service->svc_host",
                "nssi_service_init", &service->svc_host);
    }


    /* register trace groups (let someone else enable) */
    trace_register_group(TRACE_RPC_COUNTER_GNAME, &trace_counter_gid);
    trace_register_group(TRACE_RPC_INTERVAL_GNAME, &trace_interval_gid);

    /* Register standard services */
    NSSI_REGISTER_SERVER_STUB(NSSI_OP_GET_SERVICE,  rpc_get_service,  void,                   nssi_service);
    NSSI_REGISTER_SERVER_STUB(NSSI_OP_KILL_SERVICE, rpc_kill_service, nssi_kill_service_args, void);
    NSSI_REGISTER_SERVER_STUB(NSSI_OP_TRACE_RESET,  rpc_trace_reset,  nssi_trace_reset_args,  void);

    rc=NNTI_alloc(
            &transports[service->transport_id],
            service->req_size,
            2*reqs_per_queue,
            NNTI_RECV_QUEUE,
            &service->req_addr);
    if (rc != NNTI_OK) {
        log_error(rpc_debug_level, "failed registering request queue: %s",
                nnti_err_str(rc));
    }

    /* copy the service description to the local service description */
    memcpy(&local_service, service, sizeof(nssi_service));

    return rc;
}

/**
  * @brief Add operations to service.
  *
  * @param svc  @input_type  The service descriptor.
  * @param ops  @input_type  The array operations to add to the service.
  * @param len  @input_type  The number of operations to add.
  */
int nssi_service_add_op(
        const nssi_service *unused,
        const nssi_svc_op *op)
{
        int rc = NSSI_OK;

        assert(op);

        if (nthread_lock(&supported_ops_mutex)) log_warn(rpc_debug_level, "failed to get thread lock");

        if (supported_ops.find((int)op->opcode) == supported_ops.end()) {
            supported_ops[(int)op->opcode] = *op;
        }

        else {
            rc = NSSI_EEXIST;
        }

        nthread_unlock(&supported_ops_mutex);

        return rc;
}



/**
 * @brief Close down an active service.
 *
 * @ingroup rpc_server_api
 *
 * Shuts down the service and releases any associated resources.
 *
 * @param local_sd @input The local service descriptor.
 */
int nssi_service_fini(const nssi_service *service)
{
    int rc = NSSI_OK;

    nthread_lock_fini(&supported_ops_mutex);
    nthread_lock_fini(&request_args_map_mutex);
    nthread_lock_fini(&shadow_buffer_mutex);

    rc=NNTI_free((NNTI_buffer_t *)&service->req_addr);
    if (rc != NNTI_OK) {
        log_error(rpc_debug_level, "failed unregistering request queue: %s",
                nnti_err_str(rc));
    }

    time_to_die=false;

    return NSSI_OK;
}

#define NUM_QUEUES 2

/**
 * @brief Create a daemon process.
 */
void nssi_daemonize()
{
#ifndef __LIBCATAMOUNT__
    int daemon_pid = fork();

    if (daemon_pid < 0) {  /* fork error */
        log_error(rpc_debug_level, "could not fork process");
        return;
    }

    if (daemon_pid > 0) {  /* parent exits */
        exit(0);
    }

    /* child (daemon) continues */

    /* obtain a new process groupd for the daemon */
    setsid();
    umask(0);

    /* close stdin, stdout, and stderr */
    close(0), close(1); close(2);
#endif
#if 0
    i = open("/dev/null", O_RDWR);  /* open stdin */

    dup(i); /* stdout */
    dup(i); /* stderr */
#endif
}

double nssi_get_request_age(const NNTI_peer_t *caller, const int req_id)
{
    double age = -1.0;

    caller_reqid cr(caller, req_id);
    request_args_t *req_args = request_args_map[cr];

    if (req_args) {
        age = trios_get_time() - req_args->arrival_time;
    }
    else {
        age = -99.99;
    }

    return age;
}


/**
 * @brief Start the RPC server using the default request processing function.
 */
int nssi_service_start(nssi_service *svc)
{
    return nssi_service_start_wfn(svc, &nssi_process_rpc_request);
}

/**
 * @brief Start the RPC server.
 *
 * The \b nssi_service_start implements a loop that waits for incoming
 * RPC requests, then calls the process_req function pointer to
 * process those requests.
 *
 * @param service  The service descriptor.
 * @param process_req  A function pointer that takes an nssi_svc_rpc_request.
 */
int nssi_service_start_wfn(
        nssi_service *svc,
        int (*process_req)(nssi_svc_rpc_request *))
{
    int rc = NSSI_OK;
    int req_count = 0;
    double t1;
    double idle_time = 0;
    double processing_time = 0;
    trios_declare_timer(call_time);
    trios_declare_timer(loop_time);

    log_level debug_level = rpc_debug_level;

    char *req_buf;

    NNTI_work_request_t req_queue_wr;
    NNTI_status_t       status;

    progress_callback progress_cb       =NULL;
    int64_t           progress_timeout  =2000; // needs to be reasonable (2 sec)
    int64_t           progress_last_time=0;
    if (svc->progress_callback != 0) {
        progress_cb=(progress_callback)svc->progress_callback;
        progress_timeout=svc->progress_callback_timeout;
        if (progress_timeout < 100) {
            progress_timeout=100;
        }
    }

    log_debug(debug_level, "starting single-threaded rpc service");

    /* initialize indices and counters */
    req_count = 0; /* number of reqs processed */

    /* SIGINT (Ctrl-C) will get us out of this loop */
    while (!trios_exit_now()) {

        trios_start_timer(loop_time);

        log_debug(rpc_debug_level, "a");

        /* exit if the time_to_die flag is set */
        if (time_to_die) {
            rc = NSSI_OK;
            log_info(debug_level, "responding to kill request");
            goto cleanup;
        }

        /* exit if we've received our max number of reqs */
        if ((svc->max_reqs >= 0) && (req_count >= svc->max_reqs)) {
            rc = NSSI_OK;
            log_info(debug_level, "recved max_reqs=%d, exiting",svc->max_reqs);
            goto cleanup;
        }

        trace_start_interval(trace_interval_gid, 0);

        /* measure idle time */
        if (req_count > 0) {
            t1 = trios_get_time();
        }

        NNTI_create_work_request(
                &svc->req_addr,
                &req_queue_wr);
        trios_start_timer(call_time);
        rc=NNTI_wait(
                &req_queue_wr,
                progress_timeout,
                &status);
        trios_stop_timer("request queue wait", call_time);
        if (status.result == NNTI_ETIMEDOUT) {

        }
        else if (rc != NNTI_OK) {
            log_error(debug_level, "failed waiting for a request: %s",
                    nnti_err_str(rc));
        }

        if ((progress_cb) &&
            status.result == NNTI_ETIMEDOUT) {

            progress_cb(1);
            progress_last_time=trios_get_time_ms();
        }

        if (trios_exit_now()) {
            log_debug(debug_level, "time to exit");
            goto cleanup;
        }

        if (status.result == NNTI_OK) {
            req_buf = (char *)status.start+status.offset;
            log_debug(debug_level, "req_buf=%p", req_buf);

            //trace_end_interval(trace_interval_gid, TRACE_THREAD_IDLE, 0, "idle time");
            log_debug(debug_level, "after end interval");

            /* capture the idle time */
            idle_time += trios_get_time() - t1;
            log_debug(debug_level, "out of job_loop");

            /* increment the number of requests */
            req_count++;

            // this structure gets freed in the process_rpc_request function
            struct rpc_request *rpc_req = new struct rpc_request();
            rpc_req->svc           = svc;
            rpc_req->caller        = status.src;
            rpc_req->id            = req_count;
            rpc_req->arrival_time  = trios_get_time();

            // copy the short request buffer (in case of threaded servers)
            rpc_req->req_buf       = new char[status.length];  // freed in process_rpc_request
            memcpy(rpc_req->req_buf, req_buf, status.length);
            rpc_req->short_req_len = status.length;


            /* measure the processing time */
            t1 = trios_get_time();

            // The process_req function is responsible for freeing the request and
            // any memory associated with the request (this seems kind of risky).

            trios_start_timer(call_time);
            (*process_req)(rpc_req);    // call the function pointer
            trios_stop_timer("process_rpc_request", call_time);
            if (rc != NSSI_OK) {
                /* warn only... we do not exit */
                log_warn(rpc_debug_level, "main: unable to process request.");
            }

            /* measure the processing time */
            processing_time += trios_get_time() - t1;

            if ((progress_cb) &&
                (trios_get_time_ms() - progress_last_time) > progress_timeout) {

                progress_cb(0);
                progress_last_time=trios_get_time_ms();
            }
        }
        NNTI_destroy_work_request(
                &req_queue_wr);

        trios_stop_timer("service loop", loop_time);
    }


cleanup:

    debug_level = rpc_debug_level;

    /* finish any tracing */
    trace_fini();

    log_debug(debug_level, "Cleaning up...");

    /* print out stats about the server */
    log_info(debug_level, "Exiting nssi_service_start: %d "
            "reqs processed, exit_now=%d", req_count,trios_exit_now());

    FILE *fp = logger_get_file();
    fprintf(fp, "----- SERVER STATS ---------\n");
    fprintf(fp, "\tprocessed requests = %d\n", req_count);
    //fprintf(fp, "\tidle time       = %g (sec)\n", idle_time);
    fprintf(fp, "\tprocessing time = %g (sec)\n", processing_time);
    fprintf(fp, "----------------------------\n");

    return rc;
}
