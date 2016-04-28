#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "check_aec.h"

#define BUF_SIZE 1024 * 3

int check_block_sizes(struct test_state *state)
{
    int bs, status;

    for (bs = 8; bs <= 64; bs *= 2) {
        state->strm->block_size = bs;
        state->strm->rsi = (int)(state->buf_len
                                 / (bs * state->bytes_per_sample));

        status = encode_decode_large(state);
        if (status)
            return status;
    }
    return 0;
}

int check_block_sizes_short(struct test_state *state)
{
    int bs, status;
    size_t tmp;

    tmp = state->ibuf_len;
    for (bs = 8; bs <= 64; bs *= 2) {
        state->strm->block_size = bs;
        state->strm->rsi = (int)(state->buf_len
                                 / (bs * state->bytes_per_sample));
        state->ibuf_len = state->buf_len - 2 * bs + 4;
        status = encode_decode_large(state);
        if (status)
            return status;
        if (state->strm->total_out != state->buf_len) {
            printf("FAIL: Unexpected buffer length. Got %i expected %i\n",
                   (int)state->strm->total_out,
                   (int)state->buf_len);
            return 99;
        }
    }
    state->ibuf_len = tmp;
    return 0;
}

int check_rsi(struct test_state *state)
{
    int status, size;
    unsigned char *tmp;

    size = state->bytes_per_sample;

    for (tmp = state->ubuf;
         tmp < state->ubuf + state->buf_len;
         tmp += 2 * state->bytes_per_sample) {
        state->out(tmp, state->xmax, size);
        state->out(tmp + size, state->xmin, size);
    }

    printf("Checking full rsi ... ");
    status = check_block_sizes(state);
    if (status)
        return status;

    printf ("%s\n", CHECK_PASS);

    printf("Checking short rsi ... ");
    status = check_block_sizes_short(state);
    if (status)
        return status;

    printf ("%s\n", CHECK_PASS);
    return 0;
}

int main (void)
{
    int status;
    struct aec_stream strm;
    struct test_state state;

    state.buf_len = state.ibuf_len = BUF_SIZE;
    state.cbuf_len = 2 * BUF_SIZE;

    state.ubuf = (unsigned char *)malloc(state.buf_len);
    state.cbuf = (unsigned char *)malloc(state.cbuf_len);
    state.obuf = (unsigned char *)malloc(state.buf_len);

    if (!state.ubuf || !state.cbuf || !state.obuf) {
        printf("Not enough memory.\n");
        return 99;
    }

    strm.flags = AEC_DATA_PREPROCESS;
    state.strm = &strm;
    strm.bits_per_sample = 32;
    update_state(&state);

    status = check_rsi(&state);
    if (status)
        goto DESTRUCT;

DESTRUCT:
    free(state.ubuf);
    free(state.cbuf);
    free(state.obuf);

    return status;
}
