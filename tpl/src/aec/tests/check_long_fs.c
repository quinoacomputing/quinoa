#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "check_aec.h"

#define BUF_SIZE (64 * 4)

int check_long_fs(struct test_state *state)
{
    int status, size, i, bs;

    size = state->bytes_per_sample;
    bs = state->strm->block_size;

    for (i = 0; i < bs / 2; i++) {
        state->out(state->ubuf + size * i, state->xmin, size);
        state->out(state->ubuf + bs * size / 2 + size * i, 65000, size);
    }

    printf("Checking long fs ... ");

    status = state->codec(state);
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
    strm.bits_per_sample = 16;
    strm.block_size = 64;
    strm.rsi = 1;
    state.codec = encode_decode_large;
    update_state(&state);

    status = check_long_fs(&state);
    if (status)
        goto DESTRUCT;

DESTRUCT:
    free(state.ubuf);
    free(state.cbuf);
    free(state.obuf);

    return status;
}
