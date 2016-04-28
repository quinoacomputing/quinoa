/**
 * @file aec.c
 *
 * @section LICENSE
 * Copyright 2012 - 2015
 *
 * Mathis Rosenhauer, Moritz Hanke, Joerg Behrens
 * Deutsches Klimarechenzentrum GmbH
 * Bundesstr. 45a
 * 20146 Hamburg Germany
 *
 * Luis Kornblueh
 * Max-Planck-Institut fuer Meteorologie
 * Bundesstr. 53
 * 20146 Hamburg
 * Germany
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above
 *    copyright notice, this list of conditions and the following
 *    disclaimer in the documentation and/or other materials provided
 *    with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 *
 * CLI frontend for Adaptive Entropy Coding library
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libaec.h"

#define CHUNK 10485760

int get_param(unsigned int *param, int *iarg, char *argv[])
{
    if (strlen(argv[*iarg]) == 2) {
        (*iarg)++;
        if (argv[*iarg][0] == '-')
            return 1;
        else
            *param = atoi(argv[*iarg]);
    } else {
        *param = atoi(&argv[*iarg][2]);
    }
    return 0;
}

int main(int argc, char *argv[])
{
    struct aec_stream strm;
    unsigned char *in;
    unsigned char *out;
    size_t total_out;
    unsigned int chunk;
    int status;
    int input_avail, output_avail;
    char *infn, *outfn;
    FILE *infp, *outfp;
    int dflag;
    char *opt;
    int iarg;

    chunk = CHUNK;
    strm.bits_per_sample = 8;
    strm.block_size = 8;
    strm.rsi = 2;
    strm.flags = AEC_DATA_PREPROCESS;
    dflag = 0;
    iarg = 1;

    while (iarg < argc - 2) {
        opt = argv[iarg];
        if (opt[0] != '-')
            goto FAIL;
        switch (opt[1]) {
        case '3':
            strm.flags |= AEC_DATA_3BYTE;
            break;
        case 'N':
            strm.flags &= ~AEC_DATA_PREPROCESS;
            break;
        case 'b':
            if (get_param(&chunk, &iarg, argv))
                goto FAIL;
            break;
        case 'd':
            dflag = 1;
            break;
        case 'j':
            if (get_param(&strm.block_size, &iarg, argv))
                goto FAIL;
            break;
        case 'm':
            strm.flags |= AEC_DATA_MSB;
            break;
        case 'n':
            if (get_param(&strm.bits_per_sample, &iarg, argv))
                goto FAIL;
            break;
        case 'p':
            strm.flags |= AEC_PAD_RSI;
            break;
        case 'r':
            if (get_param(&strm.rsi, &iarg, argv))
                goto FAIL;
            break;
        case 's':
            strm.flags |= AEC_DATA_SIGNED;
            break;
        case 't':
            strm.flags |= AEC_RESTRICTED;
            break;
        default:
            goto FAIL;
        }
        iarg++;
    }

    if (argc - iarg < 2)
        goto FAIL;

    infn = argv[iarg];
    outfn = argv[iarg + 1];

    if (strm.bits_per_sample > 16) {
        if (strm.bits_per_sample <= 24 && strm.flags & AEC_DATA_3BYTE)
            chunk *= 3;
        else
            chunk *= 4;
    } else if (strm.bits_per_sample > 8) {
        chunk *= 2;
    }

    out = (unsigned char *)malloc(chunk);
    in = (unsigned char *)malloc(chunk);

    if (in == NULL || out == NULL)
        exit(-1);

    total_out = 0;
    strm.avail_in = 0;
    strm.avail_out = chunk;
    strm.next_out = out;

    input_avail = 1;
    output_avail = 1;

    if ((infp = fopen(infn, "rb")) == NULL) {
        fprintf(stderr, "ERROR: cannot open input file %s\n", infn);
        return 1;
    }
    if ((outfp = fopen(outfn, "wb")) == NULL) {
        fprintf(stderr, "ERROR: cannot open output file %s\n", infn);
        return 1;
    }

    if (dflag)
        status = aec_decode_init(&strm);
    else
        status = aec_encode_init(&strm);

    if (status != AEC_OK) {
        fprintf(stderr, "ERROR: initialization failed (%d)\n", status);
        return 1;
    }

    while(input_avail || output_avail) {
        if (strm.avail_in == 0 && input_avail) {
            strm.avail_in = fread(in, 1, chunk, infp);
            if (strm.avail_in != chunk)
                input_avail = 0;
            strm.next_in = in;
        }

        if (dflag)
            status = aec_decode(&strm, AEC_NO_FLUSH);
        else
            status = aec_encode(&strm, AEC_NO_FLUSH);

        if (status != AEC_OK) {
            fprintf(stderr, "ERROR: %i\n", status);
            return 1;
        }

        if (strm.total_out - total_out > 0) {
            fwrite(out, strm.total_out - total_out, 1, outfp);
            total_out = strm.total_out;
            output_avail = 1;
            strm.next_out = out;
            strm.avail_out = chunk;
        } else {
            output_avail = 0;
        }

    }

    if (dflag) {
        aec_decode_end(&strm);
    } else {
        if ((status = aec_encode(&strm, AEC_FLUSH)) != AEC_OK) {
            fprintf(stderr, "ERROR: while flushing output (%i)\n", status);
            return 1;
        }

        if (strm.total_out - total_out > 0)
            fwrite(out, strm.total_out - total_out, 1, outfp);

        aec_encode_end(&strm);
    }

    fclose(infp);
    fclose(outfp);
    free(in);
    free(out);
    return 0;

FAIL:
    fprintf(stderr, "NAME\n\taec - encode or decode files ");
    fprintf(stderr, "with Adaptive Entropy Coding\n\n");
    fprintf(stderr, "SYNOPSIS\n\taec [OPTION]... SOURCE DEST\n");
    fprintf(stderr, "\nOPTIONS\n");
    fprintf(stderr, "\t-3\n\t\t24 bit samples are stored in 3 bytes\n");
    fprintf(stderr, "\t-N\n\t\tdisable pre/post processing\n");
    fprintf(stderr, "\t-b size\n\t\tinternal buffer size in bytes\n");
    fprintf(stderr, "\t-d\n\t\tdecode SOURCE. If -d is not used: encode.\n");
    fprintf(stderr, "\t-j samples\n\t\tblock size in samples\n");
    fprintf(stderr, "\t-m\n\t\tsamples are MSB first. Default is LSB\n");
    fprintf(stderr, "\t-n bits\n\t\tbits per sample\n");
    fprintf(stderr, "\t-p\n\t\tpad RSI to byte boundary\n");
    fprintf(stderr, "\t-r blocks\n\t\treference sample interval in blocks\n");
    fprintf(stderr, "\t-s\n\t\tsamples are signed. Default is unsigned\n");
    fprintf(stderr, "\t-t\n\t\tuse restricted set of code options\n\n");
    return 1;
}
