#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "szlib.h"

#define OPTIONS_MASK (SZ_RAW_OPTION_MASK \
                      | SZ_MSB_OPTION_MASK \
                      | SZ_NN_OPTION_MASK)
#define PIXELS_PER_BLOCK (8)
#define PIXELS_PER_SCANLINE (PIXELS_PER_BLOCK*128)

int main(int argc, char *argv[])
{
    int status;
    SZ_com_t sz_param;
    unsigned char *source, *dest, *dest1;
    size_t destLen, dest1Len, sourceLen;
    FILE *fp;

    if (argc < 2) {
        fprintf(stderr, "Usage: %s file\n", argv[0]);
        return 1;
    }

    if ((fp = fopen(argv[1], "rb")) == NULL) {
        fprintf(stderr, "Can't open %s\n", argv[1]);
        return 1;
    }
    fseek(fp, 0L, SEEK_END);
    sourceLen = ftell(fp);
    fseek(fp, 0L, SEEK_SET);
    destLen = sourceLen + sourceLen / 10;

    sz_param.options_mask = OPTIONS_MASK;
    sz_param.bits_per_pixel = 64;
    sz_param.pixels_per_block = PIXELS_PER_BLOCK;
    sz_param.pixels_per_scanline = PIXELS_PER_SCANLINE;

    source = (unsigned char *)malloc(sourceLen);
    dest = (unsigned char *)malloc(destLen);
    dest1 = (unsigned char *)malloc(destLen);

    if (source == NULL || dest == NULL || dest1 == NULL)
        return 1;

    sourceLen = fread(source, 1, sourceLen, fp);

    status = SZ_BufftoBuffCompress(dest, &destLen,
                                   source, sourceLen, &sz_param);
    if (status != SZ_OK)
        return status;

    dest1Len = sourceLen;
    status = SZ_BufftoBuffDecompress(dest1, &dest1Len,
                                     dest, destLen, &sz_param);
    if (status != SZ_OK)
        return status;

    if (memcmp(source, dest1, sourceLen) != 0)
        fprintf(stderr, "File %s Buffers differ\n", argv[2]);

    free(source);
    free(dest);
    free(dest1);
    return 0;
}
