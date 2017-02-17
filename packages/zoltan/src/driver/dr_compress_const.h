/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */

#ifndef _DR_COMPRESS_CONST_H
#define _DR_COMPRESS_CONST_H

#include "zoltan.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "dr_const.h"

#if (defined ZOLTAN_GZIP)||(defined ZOLTAN_BZ2)||(defined ZOLTAN_LZMA)
#define ZOLTAN_COMPRESS
#endif

#ifdef ZOLTAN_GZIP
#include <zlib.h>
#endif /* ZOLTAN_GZIP */
#ifdef ZOLTAN_BZ2
#include <bzlib.h>
#endif /* ZOLTAN_BZ2 */
#ifdef ZOLTAN_LZMA
#include <lzma.h>
#endif /* ZOLTAN_LZMA */


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/*****************************************************************************
 *  Definitions for the File wrapper program
 *****************************************************************************/

typedef enum ZOLTAN_FILETYPE_ {
  STANDARD = 0,
  GZIP = 1,
  BZIP2 = 2,
  LZMA = 3
} ZOLTAN_FILETYPE;

#ifdef ZOLTAN_COMPRESS
typedef struct ZOLTAN_FILE_ {
  ZOLTAN_FILETYPE type;
  char * buffer;
  ssize_t size;  /* size & pos should be off_t, but because of their use */
  ssize_t pos;   /* they need to be signed (pos is initialized to -1).   */
  union {
    FILE * fileunc;
#ifdef ZOLTAN_GZIP
    gzFile filegz;
#endif /* ZOLTAN_GZIP */
#ifdef ZOLTAN_BZ2
    BZFILE * filebz;
#endif /* ZOLTAN_BZ2 */
  } strm ;
} ZOLTAN_FILE;
#else /* ZOLTAN_COMPRESS */
  typedef FILE ZOLTAN_FILE;
#endif /* ZOLTAN_COMPRESS */

ZOLTAN_FILE* ZOLTAN_FILE_open(const char *path, const char *mode, const ZOLTAN_FILETYPE type);
int ZOLTAN_FILE_printf(ZOLTAN_FILE* file, const char * format, ...);
int ZOLTAN_FILE_scanf(ZOLTAN_FILE* stream, const char * format, ... );
int ZOLTAN_FILE_puts(char *s, ZOLTAN_FILE* file);
char* ZOLTAN_FILE_gets(char * buf, int len, ZOLTAN_FILE* file);
int ZOLTAN_FILE_putc(int c, ZOLTAN_FILE* file);
int ZOLTAN_FILE_getc(ZOLTAN_FILE* file);
int ZOLTAN_FILE_ungetc(int c, ZOLTAN_FILE* file);
int ZOLTAN_FILE_flush(ZOLTAN_FILE* file);
int ZOLTAN_FILE_close(ZOLTAN_FILE* file);
void ZOLTAN_FILE_rewind(ZOLTAN_FILE* stream);

ssize_t ZOLTAN_FILE_read(char* ptr, size_t size, size_t nitems, ZOLTAN_FILE *file);

#ifndef ZOLTAN_COMPRESS
#define ZOLTAN_FILE_open(path, mode, type) fopen(path, mode)
/* #define ZOLTAN_FILE_printf(file, format ...) fprintf(file, ## format) */
#define ZOLTAN_FILE_printf fprintf
/* #define ZOLTAN_FILE_scanf(retval, stream, format ... ) (*(retval) = fscanf((stream), ## format)) */
#define ZOLTAN_FILE_scanf fscanf
#define ZOLTAN_FILE_puts(s, file) fputs(s,file)
#define ZOLTAN_FILE_gets(buf, len, file) fgets((buf), (len), (file))
#define ZOLTAN_FILE_putc(c, file) fputc((c), (file))
#define ZOLTAN_FILE_getc(file) fgetc(file)
#define ZOLTAN_FILE_ungetc(c, file) ungetc((c), (file))
#define ZOLTAN_FILE_flush(file) fflush(file)
#define ZOLTAN_FILE_close(file) fclose(file)
#define ZOLTAN_FILE_rewind(stream) rewind(stream)
#define ZOLTAN_FILE_read(ptr, size, nitems, file) fread((ptr), (size), (nitems), (file))

#else /* ZOLTAN_COMPRESS */

/*   /\*** Implemented as a macro as "vfscanf" or "vsscanf" are C99 only ***\/ */
/* #define ZOLTAN_FILE_scanf2(retval, stream, format ... )\ */
/* do { \ */
/*   char buff[1024]; \ */
/*   if (stream->type == STANDARD) { \ */
/*     *(retval) = (fscanf((stream)->strm.fileunc, ## format)); \ */
/*   } \ */
/*   else { if (ZOLTAN_FILE_gets(buff, 1024, (stream)) == NULL) \ */
/*     *(retval) = 0; \ */
/*   *(retval) = sscanf(buff, ## format); }\ */
/* } while(0) */

#endif /* ZOLTAN_COMPRESS */

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
#endif /* _DR_CONST_H */
