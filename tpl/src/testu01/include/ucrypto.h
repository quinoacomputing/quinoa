 
/*  ucrypto.h  for ANSI C  */
#ifndef UCRYPTO_H
#define UCRYPTO_H
 
#include "unif01.h"


typedef enum {
   ucrypto_OFB,                 /* Output Feedback mode */
   ucrypto_CTR,                 /* Counter mode */
   ucrypto_KTR                  /* Key counter mode */
   } ucrypto_Mode;



unif01_Gen * ucrypto_CreateAES (unsigned char *Key, int klen,
                                unsigned char *Seed, ucrypto_Mode mode,
                                int r, int s);



unif01_Gen * ucrypto_CreateSHA1 (unsigned char *Seed, int len,
                                 ucrypto_Mode mode, int r, int s);



unif01_Gen * ucrypto_CreateISAAC (int flag, unsigned int A[256]);


void ucrypto_DeleteAES (unif01_Gen * gen);
void ucrypto_DeleteSHA1 (unif01_Gen * gen);
void ucrypto_DeleteISAAC (unif01_Gen * gen);
 
#endif
 
