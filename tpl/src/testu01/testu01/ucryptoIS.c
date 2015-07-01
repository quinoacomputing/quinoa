/*************************************************************************\
 *
 * Package:        TestU01
 * File:           ucryptoIS.c
 * Environment:    ANSI C
 *
\*************************************************************************/

#include "util.h"
#include "addstr.h"

#include "ucrypto.h"
#include "unif01.h"

#include <stdio.h>
#include <string.h>


#define LEN  200                  /* Max length of strings */


static int co1 = 0;               /* Counter */


/*=========================================================================*/
/* Generator ISAAC starts here.
 *  The only changes I have made in Bob Jenkins' code is to add the keyword
 *  static before functions randinit and isaac to avoid possible name clashes,
 *  and to remove the #include directives for standard.h and rand.h, since
 *  I have included these files explicitly here. (R. Simard)
 *  I have also changed the long int's into int's.
 */

/*
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
x
x  File standard.h, from Bob Jenkins web page, for his generator ISAAC
x
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
*/

/*
------------------------------------------------------------------------------
Standard definitions and types, Bob Jenkins
------------------------------------------------------------------------------
*/
#ifndef STANDARD
# define STANDARD
# ifndef STDIO
#  include <stdio.h>
#  define STDIO
# endif
# ifndef STDDEF
#  include <stddef.h>
#  define STDDEF
# endif
typedef  unsigned long long  ub8;
#define UB8MAXVAL 0xffffffffffffffffLL
#define UB8BITS 64
typedef    signed long long  sb8;
#define SB8MAXVAL 0x7fffffffffffffffLL
typedef  unsigned  int  ub4;   /* unsigned 4-byte quantities */
#define UB4MAXVAL 0xffffffff
typedef    signed  int  sb4;
#define UB4BITS 32
#define SB4MAXVAL 0x7fffffff
typedef  unsigned short int  ub2;
#define UB2MAXVAL 0xffff
#define UB2BITS 16
typedef    signed short int  sb2;
#define SB2MAXVAL 0x7fff
typedef  unsigned       char ub1;
#define UB1MAXVAL 0xff
#define UB1BITS 8
typedef    signed       char sb1;   /* signed 1-byte quantities */
#define SB1MAXVAL 0x7f
typedef                 int  word;  /* fastest type available */

#define bis(target,mask)  ((target) |=  (mask))
#define bic(target,mask)  ((target) &= ~(mask))
#define bit(target,mask)  ((target) &   (mask))
#ifndef min
# define min(a,b) (((a)<(b)) ? (a) : (b))
#endif /* min */
#ifndef max
# define max(a,b) (((a)<(b)) ? (b) : (a))
#endif /* max */
#ifndef align
# define align(a) (((ub4)a+(sizeof(void *)-1))&(~(sizeof(void *)-1)))
#endif /* align */
#ifndef abs
# define abs(a)   (((a)>0) ? (a) : -(a))
#endif
#define TRUE  1
#define FALSE 0
#define SUCCESS 0  /* 1 on VAX */

#endif /* STANDARD */


/*
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
x
x  File rand.h, from Bob Jenkins web page, for his generator ISAAC 
x
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
*/

/*
------------------------------------------------------------------------------
rand.h: definitions for a random number generator
MODIFIED:
  960327: Creation (addition of randinit, really)
  970719: use context, not global variables, for internal state
  980324: renamed seed to flag
  980605: recommend RANDSIZL=4 for noncryptography.
------------------------------------------------------------------------------

#ifndef STANDARD
#include "standard.h"
#endif
*/

#ifndef RAND
#define RAND
#define RANDSIZL   (8)  /* I recommend 8 for crypto, 4 for simulations */
#define RANDSIZ    (1<<RANDSIZL)

/* context of random number generator */
struct randctx
{
  ub4 randcnt;
  ub4 randrsl[RANDSIZ];
  ub4 randmem[RANDSIZ];
  ub4 randa;
  ub4 randb;
  ub4 randc;
};
typedef  struct randctx  randctx;

/*
------------------------------------------------------------------------------
 If (flag==TRUE), then use the contents of randrsl[0..RANDSIZ-1] as the seed.
------------------------------------------------------------------------------
*/
static void randinit(/*_ randctx *r, word flag _*/);

static void isaac(/*_ randctx *r _*/);


/*
------------------------------------------------------------------------------
 Call rand(/o_ randctx *r _o/) to retrieve a single 32-bit random value
------------------------------------------------------------------------------
*/
#define rand(r) \
   (!(r)->randcnt-- ? \
     (isaac(r), (r)->randcnt=RANDSIZ-1, (r)->randrsl[(r)->randcnt]) : \
     (r)->randrsl[(r)->randcnt])

#endif  /* RAND */



/*
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
x
x  File rand.c, from Bob Jenkins web page, for his generator ISAAC 
x
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
*/

/*
------------------------------------------------------------------------------
rand.c: By Bob Jenkins.  My random number generator, ISAAC.
MODIFIED:
  960327: Creation (addition of randinit, really)
  970719: use context, not global variables, for internal state
  980324: added main (ifdef'ed out), also rearranged randinit()
------------------------------------------------------------------------------

#ifndef STANDARD
#include "standard.h"
#endif
#ifndef RAND
#include "rand.h"
#endif
*/

#define ind(mm,x)  (*(ub4 *)((ub1 *)(mm) + ((x) & ((RANDSIZ-1)<<2))))
#define rngstep(mix,a,b,mm,m,m2,r,x) \
{ \
  x = *m;  \
  a = (a^(mix)) + *(m2++); \
  *(m++) = y = ind(mm,x) + a + b; \
  *(r++) = b = ind(mm,y>>RANDSIZL) + x; \
}

static void     isaac(ctx)
randctx *ctx;
{
   register ub4 a,b,x,y,*m,*mm,*m2,*r,*mend;
   mm=ctx->randmem; r=ctx->randrsl;
   a = ctx->randa; b = ctx->randb + (++ctx->randc);
   for (m = mm, mend = m2 = m+(RANDSIZ/2); m<mend; )
   {
      rngstep( a<<13, a, b, mm, m, m2, r, x);
      rngstep( a>>6 , a, b, mm, m, m2, r, x);
      rngstep( a<<2 , a, b, mm, m, m2, r, x);
      rngstep( a>>16, a, b, mm, m, m2, r, x);
   }
   for (m2 = mm; m2<mend; )
   {
      rngstep( a<<13, a, b, mm, m, m2, r, x);
      rngstep( a>>6 , a, b, mm, m, m2, r, x);
      rngstep( a<<2 , a, b, mm, m, m2, r, x);
      rngstep( a>>16, a, b, mm, m, m2, r, x);
   }
   ctx->randb = b; ctx->randa = a;
}


#define mix(a,b,c,d,e,f,g,h) \
{ \
   a^=b<<11; d+=a; b+=c; \
   b^=c>>2;  e+=b; c+=d; \
   c^=d<<8;  f+=c; d+=e; \
   d^=e>>16; g+=d; e+=f; \
   e^=f<<10; h+=e; f+=g; \
   f^=g>>4;  a+=f; g+=h; \
   g^=h<<8;  b+=g; h+=a; \
   h^=a>>9;  c+=h; a+=b; \
}

/* if (flag==TRUE), then use the contents of randrsl[] to initialize mm[]. */
static void randinit(ctx, flag)
randctx *ctx;
word     flag;
{
   word i;
   ub4 a,b,c,d,e,f,g,h;
   ub4 *m,*r;
   ctx->randa = ctx->randb = ctx->randc = 0;
   m=ctx->randmem;
   r=ctx->randrsl;
   a=b=c=d=e=f=g=h=0x9e3779b9UL;  /* the golden ratio */

   for (i=0; i<4; ++i)          /* scramble it */
   {
     mix(a,b,c,d,e,f,g,h);
   }

   if (flag) 
   {
     /* initialize using the contents of r[] as the seed */
     for (i=0; i<RANDSIZ; i+=8)
     {
       a+=r[i  ]; b+=r[i+1]; c+=r[i+2]; d+=r[i+3];
       e+=r[i+4]; f+=r[i+5]; g+=r[i+6]; h+=r[i+7];
       mix(a,b,c,d,e,f,g,h);
       m[i  ]=a; m[i+1]=b; m[i+2]=c; m[i+3]=d;
       m[i+4]=e; m[i+5]=f; m[i+6]=g; m[i+7]=h;
     }
     /* do a second pass to make all of the seed affect all of m */
     for (i=0; i<RANDSIZ; i+=8)
     {
       a+=m[i  ]; b+=m[i+1]; c+=m[i+2]; d+=m[i+3];
       e+=m[i+4]; f+=m[i+5]; g+=m[i+6]; h+=m[i+7];
       mix(a,b,c,d,e,f,g,h);
       m[i  ]=a; m[i+1]=b; m[i+2]=c; m[i+3]=d;
       m[i+4]=e; m[i+5]=f; m[i+6]=g; m[i+7]=h;
     }
   }
   else
   {
     /* fill in mm[] with messy stuff */
     for (i=0; i<RANDSIZ; i+=8)
     {
       mix(a,b,c,d,e,f,g,h);
       m[i  ]=a; m[i+1]=b; m[i+2]=c; m[i+3]=d;
       m[i+4]=e; m[i+5]=f; m[i+6]=g; m[i+7]=h;
     }
   }

   isaac(ctx);            /* fill in the first set of results */
   ctx->randcnt=RANDSIZ;  /* prepare to use the first set of results */
}




/*-------------------------------------------------------------------------*/
static randctx ctx;

static unsigned long ISAAC_Bits (void *junk1, void *junk2)
{
  return rand(&ctx);
}

/*-----------------------------------------------------------------------*/

static double ISAAC_U01 (void *vpar, void *vsta)
{
   return ISAAC_Bits (vpar, vsta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/

static void WrISAAC (void *junk) {
   ub4 i;
   if (!unif01_WrLongStateFlag) {
      unif01_WrLongStateDef ();
      return;
   }
   printf ("A = {\n ");
   for (i=0; i<RANDSIZ; ++i) {
      printf (" %12u", ctx.randrsl[i]);
      if (i < RANDSIZ - 1)
         printf (",");
      if ((i % 5) == 4)
         printf ("\n ");
      };
   printf ("    }\n");
}

/*-----------------------------------------------------------------------*/

unif01_Gen * ucrypto_CreateISAAC (int flag, unsigned int A[]) {
   ub4 i;
   unif01_Gen *gen;
   size_t leng;
   char name[LEN + 1];

   util_Assert (co1 == 0,
      "ucrypto_CreateISAAC:   only 1 generator at a time can be in use");
   co1++;

   gen = util_Malloc (sizeof (unif01_Gen));

   strcpy (name, "ucrypto_CreateISAAC:");
   addstr_Int (name, "   flag = ", flag);
   if (A != NULL)
      addstr_ArrayUint (name, ",   A = ", RANDSIZ, A);
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   switch (flag) {
   case 2:
      for (i=0; i<RANDSIZ; ++i) ctx.randrsl[i] = A[i];
      break;
   case 1:
      ctx.randa = ctx.randb = ctx.randc = (ub4)0;
      for (i=0; i<RANDSIZ; ++i) ctx.randrsl[i] = A[i];
      randinit(&ctx, TRUE);
      break;
   case 0:
     /* this is Jenkins' original version without outside initialization.
        However, in his rand.c programm, Jenkins prints the numbers from
        i = 0 to RANDSIZ-1, while we use them in the reverse order. */
      ctx.randa = ctx.randb = ctx.randc = (ub4)0;
      for (i=0; i<RANDSIZ; ++i) ctx.randrsl[i] = (ub4)0;
      randinit(&ctx, TRUE);
      break;
   default:
      util_Error ("ucrypto_CreateISAAC:   flag must be in {0, 1, 2}");
   }
   ctx.randcnt = 0;

   gen->GetBits = &ISAAC_Bits;
   gen->GetU01  = &ISAAC_U01;
   gen->Write   = &WrISAAC;
   gen->param   = NULL;
   gen->state   = NULL;
   return gen;
}

/*-----------------------------------------------------------------------*/

void ucrypto_DeleteISAAC (unif01_Gen *gen)
{
   if (NULL == gen)
      return;
   gen->name = util_Free (gen->name);
   util_Free (gen);
   co1--;
}
