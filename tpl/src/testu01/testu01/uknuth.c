/*************************************************************************\
 *
 * Package:        TestU01
 * File:           uknuth.c
 * Environment:    ANSI C
 * Programmer:     Richard Simard.
 *
\*************************************************************************/

#include "util.h"
#include "addstr.h"

#include "uknuth.h"
#include "unif01.h"

#include <stdio.h>
#include <string.h>


#define LEN  200                  /* Max length of strings */


static int co1 = 0, co2 = 0, co3 = 0, co4 = 0;      /* Counters */


/*=========================  WARNING:

I HAVE CHANGED Knuth's code for the following version of rng.c. Thus it is
NOT Knuth's original code. I have renamed the variables because those of the
new version have the same name as those of the old version and we include
both versions in our file. (R. Simard)

===========================*/


/*    This program is copyright (c) 2000 by D E Knuth;
 *    you may copy it freely  AS LONG AS YOU MAKE ABSOLUTELY NO CHANGES!
 *    You could also change it, but then you must rename the file and
 *    tell people clearly that your version is not the same as mine.
 *    It is explained in Seminumerical Algorithms, 3rd edition, Section 3.6
 *    (or in the errata to the 2nd edition --- see
 *        http://www-cs-faculty.stanford.edu/~knuth/taocp.html
 *    in the changes to pages 171 and following of Volume 2).              */

/*    If you find any bugs, please report them immediately to
 *                 taocp@cs.stanford.edu
 *    (and you will be rewarded if the bug is genuine). Thanks!            */

/************ see the book for explanations and caveats! *******************/
/************ in particular, you need two's complement arithmetic **********/

/* old-style C function declarations appear here for reasons of portability */

#define KK 100                     /* the long lag */
#define LL  37                     /* the short lag */
#define MM (1L<<30)                 /* the modulus */
#define mod_diff(x,y) (((x)-(y))&(MM-1)) /* subtraction mod MM */

long ran_x1[KK];                    /* the generator state */

/* void ran_array(long aa[],int n) */
void ran_array1(aa,n)    /* put n new random numbers in aa */
  long *aa;   /* destination */
  int n;      /* array length (must be at least KK) */
{
  register int i,j;
  for (j=0;j<KK;j++) aa[j]=ran_x1[j];
  for (;j<n;j++) aa[j]=mod_diff(aa[j-KK],aa[j-LL]);
  for (i=0;i<LL;i++,j++) ran_x1[i]=mod_diff(aa[j-KK],aa[j-LL]);
  for (;i<KK;i++,j++) ran_x1[i]=mod_diff(aa[j-KK],ran_x1[i-LL]);
}

#define TT  70   /* guaranteed separation between streams */
#define is_odd(x)  ((x)&1)          /* units bit of x */
#define evenize(x) ((x)&(MM-2))   /* make x even */

/* void ran_start(long seed) */
void ran_start1(seed)    /* do this before using ran_array1 */
  long seed;            /* selector for different streams */
{
  register int t,j;
  long x[KK+KK-1];              /* the preparation buffer */
  register long ss=evenize(seed+2);
  for (j=0;j<KK;j++) {
    x[j]=ss;                      /* bootstrap the buffer */
    ss<<=1; if (ss>=MM) ss-=MM-2; /* cyclic shift 29 bits */
  }
  for (;j<KK+KK-1;j++) x[j]=0;
  x[1]++;              /* make x[1] (and only x[1]) odd */
  ss=seed&(MM-1);
  t=TT-1; while (t) {
    for (j=KK-1;j>0;j--) x[j+j]=x[j];  /* "square" */
    for (j=KK+KK-2;j>KK-LL;j-=2) x[KK+KK-1-j]=evenize(x[j]);
    for (j=KK+KK-2;j>=KK;j--) if(is_odd(x[j])) {
      x[j-(KK-LL)]=mod_diff(x[j-(KK-LL)],x[j]);
      x[j-KK]=mod_diff(x[j-KK],x[j]);
    }
    if (is_odd(ss)) {              /* "multiply by z" */
      for (j=KK;j>0;j--)  x[j]=x[j-1];
      x[0]=x[KK];            /* shift the buffer cyclically */
      if (is_odd(x[KK])) x[LL]=mod_diff(x[LL],x[KK]);
    }
    if (ss) ss>>=1; else t--;
  }
  for (j=0;j<LL;j++) ran_x1[j+KK-LL]=x[j];
  for (;j<KK;j++) ran_x1[j-LL]=x[j];
}

/* the following routines are from exercise 3.6--15 */
/* after calling ran_start1, get new randoms by, e.g., "x=ran_arr_next()" */

#define QUALITY 1009 /* recommended quality level for high-res use */
long ran_arr_buf1[QUALITY];
long ran_arr_sentinel1=-1;
long *ran_arr_ptr1=&ran_arr_sentinel1; /* the next random number, or -1 */

#define ran_arr_next1() (*ran_arr_ptr1>=0? *ran_arr_ptr1++: ran_arr_cycle1())
long ran_arr_cycle1()
{
  ran_array1(ran_arr_buf1,QUALITY);
  ran_arr_buf1[100]=-1;
  ran_arr_ptr1=ran_arr_buf1+1;
  return ran_arr_buf1[0];
}


/*-----------------------------------------------------------------------*/

static unsigned long Ran_array1_Bits (void *junk1, void *junk2)
{
   return ran_arr_next1 () << 2;
}

/*-----------------------------------------------------------------------*/

static double Ran_array1_U01 (void *vpar, void *vsta)
{
   return Ran_array1_Bits (vpar, vsta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/

static void WrRan_array1 (void *junk)
{
   int j;
   if (unif01_WrLongStateFlag) {
      printf ("ran_x1 = {\n ");
      for (j = 0; j < KK; j++) {
         printf ("%12ld", ran_x1[j]);
         if (j < KK - 1)
            printf (", ");
         if ((j % 5) == 4)
            printf ("\n ");
      };
      printf ("   }\n");
   } else
      unif01_WrLongStateDef ();
}

/*-----------------------------------------------------------------------*/

unif01_Gen *uknuth_CreateRan_array1 (long s, long A[])
{
   unif01_Gen *gen;
   size_t leng;
   char name[LEN + 1];
   int j;

   util_Assert (s <= 1073741821,
      "uknuth_CreateRan_array1:   s must be <= 1073741821");
   util_Assert (co1 == 0,
      "uknuth_CreateRan_array1:\n   only 1 such generator can be in use at a time");
   co1++;

   gen = util_Malloc (sizeof (unif01_Gen));
   strcpy (name, "uknuth_CreateRan_array1:");

   if (s < 0) {
      /* Restart with the last state A[] obtained from a previous run */
      addstr_ArrayLong (name, "   A = ", KK, A);
      for (j = 0; j < KK; j++)
         ran_x1[j] = A[j];
      *ran_arr_ptr1 = ran_arr_sentinel1;
   } else {
      /* initialize by Knuth ran_start */
      addstr_Long (name, "   s = ", s);
      ran_start1 (s);
   }
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->GetBits = &Ran_array1_Bits;
   gen->GetU01 = &Ran_array1_U01;
   gen->Write = &WrRan_array1;
   gen->param = NULL;
   gen->state = NULL;
   return gen;
}


#undef is_odd

/*=========================================================================*/

/*=========================  WARNING:

I HAVE CHANGED Knuth's code for the following version of rng-double.c. Thus
it is NOT Knuth's original code. I have renamed the variables because those
of the new version have the same name as those of the old version and we
include both versions in our file. (R. Simard)

===========================*/


/***************************** Knuth's code ********************************/

/*    This program is copyright (c) 2000 by D E Knuth;
 *    you may copy it freely  AS LONG AS YOU MAKE ABSOLUTELY NO CHANGES!
 *    You could also change it, but then you must rename the file and
 *    tell people clearly that your version is not the same as mine.
 *    It is explained in Seminumerical Algorithms, 3rd edition, Section 3.6
 *    (or in the errata to the 2nd edition --- see
 *        http://www-cs-faculty.stanford.edu/~knuth/taocp.html
 *    in the changes to pages 171 and following of Volume 2).              */

/*    If you find any bugs, please report them immediately to
 *                 taocp@cs.stanford.edu
 *    (and you will be rewarded if the bug is genuine). Thanks!            */

/************ see the book for explanations and caveats! *******************/
/************ in particular, you need two's complement arithmetic **********/

/* old-style C function declarations appear here for reasons of portability */

#define KK 100                     /* the long lag */
#define LL  37                     /* the short lag */
#define mod_sum(x,y) (((x)+(y))-(int)((x)+(y)))   /* (x+y) mod 1.0 */

double ran_u1[KK];           /* the generator state */

/* void ranf_array1(double aa[], int n) */
void ranf_array1(aa,n)    /* put n new random fractions in aa */
  double *aa;   /* destination */
  int n;      /* array length (must be at least KK) */
{
  register int i,j;
  for (j=0;j<KK;j++) aa[j]=ran_u1[j];
  for (;j<n;j++) aa[j]=mod_sum(aa[j-KK],aa[j-LL]);
  for (i=0;i<LL;i++,j++) ran_u1[i]=mod_sum(aa[j-KK],aa[j-LL]);
  for (;i<KK;i++,j++) ran_u1[i]=mod_sum(aa[j-KK],ran_u1[i-LL]);
}

#define TT  70   /* guaranteed separation between streams */
#define is_odd(s) ((s)&1)

/* void ranf_start1(long seed) */
void ranf_start1(seed)    /* do this before using ranf_array1 */
  long seed;            /* selector for different streams */
{
  register int t,s,j;
  double u[KK+KK-1],ul[KK+KK-1];
  double ulp=(1.0/(1L<<30))/(1L<<22);               /* 2 to the -52 */
  double ss=2.0*ulp*((seed&0x3fffffff)+2);

  for (j=0;j<KK;j++) {
    u[j]=ss; ul[j]=0.0;                     /* bootstrap the buffer */
    ss+=ss; if (ss>=1.0) ss-=1.0-2*ulp;  /* cyclic shift of 51 bits */
  }
  for (;j<KK+KK-1;j++) u[j]=ul[j]=0.0;
  u[1]+=ulp;ul[1]=ulp;           /* make u[1] (and only u[1]) "odd" */
  s=seed&0x3fffffff;
  t=TT-1; while (t) {
    for (j=KK-1;j>0;j--) ul[j+j]=ul[j],u[j+j]=u[j];     /* "square" */
    for (j=KK+KK-2;j>KK-LL;j-=2)
        ul[KK+KK-1-j]=0.0,u[KK+KK-1-j]=u[j]-ul[j];
    for (j=KK+KK-2;j>=KK;j--) if(ul[j]) {
      ul[j-(KK-LL)]=ulp-ul[j-(KK-LL)],
        u[j-(KK-LL)]=mod_sum(u[j-(KK-LL)],u[j]);
      ul[j-KK]=ulp-ul[j-KK],u[j-KK]=mod_sum(u[j-KK],u[j]);
    }
    if (is_odd(s)) {                             /* "multiply by z" */
      for (j=KK;j>0;j--)  ul[j]=ul[j-1],u[j]=u[j-1];
      ul[0]=ul[KK],u[0]=u[KK];       /* shift the buffer cyclically */
      if (ul[KK]) ul[LL]=ulp-ul[LL],u[LL]=mod_sum(u[LL],u[KK]);
    }
    if (s) s>>=1; else t--;
  }
  for (j=0;j<LL;j++) ran_u1[j+KK-LL]=u[j];
  for (;j<KK;j++) ran_u1[j-LL]=u[j];
}

/* the following routines are adapted from exercise 3.6--15 */
/* after calling ranf_start1, get new randoms by, e.g., "x=ranf_arr_next()" */

#define QUALITY 1009 /* recommended quality level for high-res use */
double ranf_arr_buf1[QUALITY];
double ranf_arr_sentinel1=-1.0;
double *ranf_arr_ptr1=&ranf_arr_sentinel1; /* the next random fraction, or -1 */

#define ranf_arr_next1() (*ranf_arr_ptr1>=0? *ranf_arr_ptr1++: ranf_arr_cycle1())
double ranf_arr_cycle1()
{
  ranf_array1(ranf_arr_buf1,QUALITY);
  ranf_arr_buf1[100]=-1;
  ranf_arr_ptr1=ranf_arr_buf1+1;
  return ranf_arr_buf1[0];
}


/* ========================= end of Knuth's code ========================= */



/*------------------------------- Our code --------------------------------*/

static double Ranf_array1_U01 (void *junk1, void *junk2)
{
   return ranf_arr_next1 ();
}

/*-----------------------------------------------------------------------*/

static unsigned long Ranf_array1_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (Ranf_array1_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/

static void WrRanf_array1 (void *junk)
{
   int j;
   if (unif01_WrLongStateFlag) {
      printf ("ran_u1 = {\n");
      for (j = 0; j < KK; j++) {
         printf (" %22.16f", ran_u1[j]);
         if (j < KK - 1)
            printf (",");
         if ((j % 3) == 2)
            printf ("\n");
      };
      printf ("\n     }");
   } else
      unif01_WrLongStateDef ();
}

/*-----------------------------------------------------------------------*/

unif01_Gen * uknuth_CreateRanf_array1 (long s, double B[])
{
   unif01_Gen *gen;
   size_t leng;
   char name[LEN + 1];
   int j;

   util_Assert (s <= 1073741821,
      "uknuth_CreateRanf_array1:   s must be <= 1073741821");
   util_Assert (co2 == 0,
      "uknuth_CreateRanf_array1:\n   only 1 such generator can be in use at a time");
   co2++;

   gen = util_Malloc (sizeof (unif01_Gen));
   strcpy (name, "uknuth_CreateRanf_array1:");

   if (s < 0) {
      /* Restart with the last state B[] obtained from a previous run */
      addstr_ArrayDouble (name, "   A = ", KK, B);
      for (j = 0; j < KK; j++)
         ran_u1[j] = B[j];
      *ranf_arr_ptr1 = ranf_arr_sentinel1;
   } else {
      /* initialize by Knuth ranf_start */
      addstr_Long (name, "   s = ", s);
      ranf_start1 (s);
   }
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->GetBits = &Ranf_array1_Bits;
   gen->GetU01 = &Ranf_array1_U01;
   gen->Write = &WrRanf_array1;
   gen->param = NULL;
   gen->state = NULL;
   return gen;
}

#undef is_odd



/*=========================================================================*/

/*=========================  WARNING:
 I have not made any change in the following version of Knuth's rng.c
 (R. Simard)
===========================*/



/***************************** Knuth's code ********************************/

/*    This program by D E Knuth is in the public domain and freely copyable
 *    AS LONG AS YOU MAKE ABSOLUTELY NO CHANGES!
 *    It is explained in Seminumerical Algorithms, 3rd edition, Section 3.6
 *    (or in the errata to the 2nd edition --- see
 *        http://www-cs-faculty.stanford.edu/~knuth/taocp.html
 *    in the changes to Volume 2 on pages 171 and following).              */

/*    N.B. The MODIFICATIONS introduced in the 9th printing (2002) are
      included here; there's no backwards compatibility with the original. */

/*    If you find any bugs, please report them immediately to
 *                 taocp@cs.stanford.edu
 *    (and you will be rewarded if the bug is genuine). Thanks!            */

/************ see the book for explanations and caveats! *******************/
/************ in particular, you need two's complement arithmetic **********/

#define KK 100                     /* the long lag */
#define LL  37                     /* the short lag */
#define MM (1L<<30)                 /* the modulus */
#define mod_diff(x,y) (((x)-(y))&(MM-1)) /* subtraction mod MM */

long ran_x[KK];                    /* the generator state */

#ifdef __STDC__
void ran_array(long aa[],int n)
#else
void ran_array(aa,n)    /* put n new random numbers in aa */
  long *aa;   /* destination */
  int n;      /* array length (must be at least KK) */
#endif
{
  register int i,j;
  for (j=0;j<KK;j++) aa[j]=ran_x[j];
  for (;j<n;j++) aa[j]=mod_diff(aa[j-KK],aa[j-LL]);
  for (i=0;i<LL;i++,j++) ran_x[i]=mod_diff(aa[j-KK],aa[j-LL]);
  for (;i<KK;i++,j++) ran_x[i]=mod_diff(aa[j-KK],ran_x[i-LL]);
}

/* the following routines are from exercise 3.6--15 */
/* after calling ran_start, get new randoms by, e.g., "x=ran_arr_next()" */

#define QUALITY 1009 /* recommended quality level for high-res use */
long ran_arr_buf[QUALITY];
long ran_arr_sentinel=-1;
long *ran_arr_ptr=&ran_arr_sentinel; /* the next random number, or -1 */

#define ran_arr_next() (*ran_arr_ptr>=0? *ran_arr_ptr++: ran_arr_cycle())
long ran_arr_cycle()
{
  ran_array(ran_arr_buf,QUALITY);
  ran_arr_buf[100]=-1;
  ran_arr_ptr=ran_arr_buf+1;
  return ran_arr_buf[0];
}

#define TT  70   /* guaranteed separation between streams */
#define is_odd(x)  ((x)&1)          /* units bit of x */

#ifdef __STDC__
void ran_start(long seed)
#else
void ran_start(seed)    /* do this before using ran_array */
  long seed;            /* selector for different streams */
#endif
{
  register int t,j;
  long x[KK+KK-1];              /* the preparation buffer */
  register long ss=(seed+2)&(MM-2);
  for (j=0;j<KK;j++) {
    x[j]=ss;                      /* bootstrap the buffer */
    ss<<=1; if (ss>=MM) ss-=MM-2; /* cyclic shift 29 bits */
  }
  x[1]++;              /* make x[1] (and only x[1]) odd */
  for (ss=seed&(MM-1),t=TT-1; t; ) {       
    for (j=KK-1;j>0;j--) x[j+j]=x[j], x[j+j-1]=0; /* "square" */
    for (j=KK+KK-2;j>=KK;j--)
      x[j-(KK-LL)]=mod_diff(x[j-(KK-LL)],x[j]),
      x[j-KK]=mod_diff(x[j-KK],x[j]);
    if (is_odd(ss)) {              /* "multiply by z" */
      for (j=KK;j>0;j--)  x[j]=x[j-1];
      x[0]=x[KK];            /* shift the buffer cyclically */
      x[LL]=mod_diff(x[LL],x[KK]);
    }
    if (ss) ss>>=1; else t--;
  }
  for (j=0;j<LL;j++) ran_x[j+KK-LL]=x[j];
  for (;j<KK;j++) ran_x[j-LL]=x[j];
  for (j=0;j<10;j++) ran_array(x,KK+KK-1); /* warm things up */
  ran_arr_ptr=&ran_arr_sentinel;
}


/* ========================= end of Knuth's code ========================= */


/*----------------- Our code calling Knuth's generator --------------------*/

static unsigned long Ran_array2_Bits (void *junk1, void *junk2)
{
   return ran_arr_next () << 2;
}

/*-----------------------------------------------------------------------*/

static double Ran_array2_U01 (void *vpar, void *vsta)
{
   return Ran_array2_Bits (vpar, vsta) * unif01_INV32;
}

/*-----------------------------------------------------------------------*/

static void WrRan_array2 (void *junk)
{
   int j;
   if (unif01_WrLongStateFlag) {
      printf ("ran_x = {\n ");
      for (j = 0; j < KK; j++) {
         printf ("%12ld", ran_x[j]);
         if (j < KK - 1)
            printf (", ");
         if ((j % 5) == 4)
            printf ("\n ");
      };
      printf ("   }\n");
   } else
      unif01_WrLongStateDef ();
}

/*-----------------------------------------------------------------------*/

unif01_Gen *uknuth_CreateRan_array2 (long s, long A[])
{
   unif01_Gen *gen;
   size_t leng;
   char name[LEN + 1];
   int j;

   util_Assert (s <= 1073741821,
      "uknuth_CreateRan_array2:   s must be <= 1073741821");
   util_Assert (co3 == 0,
      "uknuth_CreateRan_array2:\n   only 1 such generator can be in use at a time");
   co3++;

   gen = util_Malloc (sizeof (unif01_Gen));
   strcpy (name, "uknuth_CreateRan_array2:");

   if (s < 0) {
      /* Restart with the last state A[] obtained from a previous run */
      addstr_ArrayLong (name, "   A = ", KK, A);
      for (j = 0; j < KK; j++)
         ran_x[j] = A[j];
      *ran_arr_ptr = ran_arr_sentinel;
   } else {
      /* initialize by Knuth ran_start */
      addstr_Long (name, "   s = ", s);
      ran_start (s);
   }
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->GetBits = &Ran_array2_Bits;
   gen->GetU01 = &Ran_array2_U01;
   gen->Write = &WrRan_array2;
   gen->param = NULL;
   gen->state = NULL;
   return gen;
}


#undef is_odd


/*****************************************************************************/

/*=========================  WARNING:
 I have not made any change in the following version of Knuth's rng-double.c
 (R. Simard)
===========================*/


/***************************** Knuth's code ********************************/

/*    This program by D E Knuth is in the public domain and freely copyable
 *    AS LONG AS YOU MAKE ABSOLUTELY NO CHANGES!
 *    It is explained in Seminumerical Algorithms, 3rd edition, Section 3.6
 *    (or in the errata to the 2nd edition --- see
 *        http://www-cs-faculty.stanford.edu/~knuth/taocp.html
 *    in the changes to Volume 2 on pages 171 and following).              */

/*    N.B. The MODIFICATIONS introduced in the 9th printing (2002) are
      included here; there's no backwards compatibility with the original. */

/*    If you find any bugs, please report them immediately to
 *                 taocp@cs.stanford.edu
 *    (and you will be rewarded if the bug is genuine). Thanks!            */

/************ see the book for explanations and caveats! *******************/
/************ in particular, you need two's complement arithmetic **********/

#define KK 100                     /* the long lag */
#define LL  37                     /* the short lag */
#define mod_sum(x,y) (((x)+(y))-(int)((x)+(y)))   /* (x+y) mod 1.0 */

double ran_u[KK];           /* the generator state */

#ifdef __STDC__
void ranf_array(double aa[], int n)
#else
void ranf_array(aa,n)    /* put n new random fractions in aa */
  double *aa;   /* destination */
  int n;      /* array length (must be at least KK) */
#endif
{
  register int i,j;
  for (j=0;j<KK;j++) aa[j]=ran_u[j];
  for (;j<n;j++) aa[j]=mod_sum(aa[j-KK],aa[j-LL]);
  for (i=0;i<LL;i++,j++) ran_u[i]=mod_sum(aa[j-KK],aa[j-LL]);
  for (;i<KK;i++,j++) ran_u[i]=mod_sum(aa[j-KK],ran_u[i-LL]);
}

/* the following routines are adapted from exercise 3.6--15 */
/* after calling ranf_start, get new randoms by, e.g., "x=ranf_arr_next()" */

#define QUALITY 1009 /* recommended quality level for high-res use */
double ranf_arr_buf[QUALITY];
double ranf_arr_sentinel=-1.0;
double *ranf_arr_ptr=&ranf_arr_sentinel; /* the next random fraction, or -1 */

#define ranf_arr_next() (*ranf_arr_ptr>=0? *ranf_arr_ptr++: ranf_arr_cycle())
double ranf_arr_cycle()
{
  ranf_array(ranf_arr_buf,QUALITY);
  ranf_arr_buf[100]=-1;
  ranf_arr_ptr=ranf_arr_buf+1;
  return ranf_arr_buf[0];
}

#define TT  70   /* guaranteed separation between streams */
#define is_odd(s) ((s)&1)

#ifdef __STDC__
void ranf_start(long seed)
#else
void ranf_start(seed)    /* do this before using ranf_array */
  long seed;            /* selector for different streams */
#endif
{
  register int t,s,j;
  double u[KK+KK-1];
  double ulp=(1.0/(1L<<30))/(1L<<22);               /* 2 to the -52 */
  double ss=2.0*ulp*((seed&0x3fffffff)+2);

  for (j=0;j<KK;j++) {
    u[j]=ss;                                /* bootstrap the buffer */
    ss+=ss; if (ss>=1.0) ss-=1.0-2*ulp;  /* cyclic shift of 51 bits */
  }
  u[1]+=ulp;                     /* make u[1] (and only u[1]) "odd" */
  for (s=seed&0x3fffffff,t=TT-1; t; ) {
    for (j=KK-1;j>0;j--)
      u[j+j]=u[j],u[j+j-1]=0.0;                         /* "square" */
    for (j=KK+KK-2;j>=KK;j--) {
      u[j-(KK-LL)]=mod_sum(u[j-(KK-LL)],u[j]);
      u[j-KK]=mod_sum(u[j-KK],u[j]);
    }
    if (is_odd(s)) {                             /* "multiply by z" */
      for (j=KK;j>0;j--) u[j]=u[j-1];
      u[0]=u[KK];                    /* shift the buffer cyclically */
      u[LL]=mod_sum(u[LL],u[KK]);
    }
    if (s) s>>=1; else t--;
  }
  for (j=0;j<LL;j++) ran_u[j+KK-LL]=u[j];
  for (;j<KK;j++) ran_u[j-LL]=u[j];
  for (j=0;j<10;j++) ranf_array(u,KK+KK-1);  /* warm things up */
  ranf_arr_ptr=&ranf_arr_sentinel;
}


/* ========================= end of Knuth's code ========================= */



/*------------------------------- Our code --------------------------------*/

static double Ranf_array2_U01 (void *junk1, void *junk2)
{
   return ranf_arr_next ();
}

/*-----------------------------------------------------------------------*/

static unsigned long Ranf_array2_Bits (void *vpar, void *vsta)
{
   return (unsigned long) (Ranf_array2_U01 (vpar, vsta) * unif01_NORM32);
}

/*-----------------------------------------------------------------------*/

static void WrRanf_array2 (void *junk)
{
   int j;
   if (unif01_WrLongStateFlag) {
      printf ("ran_u = {\n");
      for (j = 0; j < KK; j++) {
         printf (" %22.16f", ran_u[j]);
         if (j < KK - 1)
            printf (",");
         if ((j % 3) == 2)
            printf ("\n");
      };
      printf ("\n     }");
   } else
      unif01_WrLongStateDef ();
}

/*-----------------------------------------------------------------------*/

unif01_Gen * uknuth_CreateRanf_array2 (long s, double B[])
{
   unif01_Gen *gen;
   size_t leng;
   char name[LEN + 1];
   int j;

   util_Assert (s <= 1073741821,
      "uknuth_CreateRanf_array2:   s must be <= 1073741821");
   util_Assert (co4 == 0,
      "uknuth_CreateRanf_array2:\n   only 1 such generator can be in use at a time");
   co4++;

   gen = util_Malloc (sizeof (unif01_Gen));
   strcpy (name, "uknuth_CreateRanf_array2:");

   if (s < 0) {
      /* Restart with the last state B[] obtained from a previous run */
      addstr_ArrayDouble (name, "   A = ", KK, B);
      for (j = 0; j < KK; j++)
         ran_u[j] = B[j];
      *ranf_arr_ptr = ranf_arr_sentinel;
   } else {
      /* initialize by Knuth ranf_start */
      addstr_Long (name, "   s = ", s);
      ranf_start (s);
   }
   leng = strlen (name);
   gen->name = util_Calloc (leng + 1, sizeof (char));
   strncpy (gen->name, name, leng);

   gen->GetBits = &Ranf_array2_Bits;
   gen->GetU01 = &Ranf_array2_U01;
   gen->Write = &WrRanf_array2;
   gen->param = NULL;
   gen->state = NULL;
   return gen;
}



/*****************************************************************************/

void uknuth_DeleteRan_array1 (unif01_Gen * gen)
{
   if (NULL == gen || co1 == 0)
      return;
   gen->name = util_Free (gen->name);
   util_Free (gen);
   co1--;
}


/*****************************************************************************/

void uknuth_DeleteRanf_array1 (unif01_Gen * gen)
{
   if (NULL == gen || co2 == 0)
      return;
   gen->name = util_Free (gen->name);
   util_Free (gen);
   co2--;
}


/*****************************************************************************/

void uknuth_DeleteRan_array2 (unif01_Gen * gen)
{
   if (NULL == gen || co3 == 0)
      return;
   gen->name = util_Free (gen->name);
   util_Free (gen);
   co3--;
}


/*****************************************************************************/

void uknuth_DeleteRanf_array2 (unif01_Gen * gen)
{
   if (NULL == gen || co4 == 0)
      return;
   gen->name = util_Free (gen->name);
   util_Free (gen);
   co4--;
}
