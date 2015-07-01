 
/*  addstr.h  for ANSI C  */
#ifndef ADDSTR_H
#define ADDSTR_H
 
#include "gdef.h"


void  addstr_Int (char *to, const char *add, int n);

void  addstr_Uint (char *to, const char *add, unsigned int n);

void  addstr_Long (char *to, const char *add, long n);

void  addstr_Ulong (char *to, const char *add, unsigned long n);

void  addstr_Double (char *to, const char *add, double x);

void  addstr_Char (char *to, const char *add, char c);

void  addstr_Bool (char *to, const char *add, int b);


#ifdef USE_LONGLONG
void  addstr_LONG (char *to, const char *add, longlong n);

void  addstr_ULONG (char *to, const char *add, ulonglong n);
#endif


void  addstr_ArrayInt (char *to, const char *add, int high, int []);

void  addstr_ArrayUint (char *to, const char *add, int high,
                        unsigned int []);

void  addstr_ArrayLong (char *to, const char *add, int high, long []);

void  addstr_ArrayUlong (char *to, const char *add, int high,
                         unsigned long []);

void  addstr_ArrayDouble (char *to, const char *add, int high, double []);

 
#endif
 

