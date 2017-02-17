/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifndef _DR_ERR_CONST_H_
#define _DR_ERR_CONST_H_

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#define MAX_ERR_MSG 1024

/* Structure to store an error message */
struct error_message
{
  int   level;
  char *err_mesg;
  int   line_no;
  char *filename;
};

typedef struct error_message ERROR_MSG;
typedef struct error_message *ERROR_MSG_PTR;


/* Macro used in the code to add an error message */
#define Gen_Error(a,b) (error_add(a, b, __FILE__, __LINE__))

/* Function prototype for error functions */
extern
void error_add(
  int   level,
  const char *message,	/* The message to add to the error list */
  const char *filename,	/* The filename in which the error occured */
  int   line		/* The line number in filename where the error
                         * was reported */
  );

extern
void error_report(int proc);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif /* _DR_ERR_CONST_H_ */
