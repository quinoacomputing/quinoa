/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include "structs.h"
#include "smalloc.h"

/* Allocate space for new orthlink, double version. */
struct orthlink *
makeorthlnk (void)
{
    struct orthlink *newlnk;

    newlnk = smalloc(sizeof(struct orthlink));
    return (newlnk);
}

/* Allocate space for new orthlink, float version. */
struct orthlink_float *
makeorthlnk_float (void)
{
    struct orthlink_float *newlnk;

    newlnk = smalloc(sizeof(struct orthlink_float));
    return (newlnk);
}
