
 
/* chrono.h for ANSI C */

#ifndef CHRONO_H
#define CHRONO_H 
#include "gdef.h"
 


typedef struct {
   unsigned long microsec;
   unsigned long second;
   } chrono_Chrono;



typedef enum {
   chrono_sec,
   chrono_min,
   chrono_hours, 
   chrono_days,
   chrono_hms
   } chrono_TimeFormat;


chrono_Chrono * chrono_Create (void);



void chrono_Delete (chrono_Chrono * C);



void chrono_Init (chrono_Chrono * C);



double chrono_Val (chrono_Chrono * C, chrono_TimeFormat Unit);



void chrono_Write (chrono_Chrono * C, chrono_TimeFormat Unit);

 
#endif
 
