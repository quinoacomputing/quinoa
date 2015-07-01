#include "gdef.h"
#include "swrite.h"
#include "bbattery.h"

int main (void)
{
   swrite_Basic = FALSE;
   bbattery_RabbitFile ("vax.bin", 1048576);
   return 0;
}
