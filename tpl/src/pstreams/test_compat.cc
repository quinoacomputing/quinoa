#define GCC_BACK_COMPAT 1
#define REDI_PSTREAMS_POPEN_USES_BIDIRECTIONAL_PIPE 1

#include "pstream_compat.h"
#include <unistd.h>

int main()
{
    {
        char c;
        redi::ipstream who("whoami");
        if (!(who >> c))
            return 1;

        redi::opstream cat("cat");
        if (!(cat << c))
            return 2;

        while (who >> c)
            cat << c;

        cat << std::endl;
    }
    // check for zombies while this process is sleeping
    sleep(10);
    return 0;
}

