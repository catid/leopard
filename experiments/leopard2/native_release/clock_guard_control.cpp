// Positive controls: each wrapped clock must terminate with status 86.
#include <cstring>
#include <ctime>
#include <sys/time.h>
int main(int argc, char** argv)
{
    if (argc != 2) return 1;
    if (std::strcmp(argv[1], "clock_gettime") == 0) {
        struct timespec value;
        return clock_gettime(CLOCK_MONOTONIC, &value);
    }
    if (std::strcmp(argv[1], "gettimeofday") == 0) {
        struct timeval value;
        return gettimeofday(&value, NULL);
    }
    if (std::strcmp(argv[1], "clock") == 0) return clock() == 0 ? 0 : 2;
    return 1;
}
