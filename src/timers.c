#include "timers.h"
#include <stdio.h>

double get_clock()
{
    struct timeval tv;
    int ret = gettimeofday(&tv,0);
    if (ret < 0)
        printf("gettimeofday error");

    return (tv.tv_sec*1.0+tv.tv_usec*1.0E-6);
}
