#ifndef MATH_HELPERS_H
#define MATH_HELPERS_H

#include <math.h>
#include <stddef.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static inline int sgn(double x) {

    if (x > 0)
        return 1;
    if (x < 0)
        return -1;
    return 0;
}

#endif
