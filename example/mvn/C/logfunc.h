#ifndef PM_LOG_FUNC
#define PM_LOG_FUNC

#include <math.h>
#include <stdint.h>

#define NDIM 4  // The number of dimensions of the domain of the objective function.

double getLogFunc   (
                    int32_t ,   // ndim
                    double []   // Point
                    );

#endif