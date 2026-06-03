#ifndef PSS_H
#define PSS_H

#include "integrator.h"
#include "model.h"
#include "types.h"

void generate_stroboscopic_map(double *q, double *p, const double *parameters,
                               size_t num_points, double time_step,
                               double *section_points, gradV_func_t dVdq,
                               gradT_func_t dTdp);

#endif
