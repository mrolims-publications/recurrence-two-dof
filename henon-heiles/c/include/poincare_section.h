#ifndef PSS_H
#define PSS_H

#include "integrator.h"
#include "models.h"
#include "types.h"

void generate_poincare_section(double *q, double *p, const double *parameters,
                               size_t num_points, double time_step,
                               size_t section_index, double section_value,
                               int crossing, double *section_points,
                               gradV_func_t dVdq, gradT_func_t dTdp);

void generate_stroboscopic_map(double *q, double *p, const double *parameters,
                               size_t num_points, double time_step,
                               double *section_points, gradV_func_t dVdq,
                               gradT_func_t dTdp);

#endif
