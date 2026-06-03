#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "types.h"
#include <stddef.h>

void yoshida4_step(double *q, double *p, const double *parameters,
                   const double time_step, double *dV, double *dT,
                   gradV_func_t dVdq, gradT_func_t dTdp);

void yoshida4_step_traj_tan(double *q, double *p, double *X,
                            const double *parameters, double time_step,
                            size_t n_dev, double *dV, double *hess_V,
                            double *dT, double *hess_T, gradV_func_t dVdq,
                            hessV_func_t hessian_V, gradT_func_t dTdp,
                            hessT_func_t hessian_T);

#endif
