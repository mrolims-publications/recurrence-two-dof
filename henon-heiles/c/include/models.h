#ifndef MODELS_H
#define MODELS_H

#include "types.h"

#define DOF 2
#define TIME_STEP 0.01

double HH_V(const double *q, const double *parameters);

void HH_grad_V(double *dV, const double *q, const double *parameters);

void HH_hess_V(double *hess_V, const double *q, const double *parameters);

double HH_T(const double *p, const double *parameters);

void HH_grad_T(double *dT, const double *p, const double *parameters);

void HH_hess_T(double *hess_T, const double *p, const double *parameters);

double PH_V(const double *q, const double *parameters);

void PH_grad_V(double *dV, const double *q, const double *parameters);

void PH_hess_V(double *hess_V, const double *q, const double *parameters);

double PH_T(const double *p, const double *parameters);

void PH_grad_T(double *dT, const double *p, const double *parameters);

void PH_hess_T(double *hess_T, const double *p, const double *parameters);

int p_from_E(double E, const double *q, double *p, const double *parameters,
             V_func_t V);

#endif
