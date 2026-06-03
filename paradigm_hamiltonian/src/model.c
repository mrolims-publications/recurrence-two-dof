#include "model.h"
#include "types.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ----------------------------------
 *
 * Paradigm Hamiltonian
 *
 * ----------------------------------- */

double PH_T(const double *p, const double *parameters) {
    (void)parameters;
    double p1 = p[0];
    double p2 = p[1];

    return 0.5 * p1 * p1 + p2;
}

void PH_grad_T(double *dT, const double *p, const double *parameters) {
    (void)parameters;
    double p1 = p[0];

    dT[0] = p1;
    dT[1] = 1.0;
}

void PH_hess_T(double *hess_T, const double *p, const double *parameters) {
    (void)p;
    (void)parameters;

    hess_T[0] = 1.0;
    hess_T[1] = 0.0;
    hess_T[2] = 0.0;
    hess_T[3] = 0.0;
}

double PH_V(const double *q, const double *parameters) {
    double q1 = q[0];
    double q2 = q[1];

    double M = parameters[0];
    double A = parameters[1];
    double k = parameters[2];

    return -M * cos(q1) - A * cos(k * (q1 - q2));
}

void PH_grad_V(double *dV, const double *q, const double *parameters) {
    double q1 = q[0];
    double q2 = q[1];

    double M = parameters[0];
    double A = parameters[1];
    double k = parameters[2];

    double s = sin(k * (q1 - q2));

    dV[0] = M * sin(q1) + k * A * s;
    dV[1] = -k * A * s;
}

void PH_hess_V(double *hess_V, const double *q, const double *parameters) {
    double q1 = q[0];
    double q2 = q[1];

    double M = parameters[0];
    double A = parameters[1];
    double k = parameters[2];

    double c = cos(k * (q1 - q2));
    double k2Ac = k * k * A * c;

    hess_V[0] = M * cos(q1) + k2Ac;
    hess_V[1] = -k2Ac;
    hess_V[2] = -k2Ac;
    hess_V[3] = k2Ac;
}

/*
 * Auxiliar function for finding valid ICs
 * */
int p_from_E(double E, const double *q, double *p, const double *parameters,
             V_func_t V) {

    double p1 = p[1];
    double arg = 2 * (E - V(q, parameters)) - p1 * p1;

    int info;
    if (arg > 0) {
        p[0] = sqrt(arg);
        info = 0;
    } else
        info = -1;

    return info;
}
