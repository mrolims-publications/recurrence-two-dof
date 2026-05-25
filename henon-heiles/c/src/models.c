#include "models.h"
#include "types.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ----------------------------------
 *
 * Hénon-Heiles Hamiltonian
 *
 * ----------------------------------- */

double HH_V(const double *q, const double *parameters) {
    (void)parameters; // explicitly mark as unused
    double q0 = q[0];
    double q1 = q[1];

    double V = (q0 * q0 + q1 * q1) / 2.0 + q0 * q0 * q1 - q1 * q1 * q1 / 3.0;

    return V;
}

void HH_grad_V(double *dV, const double *q, const double *parameters) {
    (void)parameters; // explicitly mark as unused
    double q0 = q[0];
    double q1 = q[1];
    dV[0] = q0 * (1.0 + 2.0 * q1);
    dV[1] = q1 + q0 * q0 - q1 * q1;
}

void HH_hess_V(double *hess_V, const double *q, const double *parameters) {
    (void)parameters; // explicitly mark as unused
    double q0 = q[0];
    double q1 = q[1];

    hess_V[0] = 1.0 + 2.0 * q1;
    hess_V[1] = 2.0 * q0;
    hess_V[2] = 2.0 * q0;
    hess_V[3] = 1.0 - 2.0 * q1;
}

double HH_T(const double *p, const double *parameters) {
    (void)parameters; // explicitly mark as unused
    double p0 = p[0];
    double p1 = p[1];

    return 0.5 * (p0 * p0 + p1 * p1);
}

void HH_grad_T(double *dT, const double *p, const double *parameters) {
    (void)parameters; // explicitly mark as unused
    dT[0] = p[0];
    dT[1] = p[1];
}

void HH_hess_T(double *hess_T, const double *p, const double *parameters) {
    (void)p;
    (void)parameters; // explicitly mark as unused

    hess_T[0] = 1.0;
    hess_T[1] = 0.0;
    hess_T[2] = 0.0;
    hess_T[3] = 1.0;
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
