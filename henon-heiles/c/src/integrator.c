#include "models.h"
#include "types.h"
#include <stddef.h>
#include <stdlib.h>

static const double ALPHA =
    1.0 / (2.0 - 1.25992104989487); // ~0.6756035959798289
static const double BETA =
    -1.25992104989487 / (2.0 - 1.25992104989487); // ~-0.3512071919596578

static void symplectic_step(double *q, double *p, const double *parameters,
                            const double time_step, double *dV, double *dT,
                            gradV_func_t dVdq, gradT_func_t dTdp) {

    // Half kick
    dVdq(dV, q, parameters);
    for (size_t i = 0; i < DOF; i++)
        p[i] -= 0.5 * time_step * dV[i];

    // Drift
    dTdp(dT, p, parameters);
    for (size_t i = 0; i < DOF; i++)
        q[i] += time_step * dT[i];

    // Half kick
    dVdq(dV, q, parameters);
    for (size_t i = 0; i < DOF; i++)
        p[i] -= 0.5 * time_step * dV[i];
}

void yoshida4_step(double *q, double *p, const double *parameters,
                   const double time_step, double *dV, double *dT,
                   gradV_func_t dVdq, gradT_func_t dTdp) {
    symplectic_step(q, p, parameters, ALPHA * time_step, dV, dT, dVdq, dTdp);
    symplectic_step(q, p, parameters, BETA * time_step, dV, dT, dVdq, dTdp);
    symplectic_step(q, p, parameters, ALPHA * time_step, dV, dT, dVdq, dTdp);
}
