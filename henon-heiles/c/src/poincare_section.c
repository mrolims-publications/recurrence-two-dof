#include "allocs.h"
#include "integrator.h"
#include "math_helpers.h"
#include "models.h"
#include <math_helpers.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

void generate_poincare_section(double *q, double *p, const double *parameters,
                               size_t num_points, double time_step,
                               size_t section_index, double section_value,
                               int crossing, double *section_points,
                               gradV_func_t dVdq, gradT_func_t dTdp) {

    size_t dim_pss = 2 * DOF + 1;
    double *dV = xmalloc(DOF, sizeof *dV);
    double *dT = xmalloc(DOF, sizeof *dT);
    double *q_prev = xmalloc(DOF, sizeof *q_prev);
    double *p_prev = xmalloc(DOF, sizeof *p_prev);
    double *q_cross = xmalloc(DOF, sizeof *q_cross);
    double *p_cross = xmalloc(DOF, sizeof *p_cross);

    for (size_t i = 0; i < DOF; i++) {
        q_prev[i] = q[i];
        p_prev[i] = p[i];
    }

    size_t count = 0;
    size_t n_steps = 0;
    double t_cross;

    while (count < num_points) {
        yoshida4_step(q, p, parameters, time_step, dV, dT, dVdq, dTdp);

        if ((q_prev[section_index] - section_value) *
                (q[section_index] - section_value) <
            0.0) {
            double lam = (section_value - q_prev[section_index]) /
                         (q[section_index] - q_prev[section_index]);

            for (size_t i = 0; i < DOF; i++) {
                q_cross[i] = (1 - lam) * q_prev[i] + lam * q[i];
                p_cross[i] = (1 - lam) * p_prev[i] + lam * p[i];
            }
            t_cross = n_steps * time_step + lam * time_step;

            if (crossing == 0 || sgn(p_cross[section_index]) == crossing) {
                section_points[count * dim_pss + 0] = t_cross;
                for (size_t i = 0; i < DOF; i++) {
                    section_points[count * dim_pss + 1 + i] = q_cross[i];
                    section_points[count * dim_pss + 1 + i + DOF] = p_cross[i];
                }
                count++;
            }
        }

        for (size_t i = 0; i < DOF; i++) {
            q_prev[i] = q[i];
            p_prev[i] = p[i];
        }
        n_steps++;
    }

    free(dV);
    free(dT);
    free(q_prev);
    free(p_prev);
    free(q_cross);
    free(p_cross);
}
