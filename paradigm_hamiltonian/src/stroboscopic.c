#include "allocs.h"
#include "integrator.h"
#include "math_helpers.h"
#include "model.h"
#include <math.h>
#include <math_helpers.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

void generate_stroboscopic_map(double *q, double *p, const double *parameters,
                               size_t num_points, double time_step,
                               double *section_points, gradV_func_t dVdq,
                               gradT_func_t dTdp) {
    double *dV = xmalloc(DOF, sizeof *dV);
    double *dT = xmalloc(DOF, sizeof *dT);
    double *q_prev = xmalloc(DOF, sizeof *q_prev);
    double *p_prev = xmalloc(DOF, sizeof *p_prev);
    double *q_sample = xmalloc(DOF, sizeof *q_sample);
    double *p_sample = xmalloc(DOF, sizeof *p_sample);

    for (size_t i = 0; i < DOF; i++) {
        q_prev[i] = q[i];
        p_prev[i] = p[i];
    }

    double k = parameters[2];
    double forcing_period = 2.0 * M_PI / k;

    size_t count = 0;
    size_t n_steps = 0;

    double t_prev = 0.0;
    double t_curr = 0.0;
    double next_sample_time = forcing_period;

    while (count < num_points) {
        yoshida4_step(q, p, parameters, time_step, dV, dT, dVdq, dTdp);

        t_curr = (n_steps + 1) * time_step;

        while (next_sample_time <= t_curr && count < num_points) {
            double lam = (next_sample_time - t_prev) / (t_curr - t_prev);

            for (size_t i = 0; i < DOF; i++) {
                q_sample[i] = (1.0 - lam) * q_prev[i] + lam * q[i];
                p_sample[i] = (1.0 - lam) * p_prev[i] + lam * p[i];
            }

            section_points[3 * count + 0] = next_sample_time;
            section_points[3 * count + 1] = wrap_to_pi(q_sample[0]);
            section_points[3 * count + 2] = p_sample[0];

            count++;
            next_sample_time += forcing_period;
        }

        for (size_t i = 0; i < DOF; i++) {
            q_prev[i] = q[i];
            p_prev[i] = p[i];
        }

        t_prev = t_curr;
        n_steps++;
    }

    free(dV);
    free(dT);
    free(q_prev);
    free(p_prev);
    free(q_sample);
    free(p_sample);
}
