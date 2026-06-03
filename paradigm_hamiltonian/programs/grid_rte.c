#include "allocs.h"
#include "mkdir.h"
#include "model.h"
#include "rqa.h"
#include "stroboscopic.h"
#include "vectors.h"
#include <math.h>
#include <omp.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static void print_progress(size_t done, size_t total, int width) {
    double frac = (double)done / (double)total;
    int filled = (int)(width * frac);

    printf("\r    [");
    for (int i = 0; i < width; i++) {
        putchar(i < filled ? '#' : '-');
    }
    double percentage = 100 * frac;
    printf("] %zu/%zu (%.2f %%)", done, total, percentage);
    fflush(stdout);
}

int main(int argc, char *argv[]) {
    (void)argc;
    (void)argv;

    const double M = 1.0;
    const double k = 1.0;
    const double A_vals[] = {0.092, 0.132, 0.248};
    const size_t num_A = sizeof(A_vals) / sizeof(A_vals[0]);

    const size_t n_cross = 10000;
    const size_t grid_size = 10;
    const size_t total_points = grid_size * grid_size;
    const double std_scale = 0.05;

    const size_t dim_pss = 3;
    const double q_range[2] = {-M_PI, M_PI};
    const double p_range[2] = {-M_PI, M_PI};

    const char *home = getenv("HOME");
    if (!home) {
        fprintf(stderr, "Error: HOME environment variable is not set.\n");
        return EXIT_FAILURE;
    }

    char path[512];
    snprintf(path, sizeof(path),
             "%s/Research/recurrence-two-dof/paradigm_hamiltonian/data", home);

    if (!dir_exists(path)) {
        printf("Output directory does not exist. Creating it...\n");
        if (mkdir_p(path) != 0) {
            perror("mkdir_p");
            return EXIT_FAILURE;
        }
        printf("Created directory: %s\n", path);
    } else {
        printf("Output directory exists: %s\n", path);
    }

    printf("\n=== GLOBAL PARAMETERS ===\n");
    printf("Model            : Paradigm Hamiltonian\n");
    printf("Observable       : Grid RTE\n");
    printf("n_cross          : %zu\n", n_cross);
    printf("grid_size        : %zu\n", grid_size);
    printf("Total points     : %zu\n", total_points);
    printf("std_scale        : %.5f\n", std_scale);
    printf("q range          : (%.6f, %.6f)\n", q_range[0], q_range[1]);
    printf("p range          : (%.6f, %.6f)\n", p_range[0], p_range[1]);
    printf("OpenMP threads   : %d\n", omp_get_max_threads());
    printf("Output dir       : %s\n", path);
    printf("=========================\n\n");

    size_t batch_valid = 0;
    size_t batch_invalid = 0;

    for (size_t na = 0; na < num_A; na++) {
        const double A = A_vals[na];
        const double parameters[] = {M, A, k};
        char filename[1024];
        snprintf(filename, sizeof(filename),
                 "%s/grid_rte_A=%.5f_N=%zu_grid_size=%zu_scale=%.3f.dat", path,
                 A, n_cross, grid_size, std_scale);

        FILE *fp = fopen(filename, "w");
        if (!fp) {
            perror("fopen");
            return EXIT_FAILURE;
        }

        double *rte = xmalloc(total_points, sizeof *rte);

        size_t done = 0;
        size_t n_valid = 0;
        size_t n_invalid = 0;

        printf("=== RUN %zu/%zu ===\n", na + 1, num_A);
        printf("A                : %.5f\n", A);
        printf("Output data      : %s\n", filename);
        printf("=========================\n\n");

        printf(">>> Starting grid scan for A = %.5f\n", A);

#pragma omp parallel default(none) shared(                                     \
        rte, done, n_valid, n_invalid, q_range, p_range, grid_size,            \
            total_points, A, parameters, n_cross, dim_pss, std_scale, stdout)
        {
            double *pss = xmalloc(n_cross * dim_pss, sizeof *pss);
            double *reduced_pss = xmalloc(n_cross * 2, sizeof *reduced_pss);
            uint8_t *recmat = xmalloc(n_cross * n_cross, sizeof *recmat);
            double *line_distr = xmalloc(n_cross, sizeof *line_distr);
            double mean[2];
            double std[2];

#pragma omp for
            for (size_t g = 0; g < total_points; g++) {
                const size_t i = g / grid_size;
                const size_t j = g % grid_size;

                const double Q = q_range[0] + (double)i *
                                                  (q_range[1] - q_range[0]) /
                                                  (double)(grid_size - 1);

                const double P = p_range[0] + (double)j *
                                                  (p_range[1] - p_range[0]) /
                                                  (double)(grid_size - 1);

                double q[DOF];
                double p[DOF];

                q[0] = Q;
                q[1] = 0.0;
                p[0] = P;
                p[1] = 0.0;

                generate_stroboscopic_map(q, p, parameters, n_cross, TIME_STEP,
                                          pss, PH_grad_V, PH_grad_T);

                for (size_t n = 0; n < n_cross; n++) {
                    reduced_pss[n * 2 + 0] = pss[n * dim_pss + 1];
                    reduced_pss[n * 2 + 1] = pss[n * dim_pss + 2];
                }

                const double std_max =
                    standard_deviation(n_cross, 2, reduced_pss, mean, std);
                const double eps = std_scale * std_max;

                rte[g] = recurrence_time_entropy(n_cross, 2, reduced_pss, eps,
                                                 1, recmat, line_distr);

#pragma omp atomic update
                n_valid++;

                size_t local_done;
#pragma omp atomic capture
                local_done = ++done;

                if (local_done % 10 == 0 || local_done == total_points) {
#pragma omp critical
                    {
                        print_progress(local_done, total_points, 40);
                    }
                }
            }

            free(pss);
            free(reduced_pss);
            free(recmat);
            free(line_distr);
        }

        printf("\n");

        for (size_t i = 0; i < grid_size; i++) {
            const double Q =
                q_range[0] + i * (q_range[1] - q_range[0]) / (grid_size - 1);

            for (size_t j = 0; j < grid_size; j++) {
                const double P = p_range[0] + j * (p_range[1] - p_range[0]) /
                                                  (grid_size - 1);

                fprintf(fp, "%.16f %.16f %.16e\n", Q, P,
                        rte[i * grid_size + j]);
            }
            fprintf(fp, "\n");
        }

        fclose(fp);
        free(rte);

        batch_valid += n_valid;
        batch_invalid += n_invalid;

        printf("\n=== SUMMARY ===\n");
        printf("A                : %.5f\n", A);
        printf("Total points     : %zu\n", total_points);
        printf("Valid ICs        : %zu\n", n_valid);
        printf("Invalid ICs      : %zu\n", n_invalid);
        printf("Saved data to    : %s\n", filename);
        printf("======================\n\n");
    }

    printf("=== BATCH SUMMARY ===\n");
    printf("Processed A values : %zu\n", num_A);
    printf("Total valid ICs    : %zu\n", batch_valid);
    printf("Total invalid ICs  : %zu\n", batch_invalid);
    printf("=====================\n\n");

    return EXIT_SUCCESS;
}
