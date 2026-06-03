#include "allocs.h"
#include "math_helpers.h"
#include "mkdir.h"
#include "model.h"
#include "parse_args.h"
#include "stroboscopic.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

static void print_progress(size_t done, size_t total, int width) {
    double frac = (double)done / (double)total;
    int filled = (int)(width * frac);

    printf("\r    [");
    for (int i = 0; i < width; i++) {
        putchar(i < filled ? '#' : '-');
    }
    printf("] %zu/%zu", done, total);
    fflush(stdout);
}

int main(int argc, char *argv[]) {
    (void)argc;
    (void)argv;

    const double M = 1.0;
    const double k = 1.0;
    const double A = get_double(argc, argv, 1);
    const size_t n_cross = 10000;
    const size_t num_ic = 200;
    const size_t dim_pss = 3;
    const double q_range[2] = {0, 2 * M_PI};
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

    const unsigned int initial_seed = 1312;

    printf("\n=== GLOBAL PARAMETERS ===\n");
    printf("Model            : Paradigm Hamiltonian\n");
    printf("Observable       : Stroboscopic map\n");
    printf("A                : %.5f\n", A);
    printf("n_cross          : %zu\n", n_cross);
    printf("num_ic           : %zu\n", num_ic);
    printf("dim_pss          : %zu\n", dim_pss);
    printf("q range          : (%.6f, %.6f)\n", q_range[0], q_range[1]);
    printf("p range          : (%.6f, %.6f)\n", p_range[0], p_range[1]);
    printf("Seed             : %u\n", initial_seed);
    printf("Output dir       : %s\n", path);
    printf("=========================\n\n");

    double parameters[] = {M, A, k};

    char filename[1024];
    snprintf(filename, sizeof(filename), "%s/pss_A=%.5f_n=%zu.dat", path, A,
             n_cross);

    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("fopen");
        return EXIT_FAILURE;
    }

    unsigned int seed = initial_seed;

    printf(">>> Starting Poincare section generation for A = %.5f\n", A);

    double *pss = xmalloc(n_cross * dim_pss, sizeof *pss);

    for (size_t i = 0; i < num_ic; i++) {
        double q[DOF], p[DOF];

        double rnd = (double)rand_r(&seed) / (double)RAND_MAX;
        double Q = q_range[0] + rnd * (q_range[1] - q_range[0]);

        rnd = (double)rand_r(&seed) / (double)RAND_MAX;
        double P = p_range[0] + rnd * (p_range[1] - p_range[0]);

        q[0] = Q;
        q[1] = 0.0;
        p[0] = P;
        p[1] = 0.0;

        generate_stroboscopic_map(q, p, parameters, n_cross, TIME_STEP, pss,
                                  PH_grad_V, PH_grad_T);

        for (size_t j = 0; j < n_cross; j++)
            fprintf(fp, "%.16f %.16f\n", pss[j * dim_pss + 1],
                    pss[j * dim_pss + 2]);

        if ((i + 1) % 1 == 0 || (i + 1) == num_ic) {
            print_progress(i + 1, num_ic, 40);
        }
    }

    printf("\n");

    free(pss);
    fclose(fp);

    printf("\n=== SUMMARY ===\n");
    printf("A                : %.5f\n", A);
    printf("Saved data to    : %s\n", filename);
    printf("======================\n\n");

    return EXIT_SUCCESS;
}
