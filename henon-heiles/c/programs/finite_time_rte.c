#include "allocs.h"
#include "mkdir.h"
#include "models.h"
#include "parse_args.h"
#include "poincare_section.h"
#include "rqa.h"
#include "vectors.h"
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int ensure_directory(const char *path);
static void print_progress(size_t done, size_t total, int width, FILE *fp);

int main(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr,
                "Usage: %s <energy> <save_pss> <total_crossings> "
                "<finite_crossings>\n",
                argv[0]);
        return EXIT_FAILURE;
    }

    const double energy = get_double(argc, argv, 1);
    const size_t save_pss = get_size_t(argc, argv, 2);
    const size_t total_crossings = get_size_t(argc, argv, 3);
    const size_t finite_crossings = get_size_t(argc, argv, 4);

    if (energy < 0) {
        fprintf(stderr, "Error: energy must be positive.\n");
        return EXIT_FAILURE;
    }

    if (save_pss != 0 && save_pss != 1) {
        fprintf(stderr, "Error: save_pss must be 0 or 1.\n");
        return EXIT_FAILURE;
    }

    if (finite_crossings == 0) {
        fprintf(stderr, "Error: finite_crossings must be > 0.\n");
        return EXIT_FAILURE;
    }

    if (total_crossings == 0) {
        fprintf(stderr, "Error: total_crossings must be > 0.\n");
        return EXIT_FAILURE;
    }

    if (total_crossings % finite_crossings != 0) {
        fprintf(stderr,
                "Error: total_crossings (%zu) must be divisible by "
                "finite_crossings (%zu).\n",
                total_crossings, finite_crossings);
        return EXIT_FAILURE;
    }

    const size_t num_windows = total_crossings / finite_crossings;
    const double scale = 0.05;
    const size_t dim_pss = 2 * DOF + 1;

    const char *home = getenv("HOME");
    if (!home) {
        fprintf(stderr, "Error: HOME environment variable is not set.\n");
        return EXIT_FAILURE;
    }
    char project_root[512];
    snprintf(project_root, sizeof(project_root),
             "%s/Research/recurrence-two-dof/henon-heiles", home);

    char path_to_data[512];
    size_t required_chars =
        snprintf(path_to_data, sizeof(path_to_data), "%s/data", project_root);

    if (required_chars >= sizeof(path_to_data)) {
        fprintf(stderr, "Path truncated\n");
        return EXIT_FAILURE;
    }
    if (ensure_directory(path_to_data) != 0)
        return EXIT_FAILURE;

    char data_filename[1024];
    snprintf(data_filename, sizeof(data_filename),
             "%s/finite_time_rte_E=%.5f_N=%zu_n=%zu_save_pss=%zu.dat",
             path_to_data, energy, total_crossings, finite_crossings, save_pss);
    FILE *data_file = fopen(data_filename, "w");
    if (!data_file) {
        perror("fopen");
        return EXIT_FAILURE;
    }

    char path_to_logs[512];
    required_chars =
        snprintf(path_to_logs, sizeof(path_to_logs), "%s/logs", project_root);
    if (required_chars >= sizeof(path_to_logs)) {
        fprintf(stderr, "Path truncated\n");
        return EXIT_FAILURE;
    }
    if (ensure_directory(path_to_logs) != 0)
        return EXIT_FAILURE;

    char log_filename[1024];
    snprintf(log_filename, sizeof(log_filename),
             "%s/finite_time_rte_E=%.5f_N=%zu_n=%zu_save_pss=%zu.log",
             path_to_logs, energy, total_crossings, finite_crossings, save_pss);
    FILE *log_file = fopen(log_filename, "w");
    if (!log_file) {
        perror("fopen");
        return EXIT_FAILURE;
    }

    double q[DOF], p[DOF];
    q[0] = 0.0;
    q[1] = -0.15;
    p[0] = 0.0;
    p[1] = 0.0;

    int info = p_from_E(energy, q, p, NULL, HH_V);
    if (info != 0) {
        fprintf(stderr, "Error: invalid IC for E = %.5f\n", energy);
        fclose(data_file);
        fclose(log_file);
        return EXIT_FAILURE;
    }

    double *pss = xmalloc(finite_crossings * dim_pss, sizeof *pss);
    double *reduced_pss = xmalloc(finite_crossings * 2, sizeof *reduced_pss);
    uint8_t *recmat =
        xmalloc(finite_crossings * finite_crossings, sizeof *recmat);
    double *line_distr = xmalloc(finite_crossings, sizeof *line_distr);
    double mean[2];
    double std[2];

    fprintf(log_file, "\n=== GLOBAL PARAMETERS ===\n");
    fprintf(log_file, "Model            : Henon-Heiles\n");
    fprintf(log_file, "Observable       : Finite-time RTE\n");
    fprintf(log_file, "Energy           : %.5f\n", energy);
    fprintf(log_file, "save_pss         : %zu\n", save_pss);
    fprintf(log_file, "total_crossings  : %zu\n", total_crossings);
    fprintf(log_file, "finite_crossings : %zu\n", finite_crossings);
    fprintf(log_file, "num_windows      : %zu\n", num_windows);
    fprintf(log_file, "scale            : %.5f\n", scale);
    fprintf(log_file, "Initial IC       : q=(%.6f, %.6f), p=(%.6f, %.6f)\n",
            q[0], q[1], p[0], p[1]);
    fprintf(log_file, "Output data      : %s\n", data_filename);
    fprintf(log_file, "=========================\n\n");

    fprintf(log_file, ">>> Starting finite-time RTE calculation\n");

    for (size_t n = 0; n < num_windows; n++) {
        generate_poincare_section(q, p, NULL, finite_crossings, TIME_STEP, 0,
                                  0.0, 1, pss, HH_grad_V, HH_grad_T);

        for (size_t m = 0; m < finite_crossings; m++) {
            reduced_pss[m * 2 + 0] = pss[m * dim_pss + 2];
            reduced_pss[m * 2 + 1] = pss[m * dim_pss + 4];
        }

        double std_max =
            standard_deviation(finite_crossings, 2, reduced_pss, mean, std);
        const double eps = scale * std_max;

        const double rte = recurrence_time_entropy(
            finite_crossings, 2, reduced_pss, eps, 1, recmat, line_distr);

        if (save_pss) {
            for (size_t m = 0; m < finite_crossings; m++) {
                fprintf(data_file, "%.16f %.16f %.16f\n",
                        reduced_pss[m * 2 + 0], reduced_pss[m * 2 + 1], rte);
            }
        } else {
            fprintf(data_file, "%.16f\n", rte);
        }

        if ((n + 1) % 10 == 0 || (n + 1) == num_windows) {
            print_progress(n + 1, num_windows, 40, log_file);
        }
    }

    printf("\n");

    free(pss);
    free(reduced_pss);
    free(recmat);
    free(line_distr);

    fclose(data_file);

    fprintf(log_file, "\n=== RUN SUMMARY ===\n");
    fprintf(log_file, "Processed windows : %zu\n", num_windows);
    fprintf(log_file, "Saved data to     : %s\n", data_filename);
    fprintf(log_file, "===================\n\n");
    fclose(log_file);

    return EXIT_SUCCESS;
}

int ensure_directory(const char *path) {
    if (!dir_exists(path)) {
        printf("Directory does not exist. Creating it...\n");

        if (mkdir_p(path) != 0) {
            perror("mkdir_p");
            return -1;
        }

        printf("Created directory: %s\n", path);
    } else {
        printf("Directory exists: %s\n", path);
    }

    return 0;
}

static void print_progress(size_t done, size_t total, int width, FILE *fp) {
    double frac = (double)done / (double)total;
    int filled = (int)(width * frac);

    fprintf(fp, "    [");

    for (int i = 0; i < width; i++) {
        fputc(i < filled ? '#' : '-', fp);
    }

    fprintf(fp, "] %zu/%zu\n", done, total);
    fflush(fp);
}
