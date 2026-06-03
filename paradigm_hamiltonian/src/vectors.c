#include "vectors.h"
#include <math.h>
#include <stdlib.h>

double standard_deviation(size_t N, size_t dim, const double *x, double *mean,
                          double *std) {

    // Initialize sums
    for (size_t d = 0; d < dim; d++) {
        mean[d] = 0.0;
        std[d] = 0.0;
    }

    // Compute means
    for (size_t i = 0; i < N; i++) {
        for (size_t d = 0; d < dim; d++) {
            mean[d] += x[i * dim + d];
        }
    }
    for (size_t d = 0; d < dim; d++)
        mean[d] /= N;

    // Compute variance
    for (size_t i = 0; i < N; i++) {
        for (size_t d = 0; d < dim; d++) {
            double diff = x[i * dim + d] - mean[d];
            std[d] += diff * diff;
        }
    }
    for (size_t d = 0; d < dim; d++)
        std[d] = sqrt(std[d] / N);

    double std_max = std[0];
    for (size_t d = 1; d < dim; d++)
        std_max = fmax(std_max, std[d]);

    return std_max;
}

double normalize_vector(double *v, size_t n) {
    double norm = 0.0;
    for (size_t i = 0; i < n; i++)
        norm += v[i] * v[i];
    norm = sqrt(norm);
    for (size_t i = 0; i < n; i++)
        v[i] /= norm;
    return norm;
}

void normalize_columns(double *v, size_t n, size_t k) {

    for (size_t col = 0; col < k; col++) {

        double norm = 0.0;

        /* compute column norm */

        for (size_t row = 0; row < n; row++) {

            double x = v[row * k + col];

            norm += x * x;
        }

        norm = sqrt(norm);

        /* avoid division by zero */

        if (norm == 0.0)

            continue;

        /* normalize column */

        for (size_t row = 0; row < n; row++)

            v[row * k + col] /= norm;
    }
}

void normalize_rows(double *X, size_t n_rows, size_t n_cols) {

    for (size_t row = 0; row < n_rows; row++) {

        double *r = X + row * n_cols;

        double norm = 0.0;

        for (size_t j = 0; j < n_cols; j++)

            norm += r[j] * r[j];

        norm = sqrt(norm);

        if (norm == 0.0)

            continue;

        for (size_t j = 0; j < n_cols; j++)

            r[j] /= norm;
    }
}

double sum_vector(double *v, size_t n) {
    double sum = 0.0;
    for (size_t i = 0; i < n; i++)
        sum += v[i];
    return sum;
}
