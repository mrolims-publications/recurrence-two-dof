#ifndef VECTORS_H
#define VECTORS_H

#include <stddef.h>

double standard_deviation(size_t N, size_t dim, const double *x, double *mean,
                          double *std);

double normalize_vector(double *v, size_t n);

void normalize_columns(double *v, size_t n, size_t k);

void normalize_rows(double *X, size_t n_rows, size_t n_cols);

double sum_vector(double *v, size_t n);

#endif
