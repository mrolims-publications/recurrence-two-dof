#ifndef RQA_H
#define RQA_H

#include <stddef.h>
#include <stdint.h>

/**
 * Compute the recurrence matrix of an N-point trajectory in arbitrary
 * dimension.
 *
 * The recurrence matrix R is symmetric and binary, indicating which pairs
 * of points are closer than a threshold eps according to a specified norm.
 *
 * @param N          Number of points in the trajectory.
 * @param dim        Dimension of the phase space.
 * @param x          Pointer to trajectory array of size N*dim (row-major).
 * @param eps        Recurrence threshold.
 * @param R          Output binary recurrence matrix (N x N, row-major).
 *
 * Method:
 *   - For each pair of points (i, j), compute the distance according to
 * norm_type:
 *       - Max norm: max |x_i[d] - x_j[d]|
 *   - Set R[i*N+j] = R[j*N+i] = 1 if distance <= eps, else 0.
 *
 * Notes:
 *   - R is symmetric: R[i*N+j] == R[j*N+i].
 *   - Self-recurrences (i==j) are always included.
 */
void recurrence_matrix(size_t N, size_t dim, const double *x, double eps,
                       uint8_t *R);

/**
 * Compute the distribution of vertical white lines in a binary recurrence
 * matrix.
 *
 * Only counts consecutive white cells (0) **after the first black cell (1)**
 * in each column. White lines reaching the bottom of the column are ignored.
 *
 * @param N               Size of the square recurrence matrix (N x N).
 * @param recmat          Input recurrence matrix (row-major).
 * @param min_line_length Minimum length of white lines to include.
 * @param line_dist       Output array of size N: line_dist[k] counts white
 * lines of length k+1.
 *
 * Method:
 *   - For each column:
 *       - Track the length of consecutive white cells after the first black.
 *       - If a white line ends before the bottom and length >= min_line_length,
 *         increment the corresponding entry in line_dist.
 *
 * Notes:
 *   - White lines that reach the bottom are ignored.
 *   - Only lines after the first black in each column are counted.
 */
void distr_white_vert_lines(size_t N, const uint8_t *recmat,
                            size_t min_line_length, double *line_dist);

/**
 * Compute the recurrence time entropy (RTE) of a trajectory using vertical
 * white lines.
 *
 * This is an RQA measure based on the distribution of vertical white lines
 * in the recurrence matrix.
 *
 * @param N               Number of points in the trajectory.
 * @param dim             Dimension of the phase space.
 * @param x               Pointer to trajectory array of size N*dim (row-major).
 * @param eps             Recurrence threshold.
 * @param min_line_length Minimum length of white lines to include.
 * @param recmat          Temporary recurrence matrix (N x N) for internal
 * computation.
 * @param line_distr      Temporary array of size N for white-line distribution.
 * @param norm_type       Norm to use in recurrence matrix computation:
 * NORM_MAX, NORM_L1, or NORM_L2.
 *
 * @return Recurrence time entropy (double, >= 0).
 *
 * Method:
 *   - Compute the binary recurrence matrix using recurrence_matrix().
 *   - Compute vertical white-line distribution using distr_white_vert_lines().
 *   - Compute probabilities p_k = line_distr[k] / sum(line_distr).
 *   - Compute Shannon entropy: -sum_k p_k * log(p_k).
 *
 * Notes:
 *   - Returns 0 if no white lines meet min_line_length.
 *   - relays on row-major ordering for x and recmat.
 */
double recurrence_time_entropy(size_t N, size_t dim, const double *x,
                               double eps, size_t min_line_length,
                               uint8_t *recmat, double *line_distr);

#endif
