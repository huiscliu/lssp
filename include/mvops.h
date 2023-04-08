
#ifndef LSSP_MVOPS_H
#define LSSP_MVOPS_H

#include "matrix-utils.h"
#include "vector.h"

/* y = y * beta + alpha * A * x  */
void lssp_mv_amxpby(double alpha, const lssp_mat_csr A, const lssp_vec x,  double beta, lssp_vec y);

/* z = y * beta + alpha * A * x  */
void lssp_mv_amxpbyz(double alpha, const lssp_mat_csr A, const lssp_vec x,  double beta,
        const lssp_vec y, lssp_vec z);

/* y = a * A * x  */
void lssp_mv_amxy(double a, const lssp_mat_csr A, const lssp_vec x, lssp_vec y);

/* y = A * x  */
void lssp_mv_mxy(const lssp_mat_csr A, const lssp_vec x,  lssp_vec y);

#endif
