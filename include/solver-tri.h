
#ifndef LSSP_SOLVER_TRI_H
#define LSSP_SOLVER_TRI_H

#include "type-defs.h"
#include "matrix-utils.h"

void lssp_pc_ilu_solve_lower_matrix(lssp_mat_csr L, double *x, double *rhs);
void lssp_pc_ilu_solve_upper_matrix(lssp_mat_csr L, double *x, double *rhs);

void lssp_pc_ilu_solve_lu_matrix(lssp_mat_csr L, lssp_mat_csr U, double *x, double *rhs, double *cache);
void lssp_pc_ilu_solve(LSSP_PC *pc, lssp_vec x, lssp_vec rhs);

#endif
