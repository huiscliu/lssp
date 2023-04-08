
#include "solver-tri.h"

void lssp_pc_ilu_solve_lower_matrix(lssp_mat_csr L, double *x, double *rhs)
{
    int i, j;
    int end;
    int n = L.num_rows;
    double result;

    assert(x != NULL && rhs != NULL);

    for (i = 0; i < n; i++) {
        end = L.Ap[i + 1] - 1;

        result = rhs[i];
        for (j = L.Ap[i]; j < end; j++) {
            result = result - L.Ax[j] * x[L.Aj[j]];
        }

        assert(fabs(L.Ax[end]) > 0.);
        x[i] = result / L.Ax[end];
    }
}

void lssp_pc_ilu_solve_upper_matrix(lssp_mat_csr U, double *x, double *rhs)
{
    int i, j;
    int end;
    int n = U.num_rows;
    double result;

    assert(x != NULL && rhs != NULL);

    for (i = n - 1; i >= 0; i--) {
        end = U.Ap[i];

        result = rhs[i];
        for (j = U.Ap[i + 1] - 1; j > end; j--) {
            result = result - U.Ax[j] * x[U.Aj[j]];
        }

        assert(fabs(U.Ax[end]) > 0.);
        x[i] = result / U.Ax[end];
    }
}

void lssp_pc_ilu_solve_lu_matrix(lssp_mat_csr L, lssp_mat_csr U, double *x, double *rhs, double *cache)
{
    assert(x != NULL && rhs != NULL);
    assert(cache != NULL);

    lssp_pc_ilu_solve_lower_matrix(L, cache, rhs);
    lssp_pc_ilu_solve_upper_matrix(U, x, cache);
}

void lssp_pc_ilu_solve(LSSP_PC *pc, lssp_vec x, lssp_vec rhs)
{
    lssp_pc_ilu_solve_lu_matrix(pc->L, pc->U, x.d, rhs.d, pc->cache);
}
