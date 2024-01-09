
#include "lssp.h"

static lssp_mat_csr laplacian_5pt(const int N)
{
    lssp_mat_csr csr;
    int nz = 0;
    int i, j;

    assert (N > 0);

    csr.num_rows = N * N;
    csr.num_cols = N * N;
    csr.num_nnzs = 5 * N * N - 4 * N;
    csr.Ap = lssp_malloc<int>(csr.num_rows + 1);
    csr.Aj = lssp_malloc<int>(csr.num_nnzs);
    csr.Ax = lssp_malloc<double>(csr.num_nnzs);

    csr.Ap[0] = 0;
    for(i = 0; i < csr.num_rows; i++) csr.Ap[i + 1] = 0;

    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            int indx = N * i + j;

            if (i > 0){
                csr.Aj[nz] = indx - N;
                csr.Ax[nz] = -1;
                nz++;
            }

            if (j > 0){
                csr.Aj[nz] = indx - 1;
                csr.Ax[nz] = -1;
                nz++;
            }

            csr.Aj[nz] = indx;
            csr.Ax[nz] = 4;
            nz++;

            if (j < N - 1){
                csr.Aj[nz] = indx + 1;
                csr.Ax[nz] = -1;
                nz++;
            }

            if (i < N - 1){
                csr.Aj[nz] = indx + N;
                csr.Ax[nz] = -1;
                nz++;
            }

            csr.Ap[indx + 1] = nz;
        }
    }

    return csr;
}

int main()
{
    int i, n;
    int dim = 100;
    int m = 60;
    int itr_max = 3000;

    lssp_mat_csr A;
    LSSP_SOLVER solver;
    LSSP_PC pc;
    lssp_vec xt, xg;
    lssp_vec bc;
    double tc, ta;
    int bs = 4;

    lssp_verbosity = 2;

    lssp_printf("CSR: laplacian, grid size %d\n", dim);
    A = laplacian_5pt(dim);

    n = A.num_rows;
    dim = A.num_nnzs;
    lssp_printf("CSR: rows: %d", n);
    lssp_printf(" nonzeros: %d mem (csr): %.3f Mb\n", dim, \
            ((dim + n) * sizeof(int) + dim * sizeof(double))/ 1048576.);

    /* generate vector x, b */
    xg = lssp_vec_create(n);
    bc = lssp_vec_create(n);
    xt = lssp_vec_create(n);

    for (i = 0; i < n; i++) {
        lssp_vec_set_value_by_index(xg, i, 0.);
        lssp_vec_set_value_by_index(bc, i, 1.);
    }

    /* set up solver and pc */
    lssp_solver_create(solver, LSSP_SOLVER_GMRES, pc, LSSP_PC_ILUK);

    /* settings */
    lssp_solver_set_restart(solver, m);
    lssp_solver_set_maxit(solver, itr_max);
    solver.num_blks = A.num_rows / bs;

    /* init and setup preconditioner */
    lssp_solver_assemble(solver, A, xg, bc, pc);

    /* solve */
    ta = lssp_get_time();
    lssp_solver_solve(solver, pc);
    tc = lssp_get_time() - ta;

    lssp_printf("total solver time: %g \n", tc);
    lssp_printf("solution L2 norm: %.8e residual: %.8e\n", lssp_vec_norm(xg), solver.residual);
    lssp_mv_amxpbyz(-1, A, xg, 1, bc, xt);
    lssp_printf("verification, residual: %.8e\n\n", lssp_vec_norm(xt));

    /* destroy pc and sover */
    lssp_solver_destroy(solver, pc);

    lssp_mat_destroy(A);
    lssp_vec_destroy(xg);
    lssp_vec_destroy(bc);
    lssp_vec_destroy(xt);

    return 0;
}
