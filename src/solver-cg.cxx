#include "solver-cg.h"

extern int LSSP_MAXIT;
extern double LSSP_ATOL;
extern double LSSP_RTOL;
extern double LSSP_RB;

int lssp_solver_cg(LSSP_SOLVER & solver, LSSP_PC & pc)
{
    lssp_mat_csr A = solver.A;
    int itr_max = solver.maxit;
    lssp_vec bg = solver.rhs;
    lssp_vec xg = solver.x;

    double tol_abs = solver.tol_abs;
    double tol_rel = solver.tol_rel;
    double tol_rb = solver.tol_rb;
    double b_norm;

    lssp_vec zg;
    lssp_vec rg;
    lssp_vec qg;
    lssp_vec pg;

    double rho0 = 0, rho1 = 0, beta;
    double residual;
    double tol = -1, err_rel = 0;
    int itr_out;
    int i;
    int n;
    double time;

    assert(solver.assembled);
    assert(pc.assembled);

    if (itr_max <= 0) itr_max = LSSP_MAXIT;
    if (tol_abs < 0) tol_abs = LSSP_ATOL;
    if (tol_rel < 0) tol_rel = LSSP_RTOL;

    assert(A.num_rows == A.num_cols);

    time = lssp_get_time();
    if (solver.verb >= 2) {
        lssp_printf("cg: maximal iteration: %d\n", itr_max);
        lssp_printf("cg: tolerance abs: %g\n", tol_abs);
        lssp_printf("cg: tolerance rel: %g\n", tol_rel);
        lssp_printf("cg: tolerance rbn: %g\n", tol_rb);
    }

    n = A.num_rows;
    zg = lssp_vec_create(n);
    rg = lssp_vec_create(n);
    pg = lssp_vec_create(n);
    qg = lssp_vec_create(n);

    b_norm = lssp_vec_norm(solver.rhs);
    tol_rb *= b_norm;

    lssp_mv_amxpbyz(-1, A, xg, 1, bg, rg);
    residual = lssp_vec_norm(rg);
    if (residual <= tol_abs) {
        itr_out = 0;
        goto end;
    }

    err_rel = residual;
    tol = tol_rel * err_rel;

    if (tol < tol_abs) tol = tol_abs;
    if (tol < tol_rb) tol = tol_rb;

    for (i = 0; i < n; i++) {
        zg.d[i] = 0;
    }

    for (itr_out = 0; itr_out < itr_max; itr_out++) {
        double alpha;

        pc.solve(&pc, zg, rg);
        rho1 = lssp_vec_dot(zg, rg);

        if (itr_out == 0) {
            for (i = 0; i < n; i++) {
                pg.d[i] = zg.d[i];
            }
        }
        else {
            beta = rho1 / rho0;

            for (i = 0; i < n; i++) {
                pg.d[i] = zg.d[i] + beta * pg.d[i];
            }
        }

        lssp_mv_mxy(A, pg, qg);
        alpha = lssp_vec_dot(qg, pg);

        alpha = rho1 / alpha;
        rho0 = rho1;

        for (i = 0; i < n; i++) {
            xg.d[i] = xg.d[i] + alpha * pg.d[i];
            rg.d[i] = rg.d[i] - alpha * qg.d[i];
        }

        residual = lssp_vec_norm(rg);

        if (solver.verb >= 1) {
            lssp_printf("cg: itr: %5d, abs res: %.6e, rel res: %.6e, rbn: %.6e\n",
                 itr_out, residual, (err_rel == 0 ? 0 : residual / err_rel),
                 (b_norm == 0 ? 0 : residual / b_norm));
        }

        if (residual <= tol) break;
    }

    if (itr_out < itr_max) itr_out += 1;

 end:

    solver.residual = residual;
    solver.nits = itr_out;

    lssp_vec_destroy(zg);
    lssp_vec_destroy(rg);
    lssp_vec_destroy(pg);
    lssp_vec_destroy(qg);

    time = lssp_get_time() - time;
    if (solver.verb >= 2) {
        lssp_printf("cg: total iteration: %d\n", itr_out);
        lssp_printf("cg: total time: %g\n", time);
    }

    return itr_out;
}
