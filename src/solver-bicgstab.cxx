#include "solver-bicgstab.h"

extern double LSSP_BREAKDOWN;
extern int LSSP_MAXIT;
extern double LSSP_ATOL;
extern double LSSP_RTOL;
extern double LSSP_RB;


int lssp_solver_bicgstab(LSSP_SOLVER &solver, LSSP_PC &pc)
{
    lssp_mat_csr A = solver.A;
    int itr_max = solver.maxit;
    lssp_vec bg = solver.rhs;
    lssp_vec xg = solver.x;

    double tol_abs = solver.tol_abs;
    double tol_rel = solver.tol_rel;
    double tol_rb = solver.tol_rb;
    double b_norm;

    lssp_vec rg;
    lssp_vec rh;
    lssp_vec pg;
    lssp_vec ph;
    lssp_vec sg;
    lssp_vec sh;
    lssp_vec tg;
    lssp_vec vg;

    double rho0 = 0, rho1 = 0;
    double alpha = 0, beta = 0, omega = 0;
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
        lssp_printf("bicgstab: maximal iteration: %d\n", itr_max);
        lssp_printf("bicgstab: tolerance abs: %g\n", tol_abs);
        lssp_printf("bicgstab: tolerance rel: %g\n", tol_rel);
        lssp_printf("bicgstab: tolerance rbn: %g\n", tol_rb);
    }

    n = A.num_rows;
    rg = lssp_vec_create(n);
    rh = lssp_vec_create(n);
    pg = lssp_vec_create(n);
    ph = lssp_vec_create(n);
    sg = lssp_vec_create(n);
    sh = lssp_vec_create(n);
    tg = lssp_vec_create(n);
    vg = lssp_vec_create(n);

    lssp_mv_amxpbyz(-1, A, xg, 1, bg, rg);
    for (i = 0; i < n; i++) {
        rh.d[i] = rg.d[i];
        sh.d[i] = ph.d[i] = 0.;
    }

    b_norm = lssp_vec_norm(solver.rhs);
    tol_rb *= b_norm;

    residual = err_rel = lssp_vec_norm(rg);
    if (residual <= tol_abs) {
        itr_out = 0;
        goto end;
    }

    tol = residual * tol_rel;
    if (tol < tol_abs) tol = tol_abs;
    if (tol < tol_rb) tol = tol_rb;

    for (itr_out = 0; itr_out < itr_max; itr_out++) {
        rho1 = lssp_vec_dot(rg, rh);

        if (rho1 == 0) {
            lssp_printf("bicgstab: method failed.!\n");
            break;
        }

        if (itr_out == 0) {

            for (i = 0; i < n; i++) pg.d[i] = rg.d[i];
        }
        else {
            beta = (rho1 * alpha) / (rho0 * omega);
            for (i = 0; i < n; i++) {
                pg.d[i] = rg.d[i] + beta * (pg.d[i] - omega * vg.d[i]);
            }
        }

        rho0 = rho1;

        lssp_vec_set_value(ph, 0.);
        pc.solve(&pc, ph, pg);

        lssp_mv_amxpbyz(1, A, ph, 0, pg, vg);

        alpha = rho1 / lssp_vec_dot(rh, vg);
        for (i = 0; i < n; i++) {
            sg.d[i] = rg.d[i] - alpha * vg.d[i];
        }

        if (lssp_vec_norm(sg) <= LSSP_BREAKDOWN) {
            lssp_printf("bicgstab: ||s|| is too small: %f, terminated.\n", lssp_vec_norm(sg));

            for (i = 0; i < n; i++) {
                xg.d[i] = xg.d[i] + alpha * ph.d[i];
            }

            lssp_mv_amxpbyz(-1, A, xg, 1, bg, rg);
            residual = lssp_vec_norm(rg);

            break;
        }

        lssp_vec_set_value(sh, 0.);
        pc.solve(&pc, sh, sg);

        lssp_mv_amxpbyz(1, A, sh, 0, pg, tg);

        omega = lssp_vec_dot(tg, sg) / lssp_vec_dot(tg, tg);
        for (i = 0; i < n; i++) {
            xg.d[i] = xg.d[i] + alpha * ph.d[i] + omega * sh.d[i];
            rg.d[i] = sg.d[i] - omega * tg.d[i];
        }

        residual = lssp_vec_norm(rg);

        if (solver.verb >= 1) {
            lssp_printf("bicgstab: itr: %5d, abs res: %.6e, rel res: %.6e, "
                    "rbn: %.6e\n", itr_out, residual, (err_rel == 0 ? 0 : residual / err_rel), 
                    (b_norm == 0 ? 0 : residual / b_norm));
        }

        if (residual <= tol) break;
    }

    if (itr_out < itr_max) itr_out += 1;

end:

    solver.residual = residual;
    solver.nits = itr_out;

    lssp_vec_destroy(rg);
    lssp_vec_destroy(rh);
    lssp_vec_destroy(pg);
    lssp_vec_destroy(ph);
    lssp_vec_destroy(sg);
    lssp_vec_destroy(sh);
    lssp_vec_destroy(tg);
    lssp_vec_destroy(vg);

    time = lssp_get_time() - time;
    if (solver.verb >= 2) {
        lssp_printf("bicgstab: total iteration: %d\n", itr_out);
        lssp_printf("bicgstab: total time: %g\n", time);
    }

    return itr_out;
}
