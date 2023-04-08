
#include "solver-bicrstab.h"

int lssp_solver_bicrstab(LSSP_SOLVER &solver, LSSP_PC &pc)
{
    lssp_mat_csr A = solver.A;
    lssp_vec x;
    lssp_vec r, rtld, p, s, ap, ms, map, ams, z;
    double alpha, beta, omega, rho, rho_old, tdot1, tdot2;
    double nrm2, tol, ires;
    int iter, maxiter, n;

    x = solver.x;
    n = A.num_rows;
    maxiter = solver.maxit;

    rtld = lssp_vec_create(n);
    r = lssp_vec_create(n);
    s = lssp_vec_create(n);
    ms = lssp_vec_create(n);
    ams = lssp_vec_create(n);
    p = lssp_vec_create(n);
    ap = lssp_vec_create(n);
    map = lssp_vec_create(n);
    z = lssp_vec_create(n);

    /* Initial Residual */
    lssp_mv_amxpbyz(-1, A, x, 1, solver.rhs, r);
    nrm2 = ires = lssp_vec_norm(r);

    iter = 1;
    if (nrm2 <= solver.tol_abs) {
        solver.nits = 0;
        goto end;
    }

    tol = nrm2 * solver.tol_rel;
    alpha = lssp_vec_norm(solver.rhs) * solver.tol_rb;
    if (tol < solver.tol_abs) tol = solver.tol_abs;
    if (tol < alpha) tol = alpha;

    /* init */
    lssp_vec_copy(p, r);
    lssp_mv_mxy(A, p, rtld);
    pc.solve(&pc, z, r);
    lssp_vec_copy(p, z);

    rho_old = lssp_vec_dot(rtld, z);

    for (iter = 1; iter <= maxiter; iter++) {
        lssp_mv_mxy(A, p, ap);
        pc.solve(&pc, map, ap);

        tdot1 = lssp_vec_dot(rtld, map);
        alpha = rho_old / tdot1;
        lssp_vec_axpbyz(-alpha, ap, 1, r, s);

        nrm2 = lssp_vec_norm(s);
        if (nrm2 <= tol) {
            lssp_vec_axpby(alpha, p, 1, x);
            goto end;
        }

        lssp_vec_axpbyz(-alpha, map, 1, z, ms);
        lssp_mv_mxy(A, ms, ams);
        tdot1 = lssp_vec_dot(ams, s);
        tdot2 = lssp_vec_dot(ams, ams);
        omega = tdot1 / tdot2;

        lssp_vec_axpby(alpha, p, 1, x);
        lssp_vec_axpby(omega, ms, 1, x);
        lssp_vec_axpbyz(-omega, ams, 1, s, r);

        nrm2 = lssp_vec_norm(r);
        if (solver.verb >= 2) {
            lssp_printf("bicrstab: itr: %5d, abs res: %.6e, rel res: %.6e\n",
                    iter, nrm2, nrm2 / ires);
        }

        if (tol >= nrm2) {
            goto end;
        }

        pc.solve(&pc, z, r);
        rho = lssp_vec_dot(rtld, z);

        if (rho == 0.0) {
            goto end;
        }

        beta = (rho / rho_old) * (alpha / omega);
        lssp_vec_axpby(-omega, map, 1, p);
        lssp_vec_axpby(1, z, beta, p);

        rho_old = rho;
    }

end:

    solver.nits = iter;
    solver.residual = nrm2;

    lssp_vec_destroy(rtld);
    lssp_vec_destroy(r);
    lssp_vec_destroy(s);
    lssp_vec_destroy(ms);
    lssp_vec_destroy(ams);
    lssp_vec_destroy(ap);
    lssp_vec_destroy(p);
    lssp_vec_destroy(map);
    lssp_vec_destroy(z);

    return iter;
}
