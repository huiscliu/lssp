
#include "solver-gpbicg.h"

int lssp_solver_gpbicg(LSSP_SOLVER & solver, LSSP_PC & pc)
{
    lssp_mat_csr A = solver.A;
    lssp_vec x;
    lssp_vec r, rtld, mr, p, ap, map, t, mt_old, mt, amt, y, u, w, z;
    double alpha, beta, rho, rho_old;
    double qsi, eta;
    double tmp, tdot[5];
    double bnrm2, nrm2, tol, ires;
    int iter, maxiter, n;

    x = solver.x;
    n = A.num_rows;
    maxiter = solver.maxit;

    rtld = lssp_vec_create(n);
    r = lssp_vec_create(n);
    mr = lssp_vec_create(n);
    p = lssp_vec_create(n);
    ap = lssp_vec_create(n);
    map = lssp_vec_create(n);
    t = lssp_vec_create(n);
    mt = lssp_vec_create(n);
    amt = lssp_vec_create(n);
    u = lssp_vec_create(n);
    y = lssp_vec_create(n);
    w = lssp_vec_create(n);
    z = lssp_vec_create(n);
    mt_old = lssp_vec_create(n);

    /* Initial Residual */
    lssp_mv_amxpbyz(-1, A, x, 1, solver.rhs, r);
    nrm2 = ires = bnrm2 = lssp_vec_norm(r);

    iter = 1;
    if (nrm2 <= solver.tol_abs) {
        iter = 0;
        goto end;
    }

    tol = nrm2 * solver.tol_rel;
    alpha = lssp_vec_norm(solver.rhs) * solver.tol_rb;
    if (tol < solver.tol_abs) tol = solver.tol_abs;
    if (tol < alpha) tol = alpha;

    lssp_vec_copy(rtld, r);

    pc.solve(&pc, p, r);
    rho_old = lssp_vec_dot(rtld, r);
    lssp_vec_set_value(t, 0.);
    lssp_vec_set_value(w, 0.);
    beta = 0.0;

    for (iter = 1; iter <= maxiter; iter++) {
        lssp_mv_mxy(A, p, ap);
        pc.solve(&pc, map, ap);
        tdot[0] = lssp_vec_dot(rtld, ap);

        if (tdot[0] == 0.0) {
            goto end;
        }

        alpha = rho_old / tdot[0];
        lssp_vec_axpbyz(-1, w, 1, ap, y);
        lssp_vec_axpby(1, t, alpha, y);
        lssp_vec_axpby(-1, r, 1, y);
        lssp_vec_axpbyz(-alpha, ap, 1, r, t);

        nrm2 = lssp_vec_norm(t);
        if (solver.verb >= 1) {
            lssp_printf("gpbicg: itr: %5d, abs res: %.6e, rel res: %.6e\n",
                    iter, nrm2, nrm2 / ires);
        }

        if (nrm2 <= tol) {
            lssp_vec_axpby(alpha, p, 1, x);
            goto end;
        }

        lssp_vec_axpbyz(-alpha, map, 1, mr, mt);
        lssp_mv_mxy(A, mt, amt);
        tdot[0] = lssp_vec_dot(y, y);
        tdot[1] = lssp_vec_dot(amt, t);
        tdot[2] = lssp_vec_dot(y, t);
        tdot[3] = lssp_vec_dot(amt, y);
        tdot[4] = lssp_vec_dot(amt, amt);
        if (iter == 1) {
            qsi = tdot[1] / tdot[4];
            eta = 0.0;
        }
        else {
            tmp = tdot[4] * tdot[0] - tdot[3] * tdot[3];
            qsi = (tdot[0] * tdot[1] - tdot[2] * tdot[3]) / tmp;
            eta = (tdot[4] * tdot[2] - tdot[3] * tdot[1]) / tmp;
        }

        lssp_vec_axpby(1., mt_old, beta, u);
        lssp_vec_axpby(-1, mr, 1, u);
        lssp_vec_scale(u, eta);
        lssp_vec_axpby(qsi, map, 1, u);

        lssp_vec_scale(z, eta);
        lssp_vec_axpby(qsi, mr, 1, z);
        lssp_vec_axpby(-alpha, u, 1, z);

        lssp_vec_axpby(alpha, p, 1, x);
        lssp_vec_axpby(1, z, 1., x);

        lssp_vec_axpbyz(-qsi, amt, 1, t, r);
        lssp_vec_axpby(-eta, y, 1, r);

        /* convergence check */
        nrm2 = lssp_vec_norm(r);

        if (solver.verb >= 1) {
            lssp_printf("gpbicg: itr: %5d, abs res: %.6e, rel res: %.6e\n", iter, nrm2, nrm2 / ires);
        }

        if (tol >= nrm2) {
            goto end;
        }

        pc.solve(&pc, mr, r);
        rho = lssp_vec_dot(rtld, r);
        if (rho == 0.0) {
            goto end;
        }

        beta = (rho / rho_old) * (alpha / qsi);

        lssp_vec_axpbyz(beta, ap, 1, amt, w);
        lssp_vec_axpby(-1, u, 1, p);
        lssp_vec_axpby(1., mr, beta, p);

        lssp_vec_copy(mt_old, mt);
        rho_old = rho;
    }

end:

    solver.nits = iter;
    solver.residual = nrm2;

    lssp_vec_destroy(rtld);
    lssp_vec_destroy(r);
    lssp_vec_destroy(mr);
    lssp_vec_destroy(p);
    lssp_vec_destroy(ap);
    lssp_vec_destroy(map);
    lssp_vec_destroy(t);
    lssp_vec_destroy(mt);
    lssp_vec_destroy(amt);
    lssp_vec_destroy(u);
    lssp_vec_destroy(y);
    lssp_vec_destroy(w);
    lssp_vec_destroy(z);
    lssp_vec_destroy(mt_old);

    return iter;
}
