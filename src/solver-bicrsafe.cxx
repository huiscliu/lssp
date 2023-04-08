
#include "solver-bicrsafe.h"

int lssp_solver_bicrsafe(LSSP_SOLVER &solver, LSSP_PC &pc)
{
    lssp_mat_csr A = solver.A;
    lssp_vec x;
    lssp_vec r, rtld, artld, mr, amr, p, ap, map;
    lssp_vec y, my, u, au, z;
    double alpha, beta;
    double rho, rho_old;
    double qsi, eta;
    double tmp, tdot[5];
    double nrm2, tol, ires;
    int iter, maxiter, n;

    x = solver.x;
    maxiter = solver.maxit;
    n = A.num_rows;

    rtld = lssp_vec_create(n);
    r = lssp_vec_create(n);
    mr = lssp_vec_create(n);
    amr = lssp_vec_create(n);
    p = lssp_vec_create(n);
    ap = lssp_vec_create(n);
    map = lssp_vec_create(n);
    my = lssp_vec_create(n);
    y = lssp_vec_create(n);
    u = lssp_vec_create(n);
    z = lssp_vec_create(n);
    au = lssp_vec_create(n);
    artld = lssp_vec_create(n);

    /* initial Residual */
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
    lssp_vec_copy(rtld, r);
    lssp_mv_mxy(A, rtld, artld);

    pc.solve(&pc, mr, r);
    lssp_mv_mxy(A, mr, amr);
    rho_old = lssp_vec_dot(rtld, amr);
    lssp_vec_copy(ap, amr);
    lssp_vec_copy(p, mr);
    beta = 0.0;

    for (iter = 1; iter <= maxiter; iter++) {
        pc.solve(&pc, map, ap);
        tdot[0] = lssp_vec_dot(artld, map);
        alpha = rho_old / tdot[0];

        tdot[0] = lssp_vec_dot(y, y);
        tdot[1] = lssp_vec_dot(amr, r);
        tdot[2] = lssp_vec_dot(y, r);
        tdot[3] = lssp_vec_dot(amr, y);
        tdot[4] = lssp_vec_dot(amr, amr);
        if (iter == 1) {
            qsi = tdot[1] / tdot[4];
            eta = 0.0;
        }
        else {
            tmp = tdot[4] * tdot[0] - tdot[3] * tdot[3];
            qsi = (tdot[0] * tdot[1] - tdot[2] * tdot[3]) / tmp;
            eta = (tdot[4] * tdot[2] - tdot[3] * tdot[1]) / tmp;
        }

        lssp_vec_scale(u, eta * beta);
        lssp_vec_axpby(qsi, map, 1, u);
        lssp_vec_axpby(eta, my, 1, u);
        lssp_mv_mxy(A, u, au);

        /* z = qsi*mr + eta*z - alpha*u */
        lssp_vec_scale(z, eta);
        lssp_vec_axpby(qsi, mr, 1, z);
        lssp_vec_axpby(-alpha, u, 1, z);

        lssp_vec_scale(y, eta);
        lssp_vec_axpby(qsi, amr, 1, y);
        lssp_vec_axpby(-alpha, au, 1, y);
        pc.solve(&pc, my, y);

        /* x = x + alpha*p + z */
        lssp_vec_axpby(alpha, p, 1, x);
        lssp_vec_axpby(1, z, 1, x);

        /* r = r - alpha*ap - y */
        lssp_vec_axpby(-alpha, ap, 1, r);
        lssp_vec_axpby(-1, y, 1, r);

        /* convergence check */
        nrm2 = lssp_vec_norm(r);

        if (solver.verb >= 2) {
            lssp_printf("bicrsafe: itr: %5d, abs res: %.6e, rel res: %.6e\n",
                    iter, nrm2, nrm2 / ires);
        }

        if (tol >= nrm2) goto end;

        lssp_vec_axpby(-alpha, map, 1, mr);
        lssp_vec_axpby(-1, my, 1, mr);
        lssp_mv_mxy(A, mr, amr);
        rho = lssp_vec_dot(rtld, amr);
        if (rho == 0.0) {
            goto end;
        }

        beta = (rho / rho_old) * (alpha / qsi);
        lssp_vec_axpby(-1, u, 1., p);
        lssp_vec_axpby(1, mr, beta, p);
        lssp_vec_axpby(-1, au, 1.0, ap);
        lssp_vec_axpby(1, amr, beta, ap);

        rho_old = rho;
    }

end:

    solver.nits = iter;
    solver.residual = nrm2;

    lssp_vec_destroy(rtld);
    lssp_vec_destroy(r);
    lssp_vec_destroy(mr);
    lssp_vec_destroy(amr);
    lssp_vec_destroy(map);
    lssp_vec_destroy(ap);
    lssp_vec_destroy(p);
    lssp_vec_destroy(my);
    lssp_vec_destroy(y);
    lssp_vec_destroy(u);
    lssp_vec_destroy(z);
    lssp_vec_destroy(au);
    lssp_vec_destroy(artld);

    return iter;
}
