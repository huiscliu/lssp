
#include "solver-tfqmr.h"

int lssp_solver_tfqmr(LSSP_SOLVER &solver, LSSP_PC & pc)
{
    lssp_mat_csr A = solver.A;
    lssp_vec x;
    lssp_vec r, rtld, u, p, d, t, t1, q, v;
    double tau, rho, rhoold, theta, eta, beta, alpha, w, ww, wold, s, c;
    double nrm2, tol, ires;
    int iter, maxiter;
    int m, n;

    x = solver.x;
    maxiter = solver.maxit;
    n = A.num_rows;

    r = lssp_vec_create(n);
    rtld = lssp_vec_create(n);
    u = lssp_vec_create(n);
    p = lssp_vec_create(n);
    d = lssp_vec_create(n);
    t = lssp_vec_create(n);
    t1 = lssp_vec_create(n);
    q = lssp_vec_create(n);
    v = lssp_vec_create(n);

    /* Initial Residual */
    lssp_mv_amxpbyz(-1, A, x, 1, solver.rhs, r);
    nrm2 = ires = lssp_vec_norm(r);

    iter = 1;
    if (nrm2 <= solver.tol_abs) {
        iter = 0;
        goto end;
    }

    tol = nrm2 * solver.tol_rel;
    alpha = lssp_vec_norm(solver.rhs) * solver.tol_rb;
    if (tol < solver.tol_abs) tol = solver.tol_abs;
    if (tol < alpha) tol = alpha;

    /* rtld = r
       p    = r
       u    = r
       d    = 0
       v    = Ap */
    lssp_vec_copy(rtld, r);
    lssp_vec_copy(p, r);
    lssp_vec_copy(u, r);
    lssp_vec_set_value(d, 0.);

    pc.solve(&pc, t, p);
    lssp_mv_mxy(A, t, v);

    /* rhoold = (r,rtld) */
    rhoold = lssp_vec_dot(r, rtld);
    tau = lssp_vec_norm(r);
    wold = tau;
    theta = 0.0;
    eta = 0.0;

    while (iter <= maxiter) {
        /* alpha = rho / (v,rtld) */
        s = lssp_vec_dot(v, rtld);

        /* test breakdown */
        if (fabs(s) == 0.0) {
            goto end;
        }

        alpha = rhoold / s;

        /* q = u - alpha*v */
        lssp_vec_axpbyz(-alpha, v, 1, u, q);

        /* r = r - alpha*A(u + q) */
        lssp_vec_axpbyz(1, u, 1., q, t);
        pc.solve(&pc, t1, t);
        lssp_mv_mxy(A, t1, v);
        lssp_vec_axpby(-alpha, v, 1, r);

        w = lssp_vec_norm(r);
        for (m = 0; m < 2; m++)
        {
            if (m == 0) {
                ww = sqrt(w * wold);
                lssp_vec_axpby(1, u, theta * theta * eta / alpha, d);
            }
            else {
                ww = w;
                lssp_vec_axpby(1, q, theta * theta * eta / alpha, d);
            }

            theta = ww / tau;
            c = 1.0 / sqrt(1.0 + theta * theta);
            eta = c * c * alpha;
            tau = tau * theta * c;

            pc.solve(&pc, t1, d);
            lssp_vec_axpby(eta, t1, 1, x);

            /* convergence check */
            nrm2 = tau * sqrt(1.0 + m);

            if (solver.verb >= 1) {
                lssp_printf("tfqmr: itr: %5d, abs res: %.6e, rel res: %.6e\n", iter, nrm2, nrm2 / ires);
            }

            if (tol >= nrm2) {
                goto end;
            }
        }

        rho = lssp_vec_dot(r, rtld);
        if (fabs(rho) == 0.0) {
            goto end;
        }

        beta = rho / rhoold;
        lssp_vec_axpbyz(beta, q, 1, r, u);

        lssp_vec_axpby(1, q, beta, p);
        lssp_vec_axpby(1, u, beta, p);

        pc.solve(&pc, t1, p);
        lssp_mv_mxy(A, t1, v);

        rhoold = rho;
        wold = w;
        iter++;
    }

end:
    solver.nits = iter;
    solver.residual = nrm2;

    lssp_vec_destroy(r);
    lssp_vec_destroy(rtld);
    lssp_vec_destroy(u);
    lssp_vec_destroy(p);
    lssp_vec_destroy(d);
    lssp_vec_destroy(t);
    lssp_vec_destroy(t1);
    lssp_vec_destroy(q);
    lssp_vec_destroy(v);

    return iter;
}
