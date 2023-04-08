
#include "solver-cgs.h"

int lssp_solver_cgs(LSSP_SOLVER & solver, LSSP_PC & pc)
{
    lssp_mat_csr A = solver.A;
    lssp_vec x;
    lssp_vec r, rtld, p, phat, q, qhat, u, uhat, vhat;
    double alpha, beta, rho, rho_old, tdot1;
    double nrm2, tol, ires;
    int iter, maxiter, n;

    x = solver.x;
    n = A.num_rows;
    maxiter = solver.maxit;

    r = lssp_vec_create(n);
    rtld = lssp_vec_create(n);
    p = lssp_vec_create(n);
    phat = lssp_vec_create(n);
    q = lssp_vec_create(n);
    qhat = lssp_vec_create(n);
    u = lssp_vec_create(n);
    uhat = lssp_vec_create(n);
    vhat = uhat;

    alpha = 1.0;
    rho_old = 1.0;

    /* Initial Residual */
    lssp_mv_amxpbyz(-1, A, x, 1, solver.rhs, r);

    iter = 0;
    nrm2 = ires = lssp_vec_norm(r);
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
    lssp_vec_set_value(q, 0);
    lssp_vec_set_value(p, 0);

    for (iter = 1; iter <= maxiter; iter++) {
        /* rho = <rtld,r> */
        rho = lssp_vec_dot(rtld, r);

        /* test breakdown */
        if (rho == 0.0) {
            goto end;
        }

        /* beta = (rho / rho_old) */
        beta = (rho / rho_old);

        /* u = r + beta*q */
        lssp_vec_axpbyz(beta, q, 1, r, u);

        /* p = u + beta*(q + beta*p) */
        lssp_vec_axpby(1, q, beta, p);
        lssp_vec_axpby(1, u, beta, p);

        /* phat = M^-1 * p */
        pc.solve(&pc, phat, p);

        /* v = A * phat */
        lssp_mv_mxy(A, phat, vhat);

        /* tdot1 = <rtld,vhat> */
        tdot1 = lssp_vec_dot(rtld, vhat);

        /* test breakdown */
        if (tdot1 == 0.0) {
            goto end;
        }

        /* alpha = rho / tdot1 */
        alpha = rho / tdot1;

        /* q = u - alpha*vhat */
        lssp_vec_axpbyz(-alpha, vhat, 1, u, q);

        /* phat = u + q */
        /* uhat = M^-1 * (u + q) */
        lssp_vec_axpbyz(1, u, 1, q, phat);
        pc.solve(&pc, uhat, phat);

        /* x = x + alpha*uhat */
        lssp_vec_axpby(alpha, uhat, 1, x);

        /* qhat = A * uhat */
        lssp_mv_mxy(A, uhat, qhat);

        /* r = r - alpha*qhat */
        lssp_vec_axpby(-alpha, qhat, 1, r);

        /* convergence check */
        nrm2 = lssp_vec_norm(r);

        if (solver.verb >= 1) {
            lssp_printf("cgs: itr: %5d, abs res: %.6e, rel res: %.6e\n", iter, nrm2, nrm2 / ires);
        }

        if (tol >= nrm2) {
            goto end;
        }

        rho_old = rho;
    }

end:

    lssp_vec_destroy(r);
    lssp_vec_destroy(rtld);
    lssp_vec_destroy(p);
    lssp_vec_destroy(phat);
    lssp_vec_destroy(q);
    lssp_vec_destroy(qhat);
    lssp_vec_destroy(u);
    lssp_vec_destroy(uhat);

    solver.nits = iter;
    solver.residual = nrm2;

    return iter;
}
