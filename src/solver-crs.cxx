
#include "solver-crs.h"

int lssp_solver_crs(LSSP_SOLVER & solver, LSSP_PC & pc)
{
    lssp_mat_csr A = solver.A;
    lssp_vec x;
    lssp_vec r, rtld, p, q, u, z, ap, map, uq, auq;
    double alpha, beta, rho, rho_old, tdot1;
    double nrm2, tol, ires;
    int iter, maxiter, n;

    x = solver.x;
    maxiter = solver.maxit;
    n = A.num_rows;

    r = lssp_vec_create(n);
    rtld = lssp_vec_create(n);
    p = lssp_vec_create(n);
    z = lssp_vec_create(n);
    u = z;
    uq = z;
    q = lssp_vec_create(n);
    ap = q;
    map = lssp_vec_create(n);
    auq = map;

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

    rho_old = 1.0;
    lssp_vec_set_value(q, 0.);
    lssp_vec_set_value(p, 0.);

    for (iter = 1; iter <= maxiter; iter++) {
        pc.solve(&pc, z, r);
        rho = lssp_vec_dot(rtld, z);

        /* test breakdown */
        if (rho == 0.0) {
            goto end;
        }

        beta = rho / rho_old;
        lssp_vec_axpbyz(beta, q, 1, z, u);
        lssp_vec_axpby(1, q, beta, p);
        lssp_vec_axpby(1, u, beta, p);

        lssp_mv_mxy(A, p, ap);
        pc.solve(&pc, map, ap);
        tdot1 = lssp_vec_dot(rtld, map);

        /* test breakdown */
        if (tdot1 == 0.0) {
            goto end;
        }

        alpha = rho / tdot1;
        lssp_vec_axpbyz(-alpha, map, 1, u, q);
        lssp_vec_axpbyz(1, u, 1, q, uq);

        lssp_mv_mxy(A, uq, auq);
        lssp_vec_axpby(alpha, uq, 1, x);
        lssp_vec_axpby(-alpha, auq, 1, r);

        /* convergence check */
        nrm2 = lssp_vec_norm(r);

        if (solver.verb >= 1) {
            lssp_printf("crs: itr: %5d, abs res: %.6e, rel res: %.6e\n", iter, nrm2, nrm2 / ires);
        }

        if (tol >= nrm2) {
            goto end;
        }

        rho_old = rho;
    }

end:

    solver.nits = iter;
    solver.residual = nrm2;

    lssp_vec_destroy(r);
    lssp_vec_destroy(rtld);
    lssp_vec_destroy(p);
    lssp_vec_destroy(z);
    lssp_vec_destroy(q);
    lssp_vec_destroy(map);

    return iter;
}
