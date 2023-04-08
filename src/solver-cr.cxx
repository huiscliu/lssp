
#include "solver-cr.h"

int lssp_solver_cr(LSSP_SOLVER & solver, LSSP_PC & pc)
{
    lssp_mat_csr A = solver.A;
    lssp_vec x;
    lssp_vec r, z, p, q, qtld, az;
    double alpha, beta, rho;
    double dot_rq, dot_zq;

    double nrm2, tol, ires;
    int iter, maxiter, n;

    x = solver.x;
    maxiter = solver.maxit;
    n = A.num_rows;

    r = lssp_vec_create(n);
    z = lssp_vec_create(n);
    p = lssp_vec_create(n);
    q = lssp_vec_create(n);
    qtld = lssp_vec_create(n);
    az = lssp_vec_create(n);

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
    if (tol < solver.tol_abs)
        tol = solver.tol_abs;
    if (tol < alpha)
        tol = alpha;

    /* init */
    pc.solve(&pc, p, r);
    lssp_mv_mxy(A, p, q);
    lssp_vec_copy(z, p);

    for (iter = 1; iter <= maxiter; iter++) {
        /* qtld = M^-1 * q */
        pc.solve(&pc, qtld, q);

        /* rho = <qtld,q> */
        rho = lssp_vec_dot(qtld, q);

        /* breakdown check */
        if (rho == 0.0) {
            goto end;
        }

        /* dot_rq = <r,qtld> */
        dot_rq = lssp_vec_dot(r, qtld);

        /* alpha = dot_rq / rho */
        alpha = dot_rq / rho;

        /* x = x + alpha*p */
        lssp_vec_axpby(alpha, p, 1, x);

        /* r = r - alpha*q */
        lssp_vec_axpby(-alpha, q, 1, r);

        /* convergence check */
        nrm2 = lssp_vec_norm(r);

        if (solver.verb >= 1) {
            lssp_printf("cr: itr: %5d, abs res: %.6e, rel res: %.6e\n", iter, nrm2, nrm2 / ires);
        }

        if (tol >= nrm2) {
            goto end;
        }

        /* z = z - alpha*qtld */
        lssp_vec_axpby(-alpha, qtld, 1, z);

        /* az = Az */
        lssp_mv_mxy(A, z, az);

        /* dot_zq = <az,qtld> */
        dot_zq = lssp_vec_dot(az, qtld);

        /* beta = -dot_zq / rho */
        beta = -dot_zq / rho;

        /* p = z + beta*p */
        lssp_vec_axpby(1, z, beta, p);

        /* q = az + beta*q */
        lssp_vec_axpby(1, az, beta, q);
    }

 end:

    solver.nits = iter;
    solver.residual = nrm2;

    lssp_vec_destroy(r);
    lssp_vec_destroy(z);
    lssp_vec_destroy(p);
    lssp_vec_destroy(q);
    lssp_vec_destroy(qtld);
    lssp_vec_destroy(az);

    return iter;
}
