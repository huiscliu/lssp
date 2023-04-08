
#include "solver-bicgsafe.h"

int lssp_solver_bicgsafe(LSSP_SOLVER & solver, LSSP_PC & pc)
{
    lssp_mat_csr A = solver.A;
    lssp_vec x;
    lssp_vec r, rtld, mr, amr, t, mt, p, ap;
    lssp_vec y, u, au, z;
    double alpha, beta;
    double rho, rho_old;
    double qsi, eta;
    double tmp, tdot[5];
    double nrm2, tol, ires;
    int iter, maxiter;
    int n;

    x = solver.x;
    n = A.num_rows;
    maxiter = solver.maxit;

    rtld = lssp_vec_create(n);
    r = lssp_vec_create(n);
    mr = lssp_vec_create(n);
    amr = lssp_vec_create(n);
    p = lssp_vec_create(n);
    ap = lssp_vec_create(n);
    t = lssp_vec_create(n);
    mt = lssp_vec_create(n);
    y = lssp_vec_create(n);
    u = lssp_vec_create(n);
    z = lssp_vec_create(n);
    au = lssp_vec_create(n);

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
    lssp_vec_copy(rtld, r);

    pc.solve(&pc, mr, r);
    lssp_mv_mxy(A, mr, amr);
    rho_old = lssp_vec_dot(rtld, r);
    lssp_vec_copy(ap, amr);
    lssp_vec_copy(p, mr);
    beta = 0.0;

    for (iter = 1; iter <= maxiter; iter++) {
        tdot[0] = lssp_vec_dot(rtld, ap);
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

        lssp_vec_copy(t, y);
        lssp_vec_scale(t, eta);
        lssp_vec_axpby(qsi, ap, 1, t);

        pc.solve(&pc, mt, t);

        lssp_vec_axpby(1, mt, eta * beta, u);
        lssp_mv_mxy(A, u, au);

        /* z = qsi*mr + eta*z - alpha*u */
        lssp_vec_scale(z, eta);
        lssp_vec_axpby(qsi, mr, 1, z);
        lssp_vec_axpby(-alpha, u, 1, z);

        /* y = qsi*amr + eta*y - alpha*au */
        lssp_vec_scale(y, eta);
        lssp_vec_axpby(qsi, amr, 1, y);
        lssp_vec_axpby(-alpha, au, 1, y);

        /* x = x + alpha*p + z */
        lssp_vec_axpby(alpha, p, 1, x);
        lssp_vec_axpby(1, z, 1, x);

        /* r = r - alpha*ap - y */
        lssp_vec_axpby(-alpha, ap, 1, r);
        lssp_vec_axpby(-1, y, 1.0, r);

        /* convergence check */
        nrm2 = lssp_vec_norm(r);

        if (solver.verb >= 1) {
            lssp_printf("bicgsafe: itr: %5d, abs res: %.6e, rel res: %.6e\n",
                    iter, nrm2, nrm2 / ires);
        }

        if (tol >= nrm2) {
            goto end;
        }

        /* rho = <rtld,r> */
        rho = lssp_vec_dot(rtld, r);
        if (rho == 0.0) {
            goto end;
        }

        beta = (rho / rho_old) * (alpha / qsi);
        pc.solve(&pc, mr, r);
        lssp_mv_mxy(A, mr, amr);

        lssp_vec_axpby(-1, u, 1.0, p);
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
    lssp_vec_destroy(ap);
    lssp_vec_destroy(p);
    lssp_vec_destroy(t);
    lssp_vec_destroy(mt);
    lssp_vec_destroy(y);
    lssp_vec_destroy(u);
    lssp_vec_destroy(z);
    lssp_vec_destroy(au);

    return iter;
}
