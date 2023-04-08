
#include "solver-bicgstabl.h"

int lssp_solver_bicgstabl(LSSP_SOLVER & solver, LSSP_PC & pc)
{
    lssp_mat_csr A = solver.A;
    lssp_vec x;
    lssp_vec rtld, bp, t, xp, *r, *u;

    double *tau, *gamma, *gamma1, *gamma2;
    double *sigma;
    double alpha, beta, omega, rho0, rho1;
    double nu;
    int n;

    double nrm2, tol, ires;
    int iter, maxiter;

    int l, i, j;
    int z_dim;

    A = solver.A;
    x = solver.x;
    n = A.num_rows;

    maxiter = solver.maxit;
    l = solver.bgsl;

    if (l <= 0) l = 4;

    z_dim = l + 1;
    rtld = lssp_vec_create(n);
    xp = lssp_vec_create(n);
    bp = lssp_vec_create(n);
    t = lssp_vec_create(n);

    r = lssp_malloc<lssp_vec>(l + 1);
    for (i = 0; i <= l; i++) r[i] = lssp_vec_create(n);

    u = lssp_malloc<lssp_vec>(l + 1);
    for (i = 0; i <= l; i++) u[i] = lssp_vec_create(n);

    tau = lssp_malloc<double>(z_dim * (4 + l + 1));
    gamma = &tau[z_dim * z_dim];
    gamma1 = &gamma[z_dim];
    gamma2 = &gamma1[z_dim];
    sigma = &gamma2[z_dim];

    /* copy r[0] to rtld */
    lssp_mv_amxpbyz(-1, A, x, 1, solver.rhs, r[0]);
    lssp_vec_copy(rtld, r[0]);
    lssp_vec_copy(bp, r[0]);
    lssp_vec_copy(xp, x);
    lssp_vec_set_value(u[0], 0.);

    /* Initial Residual */
    iter = 0;
    nrm2 = ires = lssp_vec_norm(r[0]);
    if (nrm2 <= solver.tol_abs) {
        solver.nits = 0;

        goto end;
    }

    tol = nrm2 * solver.tol_rel;
    alpha = lssp_vec_norm(solver.rhs) * solver.tol_rb;
    if (tol < solver.tol_abs) tol = solver.tol_abs;
    if (tol < alpha) tol = alpha;

    alpha = 0.0;
    omega = 1.0;
    rho0 = 1.0;
    while (iter <= maxiter) {
        rho0 = -omega * rho0;

        for (j = 0; j < l; j++) {
            iter++;

            /* rho1 = <rtld,r[j]> */
            rho1 = lssp_vec_dot(rtld, r[j]);

            /* test breakdown */
            if (rho1 == 0.0) {
                pc.solve(&pc, t, x);
                lssp_vec_copy(x, t);
                lssp_vec_axpby(1, xp, 1, x);

                goto end;
            }

            beta = alpha * (rho1 / rho0);
            rho0 = rho1;

            /* u[i] = r[i] - beta*u[i] (i=0,j) */
            for (i = 0; i <= j; i++) {
                lssp_vec_axpby(1, r[i], -beta, u[i]);
            }

            pc.solve(&pc, t, u[j]);
            lssp_mv_mxy(A, t, u[j + 1]);
            nu = lssp_vec_dot(rtld, u[j + 1]);

            /* test breakdown */
            if (fabs(nu) == 0.0) {
                pc.solve(&pc, t, x);
                lssp_vec_copy(x, t);
                lssp_vec_axpby(1, xp, 1, x);

                goto end;
            }

            /* alpha = rho1 / nu */
            alpha = rho1 / nu;

            /* x = x + alpha*u[0] */
            lssp_vec_axpby(alpha, u[0], 1, x);

            /* r[i] = r[i] - alpha*u[i+1] (i=0,j) */
            for (i = 0; i <= j; i++) {
                lssp_vec_axpby(-alpha, u[i + 1], 1, r[i]);
            }

            nrm2 = lssp_vec_norm(r[0]);

            if (solver.verb >= 1) {
                lssp_printf("bicgstabl: itr: %5d, abs res: %.6e, rel res: %.6e\n",
                        iter, nrm2, nrm2 / ires);
            }

            if (nrm2 <= tol) {
                pc.solve(&pc, t, x);
                lssp_vec_copy(x, t);
                lssp_vec_axpby(1, xp, 1, x);

                goto end;
            }

            pc.solve(&pc, t, r[j]);
            lssp_mv_mxy(A, t, r[j + 1]);
        }

        /* MR PART */
        for (j = 1; j <= l; j++) {
            for (i = 1; i <= j - 1; i++) {
                nu = lssp_vec_dot(r[j], r[i]);
                nu = nu / sigma[i];
                tau[i * z_dim + j] = nu;
                lssp_vec_axpby(-nu, r[i], 1, r[j]);
            }

            sigma[j] = lssp_vec_dot(r[j], r[j]);
            nu = lssp_vec_dot(r[0], r[j]);
            gamma1[j] = nu / sigma[j];
        }

        gamma[l] = gamma1[l];
        omega = gamma[l];
        for (j = l - 1; j >= 1; j--) {
            nu = 0.0;
            for (i = j + 1; i <= l; i++) {
                nu += tau[j * z_dim + i] * gamma[i];
            }
            gamma[j] = gamma1[j] - nu;
        }
        for (j = 1; j <= l - 1; j++) {
            nu = 0.0;
            for (i = j + 1; i <= l - 1; i++) {
                nu += tau[j * z_dim + i] * gamma[i + 1];
            }
            gamma2[j] = gamma[j + 1] + nu;
        }

        /* UPDATE */
        lssp_vec_axpby(gamma[1], r[0], 1, x);
        lssp_vec_axpby(-gamma1[l], r[l], 1, r[0]);
        lssp_vec_axpby(-gamma[l], u[l], 1, u[0]);

        for (j = 1; j <= l - 1; j++) {
            lssp_vec_axpby(-gamma[j], u[j], 1, u[0]);
            lssp_vec_axpby(gamma2[j], r[j], 1, x);
            lssp_vec_axpby(-gamma1[j], r[j], 1, r[0]);
        }

        nrm2 = lssp_vec_norm(r[0]);
        if (solver.verb >= 1) {
            lssp_printf("bicgstabl: itr: %5d, abs res: %.6e, rel res: %.6e\n", iter, nrm2, nrm2 / ires);
        }

        if (nrm2 < tol) {
            pc.solve(&pc, t, x);
            lssp_vec_copy(x, t);
            lssp_vec_axpby(1, xp, 1, x);

            goto end;
        }
    }

end:
    solver.nits = iter;
    solver.residual = nrm2;

    lssp_vec_destroy(rtld);
    lssp_vec_destroy(xp);
    lssp_vec_destroy(bp);
    lssp_vec_destroy(t);

    for (i = 0; i <= l; i++) {
        lssp_vec_destroy(r[i]);
        lssp_vec_destroy(u[i]);
    }
    lssp_free(r);
    lssp_free(u);

    lssp_free(tau);

    return iter;
}
