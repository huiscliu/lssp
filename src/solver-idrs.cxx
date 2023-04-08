
#include "solver-idrs.h"

static void idrs_orth(int s, lssp_vec *P)
{
    int i, j;
    double r;
    double d;

    for (j = 0; j < s; j++) {
        r = lssp_vec_norm(P[j]);
        r = 1.0 / r;

        lssp_vec_scale(P[j], r);

        for (i = j + 1; i < s; i++) {
            d = lssp_vec_dot(P[j], P[i]);
            lssp_vec_axpby(-d, P[j], 1, P[i]);
        }
    }
}

static int array_solve(int n, double *a, double *b, double *x, double *w)
{
    int i, j, k;
    double t;

    for (i = 0; i < n * n; i++)
        w[i] = a[i];

    switch (n) {
        case 1:
            x[0] = b[0] / w[0];

            break;

        case 2:
            w[0] = 1.0 / w[0];
            w[1] *= w[0];
            w[3] -= w[1] * w[2];
            w[3] = 1.0 / w[3];
            /* forward sub */
            x[0] = b[0];
            x[1] = b[1] - w[1] * x[0];
            /* backward sub */
            x[1] *= w[3];
            x[0] -= w[2] * x[1];
            x[0] *= w[0];

            break;

        default:
            for (k = 0; k < n; k++) {
                w[k + k * n] = 1.0 / w[k + k * n];
                for (i = k + 1; i < n; i++) {
                    t = w[i + k * n] * w[k + k * n];
                    for (j = k + 1; j < n; j++) {
                        w[i + j * n] -= t * w[k + j * n];
                    }
                    w[i + k * n] = t;
                }
            }

            /* forward sub */
            for (i = 0; i < n; i++) {
                x[i] = b[i];
                for (j = 0; j < i; j++) {
                    x[i] -= w[i + j * n] * x[j];
                }
            }

            /* backward sub */
            for (i = n - 1; i >= 0; i--) {
                for (j = i + 1; j < n; j++) {
                    x[i] -= w[i + j * n] * x[j];
                }
                x[i] *= w[i + i * n];
            }

            break;
    }

    return 0;
}

int lssp_solver_idrs(LSSP_SOLVER & solver, LSSP_PC & pc)
{
    lssp_mat_csr A = solver.A;
    lssp_vec x;
    lssp_vec r, t, v, av, *dX, *dR, *P;
    double om, h;
    double *M, *m, *c, *MM;
    double nrm2, tol, ires;
    int i, j, k, s, oldest;
    int iter, maxiter, n;

    x = solver.x;
    n = A.num_rows;
    maxiter = solver.maxit;
    s = solver.idrs;
    if (s <= 0) s = 4;

    r = lssp_vec_create(n);
    t = lssp_vec_create(n);
    v = lssp_vec_create(n);
    av = lssp_vec_create(n);

    dX = lssp_malloc<lssp_vec>(s);
    dR = lssp_malloc<lssp_vec>(s);
    P = lssp_malloc<lssp_vec>(s);

    for (i = 0; i < s; i++) {
        dX[i] = lssp_vec_create(n);
        dR[i] = lssp_vec_create(n);
        P[i] = lssp_vec_create(n);
    }

    m = lssp_malloc<double>(s);
    c = lssp_malloc<double>(s);
    M = lssp_malloc<double>(s * s);
    MM = lssp_malloc<double>(s * s);

    /* Initial Residual */
    iter = 0;
    lssp_mv_amxpbyz(-1, A, x, 1, solver.rhs, r);
    nrm2 = ires = lssp_vec_norm(r);

    if (nrm2 <= solver.tol_abs) {
        solver.nits = 0;
        goto end;
    }

    tol = nrm2 * solver.tol_rel;
    h = lssp_vec_norm(solver.rhs) * solver.tol_rb;
    if (tol < solver.tol_abs) tol = solver.tol_abs;
    if (tol < h) tol = h;

    /* init */
    srand(0);
    for (k = 0; k < s; k++) {
        for (i = 0; i < n; i++) {
            P[k].d[i] = (rand() * 1.) / (1. * RAND_MAX);
        }
    }

    idrs_orth(s, P);

    for (k = 0; k < s; k++) {
        pc.solve(&pc, dX[k], r);
        lssp_mv_mxy(A, dX[k], dR[k]);

        h = lssp_vec_dot(dR[k], dR[k]);
        om = lssp_vec_dot(dR[k], r);
        om = om / h;

        lssp_vec_scale(dX[k], om);
        lssp_vec_scale(dR[k], -om);

        lssp_vec_axpby(1, dX[k], 1, x);
        lssp_vec_axpby(1, dR[k], 1, r);

        /* convergence check */
        nrm2 = lssp_vec_norm(r);
        if (tol >= nrm2) {

            iter = k + 1;
            goto end;
        }

        for (i = 0; i < s; i++) {
            M[k * s + i] = lssp_vec_dot(P[i], dR[k]);
        }
    }

    iter = s;
    oldest = 0;
    for (i = 0; i < s; i++) {
        m[i] = lssp_vec_dot(P[i], r);
    }

    while (iter <= maxiter) {
        /* solve Mc=m */
        array_solve(s, M, m, c, MM);

        lssp_vec_copy(v, r);
        for (j = 0; j < s; j++) {
            lssp_vec_axpby(-c[j], dR[j], 1, v);
        }

        if ((iter % (s + 1)) == s) {
            pc.solve(&pc, av, v);
            lssp_mv_mxy(A, av, t);

            h = lssp_vec_dot(t, t);
            om = lssp_vec_dot(t, v);
            om = om / h;

            for (i = 0; i < n; i++) {
                h = om * av.d[i];
                for (j = 0; j < s; j++) {
                    h -= dX[j].d[i] * c[j];
                }

                dX[oldest].d[i] = h;
            }

            for (i = 0; i < n; i++) {
                h = -om * t.d[i];
                for (j = 0; j < s; j++) {
                    h -= dR[j].d[i] * c[j];
                }

                dR[oldest].d[i] = h;
            }
        }
        else {
            pc.solve(&pc, av, v);

            for (i = 0; i < n; i++) {
                h = om * av.d[i];
                for (j = 0; j < s; j++) {
                    h -= dX[j].d[i] * c[j];
                }

                dX[oldest].d[i] = h;
            }

            lssp_mv_mxy(A, dX[oldest], dR[oldest]);
            lssp_vec_scale(dR[oldest], -1.);
        }

        lssp_vec_axpby(1, dR[oldest], 1, r);
        lssp_vec_axpby(1, dX[oldest], 1, x);

        iter++;

        /* convergence check */
        nrm2 = lssp_vec_norm(r);

        if (solver.verb >= 1) {
            lssp_printf("idrs: itr: %5d, abs res: %.6e, rel res: %.6e\n", iter, nrm2, nrm2 / ires);
        }

        if (tol >= nrm2) {
            goto end;
        }

        for (i = 0; i < s; i++) {
            h = lssp_vec_dot(P[i], dR[oldest]);
            m[i] += h;
            M[oldest * s + i] = h;
        }

        oldest++;
        if (oldest == s) oldest = 0;
    }

end:

    solver.nits = iter;
    solver.residual = nrm2;

    lssp_vec_destroy(r);
    lssp_vec_destroy(t);
    lssp_vec_destroy(v);
    lssp_vec_destroy(av);

    for (i = 0; i < s; i++) {
        lssp_vec_destroy(dX[i]);
        lssp_vec_destroy(dR[i]);
        lssp_vec_destroy(P[i]);
    }

    lssp_free(P);
    lssp_free(dX);
    lssp_free(dR);
    lssp_free(m);
    lssp_free(c);
    lssp_free(M);
    lssp_free(MM);

    return iter;
}
