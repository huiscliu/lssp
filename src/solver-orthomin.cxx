
#include "solver-orthomin.h"

extern int LSSP_RESTART;
extern double LSSP_BREAKDOWN;
extern int LSSP_MAXIT;
extern double LSSP_ATOL;
extern double LSSP_RTOL;
extern double LSSP_RB;


int lssp_solver_orthomin(LSSP_SOLVER & solver, LSSP_PC & pc)
{
    lssp_mat_csr A = solver.A;
    int k = solver.restart;
    int itr_max = solver.maxit;
    lssp_vec rhs = solver.rhs;
    lssp_vec x = solver.x;
    double tol_abs = solver.tol_abs;
    double tol_rel = solver.tol_rel;
    double tol_rb = solver.tol_rb;
    double b_norm;
    lssp_vec z;
    lssp_vec *q;
    lssp_vec *p;
    lssp_vec r;
    lssp_vec s;
    lssp_vec sd;
    double a_j;
    double *b_j;
    double *c_j;

    int itr_out;
    int i, j = 0;
    int n;
    double tol = -1, err_rel = 0;
    double beta;
    double time;

    assert(solver.assembled);
    assert(pc.assembled);

    if (k < 0) k = LSSP_RESTART;
    if (itr_max <= 0) itr_max = LSSP_MAXIT;
    if (tol_abs < 0) tol_abs = LSSP_ATOL;
    if (tol_rel < 0) tol_rel = LSSP_RTOL;
    if (tol_rb < 0) tol_rb = LSSP_RB;

    assert(A.num_rows == A.num_cols);

    time = lssp_get_time();
    if (solver.verb >= 2) {
        lssp_printf("orthomin: k parameter m: %d\n", k);
        lssp_printf("orthomin: maximal iteration: %d\n", itr_max);
        lssp_printf("orthomin: tolerance abs: %g\n", tol_abs);
        lssp_printf("orthomin: tolerance rel: %g\n", tol_rel);
        lssp_printf("orthomin: tolerance rbn: %g\n", tol_rb);
    }

    n = A.num_rows;
    z = lssp_vec_create(n);
    r = lssp_vec_create(n);
    s = lssp_vec_create(n);
    b_j = lssp_malloc<double>(k);
    c_j = lssp_malloc<double>(k);
    q = lssp_malloc<lssp_vec >(k);
    p = lssp_malloc<lssp_vec >(k);
    sd = lssp_vec_create(n);

    for (i = 0; i < k; i++) {
        p[i] = lssp_vec_create(n);
        q[i] = lssp_vec_create(n);

        for (itr_out = 0; itr_out < n; itr_out++) q[i].d[itr_out] = 0.;
    }

    lssp_mv_amxpbyz(-1, A, x, 1, rhs, z);
    lssp_vec_set_value(r, 0.);
    pc.solve(&pc, r, z);

    lssp_vec_copy(p[0], r);
    lssp_vec_copy(sd, r);

    b_norm = lssp_vec_norm(solver.rhs);
    tol_rb *= b_norm;

    beta = lssp_vec_norm(z);
    if (beta <= tol_abs) {
        itr_out = 0;
        goto end;
    }

    err_rel = beta;
    tol = tol_rel * err_rel;
    if (tol < tol_abs) tol = tol_abs;
    if (tol < tol_rb) tol = tol_rb;

    for (itr_out = 0; itr_out < itr_max; itr_out++) {

        lssp_mv_mxy(A, sd, s);

        j = itr_out % k;
        lssp_vec_set_value(q[j], 0.);
        pc.solve(&pc, q[j], s);

        a_j = lssp_vec_dot(r, q[j]);
        c_j[j] = lssp_vec_dot(q[j], q[j]);

        if (fabs(c_j[j]) <= LSSP_BREAKDOWN) break;

        a_j = a_j / c_j[j];
        lssp_vec_axpby(a_j, p[j], 1, x);
        lssp_vec_axpby(-a_j, q[j], 1, r);
        lssp_vec_copy(sd, r);
        lssp_mv_mxy(A, r, s);

        lssp_vec_set_value(z, 0.);
        pc.solve(&pc, z, s);

        if (itr_out >= k - 1) {
            for (i = 0; i < k; i++) {
                beta = lssp_vec_dot(z, q[i]);
                b_j[i] = -beta / c_j[i];

                lssp_vec_axpby(b_j[i], p[i], 1, sd);
            }
        }
        else {
            for (i = 0; i <= itr_out; i++) {

                beta = lssp_vec_dot(z, q[i]);
                b_j[i] = -beta / c_j[i];

                lssp_vec_axpby(b_j[i], p[i], 1, sd);
            }
        }

        j = (itr_out + 1) % k;
        lssp_vec_copy(p[j], sd);
        lssp_mv_amxpbyz(-1, A, x, 1, rhs, z);
        beta = lssp_vec_norm(z);

        if (solver.verb >= 1) {
            lssp_printf("orthomin: itr: %5d, abs res: %.6e, rel res: %.6e, "
                 "rbn: %.6e\n", itr_out, beta, (err_rel == 0 ? 0 : beta / err_rel),
                 (b_norm == 0 ? 0 : beta / b_norm));
        }

        if (beta <= tol) break;
    }

    if (itr_out < itr_max) itr_out += 1;

 end:

    solver.residual = beta;
    solver.nits = itr_out;

    for (i = 0; i < k; i++) {
        lssp_vec_destroy(p[i]);
        lssp_vec_destroy(q[i]);
    }

    lssp_vec_destroy(sd);
    lssp_vec_destroy(z);
    lssp_vec_destroy(r);
    lssp_vec_destroy(s);
    lssp_free(b_j);
    lssp_free(c_j);
    lssp_free(p);
    lssp_free(q);

    time = lssp_get_time() - time;
    if (solver.verb >= 2) {
        lssp_printf("orthomin: total iteration: %d\n", itr_out);
        lssp_printf("orthomin: total time: %g\n", time);
    }

    return itr_out;
}
