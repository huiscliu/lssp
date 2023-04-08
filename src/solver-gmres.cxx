
#include "solver-gmres.h"

extern int LSSP_RESTART;
extern double LSSP_BREAKDOWN;
extern int LSSP_MAXIT;
extern double LSSP_ATOL;
extern double LSSP_RTOL;
extern double LSSP_RB;


int lssp_solver_gmres(LSSP_SOLVER &solver, LSSP_PC &pc)
{
    lssp_mat_csr A = solver.A;
    int m = solver.restart;
    int itr_max = solver.maxit;
    lssp_vec rhs = solver.rhs;
    lssp_vec x = solver.x;
    double tol_abs = solver.tol_abs;
    double tol_rel = solver.tol_rel;
    double tol_rb = solver.tol_rb;
    double b_norm;
    double rtol;

    lssp_vec wj;
    lssp_vec *v;
    lssp_vec r;

    double **Hg;
    double *ym;
    double *gg;
    lssp_vec rg;
    double *c, *s;

    int itr_inner;
    int itr_outer;
    int i, j = 0;
    int n;
    double tol = 0, err_rel = 0;
    double beta;
    int break_down = 0;
    double time;
    double gstol = 0.;

    assert(solver.assembled);
    assert(pc.assembled);

    if (m < 0) m = LSSP_RESTART;
    if (itr_max <= 0) itr_max = LSSP_MAXIT;
    if (tol_abs < 0) tol_abs = LSSP_ATOL;
    if (tol_rel < 0) tol_rel = LSSP_RTOL;
    if (tol_rb < 0) tol_rb = LSSP_RB;

    assert(A.num_rows == A.num_cols);
    assert(A.num_nnzs >= A.num_rows);

    time = lssp_get_time();
    if (solver.verb >= 2) {
        lssp_printf("gmres: restart parameter m: %d\n", m);
        lssp_printf("gmres: maximal iteration: %d\n", itr_max);
        lssp_printf("gmres: tolerance abs: %g\n", tol_abs);
        lssp_printf("gmres: tolerance rel: %g\n", tol_rel);
        lssp_printf("gmres: tolerance rbn: %g\n", tol_rb);
    }

    n = A.num_rows;
    wj = lssp_vec_create(n);
    rg = lssp_vec_create(n);
    r = lssp_vec_create(n);

    v = lssp_malloc<lssp_vec>(m);
    for (i = 0; i < m; i++) {
        v[i] = lssp_vec_create(n);
    }

    gg = lssp_malloc<double>(m + 1);
    ym = lssp_malloc<double>(m);

    Hg = (double **)malloc((m + 1) * sizeof(double *));
    for (i = 0; i <= m; i++) {
        Hg[i] = (double *)malloc(m * sizeof(double));
    }

    c = (double *) malloc(m * sizeof(*c));
    s = (double *) malloc(m * sizeof(*s));

    b_norm = lssp_vec_norm(solver.rhs);
    tol_rb *= b_norm;

    lssp_mv_amxpbyz(-1, A, x, 1, rhs, rg);
    beta = lssp_vec_norm(rg);

    if (beta <= tol_abs) {
        itr_inner = 0;
        goto end;
    }

    err_rel = beta;
    tol = tol_rel * err_rel;

    if (tol < tol_abs) tol = tol_abs;
    if (tol < tol_rb) tol = tol_rb;
    rtol = tol / beta;

    for (i = 0; i < n; i++) {
        v[0].d[i] = wj.d[i] = 0;
    }

    itr_outer = itr_inner = 0;
    while (itr_inner < itr_max) {
        int kk;
        int gs_break = 0;
        double h1, h2, gma, gs_norm = 0.;

        lssp_vec_set_value(v[0], 0.);
        pc.solve(&pc, v[0], rg);

        beta = lssp_vec_norm(v[0]);
        gg[0] = beta;
        for (kk = 1; kk <= m; kk++) {
            gg[kk] = 0;
        }

        if (itr_outer == 0) {
            gstol = rtol * beta * 0.5;
        }

        for (kk = 0; kk <= m; kk++) {
            for (i = 0; i < m; i++) Hg[kk][i] = 0;
        }

        for (kk = 0; kk < n; kk++) {
            v[0].d[kk] /= beta;
        }

        for (i = 0; i < m; i++) {
            double hij;

            itr_inner++;

            lssp_mv_mxy(A, v[i], rg);
            lssp_vec_set_value(wj, 0.);
            pc.solve(&pc, wj, rg);

            for (j = 0; j <= i; j++) {
                hij = lssp_vec_dot(wj, v[j]);
                lssp_vec_axpby(-hij, v[j], 1, wj);

                Hg[j][i] = hij;
            }

            hij = lssp_vec_norm(wj);
            Hg[i + 1][i] = hij;

            if (fabs(hij) <= LSSP_BREAKDOWN) {
                break_down = 1;
            }
            else if (i + 1 < m) {
                for (kk = 0; kk < n; kk++) {
                    v[i + 1].d[kk] = wj.d[kk] / hij;
                }
            }

            for (j = 0; j < i; j++) {

                h1 = c[j] * Hg[j][i] + s[j] * Hg[j + 1][i];

                h2 = -s[j] * Hg[j][i] + c[j] * Hg[j + 1][i];

                Hg[j][i] = h1;
                Hg[j + 1][i] = h2;
            }

            gma = sqrt(Hg[i][i] * Hg[i][i] + Hg[i + 1][i] * Hg[i + 1][i]);
            if (fabs(gma) == 0.) gma = 1e-20;

            s[i] = Hg[i + 1][i] / gma;
            c[i] = Hg[i][i] / gma;

            h1 = c[i] * Hg[i][i] + s[i] * Hg[i + 1][i];

            h2 = -s[i] * Hg[i][i] + c[i] * Hg[i + 1][i];

            Hg[i][i] = h1;
            Hg[i + 1][i] = h2;

            h1 = c[i] * gg[i] + s[i] * gg[i + 1];
            h2 = -s[i] * gg[i] + c[i] * gg[i + 1];

            gg[i] = h1;
            gg[i + 1] = h2;
            gs_norm = fabs(h2);

            if (gs_norm <= gstol) {
                gs_break = 1;
            }

            if (gs_break || break_down) {
                goto solve;
            }
        }

solve:
        if (i == m) {
            kk = m;
        }
        else {
            kk = i;
        }

        for (i = kk - 1; i >= 0; i--) {
            double upper_temp;

            ym[i] = gg[i] / Hg[i][i];
            for (j = 0; j < i; j++) {
                upper_temp = gg[j] - ym[i] * Hg[j][i];
                gg[j] = upper_temp;
            }
        }

        for (j = 0; j < n; j++) {
            double upper_temp = 0;

            for (i = 0; i < kk; i++) {
                upper_temp += v[i].d[j] * ym[i];
            }

            x.d[j] += upper_temp;
        }

        lssp_mv_amxpbyz(-1, A, x, 1, rhs, rg);
        beta = lssp_vec_norm(rg);

        if (solver.verb >= 1) {
            lssp_printf("gmres: itr: %4d / %5d, abs res: %.6e, rel res: %.6e, rbn: %.6e\n", \
                    itr_outer, itr_inner, beta, (err_rel == 0 ? 0 : beta / err_rel), \
                    (b_norm == 0 ? 0 : beta / b_norm));
        }

        if (beta <= tol) break;

        if (break_down) break;

        gstol = rtol * gs_norm / (beta / err_rel) * 0.5;

        itr_outer++;
    }

end:

    solver.residual = beta;
    solver.nits = itr_inner;

    lssp_vec_destroy(wj);
    lssp_vec_destroy(rg);
    lssp_vec_destroy(r);

    for (i = 0; i < m; i++) {
        lssp_vec_destroy(v[i]);
    }
    lssp_free(v);

    free(c);
    free(s);
    lssp_free(gg);
    lssp_free(ym);

    for (i = 0; i <= m; i++) {
        free(Hg[i]);
    }
    free(Hg);

    time = lssp_get_time() - time;

    if (solver.verb >= 2) {
        lssp_printf("gmres: total iteration: %d\n", itr_inner);
        lssp_printf("gmres: total time: %g\n", time);
    }

    return itr_inner;
}
