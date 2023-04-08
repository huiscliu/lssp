
#include "solver-qmrcgstab.h"

extern int LSSP_MAXIT;
extern double LSSP_ATOL;
extern double LSSP_RTOL;
extern double LSSP_RB;

int lssp_solver_qmrcgstab(LSSP_SOLVER &solver, LSSP_PC &pc)
{
    lssp_mat_csr A = solver.A;
    int itr_max = solver.maxit;
    lssp_vec bg = solver.rhs;
    lssp_vec xk = solver.x;
    double tol_abs = solver.tol_abs;
    double tol_rel = solver.tol_rel;
    double tol_rb = solver.tol_rb;
    double b_norm;

    lssp_vec rk;
    lssp_vec br0;
    lssp_vec pk;
    lssp_vec vk;
    lssp_vec sk;
    lssp_vec dk;
    lssp_vec tk;
    lssp_vec bdk;
    lssp_vec bxk;
    lssp_vec r;

    double rho = 1,  prho, alpha = 1, beta;
    double omega = 1, theta = 0., btheta, b_eta, eta = 0., tau, btau;
    double residual, c, err_rel, ires, rerror;
    double tol = -1;
    int itr_out;
    int n;
    double time;

    assert(solver.assembled);
    assert(pc.assembled);

    if (itr_max <= 0) itr_max = LSSP_MAXIT;
    if (tol_abs < 0) tol_abs = LSSP_ATOL;
    if (tol_rel < 0) tol_rel = LSSP_RTOL;

    assert(A.num_rows == A.num_cols);

    time = lssp_get_time();
    if (solver.verb >= 2) {
        lssp_printf("qmrcgstab: maximal iteration: %d\n", itr_max);
        lssp_printf("qmrcgstab: tolerance abs: %g\n", tol_abs);
        lssp_printf("qmrcgstab: tolerance rel: %g\n", tol_rel);
        lssp_printf("qmrcgstab: tolerance rbn: %g\n", tol_rb);
    }

    n = A.num_rows;
    rk = lssp_vec_create(n);
    r = lssp_vec_create(n);
    br0 = lssp_vec_create(n);
    pk = lssp_vec_create(n);
    vk = lssp_vec_create(n);
    sk = lssp_vec_create(n);
    dk = lssp_vec_create(n);
    tk = lssp_vec_create(n);
    bdk = lssp_vec_create(n);
    bxk = lssp_vec_create(n);

    b_norm = lssp_vec_norm(solver.rhs);
    tol_rb *= b_norm;

    lssp_mv_amxpbyz(-1, A, xk, 1, bg, tk);
    residual = err_rel = lssp_vec_norm(tk);
    if (residual <= tol_abs) {
        itr_out = 0;
        goto end;
    }

    tol = residual * tol_rel;
    if (tol < tol_abs) tol = tol_abs;
    if (tol < tol_rb) tol = tol_rb;
    tol = tol / residual;

    /* init */
    pc.solve(&pc, rk, tk);
    lssp_vec_copy(br0, rk);

    lssp_vec_set_value(pk, 0);
    lssp_vec_set_value(dk, 0);
    lssp_vec_set_value(vk, 0);
    tau = ires = lssp_vec_norm(rk);

    /* main loop */
    prho = rho;
    for (itr_out = 0; itr_out < itr_max; itr_out++) {
        rho = lssp_vec_dot(br0, rk);
        beta = rho * alpha / prho / omega;
        prho = rho;

        /* pk */
        lssp_vec_axpbyz(1, pk, -omega, vk, r);
        lssp_vec_axpbyz(beta , r, 1, rk, pk);

        /* vk = A pk */
        lssp_mv_mxy(A, pk, r);
        pc.solve(&pc, vk, r);

        /* sk */
        alpha = rho / lssp_vec_dot(br0, vk);
        lssp_vec_axpbyz(-alpha, vk, 1, rk, sk);

        /* firstt quasi-mini */
        btheta = lssp_vec_norm(sk) / tau;
        c = 1 / sqrt(1. + btheta * btheta);
        btau = tau * btheta * c;

        /* bdk */
        b_eta = c * c * alpha;
        lssp_vec_axpbyz(1., pk, theta * theta * eta / alpha, dk, bdk);

        /* bxk */
        lssp_vec_axpbyz(1, xk, b_eta, bdk, bxk);

        /* tk = A sk */
        lssp_mv_mxy(A, sk, r);
        pc.solve(&pc, tk, r);

        /* rk */
#if 1
        omega = lssp_vec_dot(sk, tk) / lssp_vec_dot(tk, tk);
#else
        omega = lssp_vec_dot(sk, sk) / lssp_vec_dot(tk, sk);
#endif
        lssp_vec_axpbyz(1., sk, -omega, tk, rk);

        /* secod quasi-mini */
        theta = lssp_vec_norm(rk) / btau;
        c = 1. / sqrt(1. + theta * theta);
        tau = btau * theta * c;
        eta = c * c * omega;

        /* dk */
        lssp_vec_axpbyz(1, sk, btheta * btheta * b_eta / omega, bdk, dk);

        /* xk */
        lssp_vec_axpbyz(1, bxk, eta, dk, xk);

        /* relative error */
        rerror = lssp_vec_norm(rk) / ires;
        if (solver.verb >= 1) {
            lssp_printf("qmrcgstab: itr: %4d, rel res: %.6e\n", itr_out, rerror);
        }

        if (rerror <= tol) {
            lssp_mv_amxpbyz(-1, A, xk, 1, bg, tk);
            residual = lssp_vec_norm(tk);
            break;
        }
    }

    if (itr_out < itr_max) itr_out += 1;

end:

    solver.residual = residual;
    solver.nits = itr_out;

    lssp_vec_destroy(rk);
    lssp_vec_destroy(r);
    lssp_vec_destroy(br0);
    lssp_vec_destroy(pk);
    lssp_vec_destroy(vk);
    lssp_vec_destroy(sk);
    lssp_vec_destroy(dk);
    lssp_vec_destroy(tk);
    lssp_vec_destroy(bdk);
    lssp_vec_destroy(bxk);

    time = lssp_get_time() - time;
    if (solver.verb >= 2) {
        lssp_printf("qmrcgstab: total iteration: %d\n", itr_out);
        lssp_printf("qmrcgstab: total time: %g\n", time);
        lssp_printf("qmrcgstab: absolute error (residual): %g\n", residual);
    }

    return itr_out;
}
