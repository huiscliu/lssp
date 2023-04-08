
#include "solver-mi20.h"

#if USE_HSL_MI20

#include <stdbool.h>
extern "C" {
#include <hsl_mi20d.h>
}

int lssp_solver_mi20amg(LSSP_SOLVER *solver)
{
    lssp_mat_csr A = solver->A;
    int n = A.num_rows;
    void *keep;
    struct mi20_control control;
    struct mi20_solve_control solve_control;
    struct mi20_info info;
    double time = lssp_get_time();

    /* pars */
    mi20_default_control(&control);
    control.f_arrays = 0;
    control.aggressive = 1;
    control.max_points = 1;
    control.st_parameter = 0.5;
    control.trunc_parameter = 1e-3;
    control.coarse_solver_its = 5;
    control.v_iterations = 1;

    mi20_default_solve_control(&solve_control);
    solve_control.abs_tol = solver->tol_abs;
    solve_control.max_its = solver->maxit;
    solve_control.init_guess = 1;

    /* setup */
    mi20_setup_csr(n, A.Ap, A.Aj, A.Ax, &keep, &control, &info);

    /* solve */
    mi20_solve(solver->rhs.d, solver->x.d, &keep, &control, &solve_control, &info);

    /* return */
    solver->nits = info.iterations;
    solver->residual = info.residual;

    /* clean up */
    mi20_finalize(&keep, &control, &info);

    time = lssp_get_time() - time;
    if (solver->verb >= 2) {
        lssp_printf("lssp: mi20amg, total time: %g\n", time);
        lssp_printf("lssp: mi20amg, nits: %d residual: %g\n", info.iterations, info.residual);
    }

    return info.iterations;
}

#endif
