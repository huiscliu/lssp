
#include "pc-mi20.h"

#if USE_HSL_MI20

#include <stdbool.h>
extern "C" {
#include <hsl_mi20d.h>
}

typedef struct MI20_AMG_
{
    void *keep;
    struct mi20_control control;
    struct mi20_info info;

} MI20_AMG;

void lssp_pc_mi20amg_create(LSSP_PC &pc)
{
    MI20_AMG *d;

    /* create */
    pc.mi20 = lssp_malloc<MI20_AMG>(1);
    d = pc.mi20;

    /* pars */
    mi20_default_control(&d->control);
    d->control.f_arrays = 0;
    d->control.aggressive = 1;
    d->control.max_points = 1;
    d->control.st_parameter = 0.5;
    d->control.trunc_parameter = 1e-3;
    d->control.coarse_solver_its = 5;
    d->control.v_iterations = 1;
}

static void lssp_pc_mi20amg_destroy(LSSP_PC *pc)
{
    MI20_AMG *d;
    if (pc == NULL) return;
    if (pc->mi20 == NULL) return;

    d = (MI20_AMG *)pc->mi20;
    mi20_finalize(&d->keep, &d->control, &d->info);
    lssp_free(d);
    pc->mi20 = NULL;
}

static void amg_solve(LSSP_PC *s, lssp_vec x, lssp_vec rhs)
{
    MI20_AMG *d;

    assert(s != NULL);
    assert(s->mi20 != NULL);

    d = (MI20_AMG *)s->mi20;

    /* solve */
    mi20_precondition(rhs.d, x.d, &d->keep, &d->control, &d->info);
}

void lssp_pc_mi20amg_assemble(LSSP_PC &pc, LSSP_SOLVER s)
{
    lssp_mat_csr A = s.A;
    MI20_AMG *d;
    int n = A.num_rows;

    /* create */
    d = pc.mi20;

    /* setup */
    mi20_setup_csr(n, A.Ap, A.Aj, A.Ax, &d->keep, &d->control, &d->info);

    pc.solve = amg_solve;
    pc.destroy = lssp_pc_mi20amg_destroy;
}

#endif
