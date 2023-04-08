
#include "pc-sxamg.h"

#if USE_SXAMG

typedef struct SX_DATA_
{
    SX_AMG_PARS pars;
    SX_AMG mg;

} SX_DATA;

void lssp_pc_sxamg_create(LSSP_PC &pc)
{
    SX_DATA *d;

    /* create */
    pc.sxamg = lssp_malloc<SX_DATA>(1);
    d = pc.sxamg;

    /* pars */
    sx_amg_pars_init(&d->pars);
    d->pars.maxit = 1;
    d->pars.verb = pc.verb;
}

static void lssp_pc_sxamg_destroy(LSSP_PC *pc)
{
    SX_DATA *d;

    if (pc == NULL) return;
    if (pc->sxamg == NULL) return;

    d = (SX_DATA *)pc->sxamg;

    sx_amg_data_destroy(&d->mg);

    lssp_free(d);
    pc->sxamg = NULL;
}

static void amg_solve(LSSP_PC *s, lssp_vec x, lssp_vec rhs)
{
    SX_DATA *d;
    SX_VEC xx, xb;
    SX_INT n, i;

    assert(s != NULL);
    assert(s->sxamg != NULL);

    d = (SX_DATA *)s->sxamg;

    /* solve */
    n = x.n;
    xx = sx_vec_create(n);
    xb = sx_vec_create(n);

    for (i = 0; i < n; i++) {
        sx_vec_set_entry(&xx, i, x.d[i]);
        sx_vec_set_entry(&xb, i, rhs.d[i]);
    }

    /* solve */
    sx_solver_amg_solve(&d->mg, &xx, &xb);

    /* write back */
    for (i = 0; i < n; i++) {
        x.d[i] = sx_vec_get_entry(&xx, i);
    }

    sx_vec_destroy(&xx);
    sx_vec_destroy(&xb);
}

void lssp_pc_sxamg_assemble(LSSP_PC &pc, LSSP_SOLVER s)
{
    lssp_mat_csr A = s.A;
    SX_DATA *d;
    int n = A.num_rows;
    SX_INT *Aj, *Ap, nnz;
    SX_FLT *Ax;
    SX_MAT XA;

    /* create */
    d = pc.sxamg;

    /* matrix */
    n = A.num_rows;
    nnz = A.num_nnzs;

    if (sizeof(int) == sizeof(SX_INT)) {
        Ap = (SX_INT *) A.Ap;
        Aj = (SX_INT *) A.Aj;
    }
    else {
        Ap = lssp_malloc<SX_INT>(n + 1);
        Aj = lssp_malloc<SX_INT>(nnz);
    }

    if (sizeof(double) == sizeof(SX_FLT)) {
        Ax = (SX_FLT *) A.Ax;
    }
    else {
        Ax = lssp_malloc<SX_FLT>(nnz);
    }

    XA = sx_mat_create(n, n, Ap, Aj, Ax);

    /* setup */
    sx_amg_setup(&d->mg, &XA, &d->pars);

    pc.solve = amg_solve;
    pc.destroy = lssp_pc_sxamg_destroy;

    /* destroy */
    sx_mat_destroy(&XA);

    if (sizeof(int) != sizeof(SX_INT)) {
        lssp_free(Ap);
        lssp_free(Aj);
    }

    if (sizeof(double) != sizeof(SX_FLT)) {
        lssp_free(Ax);
    }
}

void lssp_pc_sxamg_set_pars(LSSP_PC &pc, SX_AMG_PARS *pars)
{
    if (pars != NULL && pc.sxamg != NULL) pc.sxamg->pars = *pars;
}

#endif
