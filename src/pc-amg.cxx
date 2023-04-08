
#include "pc-amg.h"

#if USE_FASP

static short pc_type[] =
{
    PREC_AMG,
    PREC_FMG,
};

typedef struct FPAMG_DATA_
{
    AMG_param  amgpar;
    precond    *amg;
    REAL       *x;
    REAL       *b;

} FPAMG_DATA;

void lssp_pc_amg_create(LSSP_PC &pc)
{
    pc.fasp_amg = (FPAMG_DATA *)lssp_malloc<FPAMG_DATA>(1);
    fasp_param_amg_init(&pc.fasp_amg->amgpar);
}

static void lssp_pc_amg_destroy(LSSP_PC *pc)
{
    FPAMG_DATA *d;
    short precond_id = 0;

    if (pc == NULL) return;
    if (pc->fasp_amg == NULL) return;

    if (pc->type == LSSP_PC_FASP_FMG) precond_id = 1;
    d = (FPAMG_DATA *)pc->fasp_amg;
    fasp_precond_free(pc_type[precond_id], d->amg);
    lssp_free(d->x);
    lssp_free(d->b);
    lssp_free(d);
    pc->fasp_amg = NULL;
}

static void amg_solve(LSSP_PC *s, lssp_vec x, lssp_vec rhs)
{
    int n = s->A.num_rows; 
    FPAMG_DATA *d;

    assert(s != NULL);
    assert(s->fasp_amg != NULL);

    d = (FPAMG_DATA *)s->fasp_amg;

    /* init */
    if (sizeof(double) == sizeof(REAL)) {
        memcpy(d->x, x.d, n * sizeof(REAL));
        memcpy(d->b, rhs.d, n * sizeof(REAL));
    }
    else {
        for (int i = 0; i < n; i++) {
            d->x[i] = x.d[i];
            d->b[i] = rhs.d[i];
        }
    }

    /* solve */
    d->amg->fct(d->b, d->x, d->amg->data);

    /* return */
    if (sizeof(double) == sizeof(REAL)) {
        memcpy(x.d, d->x, n * sizeof(REAL));
    }
    else {
        for (int i = 0; i < n; i++) {
            x.d[i] = d->x[i];
        }
    }
}

static void lssp_amg_mat_assemble(dCSRmat *A, lssp_mat_csr Ai)
{
    int i, n;

    n = Ai.num_rows;
    A->row = n;
    A->col = n;
    A->IA  = (INT *)fasp_mem_calloc(n+1, sizeof(INT));

    for ( i = 0; i <= n; ++i ) {
        A->IA[i] = Ai.Ap[i];
    }

    INT nz = A->IA[n];
    A->nnz = nz;
    A->JA  = (INT *) fasp_mem_calloc(nz, sizeof(INT));
    A->val = (REAL *)fasp_mem_calloc(nz, sizeof(REAL));

    for ( i = 0; i < nz; ++i ) {
        A->JA[i] = Ai.Aj[i];
        A->val[i] = Ai.Ax[i];
    }
}

void lssp_pc_amg_assemble(LSSP_PC &pc, LSSP_SOLVER s)
{
    lssp_mat_csr A = s.A;
    FPAMG_DATA *d;
    dCSRmat Af;
    int n = A.num_rows;
    short precond_id = 0;  /* 0: AMG, 1: FMG */

    /* mat */
    lssp_amg_mat_assemble(&Af, A);

    /* convert */
    d = pc.fasp_amg;

    /* default par */
    if (pc.type == LSSP_PC_FASP_FMG) precond_id = 1;
    d->amg = fasp_precond_setup(pc_type[precond_id], &d->amgpar, NULL, &Af);
    fasp_dcsr_free(&Af);

    /* b and x */
    d->b = (REAL *)fasp_mem_calloc(n, sizeof(REAL));
    d->x = (REAL *)fasp_mem_calloc(n, sizeof(REAL));

    pc.solve = amg_solve;
    pc.destroy = lssp_pc_amg_destroy;
}

void lssp_pc_amg_set_pars(LSSP_PC &pc, AMG_param *amgparam)
{
    if (amgparam != NULL && pc.fasp_amg != NULL) pc.fasp_amg->amgpar = *amgparam;
}

#endif
