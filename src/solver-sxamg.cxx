
#include "solver-sxamg.h"

#if USE_SXAMG

typedef struct SXAMG_DATA_
{
    SX_AMG_PARS pars;

} SXAMG_DATA;

void lssp_solver_sxamg_create(LSSP_SOLVER &s)
{
    s.sxamg = lssp_malloc<SXAMG_DATA>(1);
    sx_amg_pars_init(&s.sxamg->pars);
}

void lssp_solver_sxamg_destroy(LSSP_SOLVER &s)
{
    if (s.sxamg != NULL) {
        lssp_free(s.sxamg);
        s.sxamg = NULL;
    }
}

int lssp_solver_sxamg(LSSP_SOLVER *solver)
{
    lssp_mat_csr A = solver->A;
    SX_INT *Aj, *Ap, nnz, n, i;
    SX_FLT *Ax;
    SX_MAT XA;
    SX_VEC b, x;
    SX_RTN rtn;

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

    /* vector */
    b = sx_vec_create(n);
    x = sx_vec_create(n);

    for (i = 0; i < n; i++) {
        sx_vec_set_entry(&x, i, solver->x.d[i]);
        sx_vec_set_entry(&b, i, solver->rhs.d[i]);
    }

    /* pars */
    solver->sxamg->pars.maxit = solver->maxit;
    solver->sxamg->pars.verb = solver->verb;
    solver->sxamg->pars.tol = solver->tol_rel;

    if (solver->log != NULL) {
        sx_set_log(solver->log);
    }

    /* solve */
    rtn = sx_solver_amg(&XA, &b, &x, &solver->sxamg->pars);

    /* return solution */
    for (i = 0; i < n; i++) {
        solver->x.d[i] = sx_vec_get_entry(&x, i);
    }

    sx_mat_destroy(&XA);
    sx_vec_destroy(&b);
    sx_vec_destroy(&x);

    if (sizeof(int) != sizeof(SX_INT)) {
        lssp_free(Ap);
        lssp_free(Aj);
    }

    if (sizeof(double) != sizeof(SX_FLT)) {
        lssp_free(Ax);
    }

    /* abs res */
    solver->residual = rtn.ares;

    return rtn.nits;
}

void lssp_solver_sxamg_set_pars(LSSP_SOLVER *solver, SX_AMG_PARS *pars)
{
    if (solver == NULL) return;
    if (solver->sxamg == NULL) return;

    if (pars != NULL) solver->sxamg->pars = *pars;
}

#endif
