
#include "solver-amg.h"

#if USE_FASP

typedef struct FASP_AMG_DATA_
{
    AMG_param pars;

} FASP_AMG_DATA;

void lssp_solver_amg_create(LSSP_SOLVER &s)
{
    s.fasp_amg = lssp_malloc<FASP_AMG_DATA>(1);
    fasp_param_amg_init(&s.fasp_amg->pars);
}

void lssp_solver_amg_destroy(LSSP_SOLVER &s)
{
    if (s.fasp_amg != NULL) {
        lssp_free(s.fasp_amg);
        s.fasp_amg = NULL;
    }
}

/* algebraic multi-grid solver */
static void lssp_solver_amg_assemble(dCSRmat *A, dvector *b, dvector *x, LSSP_SOLVER *solver)
{
    int i, n;
    
    /* mat */
    n = solver->A.num_rows;
    A->row = n;
    A->col = n;
    A->IA  = (INT *)fasp_mem_calloc(n+1, sizeof(INT));
    
    for ( i = 0; i <= n; ++i ) {
        A->IA[i] = solver->A.Ap[i];
    }
    
    INT nz = A->IA[n];
    A->nnz = nz;
    A->JA  = (INT *) fasp_mem_calloc(nz, sizeof(INT));
    A->val = (REAL *)fasp_mem_calloc(nz, sizeof(REAL));
    
    for ( i = 0; i < nz; ++i ) {
        A->JA[i] = solver->A.Aj[i];
        A->val[i] = solver->A.Ax[i];
    }
    
    /* b and x */
    b->row = n;
    b->val = (REAL *)fasp_mem_calloc(n, sizeof(REAL));

    x->row = n;
    x->val = (REAL *)fasp_mem_calloc(n, sizeof(REAL));
    
    for ( i = 0; i < n; ++i ) {
        b->val[i] = solver->rhs.d[i];
        x->val[i] = solver->x.d[i];
    }
}

int lssp_solver_amg(LSSP_SOLVER *solver)
{
    dCSRmat A;
    dvector b, x;
    int status = FASP_SUCCESS;
    int i;

    int print_level;
    double time = lssp_get_time();
    AMG_param solver_amg_pars_ = solver->fasp_amg->pars;

    // set solver parameters
    solver_amg_pars_.maxit = solver->maxit;
    solver_amg_pars_.tol = solver->tol_rel;

    print_level = PRINT_SOME;
    if (solver->verb <= 0) print_level = PRINT_NONE;

    // read A and b 
    lssp_solver_amg_assemble(&A, &b, &x, solver);

    // print out solver parameters
    if (solver->verb <= 0) print_level = PRINT_NONE;

    // solve the system //
    fasp_dvec_alloc(A.row, &x);
    fasp_dvec_set(A.row, &x, 0.0);

    if (solver->type == LSSP_SOLVER_FASP_AMG) {      // AMG as the iterative solver
        if (print_level > PRINT_NONE) fasp_param_amg_print(&solver_amg_pars_);

        fasp_solver_amg(&A, &b, &x, &solver_amg_pars_);
    }
    else if (solver->type == LSSP_SOLVER_FASP_FMG) { // Full AMG as the iterative solver 
        if (print_level > PRINT_NONE) fasp_param_amg_print(&solver_amg_pars_);
        fasp_solver_famg(&A, &b, &x, &solver_amg_pars_);
    }
    else {
        lssp_printf("lssp: fasp, wrong solver type!!!\n");
        status = ERROR_SOLVER_TYPE;
        return status;
    }

    if (status < 0) {
        lssp_printf("lssp: fasp, solver failed! Exit status = %d.\n\n", status);
    }

    for ( i = 0; i < b.row; ++i ) solver->x.d[i] = x.val[i];

    time = lssp_get_time() - time;
    if (solver->verb >= 2) {
        lssp_printf("lssp: fasp, total time: %g\n", time);
    }

    fasp_dcsr_free(&A);
    fasp_dvec_free(&b);
    fasp_dvec_free(&x);

    return status;
}

void lssp_solver_amg_reset_pars(LSSP_SOLVER &s, AMG_param par)
{
    if (s.fasp_amg != NULL) {
        s.fasp_amg->pars = par;
    }
}
#endif
