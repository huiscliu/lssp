
#include "solver-fasp.h"

#if USE_FASP

#define __STDC_VERSION__         199901L

typedef struct FASP_DATA_
{
    input_param inipar;         // parameters from input 
    itsolver_param itspar;      // parameters for itsolver
    AMG_param amgpar;           // parameters for AMG
    ILU_param ilupar;           // parameters for ILU
    Schwarz_param swzpar;       // parameters for Schwarz method

} FASP_DATA;

static short solver_type_intl[] = {
    SOLVER_CG,
    SOLVER_BiCGstab,
    SOLVER_MinRes,
    SOLVER_GMRES,
    SOLVER_VGMRES,
    SOLVER_VFGMRES,
    SOLVER_GCG,
    SOLVER_GCR,
    SOLVER_SCG,
    SOLVER_SBiCGstab,
    SOLVER_SMinRes,
    SOLVER_SGMRES,
    SOLVER_SVGMRES,
    SOLVER_SVFGMRES,
    SOLVER_SGCG,
    SOLVER_AMG,
    SOLVER_FMG,

};

static short pc_type[] = {
    PREC_NULL,
    PREC_DIAG,
    PREC_AMG,
    PREC_FMG,
    PREC_ILU,
    PREC_SCHWARZ,

};

static void lssp_solver_fasp_param_input_init(LSSP_SOLVER &s)
{
    short int solver_id = 1;
    short int precond_id = 4;
    int k;
    input_param *iniparam = &s.fasp->inipar;

    fasp_param_input_init(&s.fasp->inipar);
    fasp_param_init(&s.fasp->inipar, &s.fasp->itspar, &s.fasp->amgpar, &s.fasp->ilupar, &s.fasp->swzpar);

    k = sizeof(solver_type_intl) / sizeof(*solver_type_intl);
    if (s.fasp_solver >= 0 && s.fasp_solver < k) {
        solver_id = s.fasp_solver;
    }

    k = sizeof(pc_type) / sizeof(*pc_type);
    if (s.fasp_pc >= 0 && s.fasp_pc < k) {
        precond_id = s.fasp_pc;
    }

    /* hack default settings */
    iniparam->solver_type = solver_type_intl[solver_id];
    iniparam->precond_type = pc_type[precond_id];
    iniparam->itsolver_tol = s.tol_rel;
    iniparam->itsolver_maxit = s.maxit;
    iniparam->restart = s.restart;
}

void lssp_solver_fasp_create(LSSP_SOLVER &s)
{
    s.fasp = lssp_malloc<FASP_DATA>(1);
    lssp_solver_fasp_param_input_init(s);
}

void lssp_solver_fasp_destroy(LSSP_SOLVER &s)
{
    if (s.fasp != NULL) {
        lssp_free(s.fasp);
        s.fasp = NULL;
    }
}

static void lssp_solver_fasp_assemble(dCSRmat *A, dvector *b, dvector *x, LSSP_SOLVER *solver)
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

int lssp_solver_fasp(LSSP_SOLVER *solver)
{
    dCSRmat A;
    dvector b, x;
    int status = FASP_SUCCESS;
    int i;

    int print_level;
    int solver_type;
    int precond_type;
    double time = lssp_get_time();

    // set solver parameters
    solver->fasp->amgpar.maxit = solver->maxit;
    solver->fasp->amgpar.tol = solver->tol_rel;

    print_level = solver->fasp->inipar.print_level;
    solver_type = solver->fasp->inipar.solver_type;
    precond_type = solver->fasp->inipar.precond_type;

    // read A and b 
    lssp_solver_fasp_assemble(&A, &b, &x, solver);

    // print out solver parameters
    if (solver->verb <= 0) print_level = PRINT_NONE;

    if (print_level > PRINT_NONE) fasp_param_solver_print(&solver->fasp->itspar);

    // solve the system //
    fasp_dvec_alloc(A.row, &x);
    fasp_dvec_set(A.row, &x, 0.0);

    // preconditioned Krylov methods
    if (solver_type >= 1 && solver_type <= 20) {
        // Using no preconditioner for Krylov iterative methods
        if (precond_type == PREC_NULL) {
            status = fasp_solver_dcsr_krylov(&A, &b, &x, &solver->fasp->itspar);
        }
        // Using diag(A) as preconditioner for Krylov iterative methods
        else if (precond_type == PREC_DIAG) {
            status = fasp_solver_dcsr_krylov_diag(&A, &b, &x, &solver->fasp->itspar);
        }
        // Using AMG as preconditioner for Krylov iterative methods
        else if (precond_type == PREC_AMG || precond_type == PREC_FMG) {
            if (print_level > PRINT_NONE) fasp_param_amg_print(&solver->fasp->amgpar);
            status = fasp_solver_dcsr_krylov_amg(&A, &b, &x, &solver->fasp->itspar,
                    &solver->fasp->amgpar);
        }
        // Using ILU as preconditioner for Krylov iterative methods Q: Need to change!
        else if (precond_type == PREC_ILU) {
            if (print_level > PRINT_NONE) fasp_param_ilu_print(&solver->fasp->ilupar);
            status = fasp_solver_dcsr_krylov_ilu(&A, &b, &x, &solver->fasp->itspar,
                    &solver->fasp->ilupar);
        }
        // Using Schwarz as preconditioner for Krylov iterative methods
        else if (precond_type == PREC_SCHWARZ) {
            if (print_level > PRINT_NONE) fasp_param_Schwarz_print(&solver->fasp->swzpar);
            status = fasp_solver_dcsr_krylov_Schwarz(&A, &b, &x, &solver->fasp->itspar,
                    &solver->fasp->swzpar);
        }
        else {
            lssp_printf("lssp: fasp, wrong preconditioner type %d!!!\n", precond_type);
            status = ERROR_SOLVER_PRECTYPE;
        }

    }
    else if (solver_type == SOLVER_AMG) {       // AMG as the iterative solver
        if (print_level > PRINT_NONE) fasp_param_amg_print(&solver->fasp->amgpar);

        fasp_solver_amg(&A, &b, &x, &solver->fasp->amgpar);
    }
    else if (solver_type == SOLVER_FMG) {       // Full AMG as the iterative solver 
        if (print_level > PRINT_NONE) fasp_param_amg_print(&solver->fasp->amgpar);
        fasp_solver_famg(&A, &b, &x, &solver->fasp->amgpar);
    }
    else {
        lssp_printf("lssp: fasp, wrong solver type %d!!!\n", solver_type);
        status = ERROR_SOLVER_TYPE;
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

void lssp_solver_fasp_set_pars(LSSP_SOLVER *solver, input_param *inparam, itsolver_param *itsparam,
        AMG_param *amgparam, ILU_param *iluparam, Schwarz_param *schparam)
{
    if (solver == NULL) return;
    if (solver->fasp == NULL) return;

    if (inparam != NULL) solver->fasp->inipar = *inparam;
    if (itsparam != NULL) solver->fasp->itspar = *itsparam;
    if (amgparam != NULL) solver->fasp->amgpar = *amgparam;
    if (iluparam != NULL) solver->fasp->ilupar = *iluparam;
    if (schparam != NULL) solver->fasp->swzpar = *schparam;
}

#endif
