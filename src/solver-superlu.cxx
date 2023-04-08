
#include "solver-superlu.h"

#if USE_SUPERLU

#include <slu_ddefs.h>

typedef struct SUPERLU_DATA_
{
    superlu_options_t options;

} SUPERLU_DATA;

void lssp_solver_superlu_create(LSSP_SOLVER &s)
{
    s.superlu = lssp_malloc<SUPERLU_DATA>(1);
    set_default_options(&s.superlu->options);
}

void lssp_solver_superlu_destroy(LSSP_SOLVER &s)
{
    if (s.superlu == NULL) return;

    lssp_free(s.superlu);
    s.superlu = NULL;
}

int lssp_solver_superlu(LSSP_SOLVER *solver)
{
    SuperMatrix A, L, U, B;
    lssp_mat_csr Ai = solver->A;
    int *perm_r;
    int *perm_c, i;
    int info;
    int m = Ai.num_rows, n = Ai.num_cols, nnz = Ai.num_nnzs;
    SuperLUStat_t stat;
    superlu_options_t options;
    double time0 = lssp_get_time();
    double *BB;
    int *Aj, *Ap;
    double *b, *Ax;

    /* copy */
    Aj = lssp_copy_on<int>(Ai.Aj, nnz);
    Ap = lssp_copy_on<int>(Ai.Ap, n + 1);
    Ax = lssp_copy_on<double>(Ai.Ax, nnz);
    b = lssp_copy_on<double>(solver->rhs.d, n);

    /* matrix A */
    assert(m == n);
    dCreate_CompCol_Matrix(&A, m, n, nnz, Ax, Aj, Ap, SLU_NR, SLU_D, SLU_GE);

    /* right-hand side */
    dCreate_Dense_Matrix(&B, m, 1, b, m, SLU_DN, SLU_D, SLU_GE);

    if (!(perm_r = intMalloc(m))) ABORT("Malloc fails for perm_r[].");
    if (!(perm_c = intMalloc(n))) ABORT("Malloc fails for perm_c[].");

    /* Set the default input options. */
    options.ColPerm = COLAMD;
    StatInit(&stat);

    /* SuperLU */
    dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

    BB = (double *)((DNformat *) B.Store)->nzval;
    for (i = 0; i < n; i++) solver->x.d[i] = BB[i];

    SUPERLU_FREE(perm_r);
    SUPERLU_FREE(perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    StatFree(&stat);
    lssp_free(b);

    if (solver->verb > 0) {
        lssp_printf("superlu: solve: %lg s\n", lssp_get_time() - time0);
    }

    solver->residual = 0.0;
    return 1; 

}

#endif
