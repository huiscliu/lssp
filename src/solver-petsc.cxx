
#include "solver-petsc.h"

#if USE_PETSC

#include <petscksp.h>
#include <petscmat.h>

typedef struct PTS_DATA_
{
    PC        pc;
    KSP       ksp;
    Mat       A;
    Vec       b;
    Vec       x;
    bool      assembled;

} PTS_DATA;

static solver_petsc_setting solver_petsc_setting_func = NULL;

/* parameters */
static KSPType petsc_solver_list[] = 
{
    KSPGMRES, KSPBCGS, KSPFGMRES, KSPCG, KSPLGMRES, KSPTFQMR,
    KSPBICG, KSPMINRES, KSPLSQR
};

static PCType petsc_pc_list[] = 
{
    PCBJACOBI, PCILU, PCGAMG, PCLU
};

static int petsc_iluk_level = 1;
static int petsc_pc_right = 0;

#define NPETSolvers   (sizeof(petsc_solver_list)/sizeof(petsc_solver_list[0]))
#define NPETPCs       (sizeof(petsc_pc_list)/sizeof(petsc_pc_list[0]))

static void solver_petsc_initialize(void)
{
    int pargc = 1;
    char **pargv;

    pargv = lssp_malloc<char *>(1);
    *pargv = const_cast<char*>("petsc");
    PetscInitialize(&pargc, &pargv, NULL, NULL);
    lssp_free(pargv);
}

static void solver_petsc_finalize(void)
{
    PetscFinalize();
}

static int petsc_create(PTS_DATA *pts, LSSP_SOLVER *solver)
{
    lssp_mat_csr A = solver->A;
    PetscInt N = A.num_rows, Nlocal = A.num_rows;
    PetscInt i, prealloc, *d_nnz, *o_nnz;
    int petsc_solver_default = solver->petsc_solver;
    int petsc_pc_default = solver->petsc_pc;

    /* create data */
    pts->b = PETSC_NULL;
    pts->x = PETSC_NULL;
    pts->A = PETSC_NULL;
    pts->ksp = PETSC_NULL;
    pts->assembled = false;

    VecCreate(MPI_COMM_SELF, &pts->x);
    VecSetSizes(pts->x, Nlocal, N);
    VecSetFromOptions(pts->x);
    VecDuplicate(pts->x, &pts->b);

    /* get CSR matrix structure */
    prealloc = PETSC_DECIDE;
    d_nnz = lssp_malloc<PetscInt>(Nlocal);
    o_nnz = lssp_malloc<PetscInt>(Nlocal);

    for (i = 0; i < Nlocal; i++) {
        d_nnz[i] = A.Ap[i + 1] - A.Ap[i];
        o_nnz[i] = 0;
    }

    /* mat create */
    MatCreateAIJ(MPI_COMM_SELF, Nlocal, Nlocal, N, N, prealloc, d_nnz, prealloc, o_nnz, &pts->A);
    lssp_free(d_nnz);
    lssp_free(o_nnz);
    MatSetFromOptions(pts->A);

    /* create solver and pc */
    KSPCreate(MPI_COMM_SELF, &pts->ksp);
    KSPSetOperators(pts->ksp, pts->A, pts->A);

    KSPSetInitialGuessNonzero(pts->ksp, PETSC_TRUE);
    KSPSetFromOptions(pts->ksp);
    KSPGMRESSetRestart(pts->ksp, solver->restart);

    KSPGetPC(pts->ksp, &pts->pc);

    /* settings */
    i = sizeof(petsc_solver_list) / sizeof(*petsc_solver_list);
    if (petsc_solver_default < 0 || petsc_solver_default >= i) {
        petsc_solver_default = 1;
    }

    i = sizeof(petsc_pc_list) / sizeof(*petsc_pc_list);
    if (petsc_pc_default < 0 || petsc_pc_default >= i) {
        petsc_pc_default = 1;
    }

    KSPSetType(pts->ksp, petsc_solver_list[petsc_solver_default]);
    PCSetType(pts->pc, petsc_pc_list[petsc_pc_default]);
    PCFactorSetLevels(pts->pc, petsc_iluk_level);

    if (petsc_pc_right) KSPSetPCSide(pts->ksp, PC_RIGHT);

    /* additional settings */
    if (solver_petsc_setting_func != NULL) {
        (*solver_petsc_setting_func)(&pts->ksp, &pts->pc);
    }

    return 0;
}

static int petsc_destroy(PTS_DATA *pts)
{
    if (pts == NULL) return 0;

    if (pts->ksp != PETSC_NULL) KSPDestroy(&pts->ksp);
    if (pts->A != PETSC_NULL) MatDestroy(&pts->A);
    if (pts->b != PETSC_NULL) VecDestroy(&pts->b);
    if (pts->x != PETSC_NULL) VecDestroy(&pts->x);

    lssp_free(pts);

    return 0;
}

static int petsc_mat_add_entries(PTS_DATA *pts, int nrows, int *rows, int ncols, int *cols,
        double *values)
{
    PetscInt *prows, *pcols;

    if (nrows == 0 || ncols == 0) return 0;
    assert(pts->A != PETSC_NULL);

    /* index */
    if (sizeof(int) != sizeof(PetscInt)) {
        int i;

        prows = lssp_malloc<PetscInt>(nrows + ncols);
        pcols = prows + nrows;

        for (i = 0; i < nrows; i++) prows[i] = (PetscInt)rows[i];
        for (i = 0; i < ncols; i++) pcols[i] = (PetscInt)cols[i];
    }
    else {
        prows = (PetscInt *)rows;
        pcols = (PetscInt *)cols;
    }

    /* value */
    if (sizeof(double) == sizeof(PetscScalar)) {
        MatSetValues(pts->A, nrows, prows, ncols, pcols,
                (PetscScalar *)values, INSERT_VALUES);
    }
    else {
        PetscInt i, n = nrows * ncols;
        PetscScalar *buffer = lssp_malloc<PetscScalar>(n);

        for (i = 0; i < n; i++) buffer[i] = (PetscScalar)values[i];

        MatSetValues(pts->A, nrows, prows, ncols, pcols,
                buffer, INSERT_VALUES);
        lssp_free(buffer);
    }

    if (sizeof(int) != sizeof(PetscInt)) lssp_free(prows);

    return 0;
}

static int petsc_rhs_add_entries(PTS_DATA *pts, int n, int *ni, double *values)
{
    PetscInt *idx;

    if (n == 0) return 0;
    assert(pts->b != PETSC_NULL);

    if (sizeof(int) != sizeof(PetscInt)) {
        int i;

        idx = lssp_malloc<PetscInt>(n);
        for (i = 0; i < n; i++) idx[i] = (PetscInt)ni[i];
    }
    else {
        idx = (PetscInt *)ni;
    }

    if (sizeof(double) == sizeof(PetscScalar)) {
        VecSetValues(pts->b, n, idx, (PetscScalar*)values, INSERT_VALUES);
    }
    else {
        PetscInt i;

        PetscScalar *buffer = lssp_malloc<PetscScalar>(n);
        for (i = 0; i < n; i++) buffer[i] = (PetscScalar)values[i];

        VecSetValues(pts->b, n, idx, buffer, INSERT_VALUES);
        lssp_free(buffer);
    }

    if (sizeof(int) != sizeof(PetscInt)) lssp_free(idx);

    return 0;
}

static int petsc_assemble(PTS_DATA *pts, LSSP_SOLVER *solver)
{
    lssp_mat_csr A = solver->A;
    int N = A.num_rows;

    if (pts->assembled) return 0;

    /* matrix A */
    for (int i = 0; i < N; i++) {
        petsc_mat_add_entries(pts, 1, &i, solver->A.Ap[i + 1] - solver->A.Ap[i],
                solver->A.Aj + solver->A.Ap[i], solver->A.Ax + solver->A.Ap[i]);
    }

    /* rhs */
    for (int i = 0; i < N; i++) {
        petsc_rhs_add_entries(pts, 1, &i, solver->rhs.d + i);
    }

    /* matrix */
    MatAssemblyBegin(pts->A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(pts->A, MAT_FINAL_ASSEMBLY);

    /* rhs */
    VecAssemblyBegin(pts->b);
    VecAssemblyEnd(pts->b);

    pts->assembled = true;

    return 0;
}

void lssp_solver_petsc_create(LSSP_SOLVER &s)
{
    PTS_DATA *pts;

    solver_petsc_initialize();
    s.petsc = pts = lssp_calloc<PTS_DATA>(1);
}

void lssp_solver_petsc_destroy(LSSP_SOLVER &s)
{
    PTS_DATA *pts = s.petsc;

    if (pts == NULL) return;

    /* destroy */
    petsc_destroy(pts);
    solver_petsc_finalize();
    s.petsc = NULL;
}

int lssp_solver_petsc(LSSP_SOLVER *solver)
{
    PetscInt nits;
    PetscReal residual;
    PetscScalar *vec;
    PetscInt size;
    int i;
    double *x;
    PTS_DATA *pts = solver->petsc;

    assert(solver != NULL);
    petsc_create(pts, solver);

    assert(pts->A != PETSC_NULL);
    assert(pts->b != PETSC_NULL);
    assert(pts->x != PETSC_NULL);
    x = solver->x.d;

    /* assemble */
    petsc_assemble(pts, solver);

    /* copy initial solution */
    x = solver->x.d;
    VecGetLocalSize(pts->x, &size);
    VecGetArray(pts->x, &vec);
    for (i = 0; i < size; i++) vec[i] = x[i];
    VecRestoreArray(pts->x, &vec);

    /* reset stop standard */
    KSPSetTolerances(pts->ksp, solver->tol_rel, solver->tol_abs, PETSC_DEFAULT, solver->maxit);

    /* solve */
    KSPSolve(pts->ksp, pts->b, pts->x);
    KSPGetIterationNumber(pts->ksp, &nits);
    KSPGetResidualNorm(pts->ksp, &residual);

    /* write back solution */
    VecGetLocalSize(pts->x, &size);
    VecGetArray(pts->x, &vec);
    for (i = 0; i < size; i++) x[i] = vec[i];
    VecRestoreArray(pts->x, &vec);

    solver->residual = residual;
    solver->nits = nits;

    return nits;
}

void lssp_solver_petsc_setting(solver_petsc_setting func)
{
    solver_petsc_setting_func = func;
}

#endif
