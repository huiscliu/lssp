
#include "solver-laspack.h"

#if USE_LASPACK

#include <string.h>

extern "C" {
#include <laspack/vector.h>
#include <laspack/operats.h>
#include <laspack/qmatrix.h>

#include <laspack/itersolv.h>
#include <laspack/rtc.h>
#include <laspack/errhandl.h>
}

typedef struct LAS_DATA_
{
    Vector      *B;
    Vector      *U;
    QMatrix     *A;
    int         ksp_type;
    int         pc_type;
    int         gmres_restart;

} LAS_DATA;

static const IterProcType ksp_table[] = {
    BiCGSTABIter, GMRESIter, QMRIter,
    CGIter, CGNIter, CGSIter, BiCGIter,
};

static const PrecondProcType pc_table[] = {ILUPrecond, JacobiPrecond, SSORPrecond};

static int mat_add_entries(LAS_DATA *las_data, int nrows, int *rows, int ncols, int *cols,
        double *values)
{
    int i, j;

    for (i = 0; i < nrows; i++, rows++) {
        Q_SetLen(las_data->A, 1 + *rows, ncols);
        for (j = 0; j < ncols; j++, cols++, values++) {
            Q__SetEntry(las_data->A, 1 + *rows, j, 1 + *cols, *values);
        }
    }

    return 0;
}

static int rhs_add_entries(LAS_DATA *las_data, int n, int *ni, double *values)
{
    int i;

    for (i = 0; i < n; i++, ni++, values++)
        V__SetCmp(las_data->B, 1 + *ni, *values);

    return 0;
}

static int Init(LAS_DATA *las_data, LSSP_SOLVER *solver)
{
    int N = solver->A.num_rows;
    int k;

    /* type */
    k = sizeof(ksp_table) / sizeof(*ksp_table);
    if (solver->laspack_solver < 0 || solver->laspack_solver >= k) {
        solver->laspack_solver = 0;
    }

    las_data->ksp_type = solver->laspack_solver;

    /* type */
    k = sizeof(pc_table) / sizeof(*pc_table);
    if (solver->laspack_pc < 0 || solver->laspack_pc >= k) {
        solver->laspack_pc = 0;
    }

    las_data->pc_type = solver->laspack_pc;
    las_data->gmres_restart = solver->restart;

    /* create matrix */
    las_data->A = lssp_malloc<QMatrix>(1);
    Q_Constr(las_data->A, const_cast<char*>("LASPack matrix"), N, False, Rowws, Normal, True);

    /* create vectors */
    las_data->B = lssp_malloc<Vector>(1);
    V_Constr(las_data->B, const_cast<char*>("RHS"), N, Normal, True);

    las_data->U = lssp_malloc<Vector>(1);
    V_Constr(las_data->U, const_cast<char*>("solution"), N, Normal, True);

    /* matrix A */
    for (int i = 0; i < N; i++) {
        mat_add_entries(las_data, 1, &i, solver->A.Ap[i + 1] - solver->A.Ap[i],
                solver->A.Aj + solver->A.Ap[i], solver->A.Ax + solver->A.Ap[i]);
    }

    /* rhs */
    for (int i = 0; i < N; i++) {
        rhs_add_entries(las_data, 1, &i, solver->rhs.d + i);
    }

    return 0;
}

static void laspack_callback(int it, double rNorm, double bNorm, IterIdType id)
{
    lssp_unused(id);
    lssp_printf("laspack: itr: %5d, abs res: %le rel res: %le\n", it, rNorm,
            (double)(rNorm / (bNorm == 0.0 ? 1.0 : bNorm)));
}

void lssp_solver_laspack_create(LSSP_SOLVER &s)
{
    s.laspack = lssp_malloc<LAS_DATA>(1);
}

void lssp_solver_laspack_destroy(LSSP_SOLVER &s)
{
    LAS_DATA *las_data = s.laspack;

    if (las_data == NULL) return;

    /* destroy */
    Q_Destr(las_data->A);
    lssp_free(las_data->A);
    las_data->A = NULL;
    V_Destr(las_data->B);
    lssp_free(las_data->B);
    las_data->B = NULL;
    V_Destr(las_data->U);
    lssp_free(las_data);

    s.laspack = NULL;
}

int lssp_solver_laspack(LSSP_SOLVER *solver)
{
    int i;
    double *x;
    LAS_DATA *las_data = solver->laspack;

    assert(solver != NULL);
    x = solver->x.d;

    /* init */
    Init(las_data, solver);

    /* initial guess */
    for (i = 0; i < solver->A.num_rows; i++) V__SetCmp(las_data->U, 1 + i, x[i]);

    /* solve */
    assert(solver->tol_rel > 0. || solver->tol_rb > 0.);
    SetRTCAccuracy(solver->tol_rel > solver->tol_rb ? solver->tol_rel : solver->tol_rb);
    SetGMRESRestart(las_data->gmres_restart);
    if (solver->verb > 0) SetRTCAuxProc(laspack_callback);

    ksp_table[las_data->ksp_type](las_data->A, las_data->U, las_data->B, solver->maxit,
            pc_table[las_data->pc_type], 1.);
    solver->residual = GetLastAccuracy();
    solver->nits = GetLastNoIter();

    /* copy solution */
    assert(sizeof(*las_data->U->Cmp) == sizeof(double));
    memcpy(x, las_data->U->Cmp + 1, solver->A.num_rows * sizeof(*x));

    return solver->nits;
}

void lssp_solver_laspack_set_pars(LSSP_SOLVER *solver, int solver_type, int pc_type)
{
    if (solver == NULL) return;

    solver->laspack_solver = solver_type;
    solver->laspack_pc = pc_type;
}

#endif
