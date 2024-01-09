
#include "lssp.h"

/* constants */
int LSSP_RESTART        = 50;
int LSSP_AUG_K          = 3;
int LSSP_BGSL           = 4;
int LSSP_IDRS           = 4;
int LSSP_MAXIT          = 1000;

double LSSP_ATOL        = 1e-7;
double LSSP_RTOL        = 1e-7;
double LSSP_RB          = 1e-7;
double LSSP_BREAKDOWN   = 1e-40;

void lssp_solver_create(LSSP_SOLVER &s, LSSP_SOLVER_TYPE s_type, LSSP_PC &pc, LSSP_PC_TYPE p_type)
{
    bzero(&s, sizeof(LSSP_SOLVER));

    s.type = s_type;
    s.residual = 0;

    /* pars */
    s.tol_rel = LSSP_RTOL;
    s.tol_abs = LSSP_ATOL;
    s.tol_rb = LSSP_RB;

    s.restart = LSSP_RESTART;
    s.aug_k = LSSP_AUG_K;
    s.maxit = LSSP_MAXIT;
    s.nits = 0;
    s.bgsl = LSSP_BGSL;
    s.idrs = LSSP_IDRS;

    /* default pars */
#if USE_FASP
    s.fasp_solver = 1;
    s.fasp_pc = 4;
#endif

#if USE_LIS
    s.lis_solver = 0;
    s.lis_pc = 2;
#endif


#if USE_PETSC
    s.petsc_solver = 1;
    s.petsc_pc = 1;
#endif

#if USE_LASPACK
    s.laspack_solver = 0;
    s.laspack_pc = 0;
#endif

    lssp_mat_init(s.A);
    s.num_blks = -1;
    s.blk_size = NULL;

    s.log = NULL;
    s.verb = lssp_verbosity;

#if USE_LASPACK
    if (s_type == LSSP_SOLVER_LASPACK) {
        p_type = LSSP_PC_NON;
        lssp_solver_laspack_create(s);
    }
#endif

#if USE_SSPARSE
    if (s_type == LSSP_SOLVER_KLU) p_type = LSSP_PC_NON;

    if (s_type == LSSP_SOLVER_UMFPACK) {
        p_type = LSSP_PC_NON;
        lssp_solver_umfpack_create(s);
    }
#endif

#if USE_MUMPS
    if (s_type == LSSP_SOLVER_MUMPS) {
        p_type = LSSP_PC_NON;
        lssp_solver_mumps_create(s);
    }
#endif

#if USE_PETSC
    if (s_type == LSSP_SOLVER_PETSC) {
        p_type = LSSP_PC_NON;
        lssp_solver_petsc_create(s);
    }
#endif

#if USE_LIS
    if (s_type == LSSP_SOLVER_LIS) p_type = LSSP_PC_NON;
#endif

#if USE_QR_MUMPS
    if (s_type == LSSP_SOLVER_QR_MUMPS) p_type = LSSP_PC_NON;
#endif

    /* fasp */
#if USE_FASP
    if (s_type == LSSP_SOLVER_FASP) {
        p_type = LSSP_PC_NON;
        lssp_solver_fasp_create(s);
    }

    /* amg */
    if (s_type == LSSP_SOLVER_FASP_AMG || s_type == LSSP_SOLVER_FASP_FMG) {
        p_type = LSSP_PC_NON;
        lssp_solver_amg_create(s);
    }
#endif

#if USE_SUPERLU
    if (s_type == LSSP_SOLVER_SUPERLU) {
        p_type = LSSP_PC_NON;
        lssp_solver_superlu_create(s);
    }
#endif

#if USE_PARDISO
    if (s_type == LSSP_SOLVER_PARDISO) p_type = LSSP_PC_NON;
#endif

#if USE_HSL_MI20
    if (s_type == LSSP_SOLVER_MI20AMG) p_type = LSSP_PC_NON;
#endif

#if USE_SXAMG
    if (s_type == LSSP_SOLVER_SXAMG) {
        p_type = LSSP_PC_NON;
        lssp_solver_sxamg_create(s);
    }
#endif

    lssp_pc_create(pc, p_type);
    s.assembled = false;
}

void lssp_solver_assemble(LSSP_SOLVER &s, lssp_mat_csr &Ax, lssp_vec x, lssp_vec b, LSSP_PC &pc)
{
    lssp_mat_csr A;
    double t = 0.;

    if (Ax.num_rows <= 0) {
        lssp_error(1, "solver: wrong input matrix, number of rows "
                "should be greater than 0\n");
    }
    else {
        if (Ax.num_rows != Ax.num_cols) {
            lssp_error(1, "solver: wrong input matrix, number of rows != "
                    "number of columns\n");
        }
    }

    if (Ax.num_nnzs < Ax.num_rows) {
        lssp_error(1, "solver: wrong input matrix, singular\n");
    }

    if (s.verb > 1) {
        t = lssp_get_time();
    }

    A.num_rows = Ax.num_rows;
    A.num_cols = Ax.num_cols;
    A.num_nnzs = Ax.num_nnzs;
    A.Ap = lssp_copy_on<int>(Ax.Ap, Ax.num_rows + 1);
    A.Aj = lssp_copy_on<int>(Ax.Aj, Ax.num_nnzs);
    A.Ax = lssp_copy_on<double>(Ax.Ax, Ax.num_nnzs);

    if (!lssp_mat_csr_is_sorted(A)) lssp_mat_sort_column(A);

    s.rhs = b;
    s.x = x;

    s.A = A;

    if (s.verb > 1) {
        t = lssp_get_time() - t;

        lssp_printf("solver: assemble time: %g\n", t);
    }

    s.assembled = true;

    lssp_pc_assemble(pc, s);
}

void lssp_solver_destroy(LSSP_SOLVER &s, LSSP_PC &pc)
{
    assert(s.assembled);
    lssp_mat_destroy(s.A);

    /* fasp */
#if USE_FASP
    if (s.type == LSSP_SOLVER_FASP) {
        lssp_solver_fasp_destroy(s);
    }

    /* amg solver */
    if (s.type == LSSP_SOLVER_FASP_AMG || s.type == LSSP_SOLVER_FASP_FMG) {
        lssp_solver_amg_destroy(s);
    }
#endif

    /* laspack */
#if USE_LASPACK
    if (s.type == LSSP_SOLVER_LASPACK) {
        lssp_solver_laspack_destroy(s);
    }
#endif

    /* mumps */
#if USE_MUMPS
    if (s.type == LSSP_SOLVER_MUMPS) {
        lssp_solver_mumps_destroy(s);
    }
#endif

#if USE_PETSC
    if (s.type == LSSP_SOLVER_PETSC) {
        lssp_solver_petsc_destroy(s);
    }
#endif

#if USE_SUPERLU
    if (s.type == LSSP_SOLVER_SUPERLU) {
        lssp_solver_superlu_destroy(s);
    }
#endif

#if USE_SSPARSE
    if (s.type == LSSP_SOLVER_UMFPACK) {
        lssp_solver_umfpack_destroy(s);
    }
#endif

#if USE_SXAMG
    if (s.type == LSSP_SOLVER_SXAMG) {
        lssp_solver_sxamg_destroy(s);
    }
#endif

    lssp_pc_destroy(pc);
    s.assembled = false;
}

int lssp_solver_solve(LSSP_SOLVER &solver, LSSP_PC &pc)
{
    int it;
    LSSP_SOLVER_TYPE type = solver.type;

    assert(solver.assembled);
    assert(pc.assembled);

    switch (type) {
        case LSSP_SOLVER_GMRES:
            it = lssp_solver_gmres(solver, pc);
            return it;

        case LSSP_SOLVER_LGMRES:
            it = lssp_solver_lgmres(solver, pc);
            return it;

        case LSSP_SOLVER_RGMRES:
            it = lssp_solver_gmres_r(solver, pc);
            return it;

        case LSSP_SOLVER_RLGMRES:
            it = lssp_solver_lgmres_r(solver, pc);
            return it;

        case LSSP_SOLVER_BICGSTAB:
            it = lssp_solver_bicgstab(solver, pc);
            return it;

        case LSSP_SOLVER_BICGSTABL:
            it = lssp_solver_bicgstabl(solver, pc);
            return it;

        case LSSP_SOLVER_BICGSAFE:
            it = lssp_solver_bicgsafe(solver, pc);
            return it;

        case LSSP_SOLVER_CG:
            it = lssp_solver_cg(solver, pc);
            return it;

        case LSSP_SOLVER_CGS:
            it = lssp_solver_cgs(solver, pc);
            return it;

        case LSSP_SOLVER_GPBICG:
            it = lssp_solver_gpbicg(solver, pc);
            return it;

        case LSSP_SOLVER_GPBICR:
            it = lssp_solver_gpbicr(solver, pc);
            return it;

        case LSSP_SOLVER_CR:
            it = lssp_solver_cr(solver, pc);
            return it;

        case LSSP_SOLVER_CRS:
            it = lssp_solver_crs(solver, pc);
            return it;

        case LSSP_SOLVER_BICRSAFE:
            it = lssp_solver_bicrsafe(solver, pc);
            return it;

        case LSSP_SOLVER_BICRSTAB:
            it = lssp_solver_bicrstab(solver, pc);
            return it;

        case LSSP_SOLVER_ORTHOMIN:
            it = lssp_solver_orthomin(solver, pc);
            return it;

        case LSSP_SOLVER_QMRCGSTAB:
            it = lssp_solver_qmrcgstab(solver, pc);
            return it;

        case LSSP_SOLVER_TFQMR:
            it = lssp_solver_tfqmr(solver, pc);
            return it;

        case LSSP_SOLVER_IDRS:
            it = lssp_solver_idrs(solver, pc);
            return it;

#if USE_LASPACK
        case LSSP_SOLVER_LASPACK:
            it = lssp_solver_laspack(&solver);
            return it;
#endif

#if USE_SSPARSE
        case LSSP_SOLVER_UMFPACK:
            it = lssp_solver_umfpack(&solver);
            return it;

        case LSSP_SOLVER_KLU:
            it = lssp_solver_klu(&solver);
            return it;
#endif

#if USE_MUMPS
        case LSSP_SOLVER_MUMPS:
            it = lssp_solver_mumps(&solver);
            return it;
#endif

#if USE_PETSC
        case LSSP_SOLVER_PETSC:
            it = lssp_solver_petsc(&solver);
            return it;
#endif

#if USE_LIS
        case LSSP_SOLVER_LIS:
            it = lssp_solver_lis(&solver);
            return it;
#endif

#if USE_QR_MUMPS
        case LSSP_SOLVER_QR_MUMPS:
            it = lssp_solver_qr_mumps(&solver);
            return it;
#endif

#if USE_FASP
        case LSSP_SOLVER_FASP:
            it = lssp_solver_fasp(&solver);
            return it;

        case LSSP_SOLVER_FASP_AMG:
        case LSSP_SOLVER_FASP_FMG:
            it = lssp_solver_amg(&solver);
            return it;
#endif

#if USE_SUPERLU
        case LSSP_SOLVER_SUPERLU:
            it = lssp_solver_superlu(&solver);
            return it;
#endif

#if USE_PARDISO
        case LSSP_SOLVER_PARDISO:
            it = lssp_solver_pardiso(&solver);
            return it;
#endif

#if USE_HSL_MI20
        case LSSP_SOLVER_MI20AMG:
            it = lssp_solver_mi20amg(&solver);
            return it;
#endif

#if USE_SXAMG
        case LSSP_SOLVER_SXAMG:
            it = lssp_solver_sxamg(&solver);
            return it;
#endif

        default:
            lssp_error(1, "solver: wrong solver type!\n");
            return -1;
    }
}

void lssp_solver_reset_rhs(LSSP_SOLVER &s, lssp_vec rhs)
{
    assert(s.assembled);

    s.rhs = rhs;
}

void lssp_solver_reset_unknown(LSSP_SOLVER &s, lssp_vec x)
{
    assert(s.assembled);

    s.x = x;
}

void lssp_solver_reset_type(LSSP_SOLVER &s, LSSP_SOLVER_TYPE type)
{
    s.type = type;
}

void lssp_solver_set_rtol(LSSP_SOLVER &s, double tol)
{
    if (tol < 0) {
        lssp_warning("solver: tol is less than zero!\n");
    }
    else {
        s.tol_rel = tol;
    }
}

void lssp_solver_set_atol(LSSP_SOLVER &s, double tol)
{
    if (tol < 0) {
        lssp_warning("solver: tol is less than zero!\n");
    }
    else {
        s.tol_abs = tol;
    }
}

void lssp_solver_set_rbtol(LSSP_SOLVER &s, double tol)
{
    if (tol < 0) {
        lssp_warning("solver: tol is less than zero!\n");
    }
    else {
        s.tol_rb = tol;
    }
}

void lssp_solver_set_maxit(LSSP_SOLVER &s, int maxit)
{
    if (maxit <= 0) {
        lssp_warning("solver: maxit is less or equal to zero!\n");
    }
    else {
        s.maxit = maxit;
    }
}

void lssp_solver_set_restart(LSSP_SOLVER &s, int m)
{
    if (m <= 0) {
        lssp_warning("solver: restart is too small!\n");
    }
    else {
        s.restart = m;
    }
}

void lssp_solver_set_augk(LSSP_SOLVER &s, int k)
{
    if (k <= 0) {
        lssp_warning("solver: aug_k is less or equal to zero!\n");
    }
    else {
        s.aug_k = k;
    }
}

void lssp_solver_set_bgsl(LSSP_SOLVER &s, int k)
{
    if (k <= 0) {
        lssp_warning("solver: bgsl is less or equal to zero!\n");
    }
    else {
        s.bgsl = k;
    }
}

void lssp_solver_set_idrs(LSSP_SOLVER &s, int k)
{
    if (k <= 0) {
        lssp_warning("solver: idrs is less or equal to zero!\n");
    }
    else {
        s.idrs = k;
    }
}

void lssp_solver_reset_verbosity(LSSP_SOLVER &s, int v)
{
    s.verb = v;
}

double lssp_solver_get_residual(LSSP_SOLVER s)
{
    return s.residual;
}

int lssp_solver_get_nits(LSSP_SOLVER s)
{
    return s.nits;
}

void lssp_solver_set_log(LSSP_SOLVER &s, FILE *io)
{
    assert(io != NULL);
    s.log = io;
    lssp_set_log(io);
}
