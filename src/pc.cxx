#include "pc.h"

int lssp_pc_iluk_level_default = 1;
double lssp_pc_ilut_tol        = 1e-3;
double lssp_pc_ilut_p          = -1;
double mat_zero_diag_value       = 1e-3;
double mat_zero_diag_tol       = 1e-10;

void lssp_pc_create(LSSP_PC &pc, LSSP_PC_TYPE type)
{
    /* pars */
    pc.type = type;
    pc.iluk_level = lssp_pc_iluk_level_default;
    pc.ilut_tol = lssp_pc_ilut_tol;
    pc.ilut_p = lssp_pc_ilut_p;

    /* amg/fmg */
#if USE_FASP
    if (type == LSSP_PC_FASP_AMG || type == LSSP_PC_FASP_FMG) {
        lssp_pc_amg_create(pc);
    }
#endif

#if USE_ITSOL
    if (type == LSSP_PC_ARMS) {
        lssp_pc_arms_create(pc);
    }
#endif

#if USE_HSL_MI20
    if (type == LSSP_PC_MI20AMG) {
        lssp_pc_mi20amg_create(pc);
    }
#endif

#if USE_SXAMG
    if (type == LSSP_PC_SXAMG) {
        lssp_pc_sxamg_create(pc);
    }
#endif

    pc.cache = NULL;
    pc.solve = NULL;
    pc.destroy = NULL;

    pc.data = NULL;
    lssp_mat_init(pc.L);
    lssp_mat_init(pc.U);
    lssp_mat_init(pc.D);

    pc.verb = lssp_verbosity - 1;
    pc.log = NULL;
    pc.assembled = false;
}

void lssp_pc_destroy(LSSP_PC &pc)
{
    assert(pc.assembled);

    if (pc.destroy != NULL) {
        (*pc.destroy)(&pc);
    }

    pc.assembled = false;
}

static void lssp_pc_non_solve(LSSP_PC *pc, lssp_vec x, lssp_vec rhs)
{
    lssp_memcpy_on(x.d, rhs.d, pc->A.num_rows);
}

static void lssp_pc_non_assemble(LSSP_PC &pc, LSSP_SOLVER s)
{
    lssp_unused(s);

    pc.cache = NULL;
    pc.solve = lssp_pc_non_solve;
    pc.destroy = NULL;
}

void lssp_pc_assemble(LSSP_PC &pc, LSSP_SOLVER s)
{
    double tm = 0;

    if (!s.assembled) {
        lssp_error(1, "pc: solver hasn't been assembled, call lssp_assemble first!\n");
    }

    if (s.verb > 0) {
        tm = lssp_get_time();
    }

    pc.verb = s.verb - 1;
    pc.log = s.log;

    lssp_mat_init(pc.L);
    lssp_mat_init(pc.U);
    lssp_mat_init(pc.D);

    pc.A = s.A;
    switch (pc.type) {
        case LSSP_PC_NON:
            lssp_pc_non_assemble(pc, s);
            break;

        case LSSP_PC_ILUK:
            if (pc.verb >= 0) {
                lssp_printf("pc: type: ILUK\n");
            }

            lssp_pc_iluk_assemble(pc, s);

            break;

        case LSSP_PC_ILUT:
            if (pc.verb >= 0) {
                lssp_printf("pc: type: ILUT\n");
            }

            lssp_pc_ilut_assemble(pc, s);

            break;

#if USE_BLAS
#if USE_LAPACK
        case LSSP_PC_BILUK:
            if (pc.verb >= 0) {
                lssp_printf("pc: type: block version ILUK\n");
            }

            lssp_pc_biluk_assemble(pc, s);

            break;
#endif
#endif

#if USE_ITSOL
        case LSSP_PC_BILUT:
            if (pc.verb >= 0) {
                lssp_printf("pc: type: block version ILUT\n");
            }

            lssp_pc_bilut_assemble(pc, s);

            break;

        case LSSP_PC_ARMS:
            if (pc.verb >= 0) {
                lssp_printf("pc: type: ARMS\n");
            }

            lssp_pc_arms_assemble(pc, s);

            break;

        case LSSP_PC_VBILUK:
            if (pc.verb >= 0) {
                lssp_printf("pc: type: var block version ILUK\n");
            }

            lssp_pc_vbiluk_assemble(pc, s);

            break;

        case LSSP_PC_VBILUT:
            if (pc.verb >= 0) {
                lssp_printf("pc: type: var block version ILUT\n");
            }

            lssp_pc_vbilut_assemble(pc, s);

            break;
#endif

#if USE_FASP
        case LSSP_PC_FASP_AMG:
            if (pc.verb >= 0) {
                lssp_printf("pc: type: AMG\n");
            }

            lssp_pc_amg_assemble(pc, s);

            break;
#endif

#if USE_FASP
        case LSSP_PC_FASP_FMG:
            if (pc.verb >= 0) {
                lssp_printf("pc: type: FMG\n");
            }

            lssp_pc_amg_assemble(pc, s);

            break;
#endif

#if USE_HSL_MI20
        case LSSP_PC_MI20AMG:
            if (pc.verb >= 0) {
                lssp_printf("pc: type: MI20 AMG\n");
            }

            lssp_pc_mi20amg_assemble(pc, s);

            break;
#endif

#if USE_SXAMG
        case LSSP_PC_SXAMG:
            if (pc.verb >= 0) {
                lssp_printf("pc: type: SX-AMG\n");
            }

            lssp_pc_sxamg_assemble(pc, s);

            break;
#endif

        case LSSP_PC_USER:
            if (pc.verb >= 0) {
                lssp_printf("pc: type: user defined\n");
            }

            assert(pc.assemble != NULL);
            pc.assemble(pc, s);

            break;
        default:
            lssp_error(1, "pc: incorrect preconditioner type !\n");
            break;
    }

    if (pc.verb > 0) {
        tm = lssp_get_time() - tm;
        lssp_printf("pc: time for pc assemble: %f s\n", tm);
    }

    pc.assembled = true;
}
