
#ifndef LSSP_PC_BILUK_H
#define LSSP_PC_BILUK_H

#include "type-defs.h"
#include "matrix-utils.h"
#include "mvops.h"
#include "pc-iluk.h"

#if USE_BLAS
#if USE_LAPACK

void lssp_pc_bilu_solve(LSSP_PC *pc, lssp_vec x, lssp_vec rhs);
void lssp_pc_biluk_destroy(LSSP_PC *pc);
void lssp_pc_biluk_assemble_mat(LSSP_PC &pc, lssp_mat_bcsr A);
void lssp_pc_biluk_assemble(LSSP_PC &pc, LSSP_SOLVER s);

#endif
#endif

#endif
