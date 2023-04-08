
#ifndef LSSP_PC_BILUT_H
#define LSSP_PC_BILUT_H

#include "type-defs.h"
#include "matrix-utils.h"
#include "mvops.h"

#if USE_ITSOL

extern "C" {
#include "globheads.h"
#include "protos.h"
#include "ios.h"
}

SparMat lssp_mat_csr_to_SparMat(lssp_mat_csr A);

void lssp_pc_bilut_destroy(LSSP_PC *pc);
void lssp_pc_bilut_solve(LSSP_PC *s, lssp_vec x, lssp_vec rhs);
void lssp_pc_bilut_assemble(LSSP_PC &pc, LSSP_SOLVER s);

#endif
#endif
