
#ifndef LSSP_PC_ILUK_H
#define LSSP_PC_ILUK_H

#include "type-defs.h"
#include "matrix-utils.h"
#include "solver-tri.h"

lssp_mat_csr lssp_pc_iluk_symbolic_fac(const lssp_mat_csr A, int level, bool set_value);
void lssp_pc_iluk_assemble(LSSP_PC &pc, LSSP_SOLVER s);
void lssp_pc_iluk_destroy(LSSP_PC *pc);

/* set level */
void lssp_pc_iluk_set_level(LSSP_PC &pc, int level);

#endif
