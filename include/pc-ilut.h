
#ifndef LSSP_PC_ILUT_H
#define LSSP_PC_ILUT_H

#include "solver-tri.h"
#include "pc-iluk.h"
#include "mvops.h"

void lssp_pc_ilut_assemble(LSSP_PC &pc, LSSP_SOLVER s);
void lssp_pc_ilut_destroy(LSSP_PC *pc);

void lssp_pc_ilut_set_drop_tol(LSSP_PC &pc, double tol);
void lssp_pc_ilut_set_p(LSSP_PC &pc, int p);

#endif
