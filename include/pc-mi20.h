
#ifndef LSSP_PC_MI20_H
#define LSSP_PC_MI20_H

#include "type-defs.h"
#include "matrix-utils.h"
#include "mvops.h"

#if USE_HSL_MI20

void lssp_pc_mi20amg_create(LSSP_PC &pc);
void lssp_pc_mi20amg_assemble(LSSP_PC &pc, LSSP_SOLVER s);

#endif
#endif
