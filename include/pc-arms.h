
#ifndef LSSP_PC_ARMS_H
#define LSSP_PC_ARMS_H

#include "type-defs.h"
#include "matrix-utils.h"
#include "mvops.h"
#include "pc-bilut.h"

#if USE_ITSOL

void lssp_pc_arms_create(LSSP_PC &pc);
void lssp_pc_arms_assemble(LSSP_PC &pc, LSSP_SOLVER s);

#endif
#endif
