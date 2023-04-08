
#ifndef LSSP_PC_SXAMG_H
#define LSSP_PC_SXAMG_H

#include "type-defs.h"
#include "matrix-utils.h"
#include "mvops.h"

#if USE_SXAMG

#include "sxamg.h"

void lssp_pc_sxamg_create(LSSP_PC &pc);
void lssp_pc_sxamg_assemble(LSSP_PC &pc, LSSP_SOLVER s);

void lssp_pc_sxamg_set_pars(LSSP_PC &pc, SX_AMG_PARS *pars);

#endif
#endif
