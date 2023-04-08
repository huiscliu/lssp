
#ifndef LSSP_PC_AMG_H
#define LSSP_PC_AMG_H

#include "type-defs.h"
#include "matrix-utils.h"
#include "mvops.h"

#if USE_FASP

#define WITH_UMFPACK 0
#define WITH_PARDISO 0
#define WITH_MUMPS   0

extern "C" {
#include <fasp.h>
#include <fasp_functs.h>
}

void lssp_pc_amg_create(LSSP_PC &pc);
void lssp_pc_amg_assemble(LSSP_PC &pc, LSSP_SOLVER s);
void lssp_pc_amg_set_pars(LSSP_PC &pc, AMG_param *amgparam);

#endif
#endif
