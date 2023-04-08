
#ifndef LSSP_SOLVER_SXAMG_H
#define LSSP_SOLVER_SXAMG_H

#include "type-defs.h"
#include "matrix-utils.h"
#include "mvops.h"

#if USE_SXAMG

#include "sxamg.h"

void lssp_solver_sxamg_create(LSSP_SOLVER &s);
void lssp_solver_sxamg_destroy(LSSP_SOLVER &s);

int lssp_solver_sxamg(LSSP_SOLVER *solver);

void lssp_solver_sxamg_set_pars(LSSP_SOLVER *solver, SX_AMG_PARS *pars);

#endif

#endif
