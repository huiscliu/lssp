
#ifndef LSSP_SOLVER_LIS_H
#define LSSP_SOLVER_LIS_H

#include "type-defs.h"
#include "utils.h"

#if USE_LIS

int lssp_solver_lis(LSSP_SOLVER *solver);
void lssp_solver_lis_set_pars(LSSP_SOLVER &solver, unsigned int solver_id, unsigned int pc_id);

#endif

#endif
