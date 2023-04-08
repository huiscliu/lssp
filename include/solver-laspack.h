
#ifndef LSSP_SOLVER_LASPACK_H
#define LSSP_SOLVER_LASPACK_H

#include "type-defs.h"
#include "utils.h"

#if USE_LASPACK
void lssp_solver_laspack_create(LSSP_SOLVER &s);
void lssp_solver_laspack_destroy(LSSP_SOLVER &s);

int lssp_solver_laspack(LSSP_SOLVER *solver);

void lssp_solver_laspack_set_pars(LSSP_SOLVER *solver, int solver_type, int pc_type);
#endif

#endif
