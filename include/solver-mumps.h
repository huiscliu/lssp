
#ifndef LSSP_SOLVER_MUMPS_H
#define LSSP_SOLVER_MUMPS_H

#include "type-defs.h"
#include "utils.h"

#if USE_MUMPS
void lssp_solver_mumps_create(LSSP_SOLVER &s);
void lssp_solver_mumps_destroy(LSSP_SOLVER &s);

int lssp_solver_mumps(LSSP_SOLVER *solver);
#endif

#endif
