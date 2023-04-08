
#ifndef LSSP_SOLVER_UMFPACK_H
#define LSSP_SOLVER_UMFPACK_H

#include "type-defs.h"
#include "utils.h"

#if USE_SSPARSE
void lssp_solver_umfpack_create(LSSP_SOLVER &s);
void lssp_solver_umfpack_destroy(LSSP_SOLVER &s);

int lssp_solver_umfpack(LSSP_SOLVER *solver);
#endif

#endif
