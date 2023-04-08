
#ifndef LSSP_SOLVER_SUPERLU_H
#define LSSP_SOLVER_SUPERLU_H

#include "type-defs.h"
#include "utils.h"

#if USE_SUPERLU
void lssp_solver_superlu_create(LSSP_SOLVER &s);
void lssp_solver_superlu_destroy(LSSP_SOLVER &s);

int lssp_solver_superlu(LSSP_SOLVER *solver);
#endif

#endif
