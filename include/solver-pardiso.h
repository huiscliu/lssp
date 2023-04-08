
#ifndef LSSP_SOLVER_PARDISO_H
#define LSSP_SOLVER_PARDISO_H

#include "type-defs.h"
#include "utils.h"

#if USE_PARDISO
int lssp_solver_pardiso(LSSP_SOLVER *solver);
#endif

#endif
