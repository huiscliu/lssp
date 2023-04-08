
#ifndef LSSP_SOLVER_KLU_H
#define LSSP_SOLVER_KLU_H

#include "type-defs.h"
#include "utils.h"
#include "matrix-utils.h"

#if USE_SSPARSE
int lssp_solver_klu(LSSP_SOLVER *solver);
#endif

#endif
