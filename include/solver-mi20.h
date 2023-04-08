
#ifndef LSSP_SOLVER_MI20_H
#define LSSP_SOLVER_MI20_H

#include "type-defs.h"
#include "matrix-utils.h"
#include "mvops.h"

#if USE_HSL_MI20
int lssp_solver_mi20amg(LSSP_SOLVER *solver);
#endif

#endif
