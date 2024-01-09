
#ifndef LSSP_SOLVER_LGMRES_H
#define LSSP_SOLVER_LGMRES_H

#include "mvops.h"

int lssp_solver_lgmres(LSSP_SOLVER &solver, LSSP_PC &pc);
int lssp_solver_lgmres_r(LSSP_SOLVER &solver, LSSP_PC &pc);

#endif
