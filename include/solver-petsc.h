
#ifndef LSSP_SOLVER_PETSC_H
#define LSSP_SOLVER_PETSC_H

#include "type-defs.h"
#include "utils.h"

#if USE_PETSC

typedef void (*solver_petsc_setting)(void *ksp, void *pc);

void lssp_solver_petsc_create(LSSP_SOLVER &s);
void lssp_solver_petsc_destroy(LSSP_SOLVER &s);

int lssp_solver_petsc(LSSP_SOLVER *solver);

void lssp_solver_petsc_setting(solver_petsc_setting func);

#endif

#endif
