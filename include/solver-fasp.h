
#ifndef LSSP_SOLVER_FASP_H
#define LSSP_SOLVER_FASP_H

#include "type-defs.h"
#include "utils.h"

#if USE_FASP

#define WITH_UMFPACK             0
#define WITH_PARDISO             0
#define WITH_MUMPS               0

extern "C" {
#include <fasp.h>
#include <fasp_functs.h>
}

void lssp_solver_fasp_create(LSSP_SOLVER &s);
void lssp_solver_fasp_destroy(LSSP_SOLVER &s);

void lssp_solver_fasp_set_pars(LSSP_SOLVER *solver, input_param *inparam, itsolver_param *itsparam,
        AMG_param *amgparam, ILU_param *iluparam, Schwarz_param *schparam);

int lssp_solver_fasp(LSSP_SOLVER *solver);

#endif

#endif
