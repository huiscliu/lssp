
#ifndef LSSP_SOLVER_AMG_H
#define LSSP_SOLVER_AMG_H

#include "type-defs.h"
#include "matrix-utils.h"
#include "mvops.h"

#if USE_FASP

#define WITH_UMFPACK 0
#define WITH_PARDISO 0
#define WITH_MUMPS   0

#define __STDC_VERSION__         199901L

extern "C" {
#include <fasp.h>
#include <fasp_functs.h>
}

void lssp_solver_amg_create(LSSP_SOLVER &s);
void lssp_solver_amg_destroy(LSSP_SOLVER &s);

int lssp_solver_amg(LSSP_SOLVER *solver);
void lssp_solver_amg_reset_pars(LSSP_SOLVER &s, AMG_param par);

#endif
#endif
