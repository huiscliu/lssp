
#ifndef LSSP_LSSP_H
#define LSSP_LSSP_H

#include "pc.h"

#include "solver-gmres.h"
#include "solver-lgmres.h"

#include "solver-orthomin.h"
#include "solver-idrs.h"

#include "solver-bicgstab.h"
#include "solver-bicgstabl.h"
#include "solver-bicgsafe.h"
#include "solver-gpbicg.h"
#include "solver-cg.h"
#include "solver-cgs.h"

#include "solver-cr.h"
#include "solver-bicrstab.h"
#include "solver-bicrsafe.h"
#include "solver-crs.h"
#include "solver-gpbicr.h"

#include "solver-qmrcgstab.h"
#include "solver-tfqmr.h"

#include "solver-laspack.h"
#include "solver-mumps.h"
#include "solver-petsc.h"
#include "solver-lis.h"
#include "solver-qrmumps.h"
#include "solver-umfpack.h"
#include "solver-klu.h"
#include "solver-fasp.h"
#include "solver-superlu.h"
#include "solver-pardiso.h"
#include "solver-mi20.h"
#include "solver-amg.h"
#include "solver-sxamg.h"

/* create solver */
void lssp_solver_create(LSSP_SOLVER &s, LSSP_SOLVER_TYPE s_type, LSSP_PC &pc, LSSP_PC_TYPE p_type);

/* assemble solver */
void lssp_solver_assemble(LSSP_SOLVER &s, lssp_mat_csr &Ax, lssp_vec x, lssp_vec b, LSSP_PC &pc);

/* destroy solver */
void lssp_solver_destroy(LSSP_SOLVER &s, LSSP_PC &pc);

/* solve */
int lssp_solver_solve(LSSP_SOLVER &solver, LSSP_PC &pc);

/* reset b/rhs, where Ax = b */
void lssp_solver_reset_rhs(LSSP_SOLVER &s, lssp_vec rhs);

/* reset x, where Ax = b */
void lssp_solver_reset_unknown(LSSP_SOLVER &s, lssp_vec x);

/* reset solver type */
void lssp_solver_reset_type(LSSP_SOLVER &s, LSSP_SOLVER_TYPE type);

/* set relative tolerence */
void lssp_solver_set_rtol(LSSP_SOLVER &s, double tol);

/* set absolute tolerence */
void lssp_solver_set_atol(LSSP_SOLVER &s, double tol);

/* set relative b norm tolerence */
void lssp_solver_set_rbtol(LSSP_SOLVER &s, double tol);

/* set maximal number of iteration */
void lssp_solver_set_maxit(LSSP_SOLVER &s, int maxit);

/* set the number of restart */
void lssp_solver_set_restart(LSSP_SOLVER &s, int m);

/* set aug_k */
void lssp_solver_set_augk(LSSP_SOLVER &s, int k);

/* set bgsl */
void lssp_solver_set_bgsl(LSSP_SOLVER &s, int k);

/* set idrs */
void lssp_solver_set_idrs(LSSP_SOLVER &s, int k);

/* reset verbosity of solver */
void lssp_solver_reset_verbosity(LSSP_SOLVER &s, int v);

double lssp_solver_get_residual(LSSP_SOLVER s);
int lssp_solver_get_nits(LSSP_SOLVER s);

void lssp_solver_set_log(LSSP_SOLVER &s, FILE *io);

#endif
