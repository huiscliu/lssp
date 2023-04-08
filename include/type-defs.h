#ifndef LSSP_TYPES_H
#define LSSP_TYPES_H

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#include "config.h"

typedef struct lssp_mat_csr_
{
    int num_rows;
    int num_cols;
    int num_nnzs;
    int *Ap;
    int *Aj;
    double *Ax;

} lssp_mat_csr;

typedef struct lssp_mat_coo_
{
    int num_rows;
    int num_cols;
    int num_nnzs;
    int *Ai;
    int *Aj;
    double *Ax;

} lssp_mat_coo;

typedef struct lssp_mat_entry_
{
    int i;
    int j;
    double x;

} lssp_mat_entry;

typedef struct lssp_mat_bcsr_
{
    int num_rows;
    int num_cols;
    int num_nnzs;
    int blk_size;
    int *Ap;
    int *Aj;
    double *Ax;

} lssp_mat_bcsr;

typedef struct lssp_vec_
{
    int n;
    double *d;

} lssp_vec;

typedef enum LSSP_PC_TYPE_
{
    LSSP_PC_NON,
    LSSP_PC_ILUK,               /* ILUK */
    LSSP_PC_ILUT,               /* ILUT */

#if USE_BLAS
#if USE_LAPACK
    LSSP_PC_BILUK,              /* block-wise ILUK */
#endif
#endif

#if USE_ITSOL
    LSSP_PC_BILUT,              /* block-wise ILUT */
    LSSP_PC_VBILUT,             /* var block-wise ILUT */
    LSSP_PC_VBILUK,             /* var block-wise ILUK */
    LSSP_PC_ARMS,               /* multi-level */
#endif

#if USE_FASP
    LSSP_PC_FASP_AMG,           /* AMG from FASP */
    LSSP_PC_FASP_FMG,           /* FMG from FASP */
#endif

#if USE_HSL_MI20
    LSSP_PC_MI20AMG,            /* AMG from HSL MI20 */
#endif

#if USE_SXAMG
    LSSP_PC_SXAMG,              /* AMG from SX-AMG */
#endif

    LSSP_PC_USER,               /* user defined */

} LSSP_PC_TYPE;

struct LSSP_PC_;
struct LSSP_SOLVER_;

typedef void (*LSSP_PC_ASSEMBLE)(struct LSSP_PC_ &pc, struct LSSP_SOLVER_ s);
typedef void (*LSSP_PC_SOLVE)(struct LSSP_PC_ *s, lssp_vec x, lssp_vec rhs);
typedef void (*LSSP_PC_DESTROY)(struct LSSP_PC_ *s);

typedef struct LSSP_PC_
{
    /* ILU */
    int iluk_level;
    int ilut_p;
    double ilut_tol;
    
    /* matrix */
    lssp_mat_csr A;
    lssp_mat_csr L;
    lssp_mat_csr D;
    lssp_mat_csr U;

    /* amg */
#if USE_FASP
    struct FPAMG_DATA_ *fasp_amg;
#endif

    /* arms */
#if USE_ITSOL
    struct ARMS_DATA_  *arms_data;
#endif

    /* mi20 */
#if USE_HSL_MI20
    struct MI20_AMG_   *mi20;
#endif

#if USE_SXAMG
    struct SX_DATA_    *sxamg;
#endif

    void *data;
    double *cache;

    LSSP_PC_TYPE     type;      /* pc type */
    LSSP_PC_ASSEMBLE assemble;  /* assemble */
    LSSP_PC_SOLVE    solve;     /* solve function */
    LSSP_PC_DESTROY  destroy;   /* destroy */

    FILE *log;
    int verb;
    bool assembled;

} LSSP_PC;

/* solver type */
typedef enum LSSP_SOLVER_TYPE_
{
    LSSP_SOLVER_GMRES,     /* gmres(m) */
    LSSP_SOLVER_LGMRES,    /* lgmres(m, k) */
    LSSP_SOLVER_BICGSTAB,  /* bicgstab */
    LSSP_SOLVER_BICGSTABL, /* bicgstab(l) */
    LSSP_SOLVER_BICGSAFE,  /* bicgsafe */
    LSSP_SOLVER_CG,        /* cg */
    LSSP_SOLVER_CGS,       /* cgs */
    LSSP_SOLVER_GPBICG,    /* gpbicg */
    LSSP_SOLVER_CR,        /* cr */
    LSSP_SOLVER_CRS,       /* crs */
    LSSP_SOLVER_BICRSTAB,  /* bicrstab */
    LSSP_SOLVER_BICRSAFE,  /* bicrsafe */
    LSSP_SOLVER_GPBICR,    /* gpbicr */
    LSSP_SOLVER_QMRCGSTAB, /* qmrcgstab */
    LSSP_SOLVER_TFQMR,     /* tfqmr */
    LSSP_SOLVER_ORTHOMIN,  /* orthomin */
    LSSP_SOLVER_IDRS,      /* idrs */

#if USE_LASPACK
    LSSP_SOLVER_LASPACK,   /* laspack */
#endif

#if USE_SSPARSE
    LSSP_SOLVER_UMFPACK,   /* ssparse, umf */
    LSSP_SOLVER_KLU,       /* ssparse, klu */
#endif

#if USE_MUMPS
    LSSP_SOLVER_MUMPS,     /* mumps */
#endif

#if USE_PETSC
    LSSP_SOLVER_PETSC,     /* petsc */
#endif

#if USE_LIS
    LSSP_SOLVER_LIS,       /* lis */
#endif

#if USE_QR_MUMPS
    LSSP_SOLVER_QR_MUMPS,  /* qr-mumps */
#endif

#if USE_FASP
    LSSP_SOLVER_FASP,      /* FASP */
    LSSP_SOLVER_FASP_AMG,  /* amg from FASP */
    LSSP_SOLVER_FASP_FMG,  /* fmg from FASP */
#endif

#if USE_SUPERLU
    LSSP_SOLVER_SUPERLU,   /* superlu */
#endif

#if USE_PARDISO
    LSSP_SOLVER_PARDISO,   /* pardiso */
#endif

#if USE_HSL_MI20
    LSSP_SOLVER_MI20AMG,   /* MI20 AMG */
#endif

#if USE_SXAMG
    LSSP_SOLVER_SXAMG,     /* SX-AMG */
#endif

} LSSP_SOLVER_TYPE;

typedef struct LSSP_SOLVER_
{
    /* universal settings */
    double tol_rel;                 /* relative tolerance */
    double tol_abs;                 /* absolute tolerance */
    double tol_rb;                  /* relative b norm */
    int maxit;                      /* maximal number of iteration */
    int restart;                    /* used by gmres */
    int aug_k;                      /* used by lgmres(m, k) */
    int bgsl;                       /* l for bicgstab(l) */
    int idrs;                       /* s for idr(s) */

    /* fasp */
#if USE_FASP
    struct FASP_AMG_DATA_ *fasp_amg;
    struct FASP_DATA_     *fasp;
    int fasp_solver;
    int fasp_pc;
#endif

    /* lis */
#if USE_LIS
    int lis_solver;
    int lis_pc;
#endif

    /* petsc */
#if USE_PETSC
    struct PTS_DATA_ *petsc;
    int petsc_solver;
    int petsc_pc;
#endif

    /* laspack */
#if USE_LASPACK
    struct LAS_DATA_ *laspack;
    int laspack_solver;
    int laspack_pc;
#endif

    /* mumps */
#if USE_MUMPS
    struct MPS_DATA_ *mumps;
#endif

    /* superlu */
#if USE_SUPERLU
    struct SUPERLU_DATA_ *superlu;
#endif

    /* suitesparse */
#if USE_SSPARSE
    struct SSP_DATA_ *ssparse;
#endif

#if USE_SXAMG
    struct SXAMG_DATA_     *sxamg;
#endif

    /* data */
    lssp_mat_csr A;

    lssp_mat_bcsr Ab;
    int num_blks;
    int *blk_size;

    /* solver info */
    LSSP_SOLVER_TYPE type;
    lssp_vec rhs;
    lssp_vec x;

    /* return */
    double residual;                /* residual, return value */
    int nits;                       /* number of iteraion, return value */

    int verb;
    FILE *log;
    bool assembled;

} LSSP_SOLVER;

#endif
