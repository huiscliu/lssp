#include "pc-arms.h"

#if USE_ITSOL

#define MAX_NUM_LEV 20

typedef struct ARMS_DATA_
{
    arms ArmsSt;

    io_t arms_io;
    int arms_lfil[7];
    double droptol[7], dropcoef[7];
    int ipar[18];
    int diagscal;
    double tolind;

} ARMS_DATA;

/* sets parameters required by arms preconditioner */
static void lssp_pc_arms_pars_init(ARMS_DATA *adata, int verb)
{
    int j;

    /* default */
    adata->diagscal = 1;
    adata->tolind = 0.2;

    /* parameters */
    memset(&adata->arms_io, 0, sizeof(adata->arms_io));

    /* other pars */
    adata->arms_io.perm_type = 0; /* 0: indepent set, 1: ddPQ, 2: acoarsening-based ordering */
    adata->arms_io.Bsize = 50;
    adata->arms_io.lfil0 = 20;
    adata->arms_io.tol0 = 1e-3;
    for (j = 0; j < 17; j++) adata->ipar[j] = 0; 

    /*-------------------- */
    adata->ipar[0]  = MAX_NUM_LEV;     /* max number of levels allowed */
    adata->ipar[1]  = adata->arms_io.perm_type;   /* Indset (0) / PQ (1)    permutation   */

    /* note that these refer to completely  */
    /* different methods for reordering A   */
    /* 0 = standard ARMS independent sets   */
    /* 1 = arms with ddPQ ordering          */
    /* 2 = acoarsening-based ordering [new] */
    adata->ipar[2]  = adata->arms_io.Bsize;       /* smallest size allowed for last schur comp. */

    /* whether or not to print statistics */
    if (verb > 1) {
        adata->ipar[3]  = 1;
    }
    else {
        adata->ipar[3]  = 0;
    }

    /* interlevel methods */  
    adata->ipar[10] = 0;       /* Always do permutations - currently not used  */    
    adata->ipar[11] = 0;       /* ILUT or ILUTP - currently only ILUT is implemented */
    adata->ipar[12] = adata->diagscal;  /* diagonal row scaling before PILUT 0:no 1:yes */
    adata->ipar[13] = adata->diagscal;  /* diagonal column scaling before PILUT 0:no 1:yes */

    /* last level methods */  
    adata->ipar[14] = 1;       /* Always do permutations at last level */
    adata->ipar[15] = 0;       /* ILUTP for last level(0 = ILUT at last level) */
    adata->ipar[16] = adata->diagscal;  /* diagonal row scaling  0:no 1:yes */
    adata->ipar[17] = adata->diagscal;  /* diagonal column scaling  0:no 1:yes */

    /* set lfil */ 
    for (j = 0; j < 7; j++) adata->arms_lfil[j] = adata->arms_io.lfil0;

    /* dropcoef (droptol[k] = tol0*dropcoef[k]) ----- */
    adata->dropcoef[0] = 1.6;     /* dropcoef for L of B block */
    adata->dropcoef[1] = 1.6;     /* dropcoef for U of B block */
    adata->dropcoef[2] = 1.6;     /* dropcoef for L\ inv F */
    adata->dropcoef[3] = 1.6;     /* dropcoef for U\ inv E */
    adata->dropcoef[4] = 0.004;   /* dropcoef for forming schur comple. */
    adata->dropcoef[5] = 0.004;   /* dropcoef for last level L */
    adata->dropcoef[6] = 0.004;   /* dropcoef for last level U */
}

void lssp_pc_arms_create(LSSP_PC &pc)
{
    pc.arms_data = lssp_malloc<ARMS_DATA>(1);
    pc.arms_data->ArmsSt = lssp_malloc<armsMat>(1);

    /* pars */
    lssp_pc_arms_pars_init(pc.arms_data, 0);
}

static void lssp_pc_arms_destroy(LSSP_PC *pc)
{
    arms ArmsSt;

    if (pc == NULL) return;
    if (pc->arms_data == NULL) return;

    ArmsSt = pc->arms_data->ArmsSt;
    cleanARMS(ArmsSt);

    lssp_free(pc->arms_data);
    pc->arms_data = NULL;
}

static void arms_solve(LSSP_PC *s, lssp_vec x, lssp_vec rhs)
{
    arms ArmsSt;
    int n; 

    assert(s != NULL);
    assert(s->arms_data != NULL);

    ArmsSt = s->arms_data->ArmsSt;
    n = ArmsSt->n ;

    memcpy(x.d, rhs.d, n * sizeof(double));  
    armsol2(x.d, ArmsSt)  ;
}

void lssp_pc_arms_assemble(LSSP_PC &pc, LSSP_SOLVER s)
{
    lssp_mat_csr A = s.A;
    SparMat *As;
    arms ArmsSt = NULL;
    int n, nnz, j, lfil0, ierr;
    double tol;

    /* convert */
    As = lssp_malloc<SparMat>(1);
    *As = lssp_mat_csr_to_SparMat(A);

    n = A.num_rows;
    nnz = A.num_nnzs;
    lfil0 = pc.arms_data->arms_io.lfil0;
    tol = pc.arms_data->arms_io.tol0;
    for (j = 0; j < 7; j++) {
        pc.arms_data->arms_lfil[j] = lfil0 * ((nnz + n - 1) / n); 
        pc.arms_data->droptol[j] = tol * pc.arms_data->dropcoef[j];
    }   

    ArmsSt = pc.arms_data->ArmsSt;
    setup_arms(ArmsSt);
    ierr = arms2(As, pc.arms_data->ipar, pc.arms_data->droptol, pc.arms_data->arms_lfil,
            pc.arms_data->tolind, ArmsSt, stdout);
    if(ierr != 0) {
        lssp_printf("pc: arms, error - code %d.\n",ierr);
    }

    pc.solve = arms_solve;
    pc.destroy = lssp_pc_arms_destroy;
    cleanCS(As);
}

#endif
