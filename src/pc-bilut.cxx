#include "pc-bilut.h"

#if USE_ITSOL

extern double lssp_pc_ilut_tol;

void lssp_pc_bilut_destroy(LSSP_PC *pc)
{
    cleanVBILU((VBILUSpar *)pc->data);
}

SparMat lssp_mat_csr_to_SparMat(lssp_mat_csr A)
{
    SparMat M;
    int i, j, k;
    int n, row_end;

    M.n = n = A.num_rows;
    M.nzcount = lssp_malloc<int>(n);
    M.ja = lssp_malloc<int *>(n);
    M.ma = lssp_malloc<double *>(n);

    for (i = 0; i < n; i++) {
        M.nzcount[i] = A.Ap[i + 1] - A.Ap[i];
        M.ja[i] = lssp_malloc<int>(M.nzcount[i]);
        M.ma[i] = lssp_malloc<double>(M.nzcount[i]);
    }

    for (i = 0; i < n; i++) {
        row_end = A.Ap[i + 1];

        k = 0;
        for (j = A.Ap[i]; j < row_end; j++) {
            M.ja[i][k] = A.Aj[j];
            M.ma[i][k] = A.Ax[j];
            k++;
        }
    }

    return M;
}

void lssp_pc_bilut_solve(LSSP_PC *s, lssp_vec x, lssp_vec rhs)
{
    vblusolC(rhs.d, x.d, (VBILUSpar *)s->data);
}

void lssp_pc_bilut_assemble(LSSP_PC &pc, LSSP_SOLVER s)
{
    lssp_mat_csr A = s.A;
    SparMat *As;
    VBSparMat *vA;
    BData *w;
    int bs2;
    int ierr;
    int bcsr_block_size;

    /* convert */
    assert(s.num_blks > 0);
    assert(A.num_rows % s.num_blks == 0);

    As = lssp_malloc<SparMat>(1);
    *As = lssp_mat_csr_to_SparMat(A);
    {
        int nBlk = s.num_blks;
        int *nB;

        bcsr_block_size = A.num_rows / s.num_blks;
        assert(nBlk > 0);
        nB = lssp_malloc<int>(nBlk);
        for (int i = 0; i < nBlk; i++) nB[i] = bcsr_block_size;

        vA = lssp_malloc<VBSparMat>(1);
        csrvbsrC(1, nBlk, nB, As, vA);
        lssp_free(nB);
        cleanCS(As);
    }

    /* working array */
    w = lssp_malloc<BData>(vA->n);
    bs2 = bcsr_block_size * bcsr_block_size;
    for (int i = 0; i < vA->n; i++) w[i] = lssp_malloc<double>(bs2);

    /* check pars */
    if (pc.ilut_p <= 0) {
        pc.ilut_p = ((s.A.num_nnzs * 2) / 3 + s.A.num_rows - 1) / s.A.num_rows;
    }

    if (pc.ilut_tol < 0) {
        pc.ilut_tol = lssp_pc_ilut_tol;
    }

    pc.data = lssp_malloc<VBILUSpar>(1);
    ierr = vbilutC(vA, (VBILUSpar *)pc.data, pc.ilut_p, pc.ilut_tol, w, stdout);
    cleanVBMat(vA);

    if (ierr == -2) {
        lssp_printf("bilut: singular diagonal block\n" );
        cleanVBILU((VBILUSpar *)pc.data);
        exit(-1);
    }
    else if (ierr != 0) {
        lssp_printf("bilut: error, ierr != 0.\n" );
        exit(-1);
    }

    for (int i = 0; i < vA->n; i++) lssp_free(w[i]);
    lssp_free(w);

    pc.solve = lssp_pc_bilut_solve;
    pc.destroy = lssp_pc_bilut_destroy;
}

#endif
