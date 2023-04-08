
#include "pc-vbiluk.h"

#if USE_ITSOL

extern double lssp_pc_ilut_tol;

void lssp_pc_vbiluk_assemble(LSSP_PC &pc, LSSP_SOLVER s)
{
    lssp_mat_csr A = s.A;
    SparMat *As;
    VBSparMat *vA;
    int bs2;
    int ierr;

    /* convert */
    assert(s.num_blks > 0);
    assert(A.num_rows % s.num_blks == 0);

    As = lssp_malloc<SparMat>(1);
    *As = lssp_mat_csr_to_SparMat(A);
    {
        int nBlk = s.num_blks;

        assert(nBlk > 0);
        assert(s.blk_size != NULL);

        bs2 = 0;
        for (int i = 0; i < nBlk; i++) bs2 += s.blk_size[i];
        assert(bs2 == s.A.num_rows);

        vA = lssp_malloc<VBSparMat>(1);
        csrvbsrC(1, nBlk, s.blk_size, As, vA);
        cleanCS(As);
    }

    /* check pars */
    if (pc.ilut_p <= 0) {
        pc.ilut_p = ((s.A.num_nnzs * 2) / 3 + s.A.num_rows - 1) / s.A.num_rows;
    }

    if (pc.ilut_tol < 0) {
        pc.ilut_tol = lssp_pc_ilut_tol;
    }

    /* level */
    if (pc.iluk_level < 0) pc.iluk_level = 1;

    pc.data = lssp_malloc<VBILUSpar>(1);
    ierr = vbilukC(pc.iluk_level, vA, (VBILUSpar *)pc.data, stdout);
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

    pc.solve = lssp_pc_bilut_solve;
    pc.destroy = lssp_pc_bilut_destroy;
}

#endif
