#include "pc-biluk.h"

#if USE_BLAS
#if USE_LAPACK

#ifdef __cplusplus
extern "C" {
#endif

void FC_FUNC(dgemm,DGEMM)(char *transa, char *transb, int *l, int *m, int *n,  \
        double *alpha, double *a, int *lda, double *b, int *ldb, double *beta,\
        double *c, int *ldc);

void FC_FUNC(dgetrf,DGETRF)(int *m, int *n, double *a, int *lda, int *ipvt, int *info);
void FC_FUNC(dgetri,DGETRI)(int *n, double *a, int *lda, int *ipvt, double *work, \
        int *lwork, int *info);

#ifdef __cplusplus
}
#endif

void lssp_pc_bilu_solve(LSSP_PC *pc, lssp_vec x, lssp_vec rhs)
{
    int i, j;
    int end;
    int n = pc->L.num_rows;
    double result;
    lssp_vec y;
    lssp_vec z;

    y.n = x.n;
    y.d = pc->cache;

    z.n = x.n;
    z.d = pc->cache + n;

    for (i = 0; i < n; i++) {
        end = pc->L.Ap[i + 1] - 1;

        result = rhs.d[i];
        for (j = pc->L.Ap[i]; j < end; j++) {
            result = result - pc->L.Ax[j] * y.d[pc->L.Aj[j]];
        }

        y.d[i] = result / pc->L.Ax[end];
    }

    lssp_mv_mxy(pc->D, y, z);

    for (i = n - 1; i >= 0; i--) {
        end = pc->U.Ap[i];

        result = z.d[i];
        for (j = pc->U.Ap[i + 1] - 1; j > end; j--) {
            result = result - pc->U.Ax[j] * x.d[pc->U.Aj[j]];
        }

        x.d[i] = result / pc->U.Ax[end];
    }
}

static int lssp_pc_bilu_dense_mat_inv(double *A, int n, double *wk, int *ipiv)
{
    assert(n > 0 && A != NULL && wk != NULL && ipiv != NULL);

    if (n == 1) {
        if (A[0] == 0.0)
            return 1;
        else {
            A[0] = 1.0 / A[0];
            return 0;
        }
    }
    else {
        int lWk, info;

        FC_FUNC(dgetrf,DGETRF)(&n, &n, A, &n, ipiv, &info);

        if (info != 0) return info;

        lWk = 10 * n;
        FC_FUNC(dgetri,DGETRI)(&n, A, &n, ipiv, wk, &lWk, &info);

        return info;
    }
}

static void lssp_pc_bilu_dense_gemm(const double *A, const double *B, double alpha, double *C, \
        double beta, int n)
{
    if (n <= 0) return;

    assert(A != NULL && B != NULL && C != NULL);

    if (n == 1) {
        C[0] = A[0] * B[0] * alpha + beta * C[0];
    }
    else {
        FC_FUNC(dgemm,DGEMM)((char *)"n", (char *)"n", &n, &n, &n, &alpha, (double *)A,\
                &n, (double *)B, &n, &beta, C, &n);
    }
}

static void lssp_pc_bilu_bcsr_to_lu(const lssp_mat_bcsr A, const double *inv, lssp_mat_csr &L, \
        lssp_mat_csr &U)
{
    int cntl, cntu;
    int i, j, k;
    double *cache, *d;
    lssp_mat_coo cl, cu;
    int ridx, cidx, bs2;
    int rend;
    int nzl, nzu;

    cntl = cntu = 0;
    for (i = 0; i < A.num_rows; i++) {
        int idx;

        rend = A.Ap[i + 1];
        for (j = A.Ap[i]; j < rend; j++) {
            idx = A.Aj[j];

            if (idx < i) cntl++;
            else if (idx > i) cntu++;
        }
    }

    nzl = A.blk_size * A.blk_size * cntl + A.num_rows * A.blk_size;
    cl.num_rows = cl.num_cols = A.num_rows * A.blk_size;
    cl.num_nnzs = nzl;
    cl.Ai = lssp_malloc<int>(nzl);
    cl.Aj = lssp_malloc<int>(nzl);
    cl.Ax = lssp_malloc<double>(nzl);

    nzu = A.blk_size * A.blk_size * cntu + A.num_rows * A.blk_size;
    cu.num_rows = cu.num_cols = A.num_rows * A.blk_size;
    cu.num_nnzs = nzu;
    cu.Ai = lssp_malloc<int>(nzu);
    cu.Aj = lssp_malloc<int>(nzu);
    cu.Ax = lssp_malloc<double>(nzu);

    bs2 = A.blk_size * A.blk_size;
    nzl = nzu = 0;
    cache = lssp_malloc<double>(bs2);
    d = A.Ax;
    for (i = 0; i < A.num_rows; i++) {
        rend = A.Ap[i + 1];
        for (j = A.Ap[i]; j < rend; j++) {
            cidx = A.Aj[j] * A.blk_size;
            ridx = i * A.blk_size;

            if (A.Aj[j] < i) {
                for (k = 0; k < bs2; k++) {
                    cl.Ai[nzl] = ridx + k % A.blk_size;
                    cl.Aj[nzl] = cidx + k / A.blk_size;
                    cl.Ax[nzl] = *(d + k);

                    nzl++;
                }
            }
            else if (A.Aj[j] > i) {
                lssp_pc_bilu_dense_gemm(inv + bs2 * i, d, 1., cache, 0., A.blk_size);

                for (k = 0; k < bs2; k++) {
                    cu.Ai[nzu] = ridx + k % A.blk_size;
                    cu.Aj[nzu] = cidx + k / A.blk_size;
                    cu.Ax[nzu] = *(cache + k);

                    nzu++;
                }
            }

            d += bs2;
        }
    }

    for (i = 0; i < cl.num_rows; i++) {
        cl.Ai[nzl] = cl.Aj[nzl] = i;
        cl.Ax[nzl] = 1.;
        nzl++;

        cu.Ai[nzu] = cu.Aj[nzu] = i;
        cu.Ax[nzu] = 1.;
        nzu++;
    }

    L = lssp_mat_coo_to_csr(cl);
    U = lssp_mat_coo_to_csr(cu);

    lssp_mat_destroy(cl);
    lssp_mat_destroy(cu);
    lssp_free(cache);

    lssp_mat_sort_column(L);
    lssp_mat_sort_column(U);
}

static void lssp_pc_bilu0_fac(lssp_mat_bcsr &A, double * &inv)
{
    int i, j, k;
    int kk, n = A.num_rows;
    int *Aj;
    int row_end, m;
    double **wk, *d;
    int bs2;
    double *cache;
    double *blk;
    int *ipiv;

    bs2 = A.blk_size * A.blk_size;
    Aj = A.Aj;

    wk = lssp_malloc<double *>(n);
    inv = lssp_malloc<double>(n * bs2);
    cache = lssp_malloc<double>(10 * A.blk_size);
    ipiv = lssp_malloc<int>(A.blk_size);
    blk = lssp_malloc<double>(bs2);

    for (k = 0; k < n; k++) {
        wk[k] = NULL;
    }

    for (i = 0; i < n; i++) {
        row_end = A.Ap[i + 1];

        for (k = A.Ap[i]; k < row_end && Aj[k] < i; k++) {

            int idx = A.Aj[k];
            double *a_ik;

            d = inv + idx * bs2;
            a_ik = A.Ax + k * bs2;

            lssp_memcpy_on<double>(blk, a_ik, bs2);
            lssp_pc_bilu_dense_gemm(blk, d, 1., a_ik, 0., A.blk_size);

            m = A.Ap[idx + 1];
            for (kk = A.Ap[idx]; kk < m; kk++) {
                wk[Aj[kk]] = A.Ax + bs2 * kk;
            }

            for (j = k + 1; j < row_end; j++) {
                double *a_ij = A.Ax + j * bs2;

                if (wk[A.Aj[j]] != NULL) {
                    lssp_pc_bilu_dense_gemm(a_ik, wk[A.Aj[j]], -1., a_ij, 1., A.blk_size);
                }
            }

            for (kk = A.Ap[idx]; kk < m; kk++) {
                wk[Aj[kk]] = NULL;
            }
        }

        if (A.Aj[k] == i) {
            int rt;

            lssp_memcpy_on<double>(inv + i * bs2, A.Ax + bs2 * k, bs2);
            rt = lssp_pc_bilu_dense_mat_inv(inv + i * bs2, A.blk_size, cache, ipiv);

            if (rt) {
                lssp_error(1, "lssp: bilu(0) singular diagonal submatrix.\n");
            }
        }
        else {
            d = inv + i * bs2;
            for (i = 0; i < bs2; i++) d[i] = 0.;

            for (i = 0; i < A.blk_size; i++) d[(A.blk_size + 1) * i] = 1.;
        }
    }

    lssp_free(wk);
    lssp_free(cache);
    lssp_free(ipiv);
    lssp_free(blk);
}

static void lssp_pc_biluk_assemble_diag_mat(lssp_mat_csr &B, int num_rows, int bs, double *inv)
{
    lssp_mat_coo C;
    int i, k;
    int nz, bs2;
    int ridx;
    double *d;

    assert(bs > 0);

    C.num_rows = bs * num_rows;
    C.num_cols = bs * num_rows;
    C.num_nnzs = num_rows * bs * bs;

    C.Ai = lssp_malloc<int>(C.num_nnzs);
    C.Aj = lssp_malloc<int>(C.num_nnzs);
    C.Ax = lssp_malloc<double>(C.num_nnzs);

    bs2 = bs * bs;
    nz = 0;
    d = inv;
    for (i = 0; i < num_rows; i++) {
        ridx = i * bs;
        for (k = 0; k < bs2; k++) {
            C.Ai[nz] = ridx + k % bs;
            C.Aj[nz] = ridx + k / bs;
            C.Ax[nz] = *d;

            nz++;
            d++;
        }
    }

    B = lssp_mat_coo_to_csr(C);
    lssp_mat_destroy(C);
}

void lssp_pc_biluk_destroy(LSSP_PC *pc)
{
    assert(pc->assembled);

    lssp_mat_destroy(pc->L);
    lssp_mat_destroy(pc->U);
    lssp_mat_destroy(pc->D);
    lssp_free(pc->cache);

    pc->assembled = false;
}

static lssp_mat_bcsr lssp_pc_biluk_symbolic__(const lssp_mat_bcsr A, int level, bool sorted)
{
    lssp_mat_bcsr T;
    lssp_mat_csr M;
    lssp_mat_csr S;
    int i, j, k;
    int bs2;

    assert(level >= 0);
    assert(A.num_rows == A.num_cols && A.num_rows > 0);

    if (level == 0 && sorted) {
        T = A;
        T.Ap = lssp_copy_on<int>(A.Ap, A.num_rows + 1);
        T.Aj = lssp_copy_on<int>(A.Aj, A.num_nnzs);
        T.Ax = lssp_copy_on<double>(A.Ax, A.num_nnzs * A.blk_size * A.blk_size);

        return T;
    }

    M.num_rows = M.num_cols = A.num_rows;
    M.num_nnzs = A.num_nnzs;
    M.Ap = lssp_copy_on<int>(A.Ap, A.num_rows + 1);
    M.Aj = lssp_copy_on<int>(A.Aj, A.num_nnzs);
    M.Ax = lssp_malloc<double>(A.num_nnzs);

    if (!sorted) lssp_mat_sort_column(M);

    S = lssp_pc_iluk_symbolic_fac(M, level, false);

    T = A;
    T.Ap = lssp_copy_on<int>(S.Ap, S.num_rows + 1);
    T.Aj = lssp_copy_on<int>(S.Aj, S.num_nnzs);
    T.Ax = lssp_malloc<double>(S.num_nnzs * A.blk_size * A.blk_size);

    k = S.num_nnzs * A.blk_size * A.blk_size;
    for (i = 0; i < k; i++) T.Ax[i] = 0.;

    bs2 = A.blk_size * A.blk_size;
    for (i = 0; i < T.num_rows; i++) {
        int row_end = A.Ap[i + 1];

        for (j = A.Ap[i]; j < row_end; j++) {
            int *p;
            int pst, m;

            m = A.Aj[j];
            p = (int *)bsearch(&m, T.Aj+T.Ap[i], T.Ap[i+1]-T.Ap[i], sizeof(int), lssp_comp_int_asc);
            assert(p != NULL);
            pst = (p - T.Aj) * bs2;
            memcpy(T.Ax + pst, A.Ax + bs2 * j, bs2 * sizeof(double));
        }
    }

    lssp_mat_destroy(M);
    lssp_mat_destroy(S);

    return T;
}

void lssp_pc_biluk_assemble_mat(LSSP_PC &pc, lssp_mat_bcsr A)
{
    double *inv;
    lssp_mat_bcsr T;
    double time = lssp_get_time();
    bool A_sorted;
    int bcsr_block_size;

    A_sorted = lssp_mat_bcsr_is_sorted(A);
    T = lssp_pc_biluk_symbolic__(A, pc.iluk_level, A_sorted);

    bcsr_block_size = A.blk_size;
    assert(bcsr_block_size > 0);

    lssp_pc_bilu0_fac(T, inv);
    lssp_pc_bilu_bcsr_to_lu(T, inv, pc.L, pc.U);
    lssp_pc_biluk_assemble_diag_mat(pc.D, T.num_rows, bcsr_block_size, inv);

    pc.cache = lssp_malloc<double>(2 * pc.L.num_rows);
    pc.solve = lssp_pc_bilu_solve;

    pc.destroy = lssp_pc_biluk_destroy;

    lssp_free(inv);
    lssp_mat_destroy(T);

    time = lssp_get_time() - time;
    if (pc.verb > 0) lssp_printf("pc: BILUK assemble time: %f\n", time);
}

void lssp_pc_biluk_assemble(LSSP_PC &pc, LSSP_SOLVER s)
{
    lssp_mat_bcsr A;

    assert(s.A.num_rows == s.A.num_cols);
    assert(s.A.num_rows > 0 && s.A.num_nnzs > 0);
    assert(s.num_blks > 0);
    assert(s.A.num_rows % s.num_blks == 0);

    A = lssp_mat_csr_to_bcsr(s.A, s.A.num_rows / s.num_blks);
    lssp_pc_biluk_assemble_mat(pc, A);

    lssp_mat_destroy(A);
}

#endif
#endif
