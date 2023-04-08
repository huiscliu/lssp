#include "pc-iluk.h"

extern double mat_zero_diag_value;
extern double mat_zero_diag_tol;
extern int lssp_pc_iluk_level_default;

typedef struct {
    int nrows;
    int *nnz;
    int **Aj;

} mat_csr_ys;

typedef struct {
    mat_csr_ys L;
    mat_csr_ys U;
    double     *D;
    int        nrows;

} mat_ldu;

static void lssp_pc_iluk_symbolic__(mat_ldu &lu, mat_csr_ys A, int level)
{
    int n = A.nrows;
    int *levls = NULL, *jbuf = NULL, *iw;
    int **ulvl;
    mat_csr_ys L = lu.L, U = lu.U;
    int i, j, k, col, ip, it, jpiv;
    int incl, incu, jmin, kmin;

    assert(level >= 0);

    levls = lssp_malloc<int>(n);
    jbuf = lssp_malloc<int>(n);
    ulvl = lssp_malloc<int *>(n);
    iw = lssp_malloc<int>(n);

    for (j = 0; j < n; j++) iw[j] = -1;

    for (i = 0; i < n; i++) {
        incl = 0;
        incu = i;

        for (j = 0; j < A.nnz[i]; j++) {
            col = A.Aj[i][j];
            if (col < i) {
                jbuf[incl] = col;
                levls[incl] = 0;
                iw[col] = incl++;
            }
            else if (col > i) {
                jbuf[incu] = col;
                levls[incu] = 0;
                iw[col] = incu++;
            }
        }

        jpiv = -1;
        while (++jpiv < incl) {
            k = jbuf[jpiv];

            kmin = k;
            jmin = jpiv;
            for (j = jpiv + 1; j < incl; j++) {
                if (jbuf[j] < kmin) {
                    kmin = jbuf[j];
                    jmin = j;
                }
            }

            if (jmin != jpiv) {
                jbuf[jpiv] = kmin;
                jbuf[jmin] = k;
                iw[kmin] = jpiv;
                iw[k] = jmin;
                j = levls[jpiv];
                levls[jpiv] = levls[jmin];
                levls[jmin] = j;
                k = kmin;
            }

            for (j = 0; j < U.nnz[k]; j++) {
                col = U.Aj[k][j];
                it = ulvl[k][j] + levls[jpiv] + 1;
                if (it > level)
                    continue;
                ip = iw[col];
                if (ip == -1) {
                    if (col < i) {
                        jbuf[incl] = col;
                        levls[incl] = it;
                        iw[col] = incl++;
                    }
                    else if (col > i) {
                        jbuf[incu] = col;
                        levls[incu] = it;
                        iw[col] = incu++;
                    }
                }
                else {
                    if (levls[ip] < it) levls[ip] = it;
                }
            }
        }

        for (j = 0; j < incl; j++) iw[jbuf[j]] = -1;
        for (j = i; j < incu; j++) iw[jbuf[j]] = -1;

        L.nnz[i] = incl;
        if (incl > 0) {
            L.Aj[i] = lssp_malloc<int>(incl);
            memcpy(L.Aj[i], jbuf, sizeof(int) * incl);
        }

        k = incu - i;
        U.nnz[i] = k;
        if (k > 0) {
            U.Aj[i] = lssp_malloc<int>(k);
            memcpy(U.Aj[i], jbuf + i, sizeof(int) * k);

            ulvl[i] = lssp_malloc<int>(k);
            memcpy(ulvl[i], levls + i, k * sizeof(int));
        }
    }

    for (i = 0; i < n - 1; i++) {
        if (U.nnz[i])
            lssp_free(ulvl[i]);
    }

    lssp_free(ulvl);
    lssp_free(levls);
    lssp_free(jbuf);
    lssp_free(iw);
}

static mat_csr_ys lssp_pc_iluk_csr_to_matcsrys(lssp_mat_csr A)
{
    mat_csr_ys M;
    int i, j, k;
    int n, row_end;

    M.nrows = n = A.num_rows;
    M.nnz = lssp_malloc<int>(n);
    M.Aj = lssp_malloc<int *>(n);

    for (i = 0; i < n; i++) {
        M.nnz[i] = A.Ap[i + 1] - A.Ap[i];
        M.Aj[i] = lssp_malloc<int>(M.nnz[i]);
    }

    for (i = 0; i < n; i++) {
        row_end = A.Ap[i + 1];

        k = 0;
        for (j = A.Ap[i]; j < row_end; j++) {
            M.Aj[i][k++] = A.Aj[j];
        }
    }

    return M;
}

static void lssp_pc_iluk_assemble_mat_ldu(mat_ldu &A, int n)
{
    int i;

    assert(n > 0);

    A.nrows = n;
    A.D = lssp_malloc<double>(n);

    A.L.nrows = n;
    A.L.nnz = lssp_malloc<int>(n);
    A.L.Aj = lssp_malloc<int *>(n);

    A.U.nrows = n;
    A.U.nnz = lssp_malloc<int>(n);
    A.U.Aj = lssp_malloc<int *>(n);

    for (i = 0; i < n; i++) {
        A.L.Aj[i] = A.U.Aj[i] = NULL;
    }
}

static lssp_mat_csr lssp_pc_iluk_mat_ldu_to_mat_csr(mat_ldu A)
{
    lssp_mat_csr M;
    int i, j;
    int n, k;
    int nonzeros = 0;

    M.num_rows = M.num_cols = n = A.nrows;
    M.Ap = lssp_malloc<int>(n + 1);

    for (i = 0; i < n; i++) {
        M.Ap[i] = A.L.nnz[i] + A.U.nnz[i] + 1;
        nonzeros += M.Ap[i];
    }

    M.num_nnzs = nonzeros;

    k = M.Ap[0];
    M.Ap[0] = 0;
    for (i = 1; i <= n; i++) {
        j = M.Ap[i];
        M.Ap[i] = M.Ap[i - 1] + k;
        k = j;
    }

    M.Aj = lssp_malloc<int>(nonzeros);
    M.Ax = lssp_malloc<double>(nonzeros);

    k = 0;
    for (i = 0; i < n; i++) {

        for (j = 0; j < A.L.nnz[i]; j++) {
            M.Aj[k++] = A.L.Aj[i][j];
        }

        M.Aj[k++] = i;

        for (j = 0; j < A.U.nnz[i]; j++) {
            M.Aj[k++] = A.U.Aj[i][j];
        }
    }

    return M;
}

static void lssp_pc_iluk_csr_set_value_by_mat(lssp_mat_csr &A, lssp_mat_csr V)
{
    int i, j;
    int row_end, n;
    double *wk;

    assert(A.num_rows == V.num_rows);
    assert(A.num_cols == V.num_cols);
    assert(A.num_rows == A.num_cols);
    assert(A.num_rows > 0);
    assert(A.num_nnzs >= 0);

    if (A.num_nnzs > 0) {
        assert(A.Ap != NULL);
        assert(A.Aj != NULL);
        assert(A.Ax != NULL);

        assert(V.Ap != NULL);
        assert(V.Aj != NULL);
        assert(V.Ax != NULL);
    }

    n = A.num_rows;
    wk = lssp_malloc<double>(A.num_rows);

    for (i = 0; i < n; i++) wk[i] = 0;

    for (i = 0; i < n; i++) {

        row_end = V.Ap[i + 1];
        for (j = V.Ap[i]; j < row_end; j++) {
            wk[V.Aj[j]] = V.Ax[j];
        }

        row_end = A.Ap[i + 1];
        for (j = A.Ap[i]; j < row_end; j++) {
            A.Ax[j] = wk[A.Aj[j]];
        }

        row_end = V.Ap[i + 1];
        for (j = V.Ap[i]; j < row_end; j++) {
            wk[V.Aj[j]] = 0;
        }
    }

    lssp_free(wk);
}

lssp_mat_csr lssp_pc_iluk_symbolic_fac(const lssp_mat_csr A, int level, bool set_value)
{
    lssp_mat_csr M;
    mat_csr_ys C;
    mat_ldu I;
    int i, n = A.num_rows;

    if (level < 0) {
        lssp_warning("pc: ILUK, level may be too small: %d\n", level);
        lssp_warning("pc: ILUK, change the level to 0\n");
        level = 0;
    }

    if (level > 8) {
        lssp_warning("pc: ILUK, level may be too large: %d\n", level);
    }

    assert(A.num_rows == A.num_cols);
    assert(A.num_rows > 0);
    assert(A.num_nnzs > 0);
    assert(A.Ap != NULL);
    assert(A.Aj != NULL);
    assert(A.Ax != NULL);

    if (level == 0) {
        M = A;
        M.Ap = lssp_copy_on<int>(A.Ap, A.num_rows + 1);
        M.Aj = lssp_copy_on<int>(A.Aj, A.num_nnzs);
        M.Ax = lssp_copy_on<double>(A.Ax, A.num_nnzs);

        return M;
    }

    /* symbolic fac */
    C = lssp_pc_iluk_csr_to_matcsrys(A);
    lssp_pc_iluk_assemble_mat_ldu(I, n);
    lssp_pc_iluk_symbolic__(I, C, level);
    M = lssp_pc_iluk_mat_ldu_to_mat_csr(I);

    for (i = 0; i < n; i++) {
        lssp_free(C.Aj[i]);
    }

    lssp_free(C.Aj);
    lssp_free(C.nnz);

    for (i = 0; i < n; i++) {
        lssp_free(I.L.Aj[i]);
    }

    lssp_free(I.L.Aj);
    lssp_free(I.L.nnz);

    for (i = 0; i < n; i++) {
        lssp_free(I.U.Aj[i]);
    }

    lssp_free(I.U.Aj);
    lssp_free(I.U.nnz);
    lssp_free(I.D);

    if (set_value) lssp_pc_iluk_csr_set_value_by_mat(M, A);

    lssp_mat_sort_column(M);

    return M;
}

static void lssp_pc_ilu0_fac(lssp_mat_csr A)
{
    int i, j, k;
    int kk, n = A.num_rows;
    double *diag;
    int *C, *Ap;
    double *Ax = A.Ax;
    int row_end, m;
    double *wk;

    wk = lssp_malloc<double>(n);
    diag = lssp_malloc<double>(n);
    C = A.Aj;
    Ap = A.Ap;

    for (kk = 0; kk < n; kk++) {
        wk[kk] = 0;
    }

    diag[0] = Ax[0];
    if (fabs(diag[0]) < mat_zero_diag_tol) {
        if (diag[0] > 0) {
            diag[0] = mat_zero_diag_value;
        }
        else {
            diag[0] = - mat_zero_diag_value;
        }
    }
    diag[0] = 1. / diag[0];

    for (i = 1; i < n; i++) {
        row_end = Ap[i + 1];

        for (k = Ap[i]; C[k] < i; k++) {
            double a_ik;

            m = Ap[C[k] + 1];
            for (kk = Ap[C[k]]; kk < m; kk++) {
                wk[C[kk]] = Ax[kk];
            }

            a_ik = Ax[k] = Ax[k] * diag[C[k]];

            for (j = k + 1; j < row_end; j++) {
                if (wk[C[j]] != 0.) Ax[j] = Ax[j] - a_ik * wk[C[j]];
            }

            for (kk = Ap[C[k]]; kk < m; kk++) {
                wk[C[kk]] = 0;
            }
        }

        diag[i] = mat_zero_diag_value;
        if (A.Aj[k] == i) {
            if (fabs(A.Ax[k]) < mat_zero_diag_tol) A.Ax[k] = mat_zero_diag_value;
            diag[i] = A.Ax[k];
        }
        diag[i] = 1. / diag[i];
    }

    lssp_free(wk);
    lssp_free(diag);
}

static void lssp_pc_iluk_assemble_matrix(const lssp_mat_csr A, int blk_size, int level,
        lssp_mat_csr &L, lssp_mat_csr &U, int verb)
{
    int i, j, k, m, itr;
    int blk_num;
    int start, end;
    lssp_mat_csr B;
    lssp_mat_csr csr;
    lssp_mat_csr M;
    int nz_l, nz_u;
    double *wk;

    assert(A.num_nnzs > 0);
    assert(A.num_rows > 0);
    assert(A.num_cols > 0);
    assert(A.num_rows == A.num_cols);
    assert(blk_size > 0);

    blk_num = (A.num_rows + blk_size - 1) / blk_size;

    if (level > 0) {
        double t = lssp_get_time();

        csr = lssp_pc_iluk_symbolic_fac(A, level, true);
        t = lssp_get_time() -t;

        if (verb > 2) {
            lssp_printf("pc: time for symbolic factorization: %f\n", t);
        }

        M = lssp_mat_get_block_diag(csr, blk_size);

        lssp_mat_destroy(csr);
    }
    else {
        M = lssp_mat_get_block_diag(A, blk_size);
    }

    L.num_rows = U.num_rows = M.num_rows;
    L.num_cols = U.num_cols = M.num_cols;

    nz_l = M.num_nnzs;
    L.Ap = lssp_malloc<int>(L.num_rows + 1);
    L.Ax = lssp_malloc<double>(nz_l);
    L.Aj = lssp_malloc<int>(nz_l);

    nz_u = M.num_nnzs;
    U.Ap = lssp_malloc<int>(U.num_rows + 1);
    U.Ax = lssp_malloc<double>(nz_u);
    U.Aj = lssp_malloc<int>(nz_u);

    wk = lssp_malloc<double>(blk_size * 2);

    nz_l = nz_u = 0;
    L.Ap[0] = U.Ap[0] = 0;
    itr = 0;
    for (i = 0; i < blk_num; i++) {
        int offset;

        start = i * blk_size;
        end = (i + 1) * blk_size;
        if (end >= A.num_rows) end = A.num_rows;

        B.num_cols = B.num_rows = end - start;
        B.num_nnzs = M.Ap[end] - M.Ap[start];
        B.Ap = M.Ap + start;
        B.Ax = M.Ax + M.Ap[start];
        B.Aj = M.Aj + M.Ap[start];

        if (start > 0) {
            for (j = 0; j < B.num_nnzs; j++) {
                B.Aj[j] -= start;
            }
        }

        offset = M.Ap[start];
        if (offset > 0) {
            for (j = start; j <= end; j++) {
                M.Ap[j] -= offset;
            }
        }

        lssp_pc_ilu0_fac(B);

        if (start > 0) {
            for (j = 0; j < B.num_nnzs; j++) {
                B.Aj[j] += start;
            }
        }

        for (j = 0; j < B.num_rows; j++) {
            int r_idx = start + j;

            m = B.Ap[j + 1];
            for (k = B.Ap[j]; k < m; k++) {
                if (B.Aj[k] < r_idx) {
                    L.Aj[nz_l] = B.Aj[k];
                    L.Ax[nz_l] = B.Ax[k];
                    nz_l++;
                }
                else if (B.Aj[k] == r_idx) {
                    L.Aj[nz_l] = B.Aj[k];
                    L.Ax[nz_l] = 1;
                    nz_l++;

                    U.Aj[nz_u] = B.Aj[k];
                    U.Ax[nz_u] = B.Ax[k];
                    nz_u++;

                }
                else {
                    U.Aj[nz_u] = B.Aj[k];
                    U.Ax[nz_u] = B.Ax[k];
                    nz_u++;
                }
            }

            itr++;

            L.Ap[itr] = nz_l;
            U.Ap[itr] = nz_u;
        }

        M.Ap[end] += offset;
    }

    L.num_nnzs = L.Ap[L.num_rows];
    L.Aj = lssp_realloc<int>(L.Aj, L.num_nnzs);
    L.Ax = lssp_realloc<double>(L.Ax, L.num_nnzs);

    U.num_nnzs = U.Ap[U.num_rows];
    U.Aj = lssp_realloc<int>(U.Aj, U.num_nnzs);
    U.Ax = lssp_realloc<double>(U.Ax, U.num_nnzs);

    lssp_mat_destroy(M);
    lssp_free<double>(wk);

    if (verb > 2) {
        lssp_printf("pc: iluk, L and U filled: %d, ratio: %f\n", L.num_nnzs + U.num_nnzs, \
                (L.num_nnzs + U.num_nnzs) * 1. / A.num_nnzs);
    }
}

void lssp_pc_iluk_destroy(LSSP_PC *pc)
{
    assert(pc->assembled);

    lssp_mat_destroy(pc->L);
    lssp_mat_destroy(pc->U);

    lssp_free<double>(pc->cache);

    pc->assembled = false;
}

void lssp_pc_iluk_assemble(LSSP_PC &pc, LSSP_SOLVER s)
{
    lssp_mat_csr A;

    assert(s.A.num_rows == s.A.num_cols);
    assert(s.A.num_rows > 0 && s.A.num_nnzs > 0);

    A = lssp_mat_adjust_zero_diag(s.A, mat_zero_diag_tol);
    lssp_pc_iluk_assemble_matrix(A, A.num_rows, pc.iluk_level, pc.L, pc.U, pc.verb);
    lssp_mat_destroy(A);

    pc.cache = lssp_malloc<double>(pc.L.num_rows);
    pc.solve = lssp_pc_ilu_solve;

    pc.destroy = lssp_pc_iluk_destroy;
}

void lssp_pc_iluk_set_level(LSSP_PC &pc, int level)
{
    if (level < 0) {
        lssp_warning("pc: level is too small, set it to %d!\n", lssp_pc_iluk_level_default);
        pc.iluk_level = lssp_pc_iluk_level_default;
    }
    else {
        pc.iluk_level = level;
    }
}
