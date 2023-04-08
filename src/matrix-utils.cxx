#include "matrix-utils.h"

void lssp_mat_init(lssp_mat_csr &A)
{
    bzero(&A, sizeof(A));
}

void lssp_mat_init(lssp_mat_coo &A)
{
    bzero(&A, sizeof(A));
}

void lssp_mat_init(lssp_mat_bcsr &A)
{
    bzero(&A, sizeof(A));
}

void lssp_mat_destroy(lssp_mat_csr &A)
{
    lssp_free(A.Ap);
    lssp_free(A.Aj);
    lssp_free(A.Ax);
    lssp_mat_init(A);
}

void lssp_mat_destroy(lssp_mat_coo &A)
{
    lssp_free(A.Ai);
    lssp_free(A.Aj);
    lssp_free(A.Ax);
    lssp_mat_init(A);
}

void lssp_mat_destroy(lssp_mat_bcsr &A)
{
    lssp_free(A.Aj);
    lssp_free(A.Ap);
    lssp_free(A.Ax);
    lssp_mat_init(A);
}

lssp_mat_csr lssp_mat_create(int nrows, int ncols, int *Ap, int *Aj, double *Ax)
{
    lssp_mat_csr A;

    assert(nrows > 0);
    assert(ncols > 0);
    assert(Ap != NULL);

    lssp_mat_init(A);
    A.num_rows = nrows;
    A.num_cols = ncols;
    A.num_nnzs = Ap[nrows];

    A.Ap = lssp_copy_on<int>(Ap, nrows + 1);
    A.Aj = lssp_copy_on<int>(Aj, Ap[nrows]);
    A.Ax = lssp_copy_on<double>(Ax, Ap[nrows]);

    return A;
}

lssp_mat_bcsr lssp_mat_csr_to_bcsr(const lssp_mat_csr A, int bs)
{
    lssp_mat_bcsr B;
    int *wk, *blk, nz, bs2;
    int i, j, k;
    int Aend, rend;

    assert(A.num_rows > 0 && A.num_rows == A.num_cols && A.num_nnzs > 0);
    assert(bs > 0);

    if (A.num_rows % bs != 0) {
        lssp_error(1, "num_rows is not a multiple of block size: %d\n", bs);
    }

    lssp_mat_init(B);
    B.num_rows = B.num_cols = A.num_rows / bs;
    B.blk_size = bs;
    B.Ap = lssp_malloc<int>(B.num_rows + 1);
    B.Aj = lssp_malloc<int>(A.num_nnzs);
    wk = lssp_malloc<int>(B.num_rows);
    blk = lssp_malloc<int>(B.num_rows);

    B.Ap[0] = 0;
    Aend = 0;
    for (i = 0; i < B.num_rows; i++) {
        blk[i] = -1;
    }

    for (i = 0; i < B.num_rows; i++) {
        Aend += bs;
        nz = 0;
        for (j = i * bs; j < Aend; j++) {
            int idx;

            rend = A.Ap[j + 1];
            for (k = A.Ap[j]; k < rend; k++) {
                idx = A.Aj[k] / bs;

                if (blk[idx] == -1) {
                    blk[idx] = 1;
                    wk[nz] = idx;
                    nz++;
                }
            }
        }

        qsort(wk, nz, sizeof(int), lssp_comp_int_asc);
        lssp_memcpy_on<int>(B.Aj + B.Ap[i], wk, nz);

        for (j = 0; j < nz; j++) {
            blk[wk[j]] = -1;
        }

        B.Ap[i + 1] = B.Ap[i] + nz;
    }

    B.num_nnzs = B.Ap[B.num_rows];
    B.Aj = lssp_realloc<int>(B.Aj, B.num_nnzs);
    B.Ax = lssp_malloc<double>(B.num_nnzs * bs * bs);

    k = B.num_nnzs * bs * bs;
    for (i = 0; i < k; i++) B.Ax[i] = 0.;

    Aend = 0;
    bs2 = bs * bs;
    for (i = 0; i < B.num_rows; i++) {
        double *d = B.Ax + bs2 * B.Ap[i];
        double *p;

        rend = B.Ap[i + 1];
        nz = 0;
        for (j = B.Ap[i]; j < rend; j++) {
            blk[B.Aj[j]] = nz++;
        }

        Aend += bs;
        for (j = i * bs; j < Aend; j++) {
            int ridx, cidx;

            rend = A.Ap[j + 1];
            for (k = A.Ap[j]; k < rend; k++) {
                p = d + blk[A.Aj[k] / bs] * bs2;

                ridx = j % bs;
                cidx = A.Aj[k] % bs;

                *(p + cidx * bs + ridx) = A.Ax[k];
            }
        }

        rend = B.Ap[i + 1];
        for (j = B.Ap[i]; j < rend; j++) {
            blk[B.Aj[j]] = -1;
        }
    }

    lssp_free(blk);
    lssp_free(wk);

    return B;
}

lssp_mat_csr lssp_mat_bcsr_to_csr(const lssp_mat_bcsr A)
{
    lssp_mat_csr B;
    lssp_mat_coo C;
    int i, j, k;
    int bs2, Aend;
    int ridx, cidx;
    int nz;
    double *d;

    C.num_rows = A.blk_size * A.num_rows;
    C.num_cols = A.blk_size * A.num_cols;
    C.num_nnzs = A.num_nnzs * A.blk_size * A.blk_size;

    C.Ai = lssp_malloc<int>(C.num_nnzs);
    C.Aj = lssp_malloc<int>(C.num_nnzs);
    C.Ax = lssp_malloc<double>(C.num_nnzs);

    bs2 = A.blk_size * A.blk_size;
    nz = 0;
    d = A.Ax;
    for (i = 0; i < A.num_rows; i++) {
        Aend = A.Ap[i + 1];
        for (j = A.Ap[i]; j < Aend; j++) {
            cidx = A.Aj[j] * A.blk_size;
            ridx = i * A.blk_size;

            for (k = 0; k < bs2; k++) {
                if (fabs(*d) > 0.) {
                    C.Ai[nz] = ridx + k % A.blk_size;
                    C.Aj[nz] = cidx + k / A.blk_size;
                    C.Ax[nz] = *d;

                    nz++;
                }

                d++;
            }
        }
    }

    C.num_nnzs = nz;
    C.Ai = lssp_realloc<int>(C.Ai, nz);
    C.Aj = lssp_realloc<int>(C.Aj, nz);
    C.Ax = lssp_realloc<double>(C.Ax, nz);

    B = lssp_mat_coo_to_csr(C);
    lssp_mat_destroy(C);
    lssp_mat_sort_column(B);

    return B;
}

bool lssp_mat_bcsr_is_sorted(const lssp_mat_bcsr A)
{
    bool rtn = true;
    int i, j, row_end;

    assert(A.num_rows > 0 && A.num_cols > 0 && A.blk_size > 0);

    if (A.num_nnzs <= 0) {
        return true;
    }
    else {
        assert(A.Ap != NULL && A.Aj != NULL && A.Ax != NULL);
    }

    for (i = 0; i < A.num_rows; i++) {
        row_end = A.Ap[i + 1];

        if (row_end - A.Ap[i] > 1) {
            for (j = A.Ap[i] + 1; j < row_end; j++) {
                if (A.Aj[j - 1] > A.Aj[j]) {
                    rtn = false;
                    break;
                }
            }
        }

        if (!rtn) break;
    }

    return rtn;
}

bool lssp_mat_csr_is_sorted(const lssp_mat_csr A)
{
    bool rtn = true;
    int i, j, row_end;

    assert(A.num_rows > 0 && A.num_cols > 0);

    if (A.num_nnzs <= 0) {
        return true;
    }
    else {
        assert(A.Ap != NULL && A.Aj != NULL && A.Ax != NULL);
    }

    for (i = 0; i < A.num_rows; i++) {
        row_end = A.Ap[i + 1];

        if (row_end - A.Ap[i] > 1) {
            for (j = A.Ap[i] + 1; j < row_end; j++) {
                if (A.Aj[j - 1] > A.Aj[j]) {
                    rtn = false;
                    break;
                }
            }
        }

        if (!rtn) break;
    }

    return rtn;
}

static void lssp_mat_csr_to_coo(const int *Ap, const int *Aj, const double *Ax, const int num_rows, \
        const int num_nnzs, int *rows, int *cols, double *data)
{
    int i, j;
    int row_start, row_end;

    for(i = 0; i < num_rows; i++){
        row_start = Ap[i];
        row_end = Ap[i + 1];

        for(j = row_start; j < row_end; j++){
            rows[j] = i;
        }
    }

    for(i = 0; i < num_nnzs; i++){
        cols[i] = Aj[i];
        data[i] = Ax[i];
    }
}

lssp_mat_coo lssp_mat_csr_to_coo(const lssp_mat_csr csr)
{
    lssp_mat_coo A;

    lssp_mat_init(A);
    A.num_rows = csr.num_rows;
    A.num_cols = csr.num_cols;
    A.num_nnzs = csr.num_nnzs;

    if (csr.num_nnzs <= 0) {
        return A;
    }

    A.Ai = lssp_malloc<int>(csr.num_nnzs);
    A.Aj = lssp_malloc<int>(csr.num_nnzs);
    A.Ax = lssp_malloc<double>(csr.num_nnzs);

    lssp_mat_csr_to_coo(csr.Ap, csr.Aj, csr.Ax, A.num_rows, A.num_nnzs, A.Ai, A.Aj, A.Ax);

    return A;
}

static void lssp_mat_coo_to_csr(const int *rows, const int *cols, const double *data, \
        const int num_rows, const int num_nnzs, int *Ap, int *Aj, double *Ax)
{
    int i;
    int last, temp;
    int cumsum;

    for (i = 0; i < num_rows; i++) Ap[i] = 0;

    for (i = 0; i < num_nnzs; i++) Ap[rows[i]]++;

    for(i = 0, cumsum = 0; i < num_rows; i++){
        temp = Ap[i];

        Ap[i] = cumsum;
        cumsum += temp;
    }

    Ap[num_rows] = num_nnzs;

    for(i = 0; i < num_nnzs; i++){
        int row = rows[i];
        int dest = Ap[row];

        Aj[dest] = cols[i];
        Ax[dest] = data[i];

        Ap[row]++;
    }

    for(i = 0, last = 0; i <= num_rows; i++){
        temp = Ap[i];
        Ap[i] = last;
        last = temp;
    }
}

lssp_mat_csr lssp_mat_coo_to_csr(const lssp_mat_coo A)
{
    lssp_mat_csr csr;

    lssp_mat_init(csr);
    csr.num_rows = A.num_rows;
    csr.num_cols = A.num_cols;
    csr.num_nnzs = A.num_nnzs;

    if (A.num_nnzs <= 0) return csr;

    csr.Ap = lssp_malloc<int>(csr.num_rows + 1);
    csr.Aj = lssp_malloc<int>(csr.num_nnzs);
    csr.Ax = lssp_malloc<double>(csr.num_nnzs);

    lssp_mat_coo_to_csr(A.Ai, A.Aj, A.Ax, A.num_rows, A.num_nnzs, csr.Ap, csr.Aj, csr.Ax);
    csr.num_nnzs = csr.Ap[csr.num_rows];

    return csr;
}

static int lssp_mat_comp_row_by_col_index(const void *p, const void *n)
{
    return *(int *)p - *(int *)n;
}

void lssp_mat_sort_column(lssp_mat_csr &A)
{
    int bs;
    int *start;
    int threads;

    if (A.num_rows <= 0) {
        lssp_warning("Mat: wrong number of rows: %d %s %d\n",
                A.num_rows, __FUNCTION__, __LINE__);
        return;
    }

    if (A.num_cols <= 0) {
        lssp_warning("Mat: wrong number of columns: %d %s %d\n",
                A.num_cols, __FUNCTION__, __LINE__);
        return;
    }

    if (A.num_nnzs <= 0) {
        lssp_warning("Mat: has no non-zeros: %d %s %d\n",
                A.num_nnzs, __FUNCTION__, __LINE__);
        return;
    }

    assert(A.Ap != NULL);
    assert(A.Aj != NULL);
    assert(A.Ax != NULL);

    threads = 1;

    start = lssp_malloc<int>(threads + 1);
    bs = A.num_rows / threads;
    start[0] = 0;
    for (int i = 1; i < threads; i++) start[i] = start[i - 1] + bs;
    start[threads] = A.num_rows;

    {
        int id, i;
        lssp_mat_entry *row;
        int *p;
        int end;

        id = 0;

        row = lssp_malloc<lssp_mat_entry>(A.num_cols);
        p = lssp_malloc<int>(A.num_cols);
        end = start[id + 1];
        for (i = start[id]; i < end; i++) {
            int k, j, pt;
            bool sort_col = false;

            k = A.Ap[i + 1];
            if (k - A.Ap[i] > 1) {
                for (j = A.Ap[i] + 1; j < k; j++) {
                    if (A.Aj[j - 1] > A.Aj[j]) {
                        sort_col = true;
                        break;
                    }
                }
            }

            if (sort_col) {
                int idx;

                pt = 0;
                for (j = A.Ap[i]; j < k; j++) {
                    idx = A.Aj[j];

                    row[idx].j = A.Aj[j];
                    row[idx].x = A.Ax[j];
                    p[pt] = idx;
                    pt++;
                }

                qsort(p, pt, sizeof(int), lssp_mat_comp_row_by_col_index);

                pt = 0;
                for (j = A.Ap[i]; j < k; j++) {
                    idx = p[pt];

                    A.Aj[j] = row[idx].j;
                    A.Ax[j] = row[idx].x;
                    pt++;
                }
            }
        }

        lssp_free(row);
        lssp_free(p);
    }

    lssp_free(start);

    return;
}

lssp_mat_csr lssp_mat_adjust_zero_diag(const lssp_mat_csr A, double tol)
{
    lssp_mat_csr M = A;
    int n;
    int i, j, k;
    bool *has_diag;
    int row_end;
    double val;
    int *Aj;
    double *Ax;
    double t, t_a;

    n = A.num_rows;
    has_diag = lssp_malloc<bool>(n);

    for (i = 0; i < n; i++) {
        has_diag[i] = false;
    }

    k = 0;
    for (i = 0; i < n; i++) {
        row_end = A.Ap[i + 1];

        t = 0;
        for (j = A.Ap[i]; j < row_end; j++) {
            t_a = fabs(A.Ax[j]);

            if (t < t_a) t = t_a;
        }

        for (j = A.Ap[i]; j < row_end; j++) {
            if (A.Aj[j] == i) {
                has_diag[i] = true;
                k++;
            }
        }
    }

    M.Ap = lssp_malloc<int>(n + 1);
    M.Aj = lssp_malloc<int>(A.num_nnzs + n - k);
    M.Ax = lssp_malloc<double>(A.num_nnzs + n - k);

    for (i = 0; i < n; i++) {
        if (has_diag[i]) {
            M.Ap[i] = A.Ap[i + 1] - A.Ap[i];
        }
        else {
            M.Ap[i] = A.Ap[i + 1] - A.Ap[i] + 1;
        }
    }

    k = M.Ap[0];
    M.Ap[0] = 0;
    for (i = 1; i <= n; i++) {
        j = M.Ap[i];
        M.Ap[i] = M.Ap[i - 1] + k;
        k = j;
    }

    for (i = 0; i < n; i++) {
        Aj = M.Aj + M.Ap[i];
        Ax = M.Ax + M.Ap[i];

        row_end = A.Ap[i + 1];
        for (j = A.Ap[i]; j < row_end; j++) {
            *(Aj++) = A.Aj[j];
            *(Ax++) = A.Ax[j];
        }

        if (!has_diag[i]) {
            t = 0;
            row_end = A.Ap[i + 1];

            for (j = A.Ap[i]; j < row_end; j++) {
                t_a = fabs(A.Ax[j]);
                if (t < t_a) t = t_a;
            }

            *Aj = i;
            *Ax = 1 * tol;

            for (j = M.Ap[i + 1] - 2; j >= M.Ap[i]; j--) {
                int j1 = j + 1;

                if (M.Aj[j] > M.Aj[j1]) {
                    val = M.Ax[j1];
                    k = M.Aj[j1];

                    M.Aj[j1] = M.Aj[j];
                    M.Ax[j1] = M.Ax[j];

                    M.Aj[j] = k;
                    M.Ax[j] = val;
                }
                else {
                    break;
                }
            }
        }
    }

    lssp_free<bool>(has_diag);

    return M;
}

lssp_mat_csr lssp_mat_get_block_diag(const lssp_mat_csr A, int blk_size)
{
    int i, j, k;
    int *num_elem;
    int blk_num;
    int start, end;
    int num_nnzs;
    lssp_mat_csr M;
    int *Aj;
    double *Ax;

    assert(A.num_nnzs > 0);
    assert(A.num_rows > 0);
    assert(A.num_cols > 0);
    assert(A.num_rows == A.num_cols);
    assert(blk_size > 0);

    if (blk_size == A.num_rows) {
        M = A;
        M.Ap = lssp_copy_on<int>(A.Ap, A.num_rows + 1);
        M.Aj = lssp_copy_on<int>(A.Aj, A.num_nnzs);
        M.Ax = lssp_copy_on<double>(A.Ax, A.num_nnzs);

        return M;
    }

    blk_num = (A.num_rows + blk_size - 1) / blk_size;

    num_elem = lssp_malloc<int>(A.num_rows + 1);

    for (i = 0; i < A.num_rows; i++) {
        num_elem[i] = 0;
    }

    num_nnzs = 0;
    for (i = 0; i < blk_num; i++) {

        start = i * blk_size;
        end = (i + 1) * blk_size;
        if (end >= A.num_rows) end = A.num_rows;

        for (j = start; j < end; j++) {
            int row_end = A.Ap[j + 1];
            int col_id;

            for (k = A.Ap[j]; k < row_end; k++) {
                col_id = A.Aj[k];
                if (col_id >= start && col_id < end) num_elem[j]++;
            }

            if (num_elem[j] == 0) {
                num_elem[j] = -1;
                num_nnzs += 1;
            }
            else {
                num_nnzs += num_elem[j];
            }
        }
    }

    Aj = lssp_malloc<int>(num_nnzs);
    Ax = lssp_malloc<double>(num_nnzs);

    M.num_rows = A.num_rows;
    M.num_cols = A.num_cols;
    M.num_nnzs = num_nnzs;
    M.Ap = num_elem;
    M.Ax = Ax;
    M.Aj = Aj;

    num_nnzs = 0;
    for (i = 0; i < blk_num; i++) {
        start = i * blk_size;
        end = (i + 1) * blk_size;
        if (end >= A.num_rows) end = A.num_rows;

        for (j = start; j < end; j++) {
            int row_end = A.Ap[j + 1];
            int col_id;

            if (num_elem[j] > 0) {
                for (k = A.Ap[j]; k < row_end; k++) {
                    col_id = A.Aj[k];

                    if (col_id >= start && col_id < end) {
                        Aj[num_nnzs] = col_id;
                        Ax[num_nnzs] = A.Ax[k];
                        num_nnzs++;
                    }
                }
            }
            else {
                num_elem[j] = 1;
                Ax[num_nnzs] = 1;
                Aj[num_nnzs] = j;
                num_nnzs++;
            }
        }
    }

    k = num_elem[0];
    num_elem[0] = 0;
    for (i = 1; i <= A.num_rows; i++) {
        j = num_elem[i];
        num_elem[i] = num_elem[i - 1] + k;
        k = j;
    }

    return M;
}

static void lssp_mat_transpose_(const int *Ap, const int *Aj, const double * Ax,
        const int num_rows, const int num_cols, int *Bp, int *Bj, double * Bx)
{
    int *temp = lssp_malloc <int>(num_cols);
    int num_nnzs = Ap[num_rows];
    int i, jj;

    for (i = 0; i < num_cols; i++) {
        temp[i] = 0;
    }

    for (i = 0; i < num_nnzs; i++)
        temp[Aj[i]]++;

    Bp[0] = 0;
    for (i = 0; i < num_cols; i++) {
        Bp[i + 1] = Bp[i] + temp[i];
        temp[i] = 0;
    }

    for (i = 0; i < num_rows; i++) {
        int row_start = Ap[i];
        int row_end = Ap[i + 1];

        for (jj = row_start; jj < row_end; jj++) {
            int col;
            int offset;

            col = Aj[jj];
            offset = temp[col] + Bp[col];

            Bj[offset] = i;
            Bx[offset] = Ax[jj];

            temp[col] += 1;
        }
    }

    lssp_free<int>(temp);
}

lssp_mat_csr lssp_mat_transpose(const lssp_mat_csr A)
{
    lssp_mat_csr T;

    assert(A.num_rows > 0 && A.num_cols > 0);

    lssp_mat_init(T);
    T.num_rows = A.num_cols;
    T.num_cols = A.num_rows;
    T.num_nnzs = A.num_nnzs;

    if (A.num_nnzs <= 0) {
        return T;
    }

    assert(A.Ap != NULL && A.Aj != NULL && A.Ax != NULL);

    T.Ap = lssp_malloc <int>(T.num_rows + 1);
    T.Aj = lssp_malloc <int>(A.num_nnzs);
    T.Ax = lssp_malloc <double>(A.num_nnzs);

    lssp_mat_transpose_(A.Ap, A.Aj, A.Ax, A.num_rows, A.num_cols, T.Ap, T.Aj, T.Ax);

    return T;
}
