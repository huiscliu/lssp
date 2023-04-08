#include "pc-ilut.h"

extern double mat_zero_diag_value;
extern double mat_zero_diag_tol;
extern double lssp_pc_ilut_tol;

static int ilut_qsplit(double *a, int *ind, int n, int ncut)
{
    double tmp, abskey;
    int j, k, first, mid, last;

    first = 0;
    last = n - 1;

    if (ncut < first || ncut >= last) return 0;

loop:
    mid = first;
    abskey = fabs(a[mid]);

    for (j = first + 1; j <= last; j++) {
        if (fabs(a[j]) > abskey) {
            tmp = a[++mid];
            k = ind[mid];
            a[mid] = a[j];
            ind[mid] = ind[j];
            a[j] = tmp;
            ind[j] = k;
        }
    }

    tmp = a[mid];
    a[mid] = a[first];
    a[first] = tmp;
    k = ind[mid];
    ind[mid] = ind[first];
    ind[first] = k;

    if (mid == ncut) return 0;

    if (mid > ncut) {
        last = mid - 1;
    }
    else {
        first = mid + 1;
    }

    goto loop;
}

static void lssp_pc_ilut_fac(const lssp_mat_csr A, lssp_mat_csr &M, double ilut_tol, int p, int verb)
{

    int i,j, k;
    int kk, col, jrow, len, jwu_index, jpos;
    int n = A.num_rows;
    double t;
    double norm, mx;
    double ilut_rel_tol;

    int *C, *Ap;
    double *Ax = A.Ax;
    double *w;
    double *diag;
    int *jr;
    int *jw;
    int row_end;
    int offset = 0;
    int jwl_offset = 0;
    int jwu_offset = 0;
    int length = 8 * n;
    double a_ik;

    lssp_unused(verb);

    M = A;
    M.Ap = lssp_malloc<int>(n + 1);
    M.Aj = lssp_malloc<int>(length);
    M.Ax = lssp_malloc<double>(length);

    C = A.Aj;
    Ap = A.Ap;

    w = lssp_malloc<double>(n * 2);
    jr = lssp_malloc<int>(n * 2);
    diag = w + n;
    jw = jr + n;

    for(kk = Ap[0]; kk < Ap[1]; kk++) {
        M.Ax[kk] = Ax[kk];
        M.Aj[kk] = C[kk];
    }

    offset = Ap[1];
    M.Ap[0] = 0;
    M.Ap[1] = Ap[1];

    diag[0] = Ax[0];

    if (fabs(diag[0]) < mat_zero_diag_tol) {
        if (diag[0] > 0) {
            diag[0] = mat_zero_diag_value;
        }
        else {
            diag[0] = - mat_zero_diag_value;
        }
    }

    for (i = 0; i < n; i++) {
        jr[i] = -1;
    }

    for (i = 1; i < n; i++) {
        row_end = Ap[i + 1];

        norm = 0.0;
        for(k = Ap[i]; k < row_end; k++) {
            norm += fabs(Ax[k]);
        }
        norm /= (double)(row_end - Ap[i]);

        ilut_rel_tol = ilut_tol * norm;

        jw[i] = i;
        w[i] = 0.0;
        jr[i] = i;

        for(kk = Ap[i]; kk < row_end; kk++) {
            col = C[kk];
            if(col < i) {
                jr[col] = jwl_offset;
                jw[jwl_offset] = col;
                w[jwl_offset] = Ax[kk];
                jwl_offset++;
            }
            else if(col == i) {
                w[i] = Ax[kk];
            }
            else {
                jwu_offset++;
                jwu_index = i + jwu_offset;
                jr[col] = jwu_index;
                jw[jwu_index] = col;
                w[jwu_index] = Ax[kk];
            }
        }

        j = -1;
        len = 0;

        for (j = 0; j < jwl_offset; j++) {
            jrow = jw[j];
            jpos = j;

            for( k = j + 1; k < jwl_offset; k++) {
                if(jw[k] < jrow) {
                    jrow = jw[k];
                    jpos = k;
                }
            }

            if( jpos != j) {
                col = jw[j];
                jw[j] = jw[jpos];
                jw[jpos] = col;
                jr[jrow] = j;
                jr[col] = jpos;
                t = w[j];
                w[j] = w[jpos];
                w[jpos] = t;
            }

            jr[jrow] = -1;

            a_ik = w[j] = w[j] / diag[jrow];

            for(kk = M.Ap[jrow]; kk < M.Ap[jrow + 1]; kk++) {
                col = M.Aj[kk];

                if(col > jrow) {
                    jpos = jr[col];
                    mx = -a_ik * M.Ax[kk];

                    if(jpos == -1 && fabs(mx) < ilut_rel_tol) continue;

                    if(col < i){
                        if(jpos == -1) {
                            jw[jwl_offset] = col;
                            jr[col] = jwl_offset;
                            w[jwl_offset] = mx;
                            jwl_offset++;
                        }
                        else {
                            w[jpos] += mx;
                        }
                    }
                    else {
                        if( jpos == -1) {
                            int tmp;

                            jwu_offset++;
                            tmp = i + jwu_offset;
                            jw[tmp] = col;
                            jr[col] = tmp;
                            w[tmp] = mx;
                        }
                        else {
                            w[jpos] += mx;
                        }
                    }
                }
            }
        }

        if ((length - offset) < n ) {
            int *i_tmp;
            double *d_tmp;

            length = length + 2 * n;

            d_tmp = (double *)realloc((void *)M.Ax, length * sizeof(double));
            i_tmp = (int *)realloc((void *)M.Aj, length * sizeof(int));

            M.Ax = d_tmp;
            M.Aj = i_tmp;

            if (M.Ax == NULL || M.Aj == NULL ) {
                lssp_error(1, "pc: ilut, cannot realloc mem space: %s %d %s\n", \
                        __FILE__, __LINE__, __FUNCTION__);
            }
        }

        diag[i] = w[i];
        jr[i] = -1;

        for(j = 0; j < jwl_offset; j++) {
            jr[jw[j]] = -1;
        }

        for(j = 0; j < jwu_offset; j++) {
            jr[jw[i + j + 1]] = -1;
        }

        if (fabs(diag[i]) < mat_zero_diag_tol) {
            if (diag[i] > 0) {
                diag[i] = mat_zero_diag_value;
            }
            else {
                diag[i] = - mat_zero_diag_value;
            }
        }

        len = jwl_offset < p ? jwl_offset : p ;
        ilut_qsplit(w, jw, jwl_offset, len);

        for (kk = 0; kk < len; kk++) {
            M.Ax[offset] = w[kk];
            M.Aj[offset] = jw[kk];
            offset++;
        }

        M.Ax[offset] = diag[i];
        M.Aj[offset] = i;
        offset++;

        len = jwu_offset < p ? jwu_offset: p;
        ilut_qsplit(w + i + 1, jw + i + 1, jwu_offset, len);

        for(kk = 0; kk < len; kk++) {
            col= kk + i + 1;
            M.Ax[offset] = w[col];
            M.Aj[offset] = jw[col];
            offset++;
        }

        jwl_offset = 0;
        jwu_offset = 0;

        M.Ap[i + 1] = offset;
    }

    M.num_nnzs = offset;

    lssp_free(w);
    lssp_free(jr);
}

static void lssp_pc_ilut_assemble_matrix(const lssp_mat_csr A, int blk_size, lssp_mat_csr &L,
        lssp_mat_csr &U, double ilut_tol, int ilut_p, int verb)
{

    int i, j, k, m;
    int num_blk;
    int start, end;

    lssp_mat_csr M;
    lssp_mat_csr B;
    lssp_mat_csr T;

    int l_length, u_length;
    int l_offset, u_offset;
    int offset;

    assert(A.num_nnzs > 0);
    assert(A.num_rows > 0);
    assert(A.num_rows == A.num_cols);
    assert(blk_size > 0);
    assert(blk_size <= A.num_rows);

    num_blk = (A.num_rows + blk_size - 1) / blk_size;
    M = lssp_mat_get_block_diag(A, blk_size);

    L.num_rows = U.num_rows = M.num_rows;
    L.num_cols = U.num_cols = M.num_cols;

    l_length = 8 * L.num_rows;
    L.Ap = lssp_malloc<int>(L.num_rows + 1);
    L.Aj = lssp_malloc<int>(l_length);
    L.Ax = lssp_malloc<double>(l_length);

    u_length = 8 * U.num_rows;
    U.Ap = lssp_malloc<int>(U.num_rows + 1);
    U.Aj = lssp_malloc<int>(u_length);
    U.Ax = lssp_malloc<double>(u_length);

    L.Ap[0] = 0;
    U.Ap[0] = 0;

    l_offset = 0;
    u_offset = 0;

    for (i = 0; i < num_blk; i++) {
        start = i * blk_size;
        end = (i + 1) * blk_size;
        if (end >= A.num_rows) end = A.num_rows;

        B.num_cols = B.num_rows = end - start;
        B.num_nnzs = M.Ap[end] - M.Ap[start];
        B.Ap = M.Ap + start;
        B.Ax = M.Ax + M.Ap[start];
        B.Aj = M.Aj + M.Ap[start];

        for (j = 0; j < B.num_nnzs; j++) B.Aj[j] -= start;

        offset = M.Ap[start];
        for (j = start; j <= end; j++) {
            M.Ap[j] -= offset;
        }

        lssp_pc_ilut_fac(B, T, ilut_tol, ilut_p, verb);

        for (j = 0; j < T.num_nnzs; j++) T.Aj[j] += start;

        if ((l_length - l_offset) < T.num_nnzs) {
            l_length = l_length + T.num_nnzs + 2 * L.num_rows;

            L.Ax = (double *)realloc (L.Ax, l_length * sizeof(double));
            L.Aj = (int *)realloc (L.Aj, l_length * sizeof(int));

            if (L.Ax == NULL || L.Aj == NULL ) {
                lssp_error(1, "pc: ilut, cannot realloc mem space: %s %d %s\n", \
                        __FILE__, __LINE__, __FUNCTION__);
            }
        }

        if((u_length - u_offset) < T.num_nnzs) {
            u_length = u_length + T.num_nnzs + 2 * U.num_rows;

            U.Ax = (double *)realloc ((void *)U.Ax, (u_length) * sizeof(double));
            U.Aj = (int *)realloc ((void *)U.Aj, (u_length) * sizeof(int));

            if (U.Ax == NULL || U.Aj == NULL ) {
                lssp_error(1, "pc: ilut, cannot realloc mem space: %s %d %s\n", \
                        __FILE__, __LINE__, __FUNCTION__);
            }
        }

        for (j = 0; j < T.num_rows; j++) {
            int r_idx = start + j;

            m = T.Ap[j + 1];
            for (k = T.Ap[j]; k < m; k++) {
                if (T.Aj[k] < r_idx) {
                    L.Aj[l_offset] = T.Aj[k];
                    L.Ax[l_offset] = T.Ax[k];
                    l_offset++;
                }
                else if (T.Aj[k] == r_idx) {
                    L.Aj[l_offset] = T.Aj[k];
                    L.Ax[l_offset] = 1;
                    l_offset++;

                    U.Ax[u_offset] = T.Ax[k];
                    U.Aj[u_offset] = T.Aj[k];
                    u_offset++;
                }
                else {
                    U.Aj[u_offset] = T.Aj[k];
                    U.Ax[u_offset] = T.Ax[k];
                    u_offset++;
                }
            }
            L.Ap[r_idx+1] = l_offset;
            U.Ap[r_idx+1] = u_offset;
        }

        M.Ap[end] += offset;

        lssp_mat_destroy(T);
    }

    L.num_nnzs = L.Ap[L.num_rows];
    U.num_nnzs = U.Ap[U.num_rows];

    if (verb > 2) {
        lssp_printf("pc: ilut, L and U filled: %d, ratio: %f\n", L.num_nnzs + U.num_nnzs, \
                (L.num_nnzs + U.num_nnzs) * 1. / A.num_nnzs);
    }

    lssp_mat_destroy(M);
}

void lssp_pc_ilut_destroy(LSSP_PC *pc)
{
    assert(pc->assembled);
    lssp_pc_iluk_destroy(pc);
}

void lssp_pc_ilut_assemble(LSSP_PC &pc, LSSP_SOLVER s)
{
    lssp_mat_csr A;

    assert(s.A.num_rows == s.A.num_cols);
    assert(s.A.num_rows > 0 && s.A.num_nnzs > 0);

    if (pc.ilut_p <= 0) {
        pc.ilut_p = (s.A.num_nnzs + s.A.num_rows - 1) / s.A.num_rows;
    }

    if (pc.ilut_tol < 0) {
        pc.ilut_tol = lssp_pc_ilut_tol;
    }

    if (pc.verb > 1) {
        lssp_printf("pc: ilut, tol: %f, p: %d\n", pc.ilut_tol, pc.ilut_p);
    }

    A = lssp_mat_adjust_zero_diag(s.A, mat_zero_diag_tol);
    lssp_pc_ilut_assemble_matrix(A, A.num_rows, pc.L, pc.U, pc.ilut_tol, pc.ilut_p, pc.verb);
    lssp_mat_destroy(A);

    pc.cache = lssp_malloc<double>(pc.L.num_rows);
    pc.solve = lssp_pc_ilu_solve;

    pc.destroy = lssp_pc_ilut_destroy;
}

void lssp_pc_ilut_set_drop_tol(LSSP_PC &pc, double tol)
{
    pc.ilut_tol = fabs(tol);
}

void lssp_pc_ilut_set_p(LSSP_PC &pc, int p)
{
    pc.ilut_p = p;
}
