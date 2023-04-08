
#ifndef LSSP_MATRIX_UTILS_H
#define LSSP_MATRIX_UTILS_H

#include "type-defs.h"
#include "utils.h"

/* init matrix */
void lssp_mat_init(lssp_mat_csr &A);
void lssp_mat_init(lssp_mat_coo &A);
void lssp_mat_init(lssp_mat_bcsr &A);

/* free functions */
void lssp_mat_destroy(lssp_mat_csr &A);
void lssp_mat_destroy(lssp_mat_coo &A);
void lssp_mat_destroy(lssp_mat_bcsr &A);

/* create csr matrix */
lssp_mat_csr lssp_mat_create(int nrows, int ncols, int *Ap, int *Aj, double *Ax);

/* csr to bcsr, bs: block size */
lssp_mat_bcsr lssp_mat_csr_to_bcsr(const lssp_mat_csr A, int bs);

/* mat bcsr to csr */
lssp_mat_csr lssp_mat_bcsr_to_csr(const lssp_mat_bcsr A);  

/* csr to coo */
lssp_mat_coo lssp_mat_csr_to_coo(const lssp_mat_csr csr);

/* coo to csr */
lssp_mat_csr lssp_mat_coo_to_csr(const lssp_mat_coo A);

/* if csr mat is sorted */
bool lssp_mat_csr_is_sorted(const lssp_mat_csr A);

/* if block csr matrix is sorted */
bool lssp_mat_bcsr_is_sorted(const lssp_mat_bcsr A);

/* sort csr column index */
void lssp_mat_sort_column(lssp_mat_csr &A);

/* change zero diagonal element to some default value */
lssp_mat_csr lssp_mat_adjust_zero_diag(const lssp_mat_csr A, double tol);

/* generate block matrix */
lssp_mat_csr lssp_mat_get_block_diag(const lssp_mat_csr A, int blk_size);

/* transpose */
lssp_mat_csr lssp_mat_transpose(const lssp_mat_csr A);

#endif
