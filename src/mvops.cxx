
#include "mvops.h"

/* y = y * beta + alpha * A * x  */
static void lssp_mv_amxpby__(const int num_rows, const int *Ap, const int *Aj, \
        const double *Ax, const double *x, double alpha, double *y,  \
        double beta)  
{
    int i;

    if (Ap != NULL) {
        for (i = 0; i < num_rows; i++){
            const int row_start = Ap[i];
            const int row_end   = Ap[i + 1];
            double sum = 0.;
            int jj;

            for (jj = row_start; jj < row_end; jj++) {            
                const int j = Aj[jj];
                sum += x[j] * Ax[jj];
            }

            y[i] = sum * alpha + y[i] * beta;
        }
    }
    else {
        for (i = 0; i < num_rows; i++){
            y[i] *= beta;
        }
    }
}

void lssp_mv_amxpby(double alpha, const lssp_mat_csr A, const lssp_vec x,  double beta, lssp_vec y)
{
    assert(x.n && y.n);
    assert(x.n == A.num_cols);

    lssp_mv_amxpby__(A.num_rows, A.Ap, A.Aj, A.Ax, x.d, alpha, y.d, beta);
}

/* z = y * beta + alpha * A * x  */
static void lssp_mv_amxpbyz__(const int num_rows, const int *Ap, const int *Aj, \
        const double *Ax, const double *x, double alpha, const double *y,       \
        double beta, double *z)
{
    int i;

    if (Ap != NULL) {
        for (i = 0; i < num_rows; i++){
            const int row_start = Ap[i];
            const int row_end   = Ap[i + 1];
            double sum = 0;
            int jj;

            for (jj = row_start; jj < row_end; jj++) {
                const int j = Aj[jj];

                sum += x[j] * Ax[jj];
            }

            z[i] = y[i] * beta + alpha * sum; 
        }
    }
    else {
        for (i = 0; i < num_rows; i++){
            z[i] = y[i] * beta;
        }
    }
}

void lssp_mv_amxpbyz(double alpha, const lssp_mat_csr A, const lssp_vec x,  double beta, \
        const lssp_vec y, lssp_vec z)
{
    assert(x.n == y.n && y.n == z.n);
    assert(x.n == A.num_cols);

    lssp_mv_amxpbyz__(A.num_rows, A.Ap, A.Aj, A.Ax, x.d, alpha, y.d, beta, z.d);
}

/* y = a * A * x  */
static void lssp_mv_amxy__(const int num_rows, const int *Ap, const int *Aj, \
        const double *Ax, const double *x, double a, double *y)  
{
    int i;

    if (Ap != NULL) {
        for (i = 0; i < num_rows; i++){
            const int row_start = Ap[i];
            const int row_end   = Ap[i + 1];
            double sum = 0;
            int jj;

            for (jj = row_start; jj < row_end; jj++) {            
                const int j = Aj[jj];

                sum += x[j] * Ax[jj];
            }

            y[i] = sum * a; 
        }
    }
    else {
        for (i = 0; i < num_rows; i++){
            y[i] = 0;
        }
    }
}

void lssp_mv_amxy(double a, const lssp_mat_csr A, const lssp_vec x, lssp_vec y)
{
    assert(x.n == y.n);
    assert(x.n == A.num_cols);

    lssp_mv_amxy__(A.num_rows, A.Ap, A.Aj, A.Ax, x.d, a, y.d);
}

/* y = A * x  */
static void lssp_mv_mxy__(const int num_rows, const int *Ap, const int *Aj, \
        const double *Ax, const double *x, double *y)  
{
    int i;

    if (Ap != NULL) {
        for (i = 0; i < num_rows; i++){
            const int row_start = Ap[i];
            const int row_end   = Ap[i + 1];
            int jj;
            double sum = 0;

            for (jj = row_start; jj < row_end; jj++) {
                sum += x[Aj[jj]] * Ax[jj];
            }

            y[i] = sum; 
        }
    }
    else {
        for (i = 0; i < num_rows; i++){
            y[i] = 0;
        }
    }
}

void lssp_mv_mxy(const lssp_mat_csr A, const lssp_vec x,  lssp_vec y)
{
    assert(x.n == y.n);
    assert(x.n == A.num_cols);

    lssp_mv_mxy__(A.num_rows, A.Ap, A.Aj, A.Ax, x.d, y.d);
}
