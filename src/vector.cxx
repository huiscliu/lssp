
#include "vector.h"

lssp_vec lssp_vec_create(int n)
{
    lssp_vec v;

    assert(n >= 0);

    if (n > 0) {
        v.n = n;
        v.d = lssp_malloc<double>(n);
    }
    else {
        v.n = 0;
        v.d = NULL;
    }

    return v;
}

void lssp_vec_destroy(lssp_vec &v)
{
    if (v.n > 0) lssp_free(v.d);

    v.n = 0;
    v.d = NULL;
}

/* set value */
void lssp_vec_set_value(lssp_vec x, double val)
{
    int i, n = x.n;

    for (i = 0; i < n; i++) {
        x.d[i] = val;
    }
}

/* set value */
void lssp_vec_set_value_by_array(lssp_vec x, double *val)
{
    int n = x.n;

    assert(val != NULL);
    assert(n >= 0);
    memcpy(x.d, val, sizeof(*val) * n);
}

void lssp_vec_set_value_by_index(lssp_vec x, int i, double val)
{
    assert(i >= 0 && i < x.n);

    x.d[i] = val;
}

void lssp_vec_get_value(double *val, lssp_vec x)
{
    int n = x.n;

    assert(val != NULL);
    assert(n >= 0);
    memcpy(val, x.d, sizeof(*val) * n);
}

double lssp_vec_get_value_by_index(lssp_vec x, int i)
{
    assert(i >= 0 && i < x.n);

    return x.d[i];
}

/* copy */
void lssp_vec_copy(lssp_vec x, const lssp_vec y)
{
    int i, n = x.n;

    assert(x.n == y.n);

    for (i = 0; i < n; i++) {
        x.d[i] = y.d[i];
    }
}

/* y = beta * y + alpha * x */
void lssp_vec_axpby(double alpha, const lssp_vec x, double beta, lssp_vec y)
{
    int i, n = x.n;

    assert(x.n == y.n);

    for (i = 0; i < n; i++) {
        y.d[i] = y.d[i] * beta + x.d[i] * alpha;
    }
}

/* z = beta * y + alpha * x */
void lssp_vec_axpbyz(double alpha, const lssp_vec x, double beta, lssp_vec y, lssp_vec z)
{
    int i, n = x.n;

    assert(x.n == y.n);
    assert(z.n == y.n);

    for (i = 0; i < n; i++) {
        z.d[i] = y.d[i] * beta + x.d[i] * alpha;
    }
}

/* dot product */
double lssp_vec_dot(const lssp_vec x, const lssp_vec y)
{
    double sum = 0;
    int i, n = x.n;

    assert(x.n == y.n);
    for (i = 0; i < n; i++) sum += x.d[i] * y.d[i];

    return sum;
}

/* L2 norm */
double lssp_vec_norm(const lssp_vec x)
{
    return sqrt(lssp_vec_dot(x, x));
}

/* x = a * x */
void lssp_vec_scale(lssp_vec x, double a)
{
    int  i, n = x.n;

    for (i = 0; i < n; i++) x.d[i] *= a;
}
