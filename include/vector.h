
#ifndef LSSP_VECTOR_H
#define LSSP_VECTOR_H

#include "matrix-utils.h"

/* create and destroy vector */
lssp_vec lssp_vec_create(int n);
void lssp_vec_destroy(lssp_vec &v);

/* set value */
void lssp_vec_set_value(lssp_vec x, double val);
void lssp_vec_set_value_by_array(lssp_vec x, double *val);
void lssp_vec_set_value_by_index(lssp_vec x, int i, double val);

/* get value from vector */
void lssp_vec_get_value(double *val, lssp_vec x);
double lssp_vec_get_value_by_index(lssp_vec x, int i);

/* copy */
void lssp_vec_copy(lssp_vec des, const lssp_vec src);

/* y = beta * y + alpha * x */
void lssp_vec_axpby(double alpha, const lssp_vec x, double beta, lssp_vec y);

/* z = beta * y + alpha * x */
void lssp_vec_axpbyz(double alpha, const lssp_vec x, double beta, lssp_vec y, lssp_vec z);

/* dot product */
double lssp_vec_dot(const lssp_vec x, const lssp_vec y);

/* L2 norm */
double lssp_vec_norm(const lssp_vec x);

/* x = a * x */
void lssp_vec_scale(lssp_vec x, double a);

#endif
