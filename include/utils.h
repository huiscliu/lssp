
#ifndef LSSP_UTILS_H
#define LSSP_UTILS_H

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>
#include <math.h>

#include "config.h"

extern int lssp_verbosity;

#define lssp_unused(x) (void)(x)

/* set log */
void lssp_set_log(FILE *io);

/* compare function, ascendant, descendant order */
int lssp_comp_int_asc(const void *p, const void *n);
int lssp_comp_int_des(const void *p, const void *n);

/* timer */
double lssp_get_time();

/* memory usage */
double lssp_get_mem_usage(double *peak);

/* output to stdout and default log file */
int lssp_printf(const char *fmt, ...);

/* output error */
void lssp_error(int code, const char *fmt, ...);

/* print warning info */
void lssp_warning(const char *fmt, ...);

template <typename T> T * lssp_malloc(const int n) 
{ 
    T *ptr;

    assert(n >= 0);

    if (n == 0) return NULL;

    ptr = (T *) malloc(n * sizeof(T));

    if (ptr == NULL) {
        lssp_error(1, "lssp: failed to malloc %g MB memory: %s %d.\n", \
                n * sizeof(T) / 1048576., __FILE__, __LINE__);
    }

    return ptr;
}

template <typename T> T * lssp_calloc(const int n) 
{ 
    T *ptr;

    assert(n >= 0);

    if (n == 0) return NULL;

    ptr = (T *) calloc(n, sizeof(T));

    if (ptr == NULL) {
        lssp_error(1, "lssp: failed to calloc %g MB memory: %s %d.\n", \
                n * sizeof(T) / 1048576., __FILE__, __LINE__);
    }

    return ptr;
}

template <typename T> T * lssp_realloc(T * old, const int n) 
{ 
    T *ptr;

    assert(n >= 0);

    if (n == 0) return NULL;

    ptr = (T *) realloc((void *)old, n * sizeof(T));
    if (ptr == NULL) {
        lssp_error(1, "lssp: failed to realloc %g MB memory: %s %d.\n", \
                n * sizeof(T) / 1048576., __FILE__, __LINE__);
    }

    return ptr;
}

template <typename T> void lssp_free(T * &p) 
{ 
    if (p == NULL) return;

    free(p);

    p = NULL;
}

template<typename T> void lssp_memcpy_on(T *dst, const T *src, const int n)
{
    assert(n >= 0);

    if (n == 0) return;

    assert(dst != NULL && src != NULL);
    memcpy(dst, src, n * sizeof(T));
}

template <typename T> T * lssp_copy_on(const T *src, const int n)
{
    T *dst;

    assert(n >= 0);

    if (n == 0) return NULL;

    assert(src != NULL);

    dst = lssp_malloc<T>(n);
    lssp_memcpy_on<T>(dst, src, n);

    return dst;
}

#endif
