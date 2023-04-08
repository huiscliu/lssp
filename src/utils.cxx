
#include "utils.h"

#if HAVE_SYS_RESOURCE_H
#include <sys/resource.h>
#endif

#if TIME_WITH_SYS_TIME
#include <sys/time.h>
#include <time.h>
#else
#if HAVE_SYS_TIME_H
#include <sys/time.h>
#else
#include <time.h>
#endif
#endif

static FILE *lssp_log_handle_ = NULL;
int lssp_verbosity = 2;
static double mem_peak = 0;
static double mem_current = 0;

void lssp_set_log(FILE *io)
{
    lssp_log_handle_ = io;
}

int lssp_comp_int_asc(const void *p, const void *n)
{
    return *(int *)p - *(int *)n;
}

int lssp_comp_int_des(const void *p, const void *n)
{
    return *(int *)n - *(int *)p;
}


double lssp_get_time()
{
    struct timeval tv;

    gettimeofday(&tv, (struct timezone *)0);
    return tv.tv_sec + (double)tv.tv_usec * 1e-6;
}

#if 0 /* windows */
#include <windows.h>

double lssp_get_time()
{
    LARGE_INTEGER timer;
    static LARGE_INTEGER fre;
    static int init = 0;

    if (init != 1) {
        QueryPerformanceFrequency(&fre);
        init = 1;
    }

    QueryPerformanceCounter(&timer);

    return (double)1. * timer.QuadPart / fre.QuadPart;
}
#endif

#if HAVE_SYS_RESOURCE_H
double lssp_get_mem_usage(double *peak)
{
    struct rusage RU;

    /* getrusage() */
    getrusage(RUSAGE_SELF, &RU);
    mem_current = RU.ru_maxrss / (double)1024.;

    if (mem_current > mem_peak) mem_peak = mem_current;

    if (peak != NULL) *peak = mem_peak;

    return mem_current;
}
#else
double lssp_get_mem_usage(double *peak)
{
    mem_current = mem_peak = -1;
    if (peak != NULL) *peak = mem_peak;

    return mem_current;
}
#endif

int lssp_printf(const char *fmt, ...)
{
    va_list ap;
    int ret;

    va_start(ap, fmt);
    ret = vfprintf(stdout, fmt, ap);
    va_end(ap);

    if (lssp_log_handle_ != NULL) {
        va_start(ap, fmt);
        vfprintf(lssp_log_handle_, fmt, ap);
        va_end(ap);
    }

    fflush(stdout);
    fflush(lssp_log_handle_);

    return ret;
}

void lssp_error(int code, const char *fmt, ...)
{
    va_list ap;
    char s[2048];

    sprintf(s, "error: %s", fmt);
    va_start(ap, fmt);
    vfprintf(stdout, s, ap);
    va_end(ap);

    if (lssp_log_handle_ != NULL) {
        va_start(ap, fmt);
        vfprintf(lssp_log_handle_, fmt, ap);
        va_end(ap);
    }

    fflush(stdout);
    fflush(lssp_log_handle_);

    if (code == 0) return;
    exit(code);
}

void lssp_warning(const char *fmt, ...)
{
    va_list ap;
    char s[2028];

    sprintf(s, "warning: %s", fmt);
    va_start(ap, fmt);
    vfprintf(stdout, s, ap);
    va_end(ap);

    if (lssp_log_handle_ != NULL) {
        va_start(ap, fmt);
        vfprintf(lssp_log_handle_, fmt, ap);
        va_end(ap);
    }

    fflush(stdout);
    fflush(lssp_log_handle_);

    return;
}
