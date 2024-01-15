AC_PREREQ([2.71])

AC_INIT
m4_ifdef([AC_DISABLE_OPTION_CHECKING], [AC_DISABLE_OPTION_CHECKING()], [])

dnl This MUST precede any other macro
AC_CONFIG_AUX_DIR([utils])

AC_CONFIG_HEADERS([include/config.h])
AC_CONFIG_FILES([Makefile Makefile.inc])
LSSP_DIR="`pwd`"

AC_MSG_NOTICE(============================== compilers ===========================)

#--------------------------------------- config files and compilers
# Checks for programs.
AC_PROG_CC(gcc clang icc)
AC_PROG_CXX(g++ clang++ icpc)
AC_PROG_FC(gfortran gfortran-mp-5 g77 f77 ifort)
AC_PROG_CPP

if test "$GCC" = "yes"; then
    if $CC -V 2>&1 | grep -q Intel; then
        CFLAGS="$CFLAGS -w1"
    elif $CC -v 2>&1 | grep -q clang; then
        CFLAGS="$CFLAGS -Wall -Wextra -Warray-bounds"
    elif $CC -v 2>&1 | grep -q gcc; then
        CFLAGS="$CFLAGS -Wall -Wextra"

        for d in -Wmissing-declarations -Wmissing-prototypes -Wimplicit-function-declaration        \
            -Wunused -Wuninitialized -Wvla -Wformat -Wformat-security -Wmain -Wsequence-point       \
            -Wshadow -Wbad-function-cast -Wcomments -Wdeclaration-after-statement                   \
            -Waddress -Wpointer-sign -Wempty-body -Wnested-externs -Wold-style-declaration          \
            -Wold-style-definition -Wmissing-field-initializers -Wmaybe-uninitialized               \
            -Wno-return-local-addr -Wreturn-type -Wsequence-point -Wundef -Wunknown-pragmas         \
            -Wunused-label; do
            AC_MSG_CHECKING([whether $CC supports $d ])
            AC_LANG(C)
            CFLAGS_bak="$CFLAGS"
            CFLAGS="$CFLAGS $d"
            AC_COMPILE_IFELSE([AC_LANG_SOURCE([main(){}])], tmp=yes, tmp=no)
            if test $tmp = yes; then
                AC_MSG_RESULT(yes)
            else
                AC_MSG_RESULT(no)
                CFLAGS="$CFLAGS_bak"
            fi
        done
    fi

fi

if test "$GXX" = "yes"; then
    if $CXX -V 2>&1 | grep -q Intel; then
        CXXFLAGS="$CXXFLAGS -w1"
    elif $CXX -v 2>&1 | grep -q clang; then
        CXXFLAGS="$CXXFLAGS -Wno-deprecated -Wall -Wextra -Warray-bounds"
    elif $CXX -v 2>&1 | grep -q gcc; then
        CXXFLAGS="$CXXFLAGS -Wno-deprecated -Wall -Wextra"

        for d in -Wmissing-declarations -Wunused -Wuninitialized -Wvla -Wformat -Wformat-security      \
            -Wmain -Wsequence-point -Wshadow -Waddress -Wcomments -Wmain -Wsequence-point              \
            -Wempty-body -Wmissing-field-initializers  -Wmaybe-uninitialized -Wno-return-local-addr    \
            -Wreturn-type -Wsequence-point -Wundef -Wunknown-pragmas -Wunused-label; do
            AC_MSG_CHECKING([whether $CXX supports $d ])
            AC_LANG(C++)
            CXXFLAGS_bak="$CXXFLAGS"
            CXXFLAGS="$CXXFLAGS $d"
            AC_COMPILE_IFELSE([AC_LANG_SOURCE([main(){}])], tmp=yes, tmp=no)
            if test $tmp = yes; then
                AC_MSG_RESULT(yes)
            else
                AC_MSG_RESULT(no)
                CXXFLAGS="$CXXFLAGS_bak"
            fi
        done
    fi
fi

# Fortran
AC_FC_WRAPPERS
AC_FC_LIBRARY_LDFLAGS

LIBS="-L/usr/local/lib $LIBS -lm"

AC_ADD_LIBS(LIBS, ${FCLIBS}, append)

dnl Change relative path in compilers to full path
AC_FULLPATH(CC)
AC_FULLPATH(CXX)
AC_FULLPATH(CPP)
AC_FULLPATH(FC)

# The default linker
AC_SUBST(LINKER, "$LINKER")

AC_PREFIX_DEFAULT(/usr/local/lssp)

AC_CHECK_PROGS(AR,[ar], no)
if test $AR = "no" ; then
    AC_MSG_ERROR(* AR is nott found!);
fi

AC_SUBST(AR)
AC_PROG_RANLIB

#--------------------------------------- shared library

AC_SUBST(LIB_SUFFIX, "")
AC_SUBST(BUILD_SHARED, "")
AC_SUBST(BUILD_SHARED_LIBS, "")
AC_SUBST(BUILD_CFLAGS, "")
AC_SUBST(BUILD_CXXFLAGS, "")
dnl AC_SUBST(BUILD_FFLAGS, "")
AC_ARG_ENABLE(shared,
[  --enable-shared         build shared library liblssp.so (default)
  --disable-shared        build static library liblssp.a
  --with-build-shared     the BUILD_SHARED macro (default "gcc -shared" or "${CC} -shared")],,enable_shared="auto")

if test "$enable_shared" = "auto"; then
dnl    if test "$target_vendor" = "apple"; then
dnl	# use static link by default on MacOS, before the rpath issue is solved
dnl	enable_shared="no"
dnl    else
	enable_shared="yes"
dnl    fi
fi

if test "$enable_shared" = "yes"; then
    /bin/rm -f src/liblssp.a src/liblssp.dll
    LIB_SUFFIX=".so"

    if test -n "$with_build_shared"; then
        BUILD_SHARED="$with_build_shared"
    elif gcc -v >/dev/null 2>&1 && \
         test "$cross_compiling" != "yes" && \
         test "$target_vendor" != "apple" # MacOS: force using "${CC} -shared"
    then
        BUILD_SHARED="gcc -shared"
    else
        BUILD_SHARED="\${CC} -shared"
    fi

    # TODO: figure out a proper way to determine when to use the PIC flag
    if test "$target_cpu" != "i686" -a "$target_cpu" != "i586" -a \
            "$target_cpu" != "i486" -a "$target_cpu" != "i386"; then
        BUILD_CFLAGS="-fPIC"
        BUILD_CXXFLAGS="-fPIC"
        dnl BUILD_FFLAGS="-fPIC"
    fi

    dnl TODO: check that the above settings work!
    AC_MSG_NOTICE([*** building shared library liblssp$LIB_SUFFIX ("$BUILD_SHARED")])
else
    /bin/rm -f src/liblssp.so
    LIB_SUFFIX=".a"
    BUILD_SHARED=""
    BUILD_CFLAGS=""
    BUILD_CXXFLAGS=""
    dnl BUILD_FFLAGS=""
    AC_MSG_NOTICE([*** building static library liblssp$LIB_SUFFIX])
fi

#--------------------------------------- OpenMP
#AC_LANG(C)
#AC_SUBST(USE_OPENMP, "0")
#AC_SUBST(OMP_FLAGS, "")
#AC_ARG_ENABLE(omp,
#[  --enable-omp           enable OpenMP
#  --disable-omp          disable OpenMP (default)
#  --with-omp-flags=FLAGS compiler flags for OpenMP],,enable_omp=yes)
#AC_MSG_NOTICE([cheking for OpenMP ])
#if test "$enable_omp" != no; then
#    AC_OPENMP
#    if test -z "$with_omp_flags"; then
#        with_omp_flags=$OPENMP_CFLAGS
#    fi
#
#    CFLAGS_bak="$CFLAGS"
#    CXXFLAGS_bak="$CXXFLAGS"
#    LDFLAGS_bak="$LDFLAGS"
#    if test -n "$with_omp_flags"; then
#        AC_ADD_FLAGS(CFLAGS, $with_omp_flags)
#        AC_ADD_FLAGS(CXXFLAGS, $with_omp_flags)
#        AC_ADD_FLAGS(LDFLAGS, $with_omp_flags)
#    fi
#    
#    AC_LINK_IFELSE([AC_LANG_SOURCE([#include <omp.h>
#            main() {omp_get_max_threads;}])],,enable_omp=no)
#    AC_MSG_CHECKING([for OpenMP ])
#    AC_MSG_RESULT($enable_omp)
#    
#    if test "$enable_omp" != no; then
#        AC_DEFINE(USE_OPENMP, 1)
#    else
#        CFLAGS="$CFLAGS_bak"
#        CXXFLAGS="$CXXFLAGS_bak"
#        LDFLAGS="$LDFLAGS_bak"
#    fi
#else
#    AC_MSG_RESULT($enable_omp)
#fi

#--------------------------------------- rpath
AC_SUBST(RPATH_FLAG, "")
AC_ARG_ENABLE(rpath,
[  --enable-rpath          enable use of rpath (default)
  --disable-rpath         disable use of rpath
  --with-rpath-flag=FLAG  compiler flag for rpath (e.g., "-Wl,-rpath,")],,
enable_rpath=yes)
if test -z "$with_rpath_flag" -a "$enable_rpath" != "no"; then
    AC_MSG_CHECKING([for compiler flag for specifying rpath ])
    LIBS_bak="$LIBS"
    for with_rpath_flag in "-Wl,-rpath," "-Wl,-rpath " "-Wl,--rpath " \
			   "-rpath " "--rpath " "-Wl,-R" "-R" "-R " ; do
	LIBS="$LIBS_bak ${with_rpath_flag}/usr/lib"
	AC_LINK_IFELSE([AC_LANG_PROGRAM()], enable_rpath=yes, enable_rpath=no)
	if test "$enable_rpath" = "yes"; then break; fi
    done
    LIBS="$LIBS_bak"
    if test "$enable_rpath" = "yes"; then
	AC_MSG_RESULT("$with_rpath_flag")
	RPATH_FLAG="$with_rpath_flag"
    else
	AC_MSG_RESULT(unknown)
	RPATH_FLAG="unknown"
    fi
fi

#--------------------------------------- header files
AC_CHECK_HEADERS_ONCE([sys/time.h])
# Obsolete code to be removed.
if test $ac_cv_header_sys_time_h = yes; then
  AC_DEFINE([TIME_WITH_SYS_TIME],[1],[Define to 1 if you can safely include both <sys/time.h>
	     and <time.h>.  This macro is obsolete.])
fi
# End of obsolete code.

AC_CHECK_HEADERS([limits.h ctype.h])
AC_CHECK_HEADERS([sys/time.h sys/resource.h unistd.h sys/resource.h])

#--------------------------------------- library functions.
AC_CHECK_FUNCS([gettimeofday getrusage])

# fpe needs -lm
LIBS="$LIBS -lm"

AC_ADD_FLAGS(CPPFLAGS, -I/usr/local/include)
AC_ADD_FLAGS(LDFLAGS, -L/usr/local/lib)

#--------------------------------------- basic flags
# Note: -I`pwd`/include and -L`pwd`/src must be the first term in CPPFLAGS
# and LDFLAGS respectively, which will be changed to LSSP_INSTDIR/include by
# 'make install'
CPPFLAGS="-I${LSSP_DIR}/include $CPPFLAGS"
LDFLAGS="-L${LSSP_DIR}/src $LDFLAGS"

#--------------------------------------- Fortran
AC_FC_FUNC(dgemm,DGEMM)            dnl BLAS
AC_FC_FUNC(dgels,DGELS)            dnl LAPACK

#--------------------------------------- BLAS
AC_DEFINE(USE_BLAS, 0)
AC_LANG(C)
AC_ARG_ENABLE(blas,
[  --enable-blas         enable BLAS support (default)
  --disable-blas        disable BLAS support
  --with-blas=blas BLAS lib],,enable_blas=yes)
if test "$enable_blas" != no; then
    AC_MSG_CHECKING([whether we have BLAS ])

    if test "$with_blas" == ""; then
        LIBS_bak=$LIBS
        AC_ADD_LIBS(LIBS, -lblas)
        AC_CHECK_FORTRAN([AC_LANG_SOURCE([main(){$DGEMM();}])], blas)
        if test "$blas_ok" != "yes"; then
            enable_blas=no
            LIBS=$LIBS_bak
        fi
        with_blas=-lblas
    else
        LIBS_bak=$LIBS
        AC_ADD_LIBS(LIBS, $with_blas)
        AC_CHECK_FORTRAN([AC_LANG_SOURCE([main(){$DGEMM();}])], blas)
        if test "$blas_ok" != "yes"; then
            enable_blas=no
            LIBS=$LIBS_bak
        fi
    fi

    AC_MSG_RESULT($enable_blas)
    if test "$enable_blas" != no; then
        AC_MSG_NOTICE(* BLAS support enabled)
        AC_MSG_NOTICE(* BLAS libs = $with_blas)
        AC_DEFINE(USE_BLAS, 1)
    else
        AC_MSG_NOTICE(* BLAS is disabled)
        AC_DEFINE(USE_BLAS, 0)
    fi
else
    AC_MSG_NOTICE(* BLAS is disabled)
    AC_DEFINE(USE_BLAS, 0)
fi

#--------------------------------------- LAPACK
AC_DEFINE(USE_LAPACK, 0)
AC_LANG(C)
AC_ARG_ENABLE(lapack,
[  --enable-lapack         enable LAPACK support (default)
  --disable-lapack        disable LAPACK support
  --with-lapack=lapack LAPACK lib],,enable_lapack=yes)
if test "$enable_lapack" != no; then
    AC_MSG_CHECKING([whether we have LAPACK ])

    if test "$with_lapack" == ""; then
        LIBS_bak=$LIBS
        AC_ADD_LIBS(LIBS, -llapack)
        AC_CHECK_FORTRAN([AC_LANG_SOURCE([main(){$DGELS();}])], lapack)
        if test "$lapack_ok" != "yes"; then
            enable_lapack=no
            LIBS=$LIBS_bak
        fi
        with_lapack=-llapack
    else
        LIBS_bak=$LIBS
        AC_ADD_LIBS(LIBS, $with_lapack)
        AC_CHECK_FORTRAN([AC_LANG_SOURCE([main(){$DGELS();}])], lapack)
        if test "$lapack_ok" != "yes"; then
            enable_lapack=no
            LIBS=$LIBS_bak
        fi
    fi

    AC_MSG_RESULT($enable_lapack)
    if test "$enable_lapack" != no; then
        AC_MSG_NOTICE(* LAPACK support enabled)
        AC_MSG_NOTICE(* LAPACK libs = $with_lapack)
        AC_DEFINE(USE_LAPACK, 1)
    else
        AC_MSG_NOTICE(* LAPACK is disabled)
        AC_DEFINE(USE_LAPACK, 0)
    fi
else
    AC_MSG_NOTICE(* LAPACK is disabled)
    AC_DEFINE(USE_LAPACK, 0)
fi

#--------------------------------------- LASPACK
AC_LANG(C)
AC_ARG_ENABLE(laspack,
[  --enable-laspack      enable LASPACK support (default)
  --disable-laspack     disable LASPACK support
  --with-laspack-libdir=DIR path for LASPACK library
  --with-laspack-incdir=DIR path for LASPACK header file],,enable_laspack=yes)
if test "$enable_laspack" != no; then
    AC_MSG_CHECKING([whether we have LASPACK ])

    if test -z "$with_laspack_incdir" -a "$enable_laspack" != no; then
        for d in /usr/local/laspack/include /usr/local/laspack*/include    \
            /usr/local/include /usr/include /opt/laspack/include /opt/laspack*/include; do
            if test -r ${d}/laspack/lastypes.h; then
                with_laspack_incdir="$d"
                break
            fi
        done
        if test -z "$with_laspack_incdir"; then
            enable_laspack=no
        fi
    fi

    if test -z "$with_laspack_libdir" -a "$enable_laspack" != no; then
        for d in /usr/local/laspack/lib /usr/local/laspack*/lib  /usr/local/lib \
            /usr/local/lib64 /usr/lib /usr/lib64 /opt/laspack/lib /opt/laspack*/lib; do
            if test -r ${d}/liblaspack.*; then
                with_laspack_libdir="$d"
                break
            fi
        done
        if test -z "$with_laspack_libdir"; then
            enable_laspack=no
        fi
    fi

    AC_MSG_RESULT($enable_laspack)
    if test "$enable_laspack" != no; then
        AC_MSG_NOTICE(* LASPACK support enabled)
        AC_MSG_NOTICE(* LASPACK include dir = $with_laspack_incdir)
        AC_MSG_NOTICE(* LASPACK lib dir = $with_laspack_libdir)
        AC_MSG_NOTICE(* LASPACK libs = -llaspack)
        AC_DEFINE(USE_LASPACK, 1)
        AC_ADD_FLAGS(CPPFLAGS, -I${with_laspack_incdir})
        AC_ADD_FLAGS(LDFLAGS, -L${with_laspack_libdir})
        AC_ADD_LIBS(LIBS, -llaspack)
    else
        AC_MSG_NOTICE(* LASPACK is disabled)
        AC_DEFINE(USE_LASPACK, 0)
    fi
else
    AC_DEFINE(USE_LASPACK, 0)
    AC_MSG_NOTICE(* LASPACK is disabled)
fi

#--------------------------------------- SuiteSparse
AC_LANG(C)
AC_ARG_ENABLE(ssparse,
[  --enable-ssparse      enable SSPARSE support (default)
  --disable-ssparse     disable SSPARSE support
  --with-ssparse-libdir=DIR path for SSPARSE library
  --with-ssparse-incdir=DIR path for SSPARSE header file],,enable_ssparse=yes)
if test "$enable_ssparse" != no; then
    AC_MSG_CHECKING([whether we have SSPARSE ])

    if test -z "$with_ssparse_incdir" -a "$enable_ssparse" != no; then
        for d in /usr/local/SuiteSparse /usr/local/SuiteSparse* /usr/local /usr \
            /opt/SuiteSparse /opt/SuiteSparse*; do
            if test -r ${d}/include/umfpack.h; then
                with_ssparse_incdir="$d"
                break
            fi
        done
        if test -z "$with_ssparse_incdir"; then
            enable_ssparse=no
        fi
    fi

    if test -z "$with_ssparse_libdir" -a "$enable_ssparse" != no; then
        for d in /usr/local/SuiteSparse /usr/local/SuiteSparse* /usr/local /usr \
            /opt/SuiteSparse /opt/SuiteSparse* ; do
            if test -r ${d}/lib/libumfpack.so -o -r ${d}/lib/libumfpack.a; then
                with_ssparse_libdir="$d"
                break
            fi
        done
        if test -z "$with_ssparse_libdir"; then
            enable_ssparse=no
        fi
    fi

    AC_MSG_RESULT($enable_ssparse)
    if test "$enable_ssparse" != no; then
        AC_MSG_NOTICE(* SSPARSE support enabled)
        AC_MSG_NOTICE(* SSPARSE include dir = ${with_ssparse_incdir}/include)
        AC_MSG_NOTICE(* SSPARSE lib dir = ${with_ssparse_libdir}/lib)
        AC_MSG_NOTICE(* SSPARSE libs = -lumfpack -lamd -lklu)
        AC_DEFINE(USE_SSPARSE, 1)
        AC_ADD_FLAGS(CPPFLAGS, -I${with_ssparse_incdir}/include)
        AC_ADD_FLAGS(LDFLAGS, -L${with_ssparse_libdir}/lib)
        AC_ADD_LIBS(LIBS, -lumfpack -lamd -lklu)
    else
        AC_MSG_NOTICE(* SSPARSE is disabled)
        AC_DEFINE(USE_SSPARSE, 0)
    fi
else
    AC_DEFINE(USE_SSPARSE, 0)
    AC_MSG_NOTICE(* SSPARSE is disabled)
fi

#--------------------------------------- PETSC
AC_LANG(C)
AC_ARG_ENABLE(petsc,
[  --enable-petsc        enable PETSC solver (default)
  --disable-petsc       disable PETSC solver
  --with-petsc-incdir=DIR PETSC header files directory
  --with-petsc-libdir=DIR PETSC libraries directory],,enable_petsc=yes)

if test "$enable_petsc" != no; then
    AC_MSG_CHECKING([whether we have PETSC ])

    if test -z "$with_petsc_incdir" -a "$enable_petsc" != no; then
        for d in /usr/local/petsc-seq/include /usr/local/petsc*-seq/include \
            /opt/petsc-seq/include /opt/petsc*-seq/include; do
            if test -r ${d}/petscksp.h; then
                with_petsc_incdir="$d"
                break
            fi
        done

        if test -z "$with_petsc_incdir"; then
            enable_petsc=no
        fi
    fi

    if test -z "$with_petsc_libdir" -a "$enable_petsc" != no; then
        for d in /usr/local/petsc-seq/lib /usr/local/petsc*-seq/lib /opt/petsc-seq/lib \
            /opt/petsc*-seq/lib; do
            if test -r ${d}/libpetsc.a -o -r ${d}/libpetsc.so; then
                with_petsc_libdir="$d"
                break
            fi
        done

        if test -z "$with_petsc_libdir"; then
            enable_petsc=no
        fi
    fi

    AC_MSG_RESULT($enable_petsc)
fi

if test "$enable_petsc" != no; then
    AC_MSG_CHECKING([PETSC version ])
    tmpm=`grep "#define PETSC_VERSION_MAJOR" ${with_petsc_incdir}/petscversion.h \
        | gawk '{print $3}'`
    tmpi=`grep "#define PETSC_VERSION_MINOR" ${with_petsc_incdir}/petscversion.h \
        | gawk '{print $3}'`
    tmps=`grep "#define PETSC_VERSION_SUBMINOR" ${with_petsc_incdir}/petscversion.h \
        | gawk '{print $3}'`

    if test -n "$tmpm"; then
        AC_MSG_RESULT($tmpm.$tmpi.$tmps)
        petsc_major=$tmpm
        petsc_minor=$tmpi
    else
        AC_MSG_RESULT([unknown])
        petsc_major="0"
        petsc_minor="0"
    fi

    enable_petsc=yes;
    if (("$petsc_major" < "3")); then
        enable_petsc=no;
        AC_MSG_NOTICE([* PETSC is out of date, use 3.4.x or newer])
    elif (("$petsc_major" == "3")); then
        if (("$petsc_minor" < "4")); then
            enable_petsc=no;
            AC_MSG_NOTICE([* PETSC is out of date, use 3.4.x or newer])
        fi
    fi
fi

if test "$enable_petsc" != no; then
    if test "$enable_petsc" != no; then
        AC_MSG_NOTICE(* PETSC support enabled)
        AC_MSG_NOTICE(* PETSC include dir = $with_petsc_incdir)
        AC_MSG_NOTICE(* PETSC lib dir = $with_petsc_libdir)
        AC_MSG_NOTICE(* PETSC libs = -lpetsc -lX11)
        AC_ADD_FLAGS(CPPFLAGS, -I${with_petsc_incdir})
        AC_ADD_FLAGS(CPPFLAGS, -I${with_petsc_incdir}/petsc/mpiuni)
        AC_ADD_FLAGS(LDFLAGS, -L${with_petsc_libdir})
        AC_ADD_LIBS(LIBS, -lpetsc -lX11 -ldl)

        AC_DEFINE(USE_PETSC, 1)
        AC_DEFINE_UNQUOTED(PETSC_VERSION_MAJOR, $petsc_major)
        AC_DEFINE_UNQUOTED(PETSC_VERSION_MINOR, $petsc_minor)

        if test "$enable_rpath" = "yes"; then
            AC_ADD_LIBS(LIBS, "${with_rpath_flag}${with_petsc_libdir}")
        fi

    else
        AC_MSG_NOTICE(* PETSC support disabled)
        AC_DEFINE(USE_PETSC, 0)
    fi
else
    AC_MSG_NOTICE(* PETSC support disabled)
    AC_DEFINE(USE_PETSC, 0)
fi

#--------------------------------------- MUMPS
AC_LANG(C)
AC_ARG_ENABLE(mumps,
[  --enable-mumps        enable MUMPS solver (default)
  --disable-mumps       disable MUMPS solver
  --with-mumps-incdir=DIR MUMPS header files directory
  --with-mumps-libdir=DIR MUMPS libraries directory],,enable_mumps=yes)
if test "$enable_mumps" != no; then
    AC_MSG_CHECKING([whether we have MUMPS ])

    if test -z "$with_mumps_incdir" -a "$enable_mumps" != no; then
        for d in /usr/local/mumps-seq/include /usr/local/mumps*-seq/include /usr/local/include \
            /usr/local/include/mumps-seq/ /usr/include /usr/include/mumps-seq                  \
            /opt/mumps-seq/include /opt/mumps*-seq/include; do
            if test -r ${d}/dmumps_c.h; then
                with_mumps_incdir="$d"
                break
            fi
        done

        if test -z "$with_mumps_incdir"; then
            enable_mumps=no
        fi
    fi

    if test -z "$with_mumps_libdir" -a "$enable_mumps" != no; then
        for d in /usr/local/mumps-seq/lib /usr/local/mumps*-seq/lib /usr/local/lib    \
            /usr/local/lib64 /usr/lib /usr/lib64 /opt/mumps-seq/lib /opt/mumps*-seq/lib; do
            if test -r ${d}/libdmumps.a -o -r ${d}/libdmumps.so; then
                with_mumps_libdir="$d"
                break
            fi
        done

        if test -z "$with_mumps_libdir"; then
            enable_mumps=no
        fi
    fi

    AC_MSG_RESULT($enable_mumps)
fi

if test "$enable_mumps" != no; then
    if test "$enable_mumps" != no; then
        AC_MSG_NOTICE(* MUMPS support enabled)
        AC_MSG_NOTICE(* MUMPS include dir = $with_mumps_incdir)
        AC_MSG_NOTICE(* MUMPS lib dir = $with_mumps_libdir)
        AC_MSG_NOTICE(* MUMPS libs = -ldmumps -lmumps_common -lpord -lmpiseq -lpthread)

        AC_ADD_FLAGS(CPPFLAGS, -I${with_mumps_incdir})
        AC_ADD_FLAGS(LDFLAGS, -L${with_mumps_libdir})
        AC_ADD_LIBS(LIBS, "-ldmumps -lmumps_common -lpord -lmpiseq -lpthread")

        AC_DEFINE(USE_MUMPS, 1)
    else
        AC_MSG_NOTICE(* MUMPS support disabled)
        AC_DEFINE(USE_MUMPS, 0)
    fi
else
    AC_MSG_NOTICE(* MUMPS support disabled)
    AC_DEFINE(USE_MUMPS, 0)
fi

#--------------------------------------- ITSOL
AC_LANG(C)
AC_ARG_ENABLE(itsol,
[  --enable-itsol      enable ITSOL support (default)
  --disable-itsol     disable ITSOL support
  --with-itsol-libdir=DIR path for ITSOL library
  --with-itsol-incdir=DIR path for ITSOL header file],,enable_itsol=yes)
if test "$enable_blas" = no; then enable_itsol=no; fi
if test "$enable_itsol" != no; then
    AC_MSG_CHECKING([whether we have ITSOL ])

    if test -z "$with_itsol_incdir" -a "$enable_itsol" != no; then
        for d in /usr/local/itsol /usr/local/itsol* /usr/local /usr /opt/itsol /opt/itsol*; do
            if test -r ${d}/include/globheads.h; then
                with_itsol_incdir="$d"
                break
            fi
        done
        if test -z "$with_itsol_incdir"; then
            enable_itsol=no
        fi
    fi

    if test -z "$with_itsol_libdir" -a "$enable_itsol" != no; then
        for d in /usr/local/itsol /usr/local/itsol* /usr/local /usr /opt/itsol /opt/itsol*; do
            if test -r ${d}/lib/libitsol.so -o -r ${d}/lib/libitsol.a; then
                with_itsol_libdir="$d"
                break
            fi
        done
        if test -z "$with_itsol_libdir"; then
            enable_itsol=no
        fi
    fi

    AC_MSG_RESULT($enable_itsol)
    if test "$enable_itsol" != no; then
        AC_MSG_NOTICE(* ITSOL support enabled)
        AC_MSG_NOTICE(* ITSOL include dir = ${with_itsol_incdir}/include)
        AC_MSG_NOTICE(* ITSOL lib dir = ${with_itsol_libdir}/lib)
        AC_MSG_NOTICE(* ITSOL libs = -litsol)
        AC_DEFINE(USE_ITSOL, 1)
        AC_ADD_FLAGS(CPPFLAGS, -I${with_itsol_incdir}/include)
        AC_ADD_FLAGS(LDFLAGS, -L${with_itsol_libdir}/lib)
        AC_ADD_LIBS(LIBS, -litsol)
    else
        AC_MSG_NOTICE(* ITSOL is disabled)
        AC_DEFINE(USE_ITSOL, 0)
    fi
else
    AC_DEFINE(USE_ITSOL, 0)
    AC_MSG_NOTICE(* ITSOL is disabled)
fi

#--------------------------------------- LIS
AC_LANG(C)
AC_ARG_ENABLE(lis,
[  --enable-lis      enable LIS support (default)
  --disable-lis     disable LIS support
  --with-lis-libdir=DIR path for LIS library
  --with-lis-incdir=DIR path for LIS header file],,enable_lis=yes)
if test "$enable_lis" != no; then
    AC_MSG_CHECKING([whether we have LIS ])

    if test -z "$with_lis_incdir" -a "$enable_lis" != no; then
        for d in /usr/local/lis /usr/local/lis* /usr/local /usr /opt/lis /opt/lis*; do
            if test -r ${d}/include/lis.h; then
                with_lis_incdir="$d"
                break
            fi
        done
        if test -z "$with_lis_incdir"; then
            enable_lis=no
        fi
    fi

    if test -z "$with_lis_libdir" -a "$enable_lis" != no; then
        for d in /usr/local/lis /usr/local/lis* /usr/local /usr /opt/lis /opt/lis*; do
            if test -r ${d}/lib/liblis.so -o -r ${d}/lib/liblis.a; then
                with_lis_libdir="$d"
                break
            fi
        done
        if test -z "$with_lis_libdir"; then
            enable_lis=no
        fi
    fi

    AC_MSG_RESULT($enable_lis)
    if test "$enable_lis" != no; then
        AC_MSG_NOTICE(* LIS support enabled)
        AC_MSG_NOTICE(* LIS include dir = ${with_lis_incdir}/include)
        AC_MSG_NOTICE(* LIS lib dir = ${with_lis_libdir}/lib)
        AC_MSG_NOTICE(* LIS libs = -llis)
        AC_DEFINE(USE_LIS, 1)
        AC_ADD_FLAGS(CPPFLAGS, -I${with_lis_incdir}/include)
        AC_ADD_FLAGS(LDFLAGS, -L${with_lis_libdir}/lib)
        AC_ADD_LIBS(LIBS, -llis)
    else
        AC_MSG_NOTICE(* LIS is disabled)
        AC_DEFINE(USE_LIS, 0)
    fi
else
    AC_DEFINE(USE_LIS, 0)
    AC_MSG_NOTICE(* LIS is disabled)
fi

#--------------------------------------- QR_MUMPS
AC_LANG(C)
AC_ARG_ENABLE(qrmumps,
[  --enable-qrmumps      enable QR_MUMPS support (default)
  --disable-qrmumps     disable QR_MUMPS support
  --with-qrmumps-libdir=DIR path for QR_MUMPS library
  --with-qrmumps-incdir=DIR path for QR_MUMPS header file],,enable_qrmumps=yes)
if test "$enable_qrmumps" != no; then
    AC_MSG_CHECKING([whether we have QR_MUMPS ])

    if test -z "$with_qrmumps_incdir" -a "$enable_qrmumps" != no; then
        for d in /usr/local/qr_mumps /usr/local/qr_mumps* /usr/local /usr \
            /opt/qr_mumps /opt/qr_mumps*; do
            if test -r ${d}/include/qrm_mumps.h; then
                with_qrmumps_incdir="$d"
                break
            fi
        done
        if test -z "$with_qrmumps_incdir"; then
            enable_qrmumps=no
        fi
    fi

    if test -z "$with_qrmumps_libdir" -a "$enable_qrmumps" != no; then
        for d in /usr/local/qr_mumps /usr/local/qr_mumps* /usr/local /usr \
            /opt/qr_mumps /opt/qr_mumps*; do
            if test -r ${d}/lib/libdqrm.so -o -r ${d}/lib/libdqrm.a; then
                with_qrmumps_libdir="$d"
                break
            fi
        done
        if test -z "$with_qrmumps_libdir"; then
            enable_qrmumps=no
        fi
    fi

    AC_MSG_RESULT($enable_qrmumps)
    if test "$enable_qrmumps" != no; then
        AC_MSG_NOTICE(* QR_MUMPS support enabled)
        AC_MSG_NOTICE(* QR_MUMPS include dir = ${with_qrmumps_incdir}/include)
        AC_MSG_NOTICE(* QR_MUMPS lib dir = ${with_qrmumps_libdir}/lib)
        AC_MSG_NOTICE(* QR_MUMPS libs = -ldqrm -lqrm_common)
        AC_DEFINE(USE_QR_MUMPS, 1)
        AC_ADD_FLAGS(CPPFLAGS, -I${with_qrmumps_incdir}/include)
        AC_ADD_FLAGS(LDFLAGS, -L${with_qrmumps_libdir}/lib)
        AC_ADD_LIBS(LIBS, -ldqrm -lqrm_common)
    else
        AC_MSG_NOTICE(* QR_MUMPS is disabled)
        AC_DEFINE(USE_QR_MUMPS, 0)
    fi
else
    AC_DEFINE(USE_QR_MUMPS, 0)
    AC_MSG_NOTICE(* QR_MUMPS is disabled)
fi

#--------------------------------------- FASP
AC_LANG(C)
AC_ARG_ENABLE(fasp,
[  --enable-fasp      enable FASP support (default)
  --disable-fasp     disable FASP support
  --with-fasp-libdir=DIR path for FASP library
  --with-fasp-incdir=DIR path for FASP header file],,enable_fasp=yes)
if test "$enable_fasp" != no; then
    AC_MSG_CHECKING([whether we have FASP ])

    if test -z "$with_fasp_incdir" -a "$enable_fasp" != no; then
        for d in /usr/local/fasp /usr/local/fasp* /usr/local /usr /opt/fasp /opt/fasp*; do
            if test -r ${d}/include/fasp.h; then
                with_fasp_incdir="$d"
                break
            fi
        done
        if test -z "$with_fasp_incdir"; then
            enable_fasp=no
        fi
    fi

    if test -z "$with_fasp_libdir" -a "$enable_fasp" != no; then
        for d in /usr/local/fasp /usr/local/fasp* /usr/local /usr /opt/fasp /opt/fasp*; do
            if test -r ${d}/lib/libfasp.so -o -r ${d}/lib/libfasp.a; then
                with_fasp_libdir="$d"
                break
            fi
        done
        if test -z "$with_fasp_libdir"; then
            enable_fasp=no
        fi
    fi

    AC_MSG_RESULT($enable_fasp)
    if test "$enable_fasp" != no; then
        AC_MSG_NOTICE(* FASP support enabled)
        AC_MSG_NOTICE(* FASP include dir = ${with_fasp_incdir}/include)
        AC_MSG_NOTICE(* FASP lib dir = ${with_fasp_libdir}/lib)
        AC_MSG_NOTICE(* FASP libs = -lfasp)
        AC_DEFINE(USE_FASP, 1)
        AC_ADD_FLAGS(CPPFLAGS, -I${with_fasp_incdir}/include)
        AC_ADD_FLAGS(LDFLAGS, -L${with_fasp_libdir}/lib)
        AC_ADD_LIBS(LIBS, -lfasp)
    else
        AC_MSG_NOTICE(* FASP is disabled)
        AC_DEFINE(USE_FASP, 0)
    fi
else
    AC_DEFINE(USE_FASP, 0)
    AC_MSG_NOTICE(* FASP is disabled)
fi

#--------------------------------------- SUPERLU
AC_LANG(C)
AC_ARG_ENABLE(superlu,
[  --enable-superlu      enable SUPERLU support (default)
  --disable-superlu     disable SUPERLU support
  --with-superlu-libdir=DIR path for SUPERLU library
  --with-superlu-incdir=DIR path for SUPERLU header file],,enable_superlu=yes)
if test "$enable_superlu" != no; then
    AC_MSG_CHECKING([whether we have SUPERLU ])

    if test -z "$with_superlu_incdir" -a "$enable_superlu" != no; then
        for d in /usr/local/superlu /usr/local/superlu* /usr/local /usr \
            /opt/superlu /opt/superlu*; do
            if test -r ${d}/include/slu_ddefs.h; then
                with_superlu_incdir="$d"
                break
            fi
        done
        if test -z "$with_superlu_incdir"; then
            enable_superlu=no
        fi
    fi

    if test -z "$with_superlu_libdir" -a "$enable_superlu" != no; then
        for d in /usr/local/superlu /usr/local/superlu* /usr/local /usr \
            /opt/superlu /opt/superlu*; do
            if test -r ${d}/lib/libsuperlu.so -o -r ${d}/lib/libsuperlu.a; then
                with_superlu_libdir="$d"
                break
            fi
        done
        if test -z "$with_superlu_libdir"; then
            enable_superlu=no
        fi
    fi

    AC_MSG_RESULT($enable_superlu)
    if test "$enable_superlu" != no; then
        AC_MSG_NOTICE(* SUPERLU support enabled)
        AC_MSG_NOTICE(* SUPERLU include dir = ${with_superlu_incdir}/include)
        AC_MSG_NOTICE(* SUPERLU lib dir = ${with_superlu_libdir}/lib)
        AC_MSG_NOTICE(* SUPERLU libs = -lsuperlu)
        AC_DEFINE(USE_SUPERLU, 1)
        AC_ADD_FLAGS(CPPFLAGS, -I${with_superlu_incdir}/include)
        AC_ADD_FLAGS(LDFLAGS, -L${with_superlu_libdir}/lib)
        AC_ADD_LIBS(LIBS, -lsuperlu)
    else
        AC_MSG_NOTICE(* SUPERLU is disabled)
        AC_DEFINE(USE_SUPERLU, 0)
    fi
else
    AC_DEFINE(USE_SUPERLU, 0)
    AC_MSG_NOTICE(* SUPERLU is disabled)
fi

#--------------------------------------- PARDISO
AC_LANG(C)
AC_ARG_ENABLE(pardiso,
[  --enable-pardiso      enable PARDISO support (default)
  --disable-pardiso     disable PARDISO support
  --with-pardiso-libdir=DIR path for PARDISO library
  --with-pardiso-incdir=DIR path for PARDISO header file],,enable_pardiso=yes)
if test "$enable_pardiso" != no; then
    AC_MSG_CHECKING([whether we have PARDISO ])

    if test -z "$with_pardiso_incdir" -a "$enable_pardiso" != no; then
        for d in /usr/local/intel /opt/intel ; do
            if test -r ${d}/mkl/include/mkl_pardiso.h; then
                with_pardiso_incdir="$d/mkl/include/"
                break
            fi
        done
        if test -z "$with_pardiso_incdir"; then
            enable_pardiso=no
        fi
    fi

    if test -z "$with_pardiso_libdir" -a "$enable_pardiso" != no; then
        for d in /usr/local/pardiso /usr/local/pardiso* /usr/local/lib /usr/local/lib64 \
            /usr/lib /usr/lib64 /opt/pardiso /opt/pardiso*; do
            if test -r ${d}/libpardiso.so -o -r ${d}/libpardiso.a; then
                with_pardiso_libdir="$d"
                break
            fi
        done
        if test -z "$with_pardiso_libdir"; then
            enable_pardiso=no
        fi
    fi

    AC_MSG_RESULT($enable_pardiso)
    if test "$enable_pardiso" != no; then
        AC_MSG_NOTICE(* PARDISO support enabled)
        AC_MSG_NOTICE(* PARDISO include dir = ${with_pardiso_incdir})
        AC_MSG_NOTICE(* PARDISO lib dir = ${with_pardiso_libdir})
        AC_MSG_NOTICE(* PARDISO libs = -lpardiso)
        AC_DEFINE(USE_PARDISO, 1)
        AC_ADD_FLAGS(CPPFLAGS, -I${with_pardiso_incdir})
        AC_ADD_FLAGS(LDFLAGS, -L${with_pardiso_libdir})
        AC_ADD_LIBS(LIBS, -lpardiso)
    else
        AC_MSG_NOTICE(* PARDISO is disabled)
        AC_DEFINE(USE_PARDISO, 0)
    fi
else
    AC_DEFINE(USE_PARDISO, 0)
    AC_MSG_NOTICE(* PARDISO is disabled)
fi

#--------------------------------------- HSL_MI20
AC_LANG(C)
AC_ARG_ENABLE(hslmi20,
[  --enable-hslmi20      enable HSL_MI20 support (default)
  --disable-hslmi20     disable HSL_MI20 support
  --with-hslmi20-libdir=DIR path for HSL_MI20 library
  --with-hslmi20-incdir=DIR path for HSL_MI20 header file],,enable_hslmi20=yes)
if test "$enable_hslmi20" != no; then
    AC_MSG_CHECKING([whether we have HSL_MI20 ])

    if test -z "$with_hslmi20_incdir" -a "$enable_hslmi20" != no; then
        for d in /usr/local/hsl_mi20 /usr/local/hsl_mi20* /opt/hsl_mi20 /opt/hsl_mi20*; do
            if test -r ${d}/include/hsl_mi20d.h; then
                with_hslmi20_incdir="$d/include/"
                break
            fi
        done
        if test -z "$with_hslmi20_incdir"; then
            enable_hslmi20=no
        fi
    fi

    if test -z "$with_hslmi20_libdir" -a "$enable_hslmi20" != no; then
        for d in /usr/local/hsl_mi20 /usr/local/hsl_mi20* /opt/hsl_mi20 /opt/hsl_mi20*; do
            if test -r ${d}/lib/libhsl_mi20.so -o -r ${d}/lib/libhsl_mi20.a; then
                with_hslmi20_libdir="$d/lib"
                break
            fi
        done
        if test -z "$with_hslmi20_libdir"; then
            enable_hslmi20=no
        fi
    fi

    AC_MSG_RESULT($enable_hslmi20)
    if test "$enable_hslmi20" != no; then
        AC_MSG_NOTICE(* HSL_MI20 support enabled)
        AC_MSG_NOTICE(* HSL_MI20 include dir = ${with_hslmi20_incdir})
        AC_MSG_NOTICE(* HSL_MI20 lib dir = ${with_hslmi20_libdir})
        AC_MSG_NOTICE(* HSL_MI20 libs = -lhsl_mi20)
        AC_DEFINE(USE_HSL_MI20, 1)
        AC_ADD_FLAGS(CPPFLAGS, -I${with_hslmi20_incdir})
        AC_ADD_FLAGS(LDFLAGS, -L${with_hslmi20_libdir})
        AC_ADD_LIBS(LIBS, -lhsl_mi20)
    else
        AC_MSG_NOTICE(* HSL_MI20 is disabled)
        AC_DEFINE(USE_HSL_MI20, 0)
    fi
else
    AC_DEFINE(USE_HSL_MI20, 0)
    AC_MSG_NOTICE(* HSL_MI20 is disabled)
fi

#--------------------------------------- SXAMG
AC_LANG(C)
AC_ARG_ENABLE(sxamg,
[  --enable-sxamg      enable SXAMG support (default)
  --disable-sxamg     disable SXAMG support
  --with-sxamg-libdir=DIR path for SXAMG library
  --with-sxamg-incdir=DIR path for SXAMG header file],,enable_sxamg=yes)
if test "$enable_sxamg" != no; then
    AC_MSG_CHECKING([whether we have SXAMG ])

    if test -z "$with_sxamg_incdir" -a "$enable_sxamg" != no; then
        for d in /usr/local/sxamg /usr/local/sxamg* /opt/sxamg /opt/sxamg*; do
            if test -r ${d}/include/sxamg.h; then
                with_sxamg_incdir="$d/include/"
                break
            fi
        done
        if test -z "$with_sxamg_incdir"; then
            enable_sxamg=no
        fi
    fi

    if test -z "$with_sxamg_libdir" -a "$enable_sxamg" != no; then
        for d in /usr/local/sxamg /usr/local/sxamg* /opt/sxamg /opt/sxamg*; do
            if test -r ${d}/lib/libsxamg.so -o -r ${d}/lib/libsxamg.a; then
                with_sxamg_libdir="$d/lib"
                break
            fi
        done
        if test -z "$with_sxamg_libdir"; then
            enable_sxamg=no
        fi
    fi

    AC_MSG_RESULT($enable_sxamg)
    if test "$enable_sxamg" != no; then
        AC_MSG_NOTICE(* SXAMG support enabled)
        AC_MSG_NOTICE(* SXAMG include dir = ${with_sxamg_incdir})
        AC_MSG_NOTICE(* SXAMG lib dir = ${with_sxamg_libdir})
        AC_MSG_NOTICE(* SXAMG libs = -lsxamg)
        AC_DEFINE(USE_SXAMG, 1)
        AC_ADD_FLAGS(CPPFLAGS, -I${with_sxamg_incdir})
        AC_ADD_FLAGS(LDFLAGS, -L${with_sxamg_libdir})
        AC_ADD_LIBS(LIBS, -lsxamg)
    else
        AC_MSG_NOTICE(* SXAMG is disabled)
        AC_DEFINE(USE_SXAMG, 0)
    fi
else
    AC_DEFINE(USE_SXAMG, 0)
    AC_MSG_NOTICE(* SXAMG is disabled)
fi

#--------------------------------------- the end
if test -z "$LINKER"; then
    LINKER="$CXX"
fi

# hack ITOL and FASP, both contain "qsplit_"
if $CXX -V 2>&1 | grep -q Intel; then
    LDFLAGS="$LDFLAGS "
elif $CXX -v 2>&1 | grep -q clang; then
    LDFLAGS="$LDFLAGS "
elif $CXX -v 2>&1 | grep -q gcc; then
    LDFLAGS="$LDFLAGS -z muldefs "
fi

LIBS="-llssp $LIBS"

AC_OUTPUT
