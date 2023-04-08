# LSSP

LSSP is a **L**inear **S**olver library for **SP**arse linear system, **Ax = b**. The package is designed for Linux, Unix and Mac systems. It is also possible to compile under Windows. The code is written by a small set of C++ and C,and it is serial.

LSSP has many built-in solvers and preconditioners, and it also has interfaces to external packages, such as PETSc, MUMPS, FASP, UMFPACK (SuiteSparse), KLU (SuiteSparse), LASPACK, ITSOL, LIS, QR\_MUMPS, PARDISO, SUPERLU, and HSL MI20.

## Internal solvers
GMRES, LGMRES, BICGSTAB, BICGSTABL, BICGSAFE, CG, CGS, GPBICG, CR, CRS, BICRSTAB, BICRSAFE, GPBICR, QMRCGSTAB, TFQMR, ORTHOMIN, and IDRS

## External solvers
LASPACK, UMFPACK, KLU, MUMPS, PETSC, ITSOL, LIS, QR\_MUMPS, FASP, AMG (Algebraic Multi-Grid), SUPERLU, PARDISO, and HSL MI20 (AMG, Algebraic Multi-Grid).

PETSC, LIS, and FASP are solver collection packages, which implement their own solvers and preconditioners. For example, FASP implements CG, BiCGstab, MinRes, GMRES, VGMRES, VFGMRES, GCG, GCR, SCG, SBiCGstab, SMinRes, SGMRES, SVGMRES, SVFGMRES, SGCG, AMG and FMG.

PETSC, MUMPS and LIS are parallel solver packages and they can be (must be) compiled as serial packages. All external packages are optional.


## Internal preconditioners
ILUK, ILUT, and block-wise ILUK.

## External preconditioners
Block-wise ILUT, VBILUT, VBILUK, ARMS, FMG (FASP), AMG (FASP), and MI20 (AMG).

# How-to
## Configuration
The simplest way to configure is to run command:
```
./configure
```
This command will try to find optional packages from certain directories. Searching details are included by configure.in and some are explained below.

## Options
Run command:
```
./configure --help
```

Most optional packages are enabled by default. However, a package can be disabled when configuring, such as "--disable-itsol" to disable ITSOL. When a package needs an include path and a library path, they can be set by configure, such as --with-itsol-libdir=DIR and --with-itsol-incdir=DIR, which set library and include paths for ITSOL. If configure cannot find correct paths, users can help configure by using options.


## Compilation
After configuration, a simple **make** command can compile the package:
```
make
```

## Installation
Run command:
```
make install
```
The package will be installed to a directory. The default is /usr/local/lssp/. A different directory can be set by --prefix=DIR.

