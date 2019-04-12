************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/nobackup/vd256574/TRUST/TRUST-1.7.6/exec/TRUST_mpi_opt on a linux_opt named gre056837 with 1 processor, by vd256574 Mon Sep 24 15:37:18 2018
Using Petsc Release Version 3.6.3, Dec, 03, 2015 

                         Max       Max/Min        Avg      Total 
Time (sec):           3.256e+00      1.00000   3.256e+00
Objects:              3.500e+01      1.00000   3.500e+01
Flops:                0.000e+00      0.00000   0.000e+00  0.000e+00
Flops/sec:            0.000e+00      0.00000   0.000e+00  0.000e+00
MPI Messages:         0.000e+00      0.00000   0.000e+00  0.000e+00
MPI Message Lengths:  0.000e+00      0.00000   0.000e+00  0.000e+00
MPI Reductions:       0.000e+00      0.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 3.2181e+00  98.8%  0.0000e+00   0.0%  0.000e+00   0.0%  0.000e+00        0.0%  0.000e+00   0.0% 
 1:        KSPSolve: 3.7591e-02   1.2%  0.0000e+00   0.0%  0.000e+00   0.0%  0.000e+00        0.0%  0.000e+00   0.0% 

------------------------------------------------------------------------------------------------------------------------
See the 'Profiling' chapter of the users' manual for details on interpreting output.
Phase summary info:
   Count: number of times phase was executed
   Time and Flops: Max - maximum over all processors
                   Ratio - ratio of maximum to minimum over all processors
   Mess: number of messages sent
   Avg. len: average message length (bytes)
   Reduct: number of global reductions
   Global: entire computation
   Stage: stages of a computation. Set stages with PetscLogStagePush() and PetscLogStagePop().
      %T - percent time in this phase         %F - percent flops in this phase
      %M - percent messages in this phase     %L - percent message lengths in this phase
      %R - percent reductions in this phase
   Total Mflop/s: 10e-6 * (sum of flops over all processors)/(max time over all processors)
------------------------------------------------------------------------------------------------------------------------
Event                Count      Time (sec)     Flops                             --- Global ---  --- Stage ---   Total
                   Max Ratio  Max     Ratio   Max  Ratio  Mess   Avg len Reduct  %T %F %M %L %R  %T %F %M %L %R Mflop/s
------------------------------------------------------------------------------------------------------------------------

--- Event Stage 0: Main Stage

KSPSetUp              14 1.0 3.8147e-06 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
PCSetUp               14 1.0 5.4413e-01 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00 17  0  0  0  0  17  0  0  0  0     0
MatLUFactorSym        14 1.0 3.2720e-01 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00 10  0  0  0  0  10  0  0  0  0     0
MatLUFactorNum        14 1.0 1.5277e-01 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  5  0  0  0  0   5  0  0  0  0     0
MatAssemblyBegin      14 1.0 4.0531e-06 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatAssemblyEnd        14 1.0 4.4286e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetRowIJ           14 1.0 4.0545e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        14 1.0 6.2984e-02 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  2  0  0  0  0   2  0  0  0  0     0
VecSet                 4 1.0 7.1526e-06 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAssemblyBegin      28 1.0 9.5367e-07 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAssemblyEnd        28 1.0 4.7684e-07 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0

--- Event Stage 1: KSPSolve

KSPSolve              14 1.0 3.7560e-02 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0 100  0  0  0  0     0
PCApply               14 1.0 3.7507e-02 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0 100  0  0  0  0     0
MatSolve              14 1.0 3.7501e-02 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0 100  0  0  0  0     0
VecCopy               14 1.0 2.0027e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                14 1.0 1.7643e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

       Krylov Solver     2              2         2272     0
      Preconditioner     2              2         1968     0
              Matrix    14             14      6573224     0
           Index Set    14             14       607264     0
              Vector     2              2       173472     0
              Viewer     1              0            0     0

--- Event Stage 1: KSPSolve

========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
#PETSc Option Table entries:
-log_summary prepare_coeff_petsc.TU
-mat_mumps_icntl_14 35
#End of PETSc Option Table entries
Compiled without FORTRAN kernels
Compiled with full precision matrices (default)
sizeof(short) 2 sizeof(int) 4 sizeof(long) 8 sizeof(void*) 8 sizeof(PetscScalar) 8 sizeof(PetscInt) 4
Configure options: --prefix=/nobackup/vd256574/TRUST/TRUST-1.7.6/lib/src/LIBPETSC/petsc/linux_opt --with-single-library --with-shared-libraries=0 --with-debugging=0 --with-errorchecking=1 --COPTFLAGS="   -O3 -fPIC " --CXXOPTFLAGS="   -O3 -fPIC " --FOPTFLAGS="  -O3 -fPIC " --with-fortran=yes --with-clean=1 --with-ssl=0 --download-scalapack=../scalapack-2.0.2.tgz --download-mumps=../MUMPS_5.0.0-p1.tar.gz --download-metis=../metis-5.1.0-p1.tar.gz --download-superlu_dist=../superlu_dist_4.1.tar.gz --download-suitesparse=../SuiteSparse-4.4.3.tar.gz --download-parmetis=../parmetis-4.0.3-p2.tar.gz --download-ptscotch=../ptscotch.tar.gz --download-hypre=../hypre-2.10.0b-p2-TRUST.tar.gz --with-blas-lapack-dir=/nobackup/vd256574/TRUST/TRUST-1.7.6/lib/src/LIBLAPACK --with-mpi-dir=/nobackup/vd256574/TRUST/TRUST-1.7.6/lib/src/LIBMPI/mpich --with-x=0
-----------------------------------------
Libraries compiled on Thu Jul 26 09:35:42 2018 on gre056837 
Machine characteristics: Linux-4.9.0-6-amd64-x86_64-with-debian-9.5
Using PETSc directory: /home/vd256574/tmp-global/.tmp_TRUST_gre056837/build/LIBPETSC/petsc-3.6.3
Using PETSc arch: linux_opt
-----------------------------------------

Using C compiler: /nobackup/vd256574/TRUST/TRUST-1.7.6/lib/src/LIBMPI/mpich/bin/mpicc  -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -O3 -fPIC  ${COPTFLAGS} ${CFLAGS}
Using Fortran compiler: /nobackup/vd256574/TRUST/TRUST-1.7.6/lib/src/LIBMPI/mpich/bin/mpif90  -Wall -Wno-unused-variable -ffree-line-length-0 -Wno-unused-dummy-argument -O3 -fPIC   ${FOPTFLAGS} ${FFLAGS} 
-----------------------------------------

Using include paths: -I/home/vd256574/tmp-global/.tmp_TRUST_gre056837/build/LIBPETSC/petsc-3.6.3/linux_opt/include -I/home/vd256574/tmp-global/.tmp_TRUST_gre056837/build/LIBPETSC/petsc-3.6.3/include -I/home/vd256574/tmp-global/.tmp_TRUST_gre056837/build/LIBPETSC/petsc-3.6.3/include -I/home/vd256574/tmp-global/.tmp_TRUST_gre056837/build/LIBPETSC/petsc-3.6.3/linux_opt/include -I/nobackup/vd256574/TRUST/TRUST-1.7.6/lib/src/LIBPETSC/petsc/linux_opt/include -I/nobackup/vd256574/TRUST/TRUST-1.7.6/lib/src/LIBMPI/mpich/include
-----------------------------------------

Using C linker: /nobackup/vd256574/TRUST/TRUST-1.7.6/lib/src/LIBMPI/mpich/bin/mpicc
Using Fortran linker: /nobackup/vd256574/TRUST/TRUST-1.7.6/lib/src/LIBMPI/mpich/bin/mpif90
Using libraries: -Wl,-rpath,/home/vd256574/tmp-global/.tmp_TRUST_gre056837/build/LIBPETSC/petsc-3.6.3/linux_opt/lib -L/home/vd256574/tmp-global/.tmp_TRUST_gre056837/build/LIBPETSC/petsc-3.6.3/linux_opt/lib -lpetsc -Wl,-rpath,/nobackup/vd256574/TRUST/TRUST-1.7.6/lib/src/LIBPETSC/petsc/linux_opt/lib -L/nobackup/vd256574/TRUST/TRUST-1.7.6/lib/src/LIBPETSC/petsc/linux_opt/lib -lHYPRE -L/nobackup/vd256574/TRUST/TRUST-1.7.6/lib/src/LIBMPI/mpich/lib -L/usr/lib/gcc/x86_64-linux-gnu/6 -L/usr/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -Wl,-rpath,/nobackup/vd256574/TRUST/TRUST-1.7.6/lib/src/LIBMPI/mpich/lib -lmpicxx -lstdc++ -lsuperlu_dist_4.1 -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lumfpack -lklu -lcholmod -lbtf -lccolamd -lcolamd -lcamd -lamd -lsuitesparseconfig -Wl,-rpath,/nobackup/vd256574/TRUST/TRUST-1.7.6/lib/src/LIBLAPACK -L/nobackup/vd256574/TRUST/TRUST-1.7.6/lib/src/LIBLAPACK -llapack -lblas -lparmetis -lmetis -lpthread -lhwloc -lptesmumps -lptscotch -lptscotcherr -lscotch -lscotcherr -lm -lmpifort -lgfortran -lm -lgfortran -lm -lquadmath -lm -lmpicxx -lstdc++ -lrt -lm -lz -L/nobackup/vd256574/TRUST/TRUST-1.7.6/lib/src/LIBMPI/mpich/lib -L/usr/lib/gcc/x86_64-linux-gnu/6 -L/usr/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -ldl -Wl,-rpath,/nobackup/vd256574/TRUST/TRUST-1.7.6/lib/src/LIBMPI/mpich/lib -lmpi -lgcc_s -ldl 
-----------------------------------------

