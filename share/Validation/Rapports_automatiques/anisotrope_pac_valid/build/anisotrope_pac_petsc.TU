************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/home/vd256574/workspace/pemfc/pemfc-exec on a linux_opt named gre056837 with 1 processor, by vd256574 Mon Sep 24 10:47:54 2018
Using Petsc Release Version 3.6.3, Dec, 03, 2015 

                         Max       Max/Min        Avg      Total 
Time (sec):           1.398e+02      1.00000   1.398e+02
Objects:              3.235e+03      1.00000   3.235e+03
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
 0:      Main Stage: 1.3542e+02  96.9%  0.0000e+00   0.0%  0.000e+00   0.0%  0.000e+00        0.0%  0.000e+00   0.0% 
 1:        KSPSolve: 4.3604e+00   3.1%  0.0000e+00   0.0%  0.000e+00   0.0%  0.000e+00        0.0%  0.000e+00   0.0% 

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

KSPSetUp            1614 1.0 1.8716e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
PCSetUp             1614 1.0 5.8223e+01 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00 42  0  0  0  0  43  0  0  0  0     0
MatLUFactorSym      1614 1.0 3.6397e+01 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00 26  0  0  0  0  27  0  0  0  0     0
MatLUFactorNum      1614 1.0 1.7815e+01 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00 13  0  0  0  0  13  0  0  0  0     0
MatAssemblyBegin    1614 1.0 3.2020e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatAssemblyEnd      1614 1.0 2.1668e-01 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetRowIJ         1614 1.0 1.9773e-01 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering      1614 1.0 3.9327e+00 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  3  0  0  0  0   3  0  0  0  0     0
VecSet                 4 1.0 6.6757e-06 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAssemblyBegin    3228 1.0 1.8954e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAssemblyEnd      3228 1.0 1.0419e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0

--- Event Stage 1: KSPSolve

KSPSolve            1614 1.0 4.3566e+00 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  3  0  0  0  0 100  0  0  0  0     0
PCApply             1614 1.0 4.3509e+00 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  3  0  0  0  0 100  0  0  0  0     0
MatSolve            1614 1.0 4.3502e+00 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  3  0  0  0  0 100  0  0  0  0     0
VecCopy             1614 1.0 2.3232e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet              1614 1.0 1.7390e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

       Krylov Solver     2              2         2272     0
      Preconditioner     2              2         1968     0
              Matrix  1614           1614    757798824     0
           Index Set  1614           1614     70008864     0
              Vector     2              2       173472     0
              Viewer     1              0            0     0

--- Event Stage 1: KSPSolve

========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
#PETSc Option Table entries:
-log_summary anisotrope_pac_petsc.TU
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

