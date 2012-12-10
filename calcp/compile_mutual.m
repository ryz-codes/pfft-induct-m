% The following code compiles mutual.c with OpenMP

%%% VC++ flags
% mex mutual.c COMPFLAGS="/openmp /O2 /Wall $COMPFLAGS"

%%% GCC compiler flags
mex mutual.c -v CFLAGS="-fopenmp -Wall -O3 \$CFLAGS" LDFLAGS="-fopenmp \$LDFLAGS"

test_mutualfil_speed