SparseEigenSolver
=================

Eigenvalue solver for sparse matrices in C++ using ARPACK

To compile on a Mac type the following:

g++ -g -Wall main.cpp eigs.cpp eigs.h -o main -framework Accelerate -L[path to ARPACK] -larpack

and run with:

./main