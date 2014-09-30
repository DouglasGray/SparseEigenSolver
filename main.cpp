#include <iostream>
#include <stdlib.h>
#include "eigs.h"

int main()
{
    int n, nev, ncv, r;
    double *evals, **A;
    
    // Set matrix and problem attributes
    n = 100;
    nev = 1;
    ncv = (4*nev) > n ? n : 4*nev;
    
    if (2 > ncv - nev)
    {
        std::cout << "Incorrect number of Arnoldi basis vectors. ";
        std::cout << "Either decrease NEV or increase NCV." << std::endl;
        return 0;
    }
    
    // Generate A randomly
    A = new double*[n];
    
    for (int i=0; i< n; i++)
        A[i] = new double[n];
    
    for (int i=0; i < n; i++)
        for (int j=0; j < n; j++)
        {
            r = std::rand() % 10 + 1;
            if (r <= 2)
                A[i][j] = i + j + r;
        }
    
    // Calculate eigenvalues
    evals = new double[nev];
    try
    {
        eigs eig_solver(nev, ncv, n);
        eig_solver.eig_vals(A, evals);
        
        // If converged then output eigenvalues
        for (int i = 0; i < nev; i++)
            std::cout << "Eigenvalue " << i + 1 << ": " << evals[i] << std::endl;
        
    }
    catch (int e) {}
    
    delete evals;
    delete[] A;
    return 0;
    
    return 0;
}