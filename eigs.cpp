#include <iostream>
#include "eigs.h"

/*-- External function declarations --*/
extern "C" void dnaupd_(int *ido, char *bmat, int *n, char *which,
                        int *nev, double *tol, double *resid, int *ncv,
                        double *v, int *ldv, int *iparam, int *ipntr,
                        double *workd, double *workl, int *lworkl,
                        int *info);

extern "C" void dneupd_(int *rvec, char *All, int *select, double *dr, double *di,
                        double *z, int *ldz, double *sigmar, double* sigmai,
                        double *workev, char *bmat, int *n, char *which, int *nev,
                        double *tol, double *resid, int *ncv, double *v,
                        int *ldv, int *iparam, int *ipntr, double *workd,
                        double *workl, int *lworkl, int *ierr);

/*-- Class functions --*/

// Constructor
eigs::eigs(int _nev, int _ncv, int _n)
{
    ncv = _ncv;
    nev = _nev;
    n = _n;
    
    // dnaupd parameters
    ido = 0;
    bmat = (char *)"I";
    which = (char *)"LM";
    tol = 0.0;
    resid = new double[n];
    v = new double[n * ncv];
    ldv = n;
    iparam[0] = 1;
    iparam[2] = 500;
    iparam[6] = 1;
    workd = new double[3*n];
    workl = new double[3*ncv*(ncv+2)];
    lworkl = 3*ncv*(ncv+2);
    info = 0;
    
    // dneupd paramaters
    rvec = false;
    howmny = (char *)"A";
    select = new int[ncv];
    dr = new double[nev + 1];
    di = new double[nev + 1];
    z = new double[n * ncv];
    ldz = n;
    workev = new double[3*ncv];
}

// Destructor
eigs::~eigs()
{
    delete resid;
    delete v;
    delete workd;
    delete workl;
    delete select;
    delete dr;
    delete di;
    delete z;
    delete workev;
}

// Calculates the eigenvalues
void eigs::eig_vals(double** A, double *evals)
{
    // Repeatedly call dnaupd until either the eigenvalues converge, the maximum number of
    // iterations are reached or an error is signalled by the 'info' variable.
    do
    {
        dnaupd_(&ido, bmat, &n, which, &nev, &tol, resid,
                &ncv, v, &ldv, iparam, ipntr, workd, workl,
                &lworkl, &info);
        
        if ((ido == 1)||(ido == -1))
            av(A, n, workd+ipntr[0]-1, workd+ipntr[1]-1);
        
    } while ((ido == 1)||(ido == -1));
    
    // Check convergence
    if (info != 0)
    {
        error_string(info, "dnaupd");
    }
    
    // If successful find the eigenvalues
    dneupd_(&rvec, howmny, select, dr, di, z, &ldz,
            &sigmar, &sigmai, workev, bmat, &n, which, &nev,
            &tol, resid, &ncv, v, &ldv, iparam, ipntr,
            workd, workl, &lworkl, &info);
    
    if (info != 0)
    {
        error_string(info, "dneupd");
    }
    
    // Copy the found eigenvalues into the array
    for (int i = 0; i < nev; i++)
        evals[i] = dr[nev-1-i];
}

// Helper function to compute the matrix-vector product A*in.
// Returns the result in 'out'
void eigs::av(double **A, int n, double *in, double *out)
{
    for (int i = 0; i < n; i++)
        out[i] = 0;
    
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            out[i] += A[i][j] * in[j];
}

// Helper function which prints the errors thrown by the ARPACK functions
void eigs::error_string(int err, const char *fct)
{
    try
    {
        if (strcmp(fct, "dnaupd") == 0)
        {
            std::cout << "eigs::error_string (dnaupd): ";
            switch (err)
            {
                case 1:
                    std::cout << "Maximum number of iterations reached" << std::endl;
                    break;
                case 3:
                    std::cout << "No shifts could be applied during a cycle of the ";
                    std::cout << "Implicitly restarted Arnoldi iteration. One possibility ";
                    std::cout << "is to increase the size of NCV relative to NEV." << std::endl;
                    break;
                case -1:
                    std::cout << "Matrix order N must be positive." << std::endl;
                    throw 0;
                case -2:
                    std::cout << "Number of Arnoldi basis vectors NCV must be positive." << std::endl;
                    throw 0;
                case -3:
                    std::cout << "Require NCV-NEV >= 2 and less than or equal to N." << std::endl;
                    throw 0;
                case -4:
                    std::cout << "Max no. of iterations must be greater than zero." << std::endl;
                    throw 0;
                case -5:
                    std::cout << "WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'." << std::endl;
                    throw 0;
                case -6:
                    std::cout << "BMAT must be one of 'I' or 'G'." << std::endl;
                    throw 0;
                case -7:
                    std::cout << "Length of private work array is not sufficient." << std::endl;
                    throw 0;
                case -8:
                    std::cout << "Error return from LAPACK eigenvalue calculation." << std::endl;
                    throw 0;
                case -9:
                    std::cout << "Starting vector is zero." << std::endl;
                    throw 0;
                case -10:
                    std::cout << "IPARAM(7) must be 1,2,3,4." << std::endl;
                    throw 0;
                case -11:
                    std::cout << "IPARAM(7) = 1 and BMAT = 'G' are incompatable." << std::endl;
                    throw 0;
                case -12:
                    std::cout << "IPARAM(1) must be equal to 0 or 1." << std::endl;
                    throw 0;
                case -9999:
                    std::cout << "Could not build an Arnoldi factorization." << std::endl;
                    throw 0;
            }
        }
        else if (strcmp(fct, "dneupd") ==0)
        {
            std::cout << "eigs::error_string (dneupd): ";
            switch (err)
            {
                case 1:
                    std::cout << "The Schur form computed by LAPACK routine dlahqr ";
                    std::cout << "could not be reordered by LAPACK routine dtrsen." << std::endl;
                    throw 0;
                case -1:
                    std::cout << "Matrix order N must be positive." << std::endl;
                    throw 0;
                case -2:
                    std::cout << "Number of Arnoldi basis vectors NCV must be positive." << std::endl;
                    throw 0;
                case -3:
                    std::cout << "Require NCV-NEV >= 2 and less than or equal to N." << std::endl;
                    throw 0;
                case -5:
                    std::cout << "WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'." << std::endl;
                    throw 0;
                case -6:
                    std::cout << "BMAT must be one of 'I' or 'G'." << std::endl;
                    throw 0;
                case -7:
                    std::cout << "Length of private work array is not sufficient." << std::endl;
                    throw 0;
                case -8:
                    std::cout << "Error return from calculation of a real Schur form.";
                    std::cout << "Informational error from LAPACK routine dlahqr." << std::endl;
                    throw 0;
                case -9:
                    std::cout << "Error return from calculation of eigenvectors.";
                    std::cout << "Informational error from LAPACK routine dtrevc." << std::endl;
                    throw 0;
                case -10:
                    std::cout << "IPARAM(7) must be 1,2,3,4." << std::endl;
                    throw 0;
                case -11:
                    std::cout << "IPARAM(7) = 1 and BMAT = 'G' are incompatable." << std::endl;
                    throw 0;
                case -12:
                    std::cout << "HOWMNY = 'S' not yet implemented" << std::endl;
                    throw 0;
                case -13:
                    std::cout << "HOWMNY must be one of 'A' or 'P' if RVEC = true." << std::endl;
                    throw 0;
                case -14:
                    std::cout << "DNAUPD did not find any eigenvalues to sufficient accuracy." << std::endl;
                    throw 0;
            }
        }
        else
        {
            std::cout << "eigs::error_string: ARPACK function " << fct << " not recognised.";
            throw 0;
        }
    }
    catch (int e)
    {
        std::cout << "Terminating program." << std::endl;
        throw e;
    }
}
