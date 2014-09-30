#ifndef SparseEigenSolver_eigs_h
#define SparseEigenSolver_eigs_h

#include <string>

/*-- Class declaration --*/

class eigs
{
public:
    // Constructor & destructor
    eigs(int nev, int ncv, int n);
    ~eigs();
    
    // Computes eigenvalues
    void eig_vals(double** A, double *evals);
    
private:
    /*-- Functions --*/
    // Calculate matrix-vector product A*in and return in 'out'
    void av(double **A, int n, double *in, double *out);
    void error_string(int err, const char *fct);
    
    /*-- Parameters (dnaupd) --*/
    /* dnaupd computes the Ritz values for the specified eigenvalue problem. */
    /* These are then used by dneupd to find the actual eigenvalues */

    // Number of eigenvalues to compute
    int nev;
    
    // Order of matrix A
    int n;
    
    // Reverse communication flag. Signifies type of operation to be performed
    int ido;
    
    // Specifies the type of the matrix B that defines the semi-inner product for the operator A
    // Identity signifies a standard eigenvalue problem
    char *bmat;
    
    // Choose the largest 'nev' eigenvalues
    char *which;
    
    // Tolerance
    double tol;
    
    // Residual vector. Initialised randomly by the function dnaupd
    double *resid;
    
    // Indicates how many Arnoldi vectors are generated at each iteration (basis vectors used)
    // Must satisfy 2 <= NCV-NEV && NCV <= N
    int ncv;
    
    // Stores the final set of Arnoldi basis vectors. Represents an (n * ncv) matrix
    double *v;
    
    // Leading dimension of V, i.e first index in the dimension specified.
    int ldv;
    
    // Specify the various modes for the routines below.
    // iparam[0] = select the implicit shifts
    // iparam[2] = maximum number of iterations
    // iparam[6] = mode, i.e. standard eigenvalue problem A*x = lambda*x
    int iparam[11];
    
    // Pointers to specify various starting locations for the matrices/vectors
    // used by the Arnoldi iteration
    int ipntr[14];
    
    // Following three variables are used for the reverse communication process
    double *workd;
    
    double *workl;
    
    int lworkl;
    
    // Error flag for the dnaupd routine. Stores convergence information
    int info;
    
    /*-- Parameters (dneupd) --*/
    /* dneupd calculates the Eigenvalues using the computed Ritz values found by dnaupd */
    
    // If false do not compute the eigenvectors
    int rvec;
    
    // Following two used only if finding eigenvectors.
    char *howmny;
    int *select;
    
    // Contain the real part of converged eigenvalues
    double *dr;
    
    // Stores the imaginary part of the converged eigenvalues
    double *di;
    
    // Stores complex eigenvectors
    double *z;
    
    // Leading dimension of z
    int ldz;
    
    // Used for reverse communication
    double *workev;
    
    // Represent the real and imaginary part of the shift respectively.
    double sigmar, sigmai;
};

#endif