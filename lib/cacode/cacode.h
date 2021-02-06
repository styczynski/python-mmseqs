// Public domain code from Yi-Kuo Yu & Stephen Altschul, NCBI

#ifndef __H_INCLUDE_LAMBDA_CALCULATOR_HH
#define __H_INCLUDE_LAMBDA_CALCULATOR_HH

void ludcmp(double **a, int n, int *indx, double *d);
float **matrix(int nrl, int nrh, int ncl, int nch);
double *dvector(int nl, int nh);
double **dmatrix(int nrl, int nrh, int ncl, int nch);
int *ivector(int nl, int nh);
int **imatrix();
float **submatrix(float **a, int oldrl, int oldrh, int oldcl, int oldch, int newrl, int newcl);
void free_dvector(double *v, int nl, int nh);
void free_ivector(int *v, int nl, int nh);
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
void nrerror(const char *error_text);
double calculate_lambda(const double **mat_b, int alpha_size, double *p, double *q);
void lubksb(double **a, int n, int *indx, double b[]);

#endif