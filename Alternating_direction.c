#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#define eps 1.e-7

double f(int i, int j, int n);
void solve05 (double **A, double **A05, int n, double B, double *F, double *alpha);
void solve (double **A, double **A05, int n, double B, double *F, double *alpha);
void prog (double **A, int n, double *F, int i, int ind, double *alpha);
void print (double **A, int n, FILE *file);
void print_U (int n);
double norm(double **B, double **A, int n);
double norm1(double **A, int n);


int main (void)
{
    int n = 256;
    int i, j, ind = 0;
    double **A, **A05, *F, *alpha, **B;
    double omega = 2*n*n*sin(M_PI/n);
    double B1 = (2 + omega/n/n);
    double B2 = (2 - omega/n/n);
    double m = 1.;
    FILE *file;
    file = fopen("out.txt", "w");
    A=(double**)malloc((unsigned int)(n+1)*sizeof(double*));
    A05=(double**)malloc((unsigned int)(n+1)*sizeof(double*));
    B=(double**)malloc((unsigned int)(n+1)*sizeof(double*));
    F=(double*)malloc((unsigned int)(n)*sizeof(double));
    alpha=(double*)malloc((unsigned int)(n+1)*sizeof(double));
    alpha[1] = 0;
    for (j = 1; j < n; j++)
        alpha[j+1] = 1/(B1-alpha[j]);
    for (i = 0; i < (n+1); i++)
    {
        A[i] = (double*)malloc((unsigned int)(n+1) * sizeof(double));
        A05[i] = (double*)malloc((unsigned int)(n+1) * sizeof(double));
        B[i] = (double*)malloc((unsigned int)(n+1) * sizeof(double));
        if (i==0 || i==n)
        {
            for (j = 0; j < (n+1); j++)
            {
                A[i][j] = 0.;
                A05[i][j] = 0.;
            }
        }
        else
        {
            for (j = 0; j < (n+1); j++)
                A[i][j] = 0.;
            A05[i][0] = 0.;
            A05[i][n] = 0.;
        }
    }
    while (m > eps)
    {
        for (i = 0; i < n+1; i++)
        {
            for (j = 0; j < n+1; j++)
                B[i][j] = A[i][j];
        }
        solve05(A, A05, n, B2, F, alpha);
        solve(A, A05, n, B2, F, alpha);
        m = norm(B, A, n);
        ind++;
    }
    printf("Количество итераций = %d\n", ind);
    for (i = 0; i < n+1; i++)
    {
        for (j = 0; j < n+1; j++)
            B[i][j] = fabs(A[i][j] - sin(M_PI*i/n)*sin(M_PI*j/n));
    }
    printf("m = %lf\n", norm1(B, n));
    print(B, n, file);
    for (i = 0; i < (n+1); i++)
    {
        free(A[i]);
        free(A05[i]);
        free(B[i]);
    }
    free(A);
    free(A05);
    free(B);
    free(F);
    free(alpha);
    fclose(file);
}

void solve05 (double **A, double **A05, int n, double B, double *F, double *alpha)
{   
    for (int j = 1; j < n; j++)
    {
        for (int i = 1; i < n; i++)
            F[i]=A[i][j+1]-B*A[i][j]+A[i][j-1]+f(i,j,n)/n/n;
        prog(A05, n, F, j, 0, alpha);
    }
}

void solve (double **A, double **A05, int n, double B, double *F, double *alpha)
{
    for (int i = 1; i < n; i++)
    {
        for (int j = 1; j < n; j++)
            F[j]=A05[i+1][j]-B*A05[i][j]+A05[i-1][j]+f(i,j,n)/n/n;
        prog(A, n, F, i, 1, alpha);
    }
}

void prog (double **A, int n, double *F, int i, int ind, double *alpha)
{
    int j;
    double *beta;
    beta=(double*)malloc((unsigned int)(n+1)*sizeof(double));
    beta[1] = 0.;
    for (j = 1; j < n; j++)
        beta[j+1] = alpha[j+1] * (F[j] + beta[j]);
    if (ind == 0)
    {
        A[n][i] = 0.;
        for (j = n-1; j > 0; j--)
            A[j][i] = A[j+1][i] * alpha[j+1] + beta[j+1];
    }
    if (ind == 1)
    {
        A[i][n] = 0.;
        for (j = n-1; j > 0; j--)
            A[i][j] = A[i][j+1] * alpha[j+1] + beta[j+1];
    }
    free(beta);
}

double f(int i, int j, int n)                            //при u = sin(pi*x)*sin(pi*y)
{
    return 2*M_PI*M_PI*sin(M_PI*i/n)*sin(M_PI*j/n);
}

void print (double **A, int n, FILE *file)
{
    for (int i = 0; i < (n+1); i++)
    {
        for (int j = 0; j < (n+1); j++)
        {
            fprintf(file, "%lf", A[i][j]);
            if (j!=n)
                fprintf(file, " ");
        }
        if (i!=n)
                fprintf(file, " ");
    }
}

double norm(double **A, double **B, int n)
{
    double m;
    m = fabs(A[0][0]-B[0][0]);
    for (int i=0; i<n+1; i++)
    {
        for (int j=0; j<n+1; j++)
        {
            if (fabs(A[i][j] - B[i][j]) > m)
                m = fabs(A[i][j] - B[i][j]);
        }
    }
    return m;
}

double norm1(double **A, int n)
{
    double m;
    m = fabs(A[0][0]);
    for (int i=0; i<n+1; i++)
    {
        for (int j=0; j<n+1; j++)
        {
            if (fabs(A[i][j]) > m)
                m = fabs(A[i][j]);
        }
    }
    return m;
}