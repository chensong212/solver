/*
 *
 *  Created on: 2016.04.18
 *  Purpose: conduct the time march step
 *
 */

#include<string.h>
#include<math.h>
#include"comm.h"

/*------------------------------------------------------
 * update the conservative variables for perfect-gas
 * use Runge-kutta method
 * ------------------------------------------------*/
void updateq(int nc, double **q, double **rhs, double dt, int idt)
{
	int i, j, ii, jj, ic, ic1;
	double rho, rhoU, rhoV, rhoE, rrho, yas;
	double coef[3][3] = {{1., 3./4., 1./3.}, {0, 1./4., 2./3.},
							  {1., 1./4., 2./3.}};

	for(i=0; i<config1.ni; i++)
	{
		ii = i + config1.Ng;
		for(j=0; j<config1.nj; j++)
		{
			jj  = j + config1.Ng;

			ic  = i*config1.nj + j;
			ic1 = ii*J0 + jj;
			yas = mesh.yaks[ic1];

			rho  = q[ic][0];
			rhoU = q[ic][0]*q[ic][1];
			rhoV = q[ic][0]*q[ic][2];
			rhoE = q[ic][0]*q[ic][3];

			q[ic][0] = coef[0][idt]*qo[ic][0]           + coef[1][idt]*rho  + coef[2][idt]*dt*rhs[ic][0]/yas;
			q[ic][1] = coef[0][idt]*qo[ic][0]*qo[ic][1] + coef[1][idt]*rhoU + coef[2][idt]*dt*rhs[ic][1]/yas;
			q[ic][2] = coef[0][idt]*qo[ic][0]*qo[ic][2] + coef[1][idt]*rhoV + coef[2][idt]*dt*rhs[ic][2]/yas;
			q[ic][3] = coef[0][idt]*qo[ic][0]*qo[ic][3] + coef[1][idt]*rhoE + coef[2][idt]*dt*rhs[ic][3]/yas;

			rrho = 1./q[ic][0];
			q[ic][1] = q[ic][1]*rrho;
			q[ic][2] = q[ic][2]*rrho;
			q[ic][3] = q[ic][3]*rrho;
		}
	}
}

/*---------------------------------------------------
 * update the conservative variables for real-gas
 * use Point-implicit method
 * ------------------------------------------------*/
void updateqs(int nc, double **q, double **qs, double **rhs, double dt, int idt)
{
  // to be added...
}

/*---------------------------------------------------
 *Gaussian Elimination with pivoting to solve Ax = b.
 *Usage:     gauss(a,b,x,n)

		   a   - Matrix a[n][n]
		   b   - Right hand side vector b[n]
		   x   - Desired solution vector
		   n   - Matrix dimensions
 * ------------------------------------------------*/
void gauss(double **a, double *b, double *x, int n)
{
	int   i, j, k, m, rowx;
	double xfac, temp, amax;
	int  matrix_print_off(int nr, int nc, double **A);

	/*--------1. Do the forward reduction step.--------*/
	rowx = 0;  //Keep count of the row interchanges
	for(k=0; k<n-1; k++)
	{
		amax = fabs(a[k][k]);
		/*--------------Row pivoting--------------*/
		m = k;
		for(i=k+1; i<n; i++)
		{
			//Find the row with largest pivot
			xfac = fabs(a[i][k]);
			if(xfac > amax) {amax = xfac; m=i;}
		}
		if(m != k)
		{
			//Row interchanges i <-> k
			rowx = rowx+1;
			temp = b[k];
			b[k]  = b[m];
			b[m]  = temp;
			for(j=k; j<n; j++)
			{
				temp = a[k][j];
				a[k][j] = a[m][j];
				a[m][j] = temp;
			}
		}
		/*--------------End row pivoting--------------*/

		for(i=k+1; i<n; i++)
		{
			xfac = a[i][k]/a[k][k];
			for(j=k; j<n; j++)
			{
				a[i][j] = a[i][j] - xfac*a[k][j];
			}
			b[i] = b[i] - xfac*b[k];
		}
		/*
		printf("\n A after decomposition step %d\n\n",k);
			     matrix_print_off(n, n, a); */
	}

	/*--------2. Do the back substitution step.--------*/
	for(j=0; j<n; j++)
	{
		  k   = n-j-1;
		 x[k] = b[k];
		for(i=k+1; i<n; i++)
			x[k] = x[k] - a[k][i]*x[i];

		x[k] = x[k]/a[k][k];
	}
	//printf("\nNumber of row exchanges = %d\n",rowx);
}

int matrix_print_off (int nr, int nc, double **A)
{
	int i,j;

	if(nr<= 0) return (-1);
	if(nc<= 0) return (-2);

	for(i=0; i <nr; i++)
	{
		for(j=0; j<nc; j++)
			printf ("%9.4f  ", A[i][j]);

		printf("\n");
	}
	return (0);
}
