/*
 *
 *  Created on: 2016.04.18
 *  Purpose: Calculate thermal and chemical properties.
 *
 */
#include<string.h>
#include<math.h>
#include<mpi.h>
#include"comm.h"
#include"chemdata.h"

/*---------------------------------------------------
 * Calculate the thermal properties
 * ------------------------------------------------*/
void gettherm(int nc, double **q, double **qs, double *p,
		          double *t, double *gam1, double *rgas1, double *cv1)
{
	int ic, ns, nr, nd, ntol, icount;
	double hs, tol, told, tnew, energy, dum, tdiff,
		     tem_min, tem_max, temp, temp2, temp3, temp4,
		     cp, cps, u, v, etot, tem, rgas0;
  double c2 = 0.5, c3 = 1./3., c4 = 0.25, c5 = 0.2;
  void clean();

    /*-- Set a range for the temperature. Terminate the program
     * if un-physical conditions are encountered --*/
	 
  tem_min = 60.0;
	tem_max = 20000.0;

	if(config1.gasModel == 0)
	{
		/*-- for calorically perfect gas--*/
		rgas0 = (ru/config2.molWeight);
		for(ic=0; ic<nc; ic++)
		{
			temp = 0.5*(q[ic][1]*q[ic][1] + q[ic][2]*q[ic][2]);
			rgas1[ic] = rgas0;
			gam1[ic] = config2.gam0;
			cv1[ic]  = rgas0/(config2.gam0 - 1.);
			p[ic] = (gam1[ic] - 1)*(q[ic][3] - temp)*q[ic][0];
			t[ic] = p[ic]/q[ic][0]/rgas1[ic];
			tem = t[ic];
			if( !(tem_min<tem && tem<tem_max))
			{
				printf("T = %f (K) at at i=%d, j=%d\n", tem, ic/config1.nj, ic%config1.nj);
				printf("The temperature may be unphysical for ideal gas model \n");
#ifdef MPI_RUN
				MPI_Abort( MPI_COMM_WORLD, 99);
#else
				clean();
#endif
			}

		}
	}
	else
	{
		/*-- for none calorically perfect gas--*/
  
    // to be added...
  }
}

/*--------------------------------------------------------------
 * Calculate transport properties(viscosity, conductivity, ...
 * -------------------------------------------------------------*/
void gettrans(int nc, double **qs, double *t, double *gam1,
		          double *cv1, double *mu, double *cond, double **diff)
{
  void clean();

	if(config1.gasModel == 0)
	{
		/*---------- 1. ideal gas flow----------*/
		if(config1.visModel == 0)
			return;
		else if(config1.visModel == 1)
		{
      // constan mu and cond
    }
		else if(config1.visModel == 2)
		{
			// Sutherland's law
		}
		else
		{
			printf("unsupported viscosity model... \n");

#ifdef MPI_RUN
			MPI_Abort( MPI_COMM_WORLD, 99);
#else
			clean();
#endif

		}
	}
	else
	{
		/*---------- 2. real gas flow----------*/
	}
}

void chemsource(int nc, double **q, double **qs, double *tem, double **rhs)
{
  // to be added...
}



