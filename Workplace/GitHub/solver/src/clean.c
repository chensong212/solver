/*
 *
 *  Created on: 2016.04.18
 *  Purpose: End the job, free memory.
 *
 */

#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include"comm.h"

/*---------------------------------------------------
 * End the simulation, clean all the space
 * ------------------------------------------------*/
void clean()
{
	int nc, nc1;

	void freeU(int nlen, struct strct_U *U);
	void freeOthers();

	nc = config1.ni*config1.nj;
	nc1 = I0*J0;

	freeU(nc,  &U);
	freeU(nc1, &Ug);

	freeOthers();

	if(MyID==0)
	{
		free(mesh.x);
		free(mesh.y);
		printf("program exits! \n");
	}

#ifdef MPI_RUN
	MPI_Finalize();
#endif
	exit(0);
}

/*---------------------------------------------------
 * free memory of U vector
 * ------------------------------------------------*/
void freeU(int nlen, struct strct_U *U)
{
	int i;

	free(U->mu);
	free(U->kt);
	free(U->cv);
  free(U->rgas);
	free(U->pre);
	free(U->tem);
	free(U->gam);

	for(i=0; i<nlen; i++)
		free(U->q[i]);

	free(U->q);

	if(config1.gasModel != 0)
	{
		for(i=0; i<nlen; i++)
		{
			free(U->qs[i]);
			free(U->di[i]);
		}
	}
	free(U->qs);
	free(U->di);

}

/*---------------------------------------------------
 * free memory of other variables
 * ------------------------------------------------*/
void freeOthers()
{
	int ic, ns, nc;

	nc = config1.ni*config1.nj;

	free(mesh.xi);
	free(mesh.et);
	free(mesh.x_xi);
	free(mesh.x_et);
	free(mesh.y_xi);
	free(mesh.y_et);
	free(mesh.yaks);

	for(ic=0; ic<nc; ic++)
	{
		free(qo[ic]);
		free(rhs[ic]);
	}
	free(qo);
	free(rhs);
	if(config1.gasModel != 0)
	{
		for(ic=0; ic<nc; ic++)
			free(qso[ic]);

		for(ic=0; ic<nc; ic++)
		{
			for(ns=0; ns<config1.nspec; ns++)
			{
				free(dsdq[ic][ns]);
			}
			free(dsdq[ic]);
		}
		free(dsdq);
	}
	free(qso);

}
