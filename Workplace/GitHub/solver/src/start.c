/*
 *
 *  Created on: 2016.04.18
 *  Purpose: Start simulation, Write output information.
 *
 */

#include<string.h>
#include<mpi.h>
#include"comm.h"


void start(int argc, char *argv[])
{

	FILE *outId;

#ifdef MPI_RUN

  MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &MyID);
	MPI_Get_processor_name(processor_name, &namelen);

#else
	MyID  = 0;
	nproc = 1;
#endif

	/*------write the ouput file-------*/
	if(MyID == 0)
	{
		printf("simulation start...\n");
		outId = fopen("outInfo.dat", "w");
		fprintf(outId, "/***************************** simulation output *****************************/ \n");
		fprintf(outId, " \n");
		fclose(outId);
	}
}

