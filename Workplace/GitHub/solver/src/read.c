/*
 *
 *  Created on: 2016.04.18.
 *  Purpose: read the problem information.
 *
 */
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include"comm.h"
#include"chemdata.h"

/*---------------------------------------------------
 * Main steps of reading.
 * ------------------------------------------------*/
void read()
{
	void cheminput();
	void readconfig();
	void readic();
	void readmesh();
	void BcastData();

	neqv = 4; // solution variables without chemical terms

	if(MyID == 0)
	{
		readconfig();
		if(config1.gasModel != 0)
			cheminput();
	}

#ifdef MPI_RUN
	BcastData();
#endif

	readic();
	readmesh();
}

/*---------------------------------------------------
 * read configuration information
 * ------------------------------------------------*/
void readconfig()
{
	FILE *fp, *outId;

	fp    = fopen("config.dat", "r");
	outId = fopen("outInfo.dat", "a");
	if(fp == NULL){printf("config.dat file not found! \n");exit(0);}

	if(fscanf(fp, "%lf %lf", &config2.t0, &config2.x0) != 2)
			{printf("format error in config.dat\n");exit(0);}
	if(fscanf(fp, "%d %d %d %d %d %d %d", &config1.newrun, &config1.useDt, &config1.iStep0,
                   &config1.nStep, &config1.nRamp, &config1.Samples, &config1.ifilm) != 7)
			{printf("format error in config.dat\n");exit(0);}
	if(fscanf(fp, "%lf %lf %lf %lf",  &config2.dt0, &config2.dt1, &config2.CFL0, &config2.CFL1) != 4)
			{printf("format error in config.dat\n");exit(0);}
	if(fscanf(fp, "%d %d %d", &config1.gasModel, &config1.reacModel, &config1.visModel)!= 3)
			{printf("format error in config.dat\n");exit(0);}
	if(fscanf(fp, "%lf %lf %lf", &config2.molWeight, &config2.gam0, &config2.Pr0) != 3)
			{printf("format error in config.dat\n");exit(0);}
	if(fscanf(fp, "%lf %lf %lf %lf", &config2.p1, &config2.T1, &config2.u1, &config2.v1) != 4)
			{printf("format error in config.dat\n");exit(0);}
	if(fscanf(fp, "%lf %lf %lf %lf", &config2.p2, &config2.T2, &config2.u2, &config2.v2) != 4)
			{printf("format error in config.dat\n");exit(0);}
	fclose(fp);

	fp = fopen("gridset.dat", "r");
	if(fp == NULL){printf("gridset.dat file not found \n");exit(0);}
	if(fscanf(fp, "%d %d %d %d", &config1.ni, &config1.nj,  &config1.Ng, &config1.nblock) != 4)
			 {printf("format error in gridset.dat\n");exit(0);}
	if(fscanf(fp, "%lf %lf", &config2.Lx, &config2.Ly) != 2)
			 {printf("format error in gridset.dat\n");exit(0);}
	fclose(fp);

	fprintf(outId, "\n/---------------------------configure data---------------------------/\n");
	fprintf(outId, "\nt0=%lf, x0=%lf, Lx=%lf, Ly=%lf\n", config2.t0, config2.x0, config2.Lx, config2.Ly);
	fprintf(outId, "\ngrids point: ni=%d, nj=%d, ghost cells Ng=%d, nblock=%d\n", config1.ni, config1.nj, config1.Ng, config1.nblock);
	fprintf(outId, "\nnewrun=%d, useDt=%d, \niStep0=%d, nStep=%d, nRamp=%d, Samples=%d, ifilm=%d \n",
			            config1.newrun, config1.useDt, config1.iStep0, config1.nStep,
			            config1.nRamp, config1.Samples, config1.ifilm);
	fprintf(outId, "dt0=%lf, dt1=%lf, CFL0=%lf, CFL1=%lf, \n", config2.dt0, config2.dt1, config2.CFL0, config2.CFL1);
	fprintf(outId, "\ngasModel=%d, reacModel=%d, visModel=%d \n", config1.gasModel, config1.reacModel, config1.visModel);
	fprintf(outId, "molWeight=%lf, gam0=%lf, Pr0=%lf\n", config2.molWeight, config2.gam0, config2.Pr0);

	fprintf(outId, "\nInitial condition: \n");
	fprintf(outId, "p_1=%lf, T_1=%lf, u_1=%lf, v1=%lf \n",config2.p1, config2.T1, config2.u1, config2.v1);
	fprintf(outId, "p_2=%lf, T_2=%lf, u_2=%lf, v2=%lf \n",config2.p2, config2.T2, config2.u2, config2.v2);

	config1.thermo_base = 0;
	config1.timeOrder = 3;

	fclose(outId);
}

/*---------------------------------------------------
 * read initial condition
 * ------------------------------------------------*/
void readic()
{
	int ns, ib, nb;

	nb = 2; // No. of initial condition

	inc[0].p = config2.p1;
	inc[0].t = config2.T1;
	inc[0].u = config2.u1;
	inc[0].v = config2.v1;
	inc[1].p = config2.p2;
	inc[1].t = config2.T2;
	inc[1].u = config2.u2;
	inc[1].v = config2.v2;
	if(config1.gasModel !=0)
	{
	
    // read the species data, to be added...

  }
}

/*---------------------------------------------------
 * read mesh data
 * ------------------------------------------------*/
void readmesh()
{
	int i, j, ic,nc, mni;
	char linebuf[200], filename[20];
	FILE  *fp, *outId;

	/*---1. Grid information---*/
	I0 = config1.ni + 2*config1.Ng;
	J0 = config1.nj + 2*config1.Ng;
	nc = I0*J0;

	/*---2. Allocate memory---*/
	mesh.xi   = (double*)malloc(sizeof(double)*nc);
	mesh.et   = (double*)malloc(sizeof(double)*nc);
	mesh.x_xi = (double*)malloc(sizeof(double)*nc);
	mesh.x_et = (double*)malloc(sizeof(double)*nc);
	mesh.y_xi = (double*)malloc(sizeof(double)*nc);
	mesh.y_et = (double*)malloc(sizeof(double)*nc);
	mesh.yaks = (double*)malloc(sizeof(double)*nc);

	sprintf(filename, "set%d.dat", MyID);
	fp = fopen(filename,"r");
	if(fp == NULL)
	{
		printf("grid file set%d.dat not found! \n", MyID);

#ifdef MPI_RUN
		MPI_Abort( MPI_COMM_WORLD, 99);
#else
		exit(0);
#endif

	}

	fgets(linebuf, sizeof(linebuf), fp);
	fgets(linebuf, sizeof(linebuf), fp);
	fgets(linebuf, sizeof(linebuf), fp); // skip the title

	/*---3. Read the grids---*/
	for(j=0; j<J0; j++)
		for(i=0; i<I0; i++)
		{
			ic = i*J0 + j;
			if(fscanf(fp," %lf %lf %lf %lf %lf %lf %lf",
				  &mesh.xi[ic], &mesh.et[ic], &mesh.x_xi[ic], &mesh.x_et[ic],
				  &mesh.y_xi[ic], &mesh.y_et[ic], &mesh.yaks[ic]) != 7)
			{
				printf("format error in grid file set%d.dat! \n", MyID);

#ifdef MPI_RUN
        MPI_Abort( MPI_COMM_WORLD, 99);
#else
			  exit(0);
#endif
			}
	 }
	fclose(fp);
	
	/* After coordinate transformation,
	 * the distance between cells are equal */
	dxc = mesh.xi[1*J0] - mesh.xi[0*J0];
	dyc = mesh.et[1]    - mesh.et[0];

	/*---4. Read the physical mesh---*/

	if(MyID == 0)
	{
		outId = fopen("outInfo.dat", "a");

	  mni = nproc*config1.ni;
		nc = config1.nj*mni;
		mesh.x = (double*)malloc(sizeof(double)*nc);
		mesh.y = (double*)malloc(sizeof(double)*nc);

		fp = fopen("mesh.dat","r");
		if(fp == NULL)
		{
			printf("mesh file mesh.dat not found! \n");

#ifdef MPI_RUN
			MPI_Abort( MPI_COMM_WORLD, 99);
#else
			exit(0);
#endif

		}

		fgets(linebuf, sizeof (linebuf), fp);
		fgets(linebuf, sizeof (linebuf), fp);
		fgets(linebuf, sizeof (linebuf), fp); // skip the title
		
		for(j=0; j<config1.nj; j++)
			for(i=0; i<mni; i++)
			{
		    	ic = i*config1.nj + j;
		    	if(fscanf(fp," %lf  %lf",&mesh.x[ic], &mesh.y[ic])!= 2)
		    	{
					printf("format error in mesh file! \n");
		
#ifdef MPI_RUN
			MPI_Abort( MPI_COMM_WORLD, 99);
#else
			exit(0);
#endif

		    	};
			};

	    fclose(fp);
      fprintf(outId,"\nRead grid data complete!\n");
		  fprintf(outId,"\n/---------------------runtime information---------------------/\n");
	  	fclose(outId);
	}
}

/*---------------------------------------------------
 * read chemical and thermal data
 * ------------------------------------------------*/
void cheminput()
{

  // to be added...

}	

#ifdef MPI_RUN
/*---------------------------------------------------
 * tell other tasks the thermal-chemical data
 * ------------------------------------------------*/
void BcastData()
{
	int count1, count2, block[2];

	MPI_Aint offset[2], extent;

	/*-----------Broadcast configuration data-----------*/
	count1 = sizeof(config1)/sizeof(int);
	count2 = sizeof(config2)/sizeof(double);
	MPI_Bcast(&config1, count1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&config2, count2, MPI_DOUBLE, 0, MPI_COMM_WORLD);


	/*-----------Broadcast thermo-chemical data-----------*/
  // to be added...
}

#endif
