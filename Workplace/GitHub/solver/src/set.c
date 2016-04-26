/*
 *
 *  Created on: 2016.04.18
 *  Purpose: Conduct configuration and allocate memory.
 *
 */
#include<math.h>
#include<mpi.h>
#include"comm.h"
#include"chemdata.h"

/*---------------------------------------------------
 * set and initialize the simulation
 * ------------------------------------------------*/
void set()
{
	int nc, nc1;

	void init();
	void import();
	void setGeom();

	void allocateU(int nlen, struct strct_U *U);
	void allocateUv();
	void allocateOthers();

	if(config1.gasModel == 0)
		neqn = neqv;
	else
		neqn = neqv + config1.nspec -1;

	nc  = config1.ni*config1.nj;
	nc1 = I0*J0;

	allocateU(nc, &U);
	allocateU(nc1, &Ug);
	allocateUv();
	allocateOthers();

	if(config1.newrun == 1)
		init();
	else
		import();

	setGeom();
	

}

/*---------------------------------------------------
 * initialization
 * ------------------------------------------------*/
void init()
{
	int i, mni, i0, ir, ic, icp;
	double dis;
	void assigncells(int i1, int in, int j1, int jn, double u,
			             double v, double t, double p, double *qs);

	config2.t0 = 0;
	config1.iStep0 = 1;

	if(MyID == 0)
	{
		/* Note: only MyID==0 carries the mesh data */
		mni = nproc*config1.ni;
		dis = config2.x0*config2.Lx;
		for(i=0; i<mni; i++)
		{
			ic = i*config1.nj + 0;
			icp = (i+1)*config1.nj + 0;
			if((mesh.x[ic]<=dis) && (dis<mesh.x[icp]))
				break;
		}
	    config1.x_sh = i; 
	}

#ifdef MPI_RUN
	MPI_Bcast(&config1.x_sh, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

	// get the position of the diaphragm
    i0 = config1.x_sh/config1.ni; // the integer part
    ir = config1.x_sh%config1.ni; // the remainder part

#ifdef MPI_RUN

	if(MyID < i0)
		assigncells(0,config1.ni, 0,config1.nj, inc[0].u,inc[0].v,inc[0].t,inc[0].p, inc[0].ys);
	else if((MyID == i0) && (ir != 0))
	{
		assigncells(0,ir, 0,config1.nj, inc[0].u,inc[0].v,inc[0].t,inc[0].p, inc[0].ys);
		assigncells(ir,config1.ni, 0,config1.nj, inc[1].u,inc[1].v,inc[1].t,inc[1].p, inc[1].ys);
	} 
	else
		assigncells(0,config1.ni, 0,config1.nj, inc[1].u,inc[1].v,inc[1].t,inc[1].p, inc[1].ys);

	MPI_Barrier(MPI_COMM_WORLD);
#else

	assigncells(i0,ir, 0,config1.nj, inc[0].u,inc[0].v,inc[0].t,inc[0].p, inc[0].ys);
	assigncells(ir,config1.ni, 0,config1.nj, inc[1].u,inc[1].v,inc[1].t,inc[1].p,inc[1].ys);

#endif
}

/*-----------------------------------------------------------
 * Assign initial value to each cells
 * ---------------------------------------------------------*/
void assigncells(int i1, int in, int j1, int jn, double u,
		             double v, double t, double p, double *qs)
{
	int i, j, ic, ns;
	double ein, es, ek, rgas1, RT;
	double getenergy(double qs[], double t);
	double getes(int ns, double t);
	double getrgas(double qs[]);

	rgas1 = (ru/config2.molWeight);

	for(i=i1; i<in; i++)
		for(j=j1; j<jn; j++)
		{
			ic = i*config1.nj + j;

			U.q[ic][1] = u;
			U.q[ic][2] = v;
			U.tem[ic]  = t;
			U.pre[ic]  = p;

			ek = 0.5*(u*u + v*v);

			if(config1.gasModel == 0)
			{
				RT = rgas1*t;
				U.q[ic][0] = p/RT;
				U.q[ic][3] = ek + RT/(config2.gam0-1);
			}
			else
			{
		    // to be added...
      }
		}
}

/*-----------------------------------------------------------
 * Import the former solution.
 * ---------------------------------------------------------*/
void import()
{
	int i, j, ic, ic1, id, ns, mni, nc;
	double x, y, rho, u, v, p, T, e, dum;
	char linebuf[200], varname[200], filename[30];
	struct strct_field
	{
		double *rho, *u, *v, *p, *t, *e, *ga, **qs;
	} Uf;
	FILE *fp;
	void clean();

	mni = nproc*config1.ni;
	nc  = mni*config1.nj;

	if(MyID == 0)
	{
		Uf.rho = (double*)malloc(sizeof(double)*nc);
		Uf.u   = (double*)malloc(sizeof(double)*nc);
		Uf.v   = (double*)malloc(sizeof(double)*nc);
		Uf.p   = (double*)malloc(sizeof(double)*nc);
		Uf.t   = (double*)malloc(sizeof(double)*nc);
		Uf.e   = (double*)malloc(sizeof(double)*nc);
		if(config1.gasModel != 0)
		{
			Uf.ga   = (double*)malloc(sizeof(double)*nc);
			Uf.qs = (double**)malloc(sizeof(double*)*nc);
			for(ic=0; ic<nc; ic++)
				Uf.qs[ic] = (double*)malloc(sizeof(double)*config1.nspec);
		}
		else
			Uf.qs = NULL;

		sprintf(filename, "nstep%d_field.dat", config1.iStep0 -1);
		fp = fopen(filename,"r");
		if(fp == NULL)
		{
			printf("%s not found! \n", filename);

#ifdef MPI_RUN
			MPI_Abort( MPI_COMM_WORLD, 99);
#else
			clean();
#endif

    }
		printf("reading the field file... \n");
		fgets(linebuf,sizeof(linebuf),fp);
		fgets(varname,sizeof(linebuf),fp);
		fgets(linebuf,sizeof(linebuf),fp);

		/*----1. read the field file----*/
		for(j=0; j<config1.nj; j++)
			for(i=0; i<mni; i++)
			{
				ic = i*config1.nj + j;
				fscanf(fp," %lf %lf %lf %lf %lf %lf %lf %lf",&dum, &dum,
						&Uf.rho[ic], &Uf.u[ic], &Uf.v[ic], &Uf.p[ic], &Uf.t[ic], &Uf.e[ic]);
				if(config1.gasModel != 0)
				{
					fscanf(fp," %lf",&Uf.ga[ic]);
					for(ns=0; ns<config1.nspec; ns++)
						fscanf(fp," %le",&Uf.qs[ic][ns]);
				}
			}
		  fclose(fp);

		/*----2. write the new tcv.dat file for each processors----*/
		for(id=0; id<nproc; id++)
		{
			sprintf(filename, "tcv%d.dat", id);
			fp = fopen(filename,"w");
			fprintf(fp, "Title = \"Flow field of local processor\"\n");
			fprintf(fp, "%s", varname);
			fprintf(fp,"ZONE T='%d', I= %d, J= %d, f=point \n", MyID, config1.ni, config1.nj);

			for(j=0; j<config1.nj; j++) // without ghost cells
				for(i = 0; i<config1.ni; i++)
				{
					ic = (id*config1.ni + i)*config1.nj + j;
					ic1 = (i+config1.Ng)*J0 + j+config1.Ng;

					x   = mesh.xi[ic1];
					y   = mesh.et[ic1];
					fprintf(fp," %lf %lf %lf %lf %lf %lf %lf %lf",x, y,
							Uf.rho[ic],Uf.u[ic],Uf.v[ic],Uf.p[ic],Uf.t[ic],Uf.e[ic]);
					if(config1.gasModel != 0)
					{
					  // to be added...
          }
					fprintf(fp, "\n");
				}
			fclose(fp);
		}
	}

#ifdef MPI_RUN
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	/*----3. Each processor read the solution file----*/
  sprintf(filename, "tcv%d.dat", MyID);
  fp = fopen(filename,"r");
  if(fp == NULL)
  {
    printf("plot file tcv%d.dat not found! \n", MyID);
    clean();
  }
    
  fgets(linebuf, sizeof (linebuf), fp);
	fgets(linebuf, sizeof (linebuf), fp);
	fgets(linebuf, sizeof (linebuf), fp); // skip the title

	for(j=0; j<config1.nj; j++)
		for(i=0; i<config1.ni; i++)
		{
			ic = i*config1.nj + j;
			if(fscanf(fp," %lf %lf %lf %lf %lf %lf %lf %lf",&dum, &dum, &rho, &u, &v, &p, &T, &e) != 8)
			{
				printf("format error in tcv.dat \n");

#ifdef MPI_RUN
				MPI_Abort( MPI_COMM_WORLD, 99);
#else
				clean();
#endif

			}
			if(config1.gasModel != 0)
			{
			  // to be added...
      }

			U.q[ic][0] = rho;
			U.q[ic][1] = u;
			U.q[ic][2] = v;
			U.q[ic][3] = e;
			U.pre[ic]  = p;
			U.tem[ic]  = T;
		}
    fclose(fp);

	if(MyID == 0)
	{
		free(Uf.rho);
		free(Uf.u);
		free(Uf.v);
		free(Uf.p);
		free(Uf.t);
		free(Uf.e);
		if(config1.gasModel != 0)
		{
			for(ic=0; ic<nc; ic++)
				free(Uf.qs[ic]);
			free(Uf.qs);
		}
	}
}

/*---------------------------------------------------
 * Set the geometry related variables
 * ------------------------------------------------*/
void setGeom()
{

  //to be added...

}

/*---------------------------------------------------
 * allocate memory for U vector
 * ------------------------------------------------*/
void allocateU(int nlen, struct strct_U *U)
{
	int i;
	U->mu   = (double*)malloc(sizeof(double)*nlen);
	U->kt   = (double*)malloc(sizeof(double)*nlen);
	U->cv   = (double*)malloc(sizeof(double)*nlen);
	U->rgas = (double*)malloc(sizeof(double)*nlen);
	U->pre  = (double*)malloc(sizeof(double)*nlen);
	U->tem  = (double*)malloc(sizeof(double)*nlen);
	U->gam  = (double*)malloc(sizeof(double)*nlen);
	U->q    = (double**)malloc(sizeof(double*)*nlen);

	for(i=0; i<nlen; i++)
		U->q[i]  = (double*)malloc(sizeof(double)*neqv);

	if(config1.gasModel == 0)
	{
		U->di = NULL;
		U->qs = NULL;
	}
	else
	{
		U->qs  = (double**)malloc(sizeof(double*)*nlen);
		U->di  = (double**)malloc(sizeof(double*)*nlen);
		for(i=0; i<nlen; i++)
		{
			U->qs[i] = (double*)malloc(sizeof(double)*config1.nspec);
			U->di[i] = (double*)malloc(sizeof(double)*config1.nspec);
		}
	}

}

/*---------------------------------------------------
 * allocate memory for derivatives
 * ------------------------------------------------*/
void allocateUv()
{

  // to be added...

}

/*---------------------------------------------------
 * allocate memory for other variables
 * ------------------------------------------------*/
void allocateOthers()
{
	int ic, ns, nc;

	nc = config1.ni*config1.nj;

	/*------------Other variables-------------*/
	rhs  = (double**)malloc(sizeof(double*)*nc);
	qo   = (double**)malloc(sizeof(double*)*nc);
	for(ic=0; ic<nc; ic++)
	{
		qo[ic]   = (double*)malloc(sizeof(double)*neqv);
		rhs[ic]  = (double*)malloc(sizeof(double)*neqn);
	}
	if(config1.gasModel == 0)
		qso  = NULL;
	else
	{
		qso  = (double**)malloc(sizeof(double*)*nc);
		for(ic=0; ic<nc; ic++)
			qso[ic]  = (double*)malloc(sizeof(double)*config1.nspec);

		dsdq = (double***)malloc(sizeof(double**)*nc);
		for(ic=0; ic<nc; ic++)
		{
			dsdq[ic] = (double**)malloc(sizeof(double*)*config1.nspec);
			for(ns=0; ns<config1.nspec; ns++)
				dsdq[ic][ns] = (double*)malloc(sizeof(double)*neqn);
		}
	}
}
