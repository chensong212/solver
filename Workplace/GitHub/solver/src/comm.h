/*
 *
 *  Created on: 2016.04.18
 *  Purpose: set grobal variables for the solver
 *
 */

#ifndef COMM_H_
#define COMM_H_

#include<stdio.h>
#include<stdlib.h>

/*- constant & func. -*/
#define maxeqn 9
#define mpicell 3000 // Number of exchange cells for MPI
#define ru 8.3141
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

#define MPI_RUN

int I0, J0, neqv, neqn, nproc, MyID;
double dtc, dxc, dyc, Ttot;

struct strct_configInt
{
	int newrun, useDt, ni, nj, Ng, nblock, iStep0, nStep, timeOrder, 
		  Samples, nRamp, ifilm, gasModel, visModel, reacModel,
	    x_sh, nmix, nspec, nreac, thermo_base;
} config1;
struct strct_configDouble
{
	double t0, x0, Lx, Ly, dt0, dt1, CFL0, CFL1, molWeight, gam0,
		     Pr0, Sc0, Re0, muRef, suthC1, suthC2, condRef, diffRef, MaRef, temRef, preRef,
		     p1, T1, u1, v1, p2, T2, u2, v2;
} config2;


/*- MPI related -*/
#ifdef MPI_RUN
char processor_name[256];
int dest, ierr, nachbar_rechts, nachbar_links, namelen, NMAXproc;
struct strct_gcelltype
{
	double rho, u, v, e, p, t, mu, kt, ga, qs[6], di[6];
}mpiSend_ql[mpicell], mpiSend_qr[mpicell], mpiRecv_ql[mpicell], mpiRecv_qr[mpicell];
#endif

/*- Origin flow variables -*/
struct strct_U
{
	double **q, **qs, *pre, *tem, *mu, *kt, *cv, *rgas, *gam, **di;
} U, Ug;

/*- viscous related variables -*/
struct strct_Uv
{
	// origin variables derivatives
	double *u_xi, *u_et, *v_xi, *v_et, *T_xi, *T_et, **qs_xi, **qs_et;
} Uv;

/*- flux related variables -*/
struct strct_flux
{
	double *xix, *xiy, *etx, *ety, *yas, *rho, *du, *dv, *dt, 
         **dqs, *u, *v, *p, *e, *t, *gam, *mu, *kt, **Ds, **qs, **flux;
} U1d;

/*- initial condition variables -*/
struct strct_ic
{
	double u, v, t, p, ys[6];
} inc[2];

/*- geometry variables -*/
struct strct_metric
{
	double *x, *y, *xi, *et, *x_xi, *x_et, *y_xi, *y_et, *yaks, *qbound;
} mesh;

/*- others -*/
double **qo, **qso, **rhs, ***dsdq;

#endif /* COMM_H_ */

