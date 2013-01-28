/***************************************************************************
 *                                                                         *
 *    Simulation of a simple bacterial growth model                        *
 *    Algorithm: Finite Differences method                                 *
 *                                                                         *
 *            Jordi Garcia Ojalvo & Pau Rue, 2012/02/01                    *
 *                                                                         *
 ***************************************************************************/

#include "sptk.h"
double *xrn;

int main(int argc, char *argv[])
{
  double aux;
  int i,j,i2,j2,i3;
  int ij,it,seed,L,N,Lh,Mh,nst,nme,ict,npat;
  /* Integration parameters */
  double tt,t,st,dx,dt,dth,r2;
  double cx,cy;
  /* Initial conditions */
  double umin,umax,rcol,rring;
  int *n1,*n2,*n3,*n4;
  /* Model parameters */
  double agr0,agrcl,agrcp;
  double x0;
  double xrn_avg, xrn2_avg, xrn_std;  
  double am,km,lbd,nh,ikm,akm,gd,kd,ph;
  double aw,lbdw;
  double as,ks,lbds,ns,iks,aks;
  double *Du,fD,fD1,diffpref;
  double epsK,Kmax0, Kmin, lbdK, kk, nk,ikk; 
  double Kmax,Kmaxi, u_avg,ui_avg, reac_Kmax, reac_Kmaxi;
  /* Integrator variables */
  double K, Ki;
  double *u,*ui,*agr,*mtx,*mtxi,*w,*wi,*s,*si;
  double *reac_u,*reac_ui;
  double *reac_m,*reac_mi;
  double *reac_w,*reac_wi;
  double *reac_s,*reac_si;
  double *diff_u,*diff_ui;
  double *diff_u2,*diff_u2i;
  double mn,mx;
  /* GSL rng stuff */
  const gsl_rng_type * type = gsl_rng_default;
  gsl_rng * rng;
  /* File handling */
  char name1[60],name3[60];
  FILE *input,*outdat;
  /* Fourier space */
  double ki,kj,dk,*k2vec;
  /* FFTW stuff */
  fftw_plan plan_u1f,plan_u1b;
  fftw_complex *xrn_ft;
  double complex *Cu1;
  int flag_save_u,flag_save_m,flag_save_w,flag_save_s;

/**********  INPUT DATA  **********/
  strcpy(name1,argv[1]);

  sprintf(name3,"parfiles/%s.par",name1);
  input = fopen(name3,"r");
  fscanf(input,"%lg",&agr0);	  next_line(input);
  fscanf(input,"%lg",&agrcl);	  next_line(input);
  fscanf(input,"%lg",&agrcp);	  next_line(input);
  fscanf(input,"%lg",&Kmax0);     next_line(input);
  fscanf(input,"%lg",&Kmin);     next_line(input);
  fscanf(input,"%lg",&lbdK);     next_line(input);
  fscanf(input,"%lg",&epsK);     next_line(input);
  fscanf(input,"%lg",&kk);     next_line(input);
  fscanf(input,"%lg",&nk);     next_line(input);
  fscanf(input,"%lg",&am);	  next_line(input);
  fscanf(input,"%lg",&km);     next_line(input);
  fscanf(input,"%lg",&lbd);    next_line(input);
  fscanf(input,"%lg",&nh);	  next_line(input);
  fscanf(input,"%lg",&aw);	  next_line(input);
  fscanf(input,"%lg",&lbdw);    next_line(input);
  fscanf(input,"%lg",&as);	  next_line(input);
  fscanf(input,"%lg",&ks);     next_line(input);
  fscanf(input,"%lg",&lbds);    next_line(input);
  fscanf(input,"%lg",&ns);	  next_line(input);
  fscanf(input,"%lg",&gd);	  next_line(input);
  fscanf(input,"%lg",&kd);	  next_line(input);
  fscanf(input,"%lg",&ph);	  next_line(input);
  fscanf(input,"%lg",&umin);  next_line(input);
  fscanf(input,"%lg",&umax);  next_line(input);
  fscanf(input,"%i", &ict);   next_line(input);
  fscanf(input,"%lg",&rcol);  next_line(input);
  fscanf(input,"%lg",&rring);  next_line(input);
  fscanf(input,"%i",&L);      next_line(input);
  fscanf(input,"%lg",&tt);    next_line(input);
  fscanf(input,"%lg",&st);    next_line(input);
  fscanf(input,"%lg",&dx);    next_line(input);
  fscanf(input,"%lg",&dt);    next_line(input);
  fscanf(input,"%i",&seed);   next_line(input);
  fscanf(input,"%i",&flag_save_u);   next_line(input);
  fscanf(input,"%i",&flag_save_m);   next_line(input);
  fscanf(input,"%i",&flag_save_w);   next_line(input);
  fscanf(input,"%i",&flag_save_s);   next_line(input);
  fclose(input); 

  sprintf(name3,"output/%s_gr.prof",name1);
  outdat = fopen(name3,"w");
  fwrite(&L,sizeof(int),1,outdat);
  fclose(outdat);

  if (flag_save_u) {
    sprintf(name3,"output/%s_u.pat",name1);
    outdat = fopen(name3,"w");
    fwrite(&L,sizeof(int),1,outdat);
    npat = (int)(tt/st+0.5);
    fwrite(&npat,sizeof(int),1,outdat);
    fclose(outdat);
  }
  if (flag_save_m) {
    sprintf(name3,"output/%s_m.pat",name1);
    outdat = fopen(name3,"w");
    fwrite(&L,sizeof(int),1,outdat);
    npat = (int)(tt/st+0.5);
    fwrite(&npat,sizeof(int),1,outdat);
    fclose(outdat);
  }
  if (flag_save_w) {
    sprintf(name3,"output/%s_w.pat",name1);
    outdat = fopen(name3,"w");
    fwrite(&L,sizeof(int),1,outdat);
    npat = (int)(tt/st+0.5);
    fwrite(&npat,sizeof(int),1,outdat);
    fclose(outdat);
  }
  if (flag_save_s) {  
    sprintf(name3,"output/%s_s.pat",name1);
    outdat = fopen(name3,"w");
    fwrite(&L,sizeof(int),1,outdat);
    npat = (int)(tt/st+0.5);
    fwrite(&npat,sizeof(int),1,outdat);
    fclose(outdat);
  }
  
/**********  SIMULATION CONSTANTS  **********/
  N = L*L;
  Lh = L/2;
  Mh = L*(Lh+1);
  nst = (int)(st/dt+0.5);
  nme = (int)(tt/dt+0.5);
  dth = dt*0.5;
  fD = 1/(dx*dx);
  fD1 = 1/(4*dx*dx);
  dk=2.0*M_PI/((double)L * dx);
  akm = am/pow(km,nh);
  ikm = 1/pow(km,nh);
  aks = as/pow(ks,ns);
  iks = 1/pow(ks,ns);
  ikk = 1/pow(kk,nk);
  mn = 0;
  mx = 0;
  
  /* Memory allocation */
  u = (double *)calloc(N,sizeof(double));
  ui = (double *)calloc(N,sizeof(double));
  reac_u = (double *)calloc(N,sizeof(double));
  reac_ui = (double *)calloc(N,sizeof(double));
  mtx = (double *)calloc(N,sizeof(double));
  mtxi = (double *)calloc(N,sizeof(double));
  reac_m = (double *)calloc(N,sizeof(double));
  reac_mi = (double *)calloc(N,sizeof(double));
  w = (double *)calloc(N,sizeof(double));
  wi = (double *)calloc(N,sizeof(double));
  reac_w = (double *)calloc(N,sizeof(double));
  reac_wi = (double *)calloc(N,sizeof(double));

  s = (double *)calloc(N,sizeof(double));
  si = (double *)calloc(N,sizeof(double));
  reac_s = (double *)calloc(N,sizeof(double));
  reac_si = (double *)calloc(N,sizeof(double));
  diff_u = (double *)calloc(N,sizeof(double));
  diff_ui = (double *)calloc(N,sizeof(double));
  diff_u2 = (double *)calloc(N,sizeof(double));
  diff_u2i = (double *)calloc(N,sizeof(double));
  agr=(double *)calloc(N,sizeof(double));
  k2vec=(double *)calloc(N,sizeof(double));
  Du=(double *)calloc(N,sizeof(double));
  xrn=(double *)calloc(N,sizeof(double));
  xrn_ft =(fftw_complex *) fftw_malloc(L*(L/2+1)*sizeof(fftw_complex));
  Cu1 = (double complex *)calloc(Mh,sizeof(double complex));

  n1 = (int *)calloc(N,sizeof(int));
  n2 = (int *)calloc(N,sizeof(int));
  n3 = (int *)calloc(N,sizeof(int));
  n4 = (int *)calloc(N,sizeof(int));

  for (i=0;i<L;i++)
    for (j=0;j<L;j++) {
      ij=ncord(L,i,j,0,0);
      n1[ij]=ncord(L,i,j,1,0);
      n2[ij]=ncord(L,i,j,0,1);
      n3[ij]=ncord(L,i,j,-1,0);
      n4[ij]=ncord(L,i,j,0,-1);
    }

/* FFT Plans */
  plan_u1f = fftw_plan_dft_r2c_2d(L,L,xrn,xrn_ft,FFTW_MEASURE);
  plan_u1b = fftw_plan_dft_c2r_2d(L,L,xrn_ft,xrn, FFTW_MEASURE);

/**********  RNG initialisation  **********/
  gsl_rng_env_setup();
  rng = gsl_rng_alloc (type);
  gsl_rng_set (rng, seed);

/**********  Growth rate heterogeneity  **********/
  for(ij=0;ij<N;ij++)
    xrn[ij]=gsl_ran_gaussian(rng,1.0);

  for(i=0;i<L;i++) {
    if(i<Lh) ki=dk*(double)i;
    else ki=dk*(double)(i-L);
    for(j=0;j<Lh+1;j++) {
      ij=j+(Lh+1)*i;
      kj=dk*(double)j;
      k2vec[ij]=ki*ki+kj*kj;
    }
  }
  for (ij=0;ij<Mh;ij++) {
    Cu1[ij]=exp(-pow(k2vec[ij]/(2*agrcl*agrcl),agrcp));
  }
  fftw_execute(plan_u1f);
  for (ij=0;ij<Mh;ij++) {
    xrn_ft[ij]*=Cu1[ij];
  }
  fftw_execute(plan_u1b);
  for (ij=0;ij<N;ij++) {
    xrn[ij]=xrn[ij]/(double)N;
  }

  // normalize xrn
  xrn_avg=0;
  xrn2_avg=0;

  for (ij=0;ij<N;ij++) {
    xrn_avg+=(xrn[ij]-xrn_avg)/(((double) ij+1));
    xrn2_avg+=(xrn[ij]*xrn[ij]-xrn2_avg)/(((double) ij+1));
  }
  xrn_std=sqrt(xrn2_avg-xrn_avg*xrn_avg);
  //printf("%g %g\n",xrn_avg, xrn_std);

  for (ij=0;ij<N;ij++) {
    xrn[ij]=(xrn[ij]-xrn_avg)/xrn_std;
  }
  // now xrn is approximately distributed as a standard normal distrib.


  x0 = log(agr0);
  for (ij=0;ij<N;ij++) {
	agr[ij] = agr0*(0.25+0.75/(1.0+exp(-1 -xrn[ij])));
  }

  sprintf(name3,"output/%s_gr.prof",name1);
  outdat = fopen(name3,"a");
  t = 0;
  fwrite(agr,sizeof(double),N,outdat);
  fclose(outdat);

/**********  INITIAL CONDITIONS  **********/
  Kmax=Kmax0;
  Kmaxi=Kmax0;

  for (ij=0;ij<N;ij++) {
     u[ij]=0.0;
     mtx[ij]=0.0;
     w[ij]=0.0;
	 //K[ij]=Kmax;
     s[ij]=0.0;
  }

/**** I.C.: random spatially uncorrelated noise ****/
  if (ict==0) {
    for (i=0;i<L;i++)
      for (j=0;j<L;j++) {
        ij=ncord(L,i,j,0,0);
        u[ij]=fabs(gsl_ran_flat(rng,umin,umax));
      }
  }

/**** I.C.: random punctual perturbation in the center ****/
  if (ict==1) {
	i=L/2;
	j=L/2;
	for (i2=0;i2<L;i2++)
		for (j2=0;j2<L;j2++) {
			r2 = sqrt(pow(i2-i,2.0)+pow(j2-j,2.0));
			if (r2<5) {
				ij=ncord(L,i2,j2,0,0);
      			u[ij]+=umax;
      			}
		}
  }


/**** I.C.: random coffee ring ****/
  if (ict>1) {
	for(i3=0; i3<ict;i3++) {
	  aux=gsl_ran_flat(rng,0,2*M_PI);
	  //if(gsl_ran_bernoulli(rng,0.025)) {
	  if(gsl_ran_bernoulli(rng,0.05)) {
		do {
		  cx=gsl_ran_flat(rng,0,L);
		  cy=gsl_ran_flat(rng,0,L);
		  r2=sqrt(pow(cx-Lh,2.0)+pow(cy-Lh,2.0));
		} while(r2>rring);
	  }
	  else {
	    r2=gsl_ran_exponential(rng,0.005);
	    if(r2>1) continue;
	    cx=Lh+rring*(1-r2)*cos(aux);
	    cy=Lh+rring*(1-r2)*sin(aux);
	  }
	  i=(int) cx;
	  j=(int) cy;
	  for (i2=(i-10);i2<(i+10);i2++)
		for (j2=(j-10);j2<(j+10);j2++) {
	  	  r2 = sqrt(pow(i2-cx,2.0)+pow(j2-cy,2.0));
		  if (r2<rcol) {
			ij=ncord(L,i2,j2,0,0);
			u[ij]+=umax;
		  }
		}
    }
  }

/**********  SIMULATION  **********/
  for (it=0;it<=nme;it++) {
	if (it%nst==0) {
      t=dt*(double)it;
	  printf("%g ",t);
      if (flag_save_u) {	  
        sprintf(name3,"output/%s_u.pat",name1);
        outdat = fopen(name3,"a");
        fwrite(&t,sizeof(double),1,outdat); 
        fwrite(u,sizeof(double),N,outdat); 
        fclose(outdat);
	  }
      if (flag_save_m) {	  
        sprintf(name3,"output/%s_m.pat",name1);
        outdat = fopen(name3,"a");
        fwrite(&t,sizeof(double),1,outdat); 
        fwrite(mtx,sizeof(double),N,outdat); 
        fclose(outdat);
	  }
      if (flag_save_w) {
        sprintf(name3,"output/%s_w.pat",name1);
        outdat = fopen(name3,"a");
        fwrite(&t,sizeof(double),1,outdat); 
        fwrite(w,sizeof(double),N,outdat); 
        fclose(outdat);
	  }
      if (flag_save_s) {
        sprintf(name3,"output/%s_s.pat",name1);
        outdat = fopen(name3,"a");
        fwrite(&t,sizeof(double),1,outdat); 
        fwrite(s,sizeof(double),N,outdat); 
        fclose(outdat);
	  }
    }

    /* Heun method in time - FD in space */
 
    u_avg=0;
    for (ij=0;ij<N;ij++)
      u_avg+=u[ij];
    u_avg/=(double)N;
    reac_Kmax=-lbdK*u_avg;
    Kmaxi=Kmax+dt*reac_Kmax;
    for (ij=0;ij<N;ij++) {
      Du[ij] = gd*exp(-mtx[ij]/kd);	
      K=(Kmax+Kmin*ikk*pow(w[ij],nk))/(1+ikk*pow(w[ij],nk));
	  reac_u[ij]=agr[ij]*u[ij]*(1-u[ij]/K);
	  diff_u[ij]=fD*Du[ij]*(u[n1[ij]]+u[n2[ij]]+u[n3[ij]]+u[n4[ij]]-4*u[ij]);
	  diffpref = -fD1*Du[ij]/kd;	  
	  diff_u2[ij]= diffpref*((u[n1[ij]]-u[n3[ij]])*(mtx[n1[ij]]-mtx[n3[ij]])+(u[n2[ij]]-u[n4[ij]])*(mtx[n2[ij]]-mtx[n4[ij]]));
      ui[ij]=u[ij]+dt*(reac_u[ij]+diff_u[ij]+diff_u2[ij]);
	  reac_m[ij]=akm*pow(u[ij],nh)/(1+ikm*pow(u[ij],nh))-lbd*mtx[ij];
      mtxi[ij]=mtx[ij]+dt*reac_m[ij];
      reac_w[ij]=aw*u[ij]-lbdw*w[ij];
      wi[ij]=w[ij]+dt*reac_w[ij];
	  
      reac_s[ij]=aks*pow(u[ij],ns)/(1+iks*pow(u[ij],ns))-lbds*s[ij];
      si[ij]=s[ij]+dt*reac_s[ij];
    }
	if (it%nst==0) {
		printf(" (%g, %g) ",u_avg, Kmax);
		mn=1e10;
		mx=0;
		for (ij=0;ij<N;ij++) {
			if (mn>u[ij]) mn=u[ij];
	  		if (mx<u[ij]) mx=u[ij];
		}
		printf("%g %g %g | ", mn,mx,u[0]);
		mn=1e10;
		mx=0;
		for (ij=0;ij<N;ij++) {
			if (mn>w[ij]) mn=w[ij];
	  		if (mx<w[ij]) mx=w[ij];
		}
		printf("%g %g\n", mn,mx);


	}

	ui_avg=0;
    for (ij=0;ij<N;ij++)
      ui_avg+=ui[ij];
    ui_avg/=(double)N;
    reac_Kmaxi=-lbdK*ui_avg;
    Kmax=Kmax+dth*(reac_Kmax+reac_Kmaxi);
    for (ij=0;ij<N;ij++) {
      Du[ij] = gd*exp(-mtxi[ij]/kd);	  	 
      Ki=(Kmaxi+Kmin*ikk*pow(wi[ij],nk))/(1+ikk*pow(wi[ij],nk));	  
      reac_ui[ij]=agr[ij]*ui[ij]*(1-ui[ij]/Ki);
	  diff_ui[ij]=fD*Du[ij]*(ui[n1[ij]]+ui[n2[ij]]+ui[n3[ij]]+ui[n4[ij]]-4*ui[ij]);
	  diffpref = -fD1*Du[ij]/kd;	  
	  diff_u2i[ij]= diffpref*((ui[n1[ij]]-ui[n3[ij]])*(mtxi[n1[ij]]-mtxi[n3[ij]])+(ui[n2[ij]]-ui[n4[ij]])*(mtxi[n2[ij]]-mtxi[n4[ij]]));
      aux=u[ij]+dth*(reac_u[ij]+reac_ui[ij]+diff_u[ij]+diff_ui[ij]+diff_u2[ij]+diff_u2i[ij]);
	  u[ij]=aux;
      reac_mi[ij]=akm*pow(ui[ij],nh)/(1+ikm*pow(ui[ij],nh))-lbd*mtxi[ij];
      mtx[ij]=mtx[ij]+dth*(reac_m[ij]+reac_mi[ij]);
      reac_wi[ij]=aw*ui[ij]-lbdw*wi[ij];
      w[ij]=w[ij]+dth*(reac_w[ij]+reac_wi[ij]);
      reac_si[ij]=aks*pow(ui[ij],ns)/(1+iks*pow(ui[ij],ns))-lbds*si[ij];
      s[ij]=s[ij]+dth*(reac_s[ij]+reac_si[ij]);
    }
  }

  /**** writing final configuration to a file ****/

  if (flag_save_u) {  
    sprintf(name3,"output/%s_u.fc",name1);
    outdat = fopen(name3,"w");
    fwrite(&L,sizeof(int),1,outdat);
    fwrite(u,sizeof(double),N,outdat);
    fclose(outdat);
  }
  if (flag_save_m) {
    sprintf(name3,"output/%s_m.fc",name1);
    outdat = fopen(name3,"w");
    fwrite(&L,sizeof(int),1,outdat);
    fwrite(mtx,sizeof(double),N,outdat);
    fclose(outdat);
  }
  if (flag_save_w) {
    sprintf(name3,"output/%s_w.fc",name1);
    outdat = fopen(name3,"w");
    fwrite(&L,sizeof(int),1,outdat);
    fwrite(w,sizeof(double),N,outdat);
    fclose(outdat);
  }
  if (flag_save_s) {
    sprintf(name3,"output/%s_s.fc",name1);
    outdat = fopen(name3,"w");
    fwrite(&L,sizeof(int),1,outdat);
    fwrite(s,sizeof(double),N,outdat);
    fclose(outdat);
  }
  printf("\nQuitting.\n\n");

  return 0;
}
