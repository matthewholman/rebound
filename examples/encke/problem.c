#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "rebound.h"
#include "universal.h"
#include "integrator_encke.h"

static struct reb_dpconst7 dpcast(struct reb_dp7 dp){
    struct reb_dpconst7 dpc = {
        .p0 = dp.p0, 
        .p1 = dp.p1, 
        .p2 = dp.p2, 
        .p3 = dp.p3, 
        .p4 = dp.p4, 
        .p5 = dp.p5, 
        .p6 = dp.p6, 
    };
    return dpc;
}

static inline void add_cs(double* p, double* csp, double inp){
  double y = inp - *csp;

    double t = *p + y;
    *csp = (t - *p) - y;
    *p = t;
}

void enckeias15_0(double *apert, double *atc, double *xtrue, double *m, int n);

void insolarencke(struct reb_simulation* r, double epsilon, int is);

void consqp(struct reb_simulation* r, double *Ep, double L[NDIM], double *Lmp);

int main(int argc, char **argv)
{
    /*begin initialization*/
	//This is a test comment DH

    struct reb_simulation* r = reb_create_simulation();

    double tmax;

    int nsteps=1000,flgt=0,flg=0;    

    int iseed=1; // = 1000;

    double tv[NMAX2],diff,tticks,tav=0.;

    int ccount = 0,rcount = 0;
    double itav1=0.,itav2=0.;
    int ct1=0,ct2=0;

    double E0, E;
    double L0[NDIM], L0m, Lm, L[NDIM];
    clock_t timek;

    // command line input
    r->dt = atof(argv[1]);
    tmax = atof(argv[2]);
    double epsilon = atof(argv[3]);

    /*log spaced outputs*/
    /*diff = log10(tmax) - 0.;
    tticks = diff/((float)nsteps);
    tvl[0] = 0;
    tv[0] = 1; */
    
    diff = tmax;
   tticks = diff/((float)nsteps);
    tv[0] = 0.;

    // Set up output times.
    for(int i=1; i<(nsteps+1); i++){    
      /*tvl[i] = tvl[i-1] + tticks;
      tv[i] = pow(10,tvl[i]); */
      tv[i] = tv[i-1] + tticks;
    }

    for(int is=0; is<iseed; is++){        
      printf("seed=%d\n",is);
      timek = clock();
      printf("%g done\n",(float)is/(float)iseed);

      printf("here\n");
      fflush(stdout);

      // Set up initial conditions
      insolarencke(r, epsilon, is);

      // Calculate the initial energy and angular momentum
      consqp(r, &E0, L0, &L0m);

      int ic = 0;
      int i = 0;
      double hav = 0.;
      
      //t = 0.;
      r->ri_encke.t0 = 0.;

      // Open output file 
      char filen[50];
      FILE *fcons;      
      sprintf(filen,"data/fcons%d.txt",is);
      fcons = fopen(filen,"w");

      fprintf(fcons,"%0.16le %0.16le\n",r->t,E0);
      fflush(fcons);

      r->ri_encke.ccount=0;	
      
      while(r->t < tmax){
        i++;

	r->ri_encke.iter=0;

	reb_integrator_encke_step(r);

	if(rcount == 0){
            itav1 = itav1 + (float)(r->ri_encke.iter);
            ct1++;
	}else{
            itav2 = itav2 + (float)(r->ri_encke.iter);
            ct2++;
	}

        if((r->t) > (r->ri_encke.t0)){
	    // Updates the reference trajectory
	    update_two_bod(r,&flg,&rcount);

          flgt += flg;
        }
	
        hav += (r->dt);

    	
        if((r->t) > tv[ic]){
          ic++;
          consqp(r, &E,L,&Lm);
          fprintf(fcons,"%0.16le %0.16le\n",(r->t),E);
        }
        if(ic == (nsteps+1)) break; //done with integration
      }
      printf("it1=%g it2=%g\n",itav1/((float)ct1),itav2/((float)ct2));

      hav /= (double)i;
      consqp(r,&E,L,&Lm);
      timek = clock()-timek;
      double time_taken = ((double)timek)/CLOCKS_PER_SEC; // in seconds 
      tav += time_taken;

      printf("ccount=%d rcount=%d\n",ccount,rcount);
      printf("%g done, dE/E=%g, timet=%e hav=%g tdone=%g\n",(float)is/(float)iseed,(E0-E)/E0,time_taken,hav,tv[ic]);
      fclose(fcons);
    }

    tav /= (float)(iseed+1);
    printf("av compute time=%e\n",tav);

}


void insolarencke(struct reb_simulation* r, double epsilon, int is)
{

    /*Below data from Hairer pg. 13-14.  Convert time to years.
      Distance is in AU, mass in solar masses */

    int nprog = 5;

    r->N = nprog - 1;
    int N = r->N;
    int n3 = N*3;

    reb_integrator_encke_alloc(r);

    r->ri_encke.epsilon = epsilon;

    r->ri_encke.ccount = 0;
    r->ri_encke.iter = 0;    
    
    double* restrict const m = r->ri_encke.m;

    double* restrict const at = r->ri_encke.at;

    double* restrict const csx = r->ri_encke.csx;
    double* restrict const csv = r->ri_encke.csv;        
    double* restrict const csa0 = r->ri_encke.csa0;

    orbstruct* restrict const orb = r->ri_encke.orb;

    const struct reb_dpconst7 e  = dpcast(r->ri_encke.e);
    const struct reb_dpconst7 b  = dpcast(r->ri_encke.b);
    const struct reb_dpconst7 er = dpcast(r->ri_encke.er);
    const struct reb_dpconst7 br = dpcast(r->ri_encke.br);

    double mtemp[5] =
	{
	    1.00000597682, // Sun + inner planets
	    1. / 1047.355, // Jupiter
	    1. / 3501.6,   // Saturn
	    1. / 22869.,   // Uranus
	    1. / 19314.,   // Neptune
	};

    double xtemp[5][3] = 
	{
	    {-4.06428567034226e-3, -6.08813756435987e-3, -1.66162304225834e-6}, // Sun
	    {+3.40546614227466e+0, +3.62978190075864e+0, +3.42386261766577e-2}, // Jupiter
	    {+6.60801554403466e+0, +6.38084674585064e+0, -1.36145963724542e-1}, // Saturn
	    {+1.11636331405597e+1, +1.60373479057256e+1, +3.61783279369958e-1}, // Uranus
	    {-3.01777243405203e+1, +1.91155314998064e+0, -1.53887595621042e-1}, // Neptune
	};
    double vtemp[5][3] =
	{
	    {+6.69048890636161e-6, -6.33922479583593e-6, -3.13202145590767e-9}, // Sun
	    {-5.59797969310664e-3, +5.51815399480116e-3, -2.66711392865591e-6}, // Jupiter
	    {-4.17354020307064e-3, +3.99723751748116e-3, +1.67206320571441e-5}, // Saturn
	    {-3.25884806151064e-3, +2.06438412905916e-3, -2.17699042180559e-5}, // Uranus
	    {-2.17471785045538e-4, -3.11361111025884e-3, +3.58344705491441e-5}, // Neptune
	};

    double* xhel = malloc(sizeof(double)*n3);
    double* vhel = malloc(sizeof(double)*n3);        

    r->ri_encke.Msun = mtemp[0];

    double Msun = mtemp[0];    
    
    //Jupiter perturbing
    xtemp[1][0] = xtemp[1][0]+(float)is * 1e-14;

    for(int i=0; i<nprog; i++){    
	for(int k=0; k<NDIM; k++){                  
	    vtemp[i][k] = vtemp[i][k]*YEAR;
	}
    }
    
    /*force center of mass to be stationary and centered (not needed) */

    double xcm[NDIM],vcm[NDIM];
    
    double mt = 0.0;
    for(int i=0; i<nprog; i++){
	mt += mtemp[i];
    }
    for(int i=0; i<NDIM; i++){                    
	xcm[i] = 0.0;
	vcm[i] = 0.0;
    }
    for(int i=0; i<nprog; i++){                    
	for(int k=0; k<NDIM; k++){                      
	    xcm[k] += mtemp[i]*xtemp[i][k];
	    vcm[k] += mtemp[i]*vtemp[i][k];
	}
    }

    for(int i=0; i<NDIM; i++){                        
	vcm[i] /= mt;
	xcm[i] /= mt;
    }
    for(int i=0; i<nprog; i++){                          
	for(int k=0; k<NDIM; k++){                          
	    xtemp[i][k] -= xcm[k];
	    vtemp[i][k] -= vcm[k];
	}
    }

    for(int i=0; i<N; i++){
	m[i] = mtemp[i+1]; //Solar mass
    }

    for(int i=0; i<N; i++){
	for(int k=0; k<NDIM; k++){                                
	    xhel[3*i+k] = xtemp[i+1][k] - xtemp[0][k];
	    vhel[3*i+k] = vtemp[i+1][k] - vtemp[0][k];
	    orb[i].csxh[k] = 0.0;
	    orb[i].csvh[k] = 0.0;	    
	    at[3*i+k] = 0.;
	}
    }

    // Initializing the IAS15 constants
    for(int k=0; k<n3; k++){                              
	e.p0[k] = 0.; e.p1[k] = 0.; e.p2[k] = 0.; e.p3[k] = 0.; e.p4[k] = 0.; e.p5[k] = 0.; e.p6[k] = 0.;
	b.p0[k] = 0.; b.p1[k] = 0.; b.p2[k] = 0.; b.p3[k] = 0.; b.p4[k] = 0.; b.p5[k] = 0.; b.p6[k] = 0.;
	er.p0[k] = 0.; er.p1[k] = 0.; er.p2[k] = 0.; er.p3[k] = 0.; er.p4[k] = 0.; er.p5[k] = 0.;
	br.p6[k] = 0.; br.p0[k] = 0.; br.p1[k] = 0.; br.p2[k] = 0.; br.p3[k] = 0.; br.p4[k] = 0.; br.p5[k] = 0.; br.p6[k] = 0.;

	csa0[k] = 0.; csx[k] = 0.; csv[k] = 0;
    }

    // Initializing the orbit constants    
    for(int i=0; i<N; i++){                                  
	getorb(Msun, m[i], &xhel[3*i], &vhel[3*i], 1, &(orb[i]));	
    }

    // Calculate perturbed acceleration
    enckeias15_0(at,csa0,xhel,m,N);

    r->dt_last_done = 0.;
    r->N = 0;    

    // Add and initialize particles    
    for(int i=0; i<N; i++){

	struct reb_particle tp = {0};

	// set the mass.
	tp.m  = m[i];
	
	// perturbations are zero at first.
	tp.x  = 0.0;
	tp.y  = 0.0;
	tp.z  = 0.0;
	tp.vx = 0.0;
	tp.vy = 0.0;
	tp.vz = 0.0;
	
	// perturbed accelerations are not zero.
	tp.ax = at[3*i+0];		
	tp.ay = at[3*i+1];		
	tp.az = at[3*i+2];		

	reb_add(r, tp);
    }

    free(xhel);
    free(vhel);
    

}


// Calculate energy and angular momentum.
void consqp(struct reb_simulation* r,
	    double *Ep, double L[NDIM], double *Lmp)
{

    int n = r->N;

    int nprog = n + 1;

    double* restrict const m = r->ri_encke.m;
    orbstruct* restrict const orb = r->ri_encke.orb;    

    double xcm[NDIM],vcm[NDIM],sv[NDIM],sx[NDIM];
    double M,x[NMAX][NDIM],v[NMAX][NDIM];

    double xpert[NMAX][NDIM],vpert[NMAX][NDIM];
    double xhel[NMAX][NDIM],vhel[NMAX][NDIM];

    double E=0., Lm=0.;
    double diff, rij;

    double mtemp[NMAX];

    int* map = r->ri_encke.map;
    struct reb_particle* const particles = r->particles;    

    for(int i=0; i<n; i++){
        int mi = map[i];
	xpert[i][0]=particles[mi].x;
	xpert[i][1]=particles[mi].y;
	xpert[i][2]=particles[mi].z;			
	vpert[i][0]=particles[mi].vx;
	vpert[i][1]=particles[mi].vy;
	vpert[i][2]=particles[mi].vz;

	xhel[i][0] = orb[mi].x;
	xhel[i][1] = orb[mi].y;
	xhel[i][2] = orb[mi].z;

	vhel[i][0] = orb[mi].vx;
	vhel[i][1] = orb[mi].vy;
	vhel[i][2] = orb[mi].vz;
    }
    
    for(int i=0; i<n; i++){    
	for(int k=0; k<NDIM; k++){		
	    xhel[i][k] += xpert[i][k];
	    vhel[i][k] += vpert[i][k];
	}
    }

    M = r->ri_encke.Msun;
    E = 0.;
    Lm = 0.;

    mtemp[0] = M;
    for(int i=0; i<n; i++){
	mtemp[i+1] = m[i];
	M += m[i];
    }

    for(int i=0; i<NDIM; i++){        
	xcm[i] = 0;
	vcm[i] = 0;
	sx[i] = 0;
	sv[i] = 0;
	L[i] = 0;
    }

    for(int i=0; i<n; i++){        
	for(int k=0; k<NDIM; k++){        	
	    sx[k] -= m[i]*xhel[i][k];
	    sv[k] -= m[i]*vhel[i][k];
	} 
    }

    for(int i=0; i<NDIM; i++){            
	x[0][i] = sx[i]/M;
	v[0][i] = sv[i]/M;
    }

    for(int i=0; i<n; i++){            
	for(int k=0; k<NDIM; k++){        		
	    x[i+1][k] = xhel[i][k] + x[0][k];
	    v[i+1][k] = vhel[i][k] + v[0][k];
	}
    }

    for(int i=0; i<nprog; i++){            
	L[0] += mtemp[i]*(x[i][1]*v[i][2] - x[i][2]*v[i][1]);
	L[1] += mtemp[i]*(x[i][2]*v[i][0] - x[i][0]*v[i][2]);
	L[2] += mtemp[i]*(x[i][0]*v[i][1] - x[i][1]*v[i][0]);
	for(int k=0; k<NDIM; k++){        			
	    E += 1.0/2.0*mtemp[i]*v[i][k]*v[i][k];
	}
	for(int j=i+1; j<nprog; j++){
	    rij = 0;
	    for(int k=0; k<NDIM; k++){        			
		diff = x[i][k]-x[j][k];
		rij += diff*diff;
	    }
	    rij = sqrt(rij);
	    E -= GNEWT*mtemp[i]*mtemp[j]/rij;
        }
    }
    Lm = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]);

    *Lmp = Lm;
    *Ep = E;

}



// I believe this only gets called for the initial conditions
// This is calculating the mutual accelations among the non-sun masses, along with
// the associated compensated sums.
//
void enckeias15_0(double *apert, double *atc, double *xtrue, double *m, int n)
{
    double factor1[NMAX][NDIM];
    double r2,rij[NDIM],r3;

    int n3 = 3*n;
    
    for(int i=0; i<n3; i++){
	apert[i] = 0.;
	atc[i] = 0.;
    }

    // start at 0.						    	    	        
    for(int i=0; i<n; i++){  
	double xtrue2 = 0.;
	for(int k=0; k<NDIM; k++){          
	    xtrue2 += xtrue[3*i+k] * xtrue[3*i+k];
	}
	double xtrue3 = xtrue2 * sqrt(xtrue2);
	double tmp = GNEWT*m[i]/xtrue3;
	for(int k=0; k<NDIM; k++){                
	    factor1[i][k] = xtrue[3*i+k]*tmp;
	}
    }

    // start at 0.						    	    	        
    for(int i=0; i<n; i++){    
	for(int j=i+1;j< n; j++){
	    r2 = 0.;
	    for(int k=0; k<NDIM; k++){          	  
		rij[k] = xtrue[3*i+k] - xtrue[3*j+k];
		r2 += rij[k] * rij[k];
	    }
	    r3 = r2*sqrt(r2);
	    for(int k=0; k<NDIM; k++){          	  	  
		double factor = GNEWT*rij[k]/r3;
		add_cs(&(apert[3*i+k]), &(atc[3*i+k]), -factor*m[j]);
		add_cs(&(apert[3*i+k]), &(atc[3*i+k]), -factor1[j][k]);

		add_cs(&(apert[3*j+k]), &(atc[3*j+k]), factor*m[i]);
		add_cs(&(apert[3*j+k]), &(atc[3*j+k]), -factor1[i][k]);
	    }
	}
    }
}
