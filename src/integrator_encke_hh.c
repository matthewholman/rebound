/**
 * @file    integrator_encke_hh.c
 * @brief   EnckeHH integrator.
 * @author  David Hernandez <dmhernandez@cfa.harvard.edu>
 * @author  Matthew Holman <mholmancfa.harvard.edu>
 * @details This file implements the EnckeHH integration scheme.  
 * This routine implements Encke's method, in which the equations of
 * motion represent the difference between the actual trajectory and
 * a nearby, Keplerian reference trajectory.  EnckeHH builds upon the
 * IAS15 integrator (Rein & Spiegel 2014), using the same fifteenth
 * order integrator and care to control roundoff error, but Encke's
 * equations of motion. 
 *
 * As noted in integrator_ias15.c:
 * For more details see Rein & Spiegel 2014. Also see Everhart, 1985,
 * ASSL Vol. 115, IAU Colloq. 83, Dynamics of Comets, Their Origin 
 * and Evolution, 185 for the original implementation by Everhart.
 * Part of this code is based a function from the ORSE package.
 * See orsa.sourceforge.net for more details on their implementation.
 * 
 * @section     LICENSE
 * Copyright (c) 2020      David Hernandez, Matthew Holman.
 * Copyright (c) 2011-2012 Hanno Rein, Dave Spiegel.
 * Copyright (c) 2002-2004 Pasquale Tricarico.
 *
 * This file was developed within a branch of rebound, but the authors
 * of rebound are not responsible for the content of this file.
 *
 * EnckeHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "rebound.h"
#include "particle.h"
#include "gravity.h"
#include "tools.h"
#include "integrator.h"
#include "integrator_encke_hh.h"
#include "integrator_ias15.h"


#define EPSILONRESET 0.01
#include "universal.h"

/////////////////////////
//   Constants 

static const double safety_factor = 0.25;
static const double safety_factor2 = 0.75;

// Gauss Radau spacings
// The final value 1.0 is included for convenience.
static const double h[9]    = {0.0, 0.0562625605369221464656521910318, 0.180240691736892364987579942780, 0.352624717113169637373907769648, 0.547153626330555383001448554766, 0.734210177215410531523210605558, 0.885320946839095768090359771030, 0.977520613561287501891174488626, 1.0};
// Other constants
static const double rr[28] = {0.0562625605369221464656522, 0.1802406917368923649875799, 0.1239781311999702185219278, 0.3526247171131696373739078, 0.2963621565762474909082556, 0.1723840253762772723863278, 0.5471536263305553830014486, 0.4908910657936332365357964, 0.3669129345936630180138686, 0.1945289092173857456275408, 0.7342101772154105315232106, 0.6779476166784883850575584, 0.5539694854785181665356307, 0.3815854601022408941493028, 0.1870565508848551485217621, 0.8853209468390957680903598, 0.8290583863021736216247076, 0.7050802551022034031027798, 0.5326962297259261307164520, 0.3381673205085403850889112, 0.1511107696236852365671492, 0.9775206135612875018911745, 0.9212580530243653554255223, 0.7972799218243951369035945, 0.6248958964481178645172667, 0.4303669872307321188897259, 0.2433104363458769703679639, 0.0921996667221917338008147};
static const double c[21] = {-0.0562625605369221464656522, 0.0101408028300636299864818, -0.2365032522738145114532321, -0.0035758977292516175949345, 0.0935376952594620658957485, -0.5891279693869841488271399, 0.0019565654099472210769006, -0.0547553868890686864408084, 0.4158812000823068616886219, -1.1362815957175395318285885, -0.0014365302363708915424460, 0.0421585277212687077072973, -0.3600995965020568122897665, 1.2501507118406910258505441, -1.8704917729329500633517991, 0.0012717903090268677492943, -0.0387603579159067703699046, 0.3609622434528459832253398, -1.4668842084004269643701553, 2.9061362593084293014237913, -2.7558127197720458314421588};
static const double d[21] = {0.0562625605369221464656522, 0.0031654757181708292499905, 0.2365032522738145114532321, 0.0001780977692217433881125, 0.0457929855060279188954539, 0.5891279693869841488271399, 0.0000100202365223291272096, 0.0084318571535257015445000, 0.2535340690545692665214616, 1.1362815957175395318285885, 0.0000005637641639318207610, 0.0015297840025004658189490, 0.0978342365324440053653648, 0.8752546646840910912297246, 1.8704917729329500633517991, 0.0000000317188154017613665, 0.0002762930909826476593130, 0.0360285539837364596003871, 0.5767330002770787313544596, 2.2485887607691597933926895, 2.7558127197720458314421588};

// From integrator_ias15.c
static void copybuffers(const struct reb_dpconst7 _a, const struct reb_dpconst7 _b, int N3);

// Machine independent implementation of pow(*,1./7.)
// From integrator_ias15.c
static double sqrt7(double a){
    double x = 1.;
    for (int k=0; k<20;k++){  // A smaller number should be ok too.
        double x6 = x*x*x*x*x*x;
        x += (a/x6-x)/7.;
    }
    return x;
}

// From integrator_ias15.c
static inline void add_cs(double* p, double* csp, double inp){
  double y = inp - *csp;

    double t = *p + y;
    *csp = (t - *p) - y;
    *p = t;
}

// From integrator_ias15.c
static void free_dp7(struct reb_dp7* dp7){
    free(dp7->p0);
    free(dp7->p1);
    free(dp7->p2);
    free(dp7->p3);
    free(dp7->p4);
    free(dp7->p5);
    free(dp7->p6);
    dp7->p0 = NULL;
    dp7->p1 = NULL;
    dp7->p2 = NULL;
    dp7->p3 = NULL;
    dp7->p4 = NULL;
    dp7->p5 = NULL;
    dp7->p6 = NULL;
}

// From integrator_ias15.c
static void clear_dp7(struct reb_dp7* const dp7, const int N3){
    for (int k=0;k<N3;k++){
        dp7->p0[k] = 0.;
        dp7->p1[k] = 0.;
        dp7->p2[k] = 0.;
        dp7->p3[k] = 0.;
        dp7->p4[k] = 0.;
        dp7->p5[k] = 0.;
        dp7->p6[k] = 0.;
    }
}

// From integrator_ias15.c
static void realloc_dp7(struct reb_dp7* const dp7, const int N3){
    dp7->p0 = realloc(dp7->p0,sizeof(double)*N3);
    dp7->p1 = realloc(dp7->p1,sizeof(double)*N3);
    dp7->p2 = realloc(dp7->p2,sizeof(double)*N3);
    dp7->p3 = realloc(dp7->p3,sizeof(double)*N3);
    dp7->p4 = realloc(dp7->p4,sizeof(double)*N3);
    dp7->p5 = realloc(dp7->p5,sizeof(double)*N3);
    dp7->p6 = realloc(dp7->p6,sizeof(double)*N3);
    clear_dp7(dp7,N3);
}

// From integrator_ias15.c
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

void reb_integrator_encke_hh_reset(struct reb_simulation* r){
    r->ri_encke_hh.allocatedN  = 0;
    r->ri_encke_hh.map_allocated_N  = 0;
    free_dp7(&(r->ri_encke_hh.g));
    free_dp7(&(r->ri_encke_hh.e));
    free_dp7(&(r->ri_encke_hh.b));
    free_dp7(&(r->ri_encke_hh.csb));
    free_dp7(&(r->ri_encke_hh.er));
    free_dp7(&(r->ri_encke_hh.br));

    free(r->ri_encke_hh.at);
    r->ri_encke_hh.at =  NULL;

    free(r->ri_encke_hh.x0);
    r->ri_encke_hh.x0 =  NULL;
    free(r->ri_encke_hh.v0);
    r->ri_encke_hh.v0 =  NULL;
    free(r->ri_encke_hh.a0);
    r->ri_encke_hh.a0 =  NULL;

    free(r->ri_encke_hh.xt);
    r->ri_encke_hh.xt =  NULL;

    free(r->ri_encke_hh.atr);
    r->ri_encke_hh.atr =  NULL;
    
    free(r->ri_encke_hh.csx);
    r->ri_encke_hh.csx=  NULL;
    free(r->ri_encke_hh.csv);
    r->ri_encke_hh.csv=  NULL;
    free(r->ri_encke_hh.csa0);
    r->ri_encke_hh.csa0 =  NULL;

    for(int i=0; i<NSTEP; i++){
	free(r->ri_encke_hh.xkep[i]);
    }
    free(r->ri_encke_hh.vkep);
    r->ri_encke_hh.vkep = NULL;
    for(int i=0; i<NSTEP; i++){
	free(r->ri_encke_hh.vkep[i]);	
    }
    free(r->ri_encke_hh.vkep);
    r->ri_encke_hh.vkep = NULL;

    free(r->ri_encke_hh.csxhtemp);
    for(int i=0; i<NSTEP; i++){
	free(r->ri_encke_hh.csxhtemp[i]);
    }
    free(r->ri_encke_hh.csxhtemp);    
    r->ri_encke_hh.csxhtemp = NULL;

    for(int i=0; i<NSTEP; i++){
	free(r->ri_encke_hh.csvhtemp[i]);	
    }
    free(r->ri_encke_hh.csvhtemp);
    r->ri_encke_hh.csvhtemp = NULL;

    for(int i=0; i<NSTEP; i++){
	free(r->ri_encke_hh.xtwobod2[i]);	
    }
    free(r->ri_encke_hh.xtwobod2);
    r->ri_encke_hh.xtwobod2 = NULL;

    for(int i=0; i<NSTEP; i++){
	free(r->ri_encke_hh.xtwobod3[i]);	
    }
    free(r->ri_encke_hh.xtwobod3);
    r->ri_encke_hh.xtwobod3 = NULL;
    
    for(int i=0; i<NSTEP; i++){
	free(r->ri_encke_hh.factor1keep[i]);	
    }
    free(r->ri_encke_hh.factor1keep);
    r->ri_encke_hh.factor1keep = NULL;
    
    free(r->ri_encke_hh.gravitycsvr);
    r->ri_encke_hh.gravitycsvr =  NULL;

    free(r->ri_encke_hh.m);    
    r->ri_encke_hh.m =  NULL;
    
    free(r->ri_encke_hh.map);
    r->ri_encke_hh.map =  NULL;
}

void reb_integrator_encke_hh_alloc(struct reb_simulation* r){
    
    int N3 = 3*r->N;
    int N = r->N;    

    if (N3 > r->ri_encke_hh.allocatedN) {
        realloc_dp7(&(r->ri_encke_hh.g),N3);
        realloc_dp7(&(r->ri_encke_hh.b),N3);
        realloc_dp7(&(r->ri_encke_hh.csb),N3);
        realloc_dp7(&(r->ri_encke_hh.e),N3);
        realloc_dp7(&(r->ri_encke_hh.br),N3);
        realloc_dp7(&(r->ri_encke_hh.er),N3);
        r->ri_encke_hh.m = realloc(r->ri_encke_hh.m,sizeof(double)*N);	
        r->ri_encke_hh.at = realloc(r->ri_encke_hh.at,sizeof(double)*N3);
        r->ri_encke_hh.x0 = realloc(r->ri_encke_hh.x0,sizeof(double)*N3);
        r->ri_encke_hh.v0 = realloc(r->ri_encke_hh.v0,sizeof(double)*N3);
        r->ri_encke_hh.a0 = realloc(r->ri_encke_hh.a0,sizeof(double)*N3);

        r->ri_encke_hh.xt = realloc(r->ri_encke_hh.xt,sizeof(double)*N3);
	
        //r->ri_encke_hh.atr = realloc(r->ri_encke_hh.atr,sizeof(double)*N3);	
        r->ri_encke_hh.csx= realloc(r->ri_encke_hh.csx,sizeof(double)*N3);
        r->ri_encke_hh.csv= realloc(r->ri_encke_hh.csv,sizeof(double)*N3);
	r->ri_encke_hh.csa0 = realloc(r->ri_encke_hh.csa0,sizeof(double)*N3);

	r->ri_encke_hh.orb = realloc(r->ri_encke_hh.orb,sizeof(orbstruct)*N);

	r->ri_encke_hh.gravitycsvr = realloc(r->ri_encke_hh.gravitycsvr,sizeof(double)*N3);	

        double* restrict const csx = r->ri_encke_hh.csx; 
        double* restrict const csv = r->ri_encke_hh.csv; 

        orbstruct* restrict const orb = r->ri_encke_hh.orb;

	r->ri_encke_hh.xkep = realloc(r->ri_encke_hh.xkep,sizeof(double *)*NSTEP);
	for (int i=0; i<NSTEP; i++) {
	    r->ri_encke_hh.xkep[i] = realloc(r->ri_encke_hh.xkep[i],sizeof(double)*N3);
	}

	r->ri_encke_hh.vkep = realloc(r->ri_encke_hh.vkep,sizeof(double *)*NSTEP);
	for (int i=0; i<NSTEP; i++) {
	    r->ri_encke_hh.vkep[i] = realloc(r->ri_encke_hh.vkep[i],sizeof(double)*N3);
	}

	r->ri_encke_hh.csxhtemp = realloc(r->ri_encke_hh.csxhtemp,sizeof(double *)*NSTEP);
	for (int i=0; i<NSTEP; i++) {
	    r->ri_encke_hh.csxhtemp[i] = realloc(r->ri_encke_hh.csxhtemp[i],sizeof(double)*N3);
	}

	r->ri_encke_hh.csvhtemp = realloc(r->ri_encke_hh.csvhtemp,sizeof(double *)*NSTEP);
	for (int i=0; i<NSTEP; i++) {
	    r->ri_encke_hh.csvhtemp[i] = realloc(r->ri_encke_hh.csvhtemp[i],sizeof(double)*N3);
	}

	r->ri_encke_hh.xtwobod2 = realloc(r->ri_encke_hh.xtwobod2,sizeof(double *)*NSTEP);
	for (int i=0; i<NSTEP; i++) {
	    r->ri_encke_hh.xtwobod2[i] = realloc(r->ri_encke_hh.xtwobod2[i],sizeof(double)*N);
	}

	r->ri_encke_hh.xtwobod3 = realloc(r->ri_encke_hh.xtwobod3,sizeof(double *)*NSTEP);
	for (int i=0; i<NSTEP; i++) {
	    r->ri_encke_hh.xtwobod3[i] = realloc(r->ri_encke_hh.xtwobod3[i],sizeof(double)*N);
	}

	r->ri_encke_hh.factor1keep = realloc(r->ri_encke_hh.factor1keep,sizeof(double *)*NSTEP);
	for (int i=0; i<NSTEP; i++) {
	    r->ri_encke_hh.factor1keep[i] = realloc(r->ri_encke_hh.factor1keep[i],sizeof(double)*N3);
	}
	
        for (int i=0;i<N3;i++){
            // Initialize compensated summation coefficients
            csx[i] = 0.;
            csv[i] = 0.;
	}
	for (int i=0;i<N;i++){
	    for(int k=0; k<NDIM; k++){
		orb[i].csxh[k] = 0.;
		orb[i].csvh[k] = 0.;
	    }
        }
        r->ri_encke_hh.allocatedN = N3;
    }
    if (N3/3 > r->ri_encke_hh.map_allocated_N){
        r->ri_encke_hh.map = realloc(r->ri_encke_hh.map,sizeof(int)*(N3/3));
        for (int i=0;i<N3/3;i++){
            r->ri_encke_hh.map[i] = i;
        }
        r->ri_encke_hh.map_allocated_N = N3/3;
    }

}

// From integrator_ias15.c
void predict_next_step(double ratio, int N3,  const struct reb_dpconst7 _e, const struct reb_dpconst7 _b, const struct reb_dpconst7 e, const struct reb_dpconst7 b);

// Forward declaration
void predict_kepler(orbstruct* restrict const orb, 
		    double dt, double del0,
		    double **xkep,
		    double **vkep,
		    double **csxhtemp, double **csvhtemp,
		    double **xtwobod3, double **xtwobod2,
		    int n,
		    int si);

// This is the main routine.
int reb_integrator_encke_hh_step(struct reb_simulation* r)
{

    reb_integrator_encke_hh_alloc(r);

    struct reb_particle* const particles = r->particles;
    // map is not fully implemented yet
    int* map; // this map allow for integrating only a selection of particles
    int N = r->N;
    map = r->ri_encke_hh.map; // identity map
    
    double s[9];

    // Should this be moved into the simulation structure?
    static double gravitycsv[NMAX];

    double t = r->t;
    
    double dt_last_done = r->dt_last_done;
    double dt = r->dt;

    double const Msun = r->ri_encke_hh.Msun;
    double* restrict const m = r->ri_encke_hh.m;

    double const epsilon = r->ri_encke_hh.epsilon;

    // compensated sums for the perturbations
    double* restrict const csx = r->ri_encke_hh.csx; 
    double* restrict const csv = r->ri_encke_hh.csv; 
    double* restrict const csa0 = r->ri_encke_hh.csa0;

    double* restrict const x0 = r->ri_encke_hh.x0; 
    double* restrict const v0 = r->ri_encke_hh.v0; 
    double* restrict const a0 = r->ri_encke_hh.a0;

    double* restrict const at = r->ri_encke_hh.at; 
    double* restrict const atr = r->ri_encke_hh.atr;

    double* restrict const xt = r->ri_encke_hh.xt;

    double* restrict const gravitycsvr = r->ri_encke_hh.gravitycsvr;        
    struct reb_vec3d* gravity_cs = r->gravity_cs; 
    const struct reb_dpconst7 g  = dpcast(r->ri_encke_hh.g);
    const struct reb_dpconst7 e  = dpcast(r->ri_encke_hh.e);
    const struct reb_dpconst7 b  = dpcast(r->ri_encke_hh.b);
    const struct reb_dpconst7 csb = dpcast(r->ri_encke_hh.csb);
    const struct reb_dpconst7 er = dpcast(r->ri_encke_hh.er);
    const struct reb_dpconst7 br = dpcast(r->ri_encke_hh.br);

    double** restrict const xkep = r->ri_encke_hh.xkep;
    double** restrict const vkep = r->ri_encke_hh.vkep;        

    double** restrict const csxhtemp = r->ri_encke_hh.csxhtemp;
    double** restrict const csvhtemp = r->ri_encke_hh.csvhtemp;    

    double** restrict const xtwobod2 = r->ri_encke_hh.xtwobod2;
    double** restrict const xtwobod3 = r->ri_encke_hh.xtwobod3;

    double** restrict const factor1keep = r->ri_encke_hh.factor1keep;        

    int iterations;
    int const si=8; 

    double t_beginning = t;
    double dtinit = dt;
    double dt_new;

    double sign();

    int N3 = 3*N;
    double min_dt = 0;

    double del0 = t_beginning-(r->ri_encke_hh.t0);

    // Clear arrays.
    for(int i=0; i<N3; i++){
	csb.p0[i] = 0.; csb.p1[i] = 0.; csb.p2[i] = 0.; csb.p3[i] = 0.; csb.p4[i] = 0.; csb.p5[i] = 0.; csb.p6[i] = 0.;
    }

    // This is the point where IAS15 loads the particle positions, velocities, and accelerations.

    for(int k=0;k<N;k++) {
        int mk = map[k];
        m[k]      = particles[mk].m;	
        x0[3*k]   = particles[mk].x;
        x0[3*k+1] = particles[mk].y;
        x0[3*k+2] = particles[mk].z;
        v0[3*k]   = particles[mk].vx;
        v0[3*k+1] = particles[mk].vy;
        v0[3*k+2] = particles[mk].vz;
        a0[3*k]   = particles[mk].ax;
        a0[3*k+1] = particles[mk].ay; 
        a0[3*k+2] = particles[mk].az;
    }

    /*
    for(int i=0; i<N3; i++){    
	x0[i] = x[i];
	v0[i] = v[i];
	a0[i] = at[i];
    }
    */

    // We know what time step we are going to try to take.
    // Get the positions and the velocities at the substeps and at the end
    // of the step.
    orbstruct* restrict const orb = r->ri_encke_hh.orb;

    predict_kepler(orb,dtinit,del0,xkep,vkep,csxhtemp,csvhtemp,xtwobod3,xtwobod2,N,si);

    // Initializing the g values
    for(int k=0;k<N3;++k) {
	g.p0[k] = b.p6[k]*d[15] + b.p5[k]*d[10] + b.p4[k]*d[6] + b.p3[k]*d[3]  + b.p2[k]*d[1]  + b.p1[k]*d[0]  + b.p0[k];
	g.p1[k] = b.p6[k]*d[16] + b.p5[k]*d[11] + b.p4[k]*d[7] + b.p3[k]*d[4]  + b.p2[k]*d[2]  + b.p1[k];
	g.p2[k] = b.p6[k]*d[17] + b.p5[k]*d[12] + b.p4[k]*d[8] + b.p3[k]*d[5]  + b.p2[k];
	g.p3[k] = b.p6[k]*d[18] + b.p5[k]*d[13] + b.p4[k]*d[9] + b.p3[k];
	g.p4[k] = b.p6[k]*d[19] + b.p5[k]*d[14] + b.p4[k];
	g.p5[k] = b.p6[k]*d[20] + b.p5[k];
	g.p6[k] = b.p6[k];
    }
    
    double predictor_corrector_error = 1e300;
    double predictor_corrector_error_last = 2;
    iterations = 0;
  
    while(1){

        if(predictor_corrector_error<1e-16){
            break;
        }
        if(iterations > 2 && predictor_corrector_error_last <= predictor_corrector_error){
	    r->ri_encke_hh.ccount +=1;
            break;
        }
        if (iterations>=12){
            printf("Too many iterations.  Reduce timestep\n");
	    fflush(stdout);
            exit(0);                            
        }
        predictor_corrector_error_last = predictor_corrector_error;
        predictor_corrector_error = 0.;
        iterations++;
        r->ri_encke_hh.iter = iterations;

	for(int n=1;n<si;n++) {                   // Loop over interval using Gauss-Radau spacings

            
            s[0] = dtinit * h[n];
            s[1] = s[0] * s[0] / 2.;
            s[2] = s[1] * h[n] / 3.;

            s[3] = s[2] * h[n] / 2.;
            s[4] = 3. * s[3] * h[n] / 5.;
            s[5] = 2. * s[4] * h[n] / 3.;
            s[6] = 5. * s[5] * h[n] / 7.;
            s[7] = 3. * s[6] * h[n] / 4.;
            s[8] = 7. * s[7] * h[n] / 9.;
            

	    for(int i=0; i<N; i++){  	    
	        // Predict perturbation positions at interval n using b values

		// These statements could be eliminated but leaving them for clarity.
                int k0 = 3*i;
                int k1 = 3*i+1;
                int k2 = 3*i+2;

		// The order of the sums below, which proceed left to right, matters.


                xt[k0] = -csx[k0] + (s[8]*b.p6[k0] + s[7]*b.p5[k0] + s[6]*b.p4[k0] + s[5]*b.p3[k0] + s[4]*b.p2[k0] + s[3]*b.p1[k0] + s[2]*b.p0[k0] + s[1]*a0[k0] + s[0]*v0[k0] ) + x0[k0];		

                xt[k1] = -csx[k1] + (s[8]*b.p6[k1] + s[7]*b.p5[k1] + s[6]*b.p4[k1] + s[5]*b.p3[k1] + s[4]*b.p2[k1] + s[3]*b.p1[k1] + s[2]*b.p0[k1] + s[1]*a0[k1] + s[0]*v0[k1] ) + x0[k1];
		
                xt[k2] = -csx[k2] + (s[8]*b.p6[k2] + s[7]*b.p5[k2] + s[6]*b.p4[k2] + s[5]*b.p3[k2] + s[4]*b.p2[k2] + s[3]*b.p1[k2] + s[2]*b.p0[k2] + s[1]*a0[k2] + s[0]*v0[k2] ) + x0[k2];


/*
also comment out s vectors above if use this

                    xt[k0] = -csx[k0] + ((((((((b.p6[k0]*7.*h[n]/9. + b.p5[k0])*3.*h[n]/4. + b.p4[k0])*5.*h[n]/7. + b.p3[k0])*2.*h[n]/3. + b.p2[k0])*3.*h[n]/5. + b.p1[k0])*h[n]/2. + b.p0[k0])*h[n]/3. + a0[k0])*dtinit*h[n]/2. + v0[k0])*dtinit*h[n] + x0[k0];

                    xt[k1] = -csx[k1] + ((((((((b.p6[k1]*7.*h[n]/9. + b.p5[k1])*3.*h[n]/4. + b.p4[k1])*5.*h[n]/7. + b.p3[k1])*2.*h[n]/3. + b.p2[k1])*3.*h[n]/5. + b.p1[k1])*h[n]/2. + b.p0[k1])*h[n]/3. + a0[k1])*dtinit*h[n]/2. + v0[k1])*dtinit*h[n] + x0[k1];

                    xt[k2] = -csx[k2] + ((((((((b.p6[k2]*7.*h[n]/9. + b.p5[k2])*3.*h[n]/4. + b.p4[k2])*5.*h[n]/7. + b.p3[k2])*2.*h[n]/3. + b.p2[k2])*3.*h[n]/5. + b.p1[k2])*h[n]/2. + b.p0[k2])*h[n]/3. + a0[k2])*dtinit*h[n]/2. + v0[k2])*dtinit*h[n] + x0[k2];
  */                  
            }

	    // Calculate acceleration based on the current positions, 
	    // for substep n.
	    encke_hh_acceleration(at,gravitycsv,
			    NULL, NULL,
			    xt,Msun,m,
			    xkep[n],vkep[n],
			    xtwobod2[n],xtwobod3[n],
			    csxhtemp[n],csvhtemp[n],
			    factor1keep[n],
			    N);
	    
            switch (n) {                            // Improve b and g values
	    case 1:
		for(int k=0;k<N3;++k) {
		    double tmp = g.p0[k];
		    double gk = at[k];
		    double gk_cs = gravitycsv[k];
		    add_cs(&gk, &gk_cs, -a0[k]);
		    add_cs(&gk, &gk_cs, csa0[k]);
		    g.p0[k]  = gk/rr[0];
		    add_cs(&(b.p0[k]), &(csb.p0[k]), g.p0[k]-tmp);
		} 
		break;
	    case 2: 
		for(int k=0;k<N3;++k) {
		    double tmp = g.p1[k];
		    double gk = at[k];
		    double gk_cs = gravitycsv[k];
		    add_cs(&gk, &gk_cs, -a0[k]);
		    add_cs(&gk, &gk_cs, csa0[k]);
		    g.p1[k] = (gk/rr[1] - g.p0[k])/rr[2];
		    tmp = g.p1[k] - tmp;
		    add_cs(&(b.p0[k]), &(csb.p0[k]), tmp * c[0]);
		    add_cs(&(b.p1[k]), &(csb.p1[k]), tmp);
		} 
		break;
	    case 3:
		for(int k=0;k<N3;++k) {
		    double tmp = g.p2[k];
		    double gk = at[k];
		    double gk_cs = gravitycsv[k];
		    add_cs(&gk, &gk_cs, -a0[k]);
		    add_cs(&gk, &gk_cs, csa0[k]);
		    g.p2[k] = ((gk/rr[3] - g.p0[k])/rr[4] - g.p1[k])/rr[5];
		    tmp = g.p2[k] - tmp;
		    add_cs(&(b.p0[k]), &(csb.p0[k]), tmp * c[1]);
		    add_cs(&(b.p1[k]), &(csb.p1[k]), tmp * c[2]);
		    add_cs(&(b.p2[k]), &(csb.p2[k]), tmp);
		} break;
	    case 4:
		for(int k=0;k<N3;++k) {
		    double tmp = g.p3[k];
		    double gk = at[k];
		    double gk_cs = gravitycsv[k];
		    add_cs(&gk, &gk_cs, -a0[k]);
		    add_cs(&gk, &gk_cs, csa0[k]);
		    g.p3[k] = (((gk/rr[6] - g.p0[k])/rr[7] - g.p1[k])/rr[8] - g.p2[k])/rr[9];
		    tmp = g.p3[k] - tmp;
		    add_cs(&(b.p0[k]), &(csb.p0[k]), tmp * c[3]);
		    add_cs(&(b.p1[k]), &(csb.p1[k]), tmp * c[4]);
		    add_cs(&(b.p2[k]), &(csb.p2[k]), tmp * c[5]);
		    add_cs(&(b.p3[k]), &(csb.p3[k]), tmp);
		} break;
	    case 5:
		for(int k=0;k<N3;++k) {
		    double tmp = g.p4[k];
		    double gk = at[k];
		    double gk_cs = gravitycsv[k];
		    add_cs(&gk, &gk_cs, -a0[k]);
		    add_cs(&gk, &gk_cs, csa0[k]);
		    g.p4[k] = ((((gk/rr[10] - g.p0[k])/rr[11] - g.p1[k])/rr[12] - g.p2[k])/rr[13] - g.p3[k])/rr[14];
		    tmp = g.p4[k] - tmp;
		    add_cs(&(b.p0[k]), &(csb.p0[k]), tmp * c[6]);
		    add_cs(&(b.p1[k]), &(csb.p1[k]), tmp * c[7]);
		    add_cs(&(b.p2[k]), &(csb.p2[k]), tmp * c[8]);
		    add_cs(&(b.p3[k]), &(csb.p3[k]), tmp * c[9]);
		    add_cs(&(b.p4[k]), &(csb.p4[k]), tmp);
		} break;
	    case 6:
		for(int k=0;k<N3;++k) {
		    double tmp = g.p5[k];
		    double gk = at[k];
		    double gk_cs = gravitycsv[k];
		    add_cs(&gk, &gk_cs, -a0[k]);
		    add_cs(&gk, &gk_cs, csa0[k]);
		    g.p5[k] = (((((gk/rr[15] - g.p0[k])/rr[16] - g.p1[k])/rr[17] - g.p2[k])/rr[18] - g.p3[k])/rr[19] - g.p4[k])/rr[20];
		    tmp = g.p5[k] - tmp;
		    add_cs(&(b.p0[k]), &(csb.p0[k]), tmp * c[10]);
		    add_cs(&(b.p1[k]), &(csb.p1[k]), tmp * c[11]);
		    add_cs(&(b.p2[k]), &(csb.p2[k]), tmp * c[12]);
		    add_cs(&(b.p3[k]), &(csb.p3[k]), tmp * c[13]);
		    add_cs(&(b.p4[k]), &(csb.p4[k]), tmp * c[14]);
		    add_cs(&(b.p5[k]), &(csb.p5[k]), tmp);
		} break;
	    case 7:
                {
                    double maxak = 0.0;
                    double maxb6ktmp = 0.0;
                    for(int k=0;k<N3;++k) {
                        double tmp = g.p6[k];
                        double gk = at[k];
                        double gk_cs = gravitycsv[k];
                        add_cs(&gk, &gk_cs, -a0[k]);
                        add_cs(&gk, &gk_cs, csa0[k]);
                        g.p6[k] = ((((((gk/rr[21] - g.p0[k])/rr[22] - g.p1[k])/rr[23] - g.p2[k])/rr[24] - g.p3[k])/rr[25] - g.p4[k])/rr[26] - g.p5[k])/rr[27];
                        tmp = g.p6[k] - tmp;    
                        add_cs(&(b.p0[k]), &(csb.p0[k]), tmp * c[15]);
                        add_cs(&(b.p1[k]), &(csb.p1[k]), tmp * c[16]);
                        add_cs(&(b.p2[k]), &(csb.p2[k]), tmp * c[17]);
                        add_cs(&(b.p3[k]), &(csb.p3[k]), tmp * c[18]);
                        add_cs(&(b.p4[k]), &(csb.p4[k]), tmp * c[19]);
                        add_cs(&(b.p5[k]), &(csb.p5[k]), tmp * c[20]);
                        add_cs(&(b.p6[k]), &(csb.p6[k]), tmp);

                        double ak  = fabs(at[k]);
                        if (isnormal(ak) && ak>maxak){
                            maxak = ak;
                        }
                        double b6ktmp = fabs(tmp);  // change of b6ktmp coefficient
                        if (isnormal(b6ktmp) && b6ktmp>maxb6ktmp){
                            maxb6ktmp = b6ktmp;
                        }
                    }
                    predictor_corrector_error = maxb6ktmp/maxak;
                    break;
                }
            }
        }
    }
    
    /* ----------------------------------------
       -----------------------------------*/
    /*BEGIN TIMESTEP ROUTINE */
    
    // Set time back to initial value (will be updated below) 
    t = t_beginning;
    // Find new timestep
    
    double dt_done = dtinit;

    if ((epsilon > 0.) && (t != 0.0)){

        // if error estimate is available increase by more educated guess
	adjust_step(xkep,xtwobod3,at,dt_done,&dt_new,epsilon,si,N,factor1keep);
        
        if (fabs(dt_new)< min_dt) dt_new = fabs(min_dt) * sign(dt_new);
        if (fabs(dt_new/dt_done) < safety_factor) { // New timestep is significantly smaller.

            // Reset particles
            for(int k=0;k<N;++k) {
                int mk = map[k];
                particles[mk].x = x0[3*k+0]; // Set inital position
                particles[mk].y = x0[3*k+1];
                particles[mk].z = x0[3*k+2];

                particles[mk].vx = v0[3*k+0];    // Set inital velocity
                particles[mk].vy = v0[3*k+1];
                particles[mk].vz = v0[3*k+2];
                
                particles[mk].ax = a0[3*k+0];    // Set inital acceleration
                particles[mk].ay = a0[3*k+1];
                particles[mk].az = a0[3*k+2];
            }

            dt = dt_new;
            if (dt_last_done!=0.){       // Do not predict next e/b values if this is the first time step.
                double ratio = dt/ dt_last_done;
                predict_next_step(ratio, N3, er, br, e, b);
            }
            return 0; // Step rejected. Do again. 
        }       
        if (fabs(dt_new/dt_done) > 1.0) {   // New timestep is larger.
            if (dt_new/dt_done > 1./safety_factor2) dt_new = dt_done /safety_factor2; // Don't increase the timestep by too much compared to the last one.
        }
        dt = dt_new;
    } else{

	dt_new = dt_done;
	dt = dt_new;
	
    }

    // Find new position and velocity values at end of the sequence
    double dt_done2 = dt_done * dt_done;

    //END TIME STEP STUFF

    /*************************************/

    for(int k=0;k<N3;++k) {

        
	{
	    add_cs(&(x0[k]), &(csx[k]), b.p6[k]/72.*dt_done2);
	    add_cs(&(x0[k]), &(csx[k]), b.p5[k]/56.*dt_done2);
	    add_cs(&(x0[k]), &(csx[k]), b.p4[k]/42.*dt_done2);
	    add_cs(&(x0[k]), &(csx[k]), b.p3[k]/30.*dt_done2);
	    add_cs(&(x0[k]), &(csx[k]), b.p2[k]/20.*dt_done2);
	    add_cs(&(x0[k]), &(csx[k]), b.p1[k]/12.*dt_done2);
	    add_cs(&(x0[k]), &(csx[k]), b.p0[k]/6.*dt_done2);
	    add_cs(&(x0[k]), &(csx[k]), a0[k]/2.*dt_done2);
	    add_cs(&(x0[k]), &(csx[k]), v0[k]*dt_done);
	}
	{
	    add_cs(&(v0[k]), &(csv[k]), b.p6[k]/8.*dt_done);
	    add_cs(&(v0[k]), &(csv[k]), b.p5[k]/7.*dt_done);
	    add_cs(&(v0[k]), &(csv[k]), b.p4[k]/6.*dt_done);
	    add_cs(&(v0[k]), &(csv[k]), b.p3[k]/5.*dt_done);
	    add_cs(&(v0[k]), &(csv[k]), b.p2[k]/4.*dt_done);
	    add_cs(&(v0[k]), &(csv[k]), b.p1[k]/3.*dt_done);
	    add_cs(&(v0[k]), &(csv[k]), b.p0[k]/2.*dt_done);
	    add_cs(&(v0[k]), &(csv[k]), a0[k]*dt_done);
	}
    
    /*


            add_cs(&(x0[k]), &(csx[k]), b.p6[k]/72.*dt_done*dt_done);
            add_cs(&(x0[k]), &(csx[k]), b.p5[k]/56.*dt_done*dt_done);
            add_cs(&(x0[k]), &(csx[k]), b.p4[k]/42.*dt_done*dt_done);
            add_cs(&(x0[k]), &(csx[k]), b.p3[k]/30.*dt_done*dt_done);
            add_cs(&(x0[k]), &(csx[k]), b.p2[k]/20.*dt_done*dt_done);
            add_cs(&(x0[k]), &(csx[k]), b.p1[k]/12.*dt_done*dt_done);
            add_cs(&(x0[k]), &(csx[k]), b.p0[k]/6.*dt_done*dt_done);
            add_cs(&(x0[k]), &(csx[k]), a0[k]/2.*dt_done*dt_done);
            add_cs(&(x0[k]), &(csx[k]), v0[k]*dt_done);

            add_cs(&(v0[k]), &(csv[k]), b.p6[k]/8.*dt_done);
            add_cs(&(v0[k]), &(csv[k]), b.p5[k]/7.*dt_done);
            add_cs(&(v0[k]), &(csv[k]), b.p4[k]/6.*dt_done);
            add_cs(&(v0[k]), &(csv[k]), b.p3[k]/5.*dt_done);
            add_cs(&(v0[k]), &(csv[k]), b.p2[k]/4.*dt_done);
            add_cs(&(v0[k]), &(csv[k]), b.p1[k]/3.*dt_done);
            add_cs(&(v0[k]), &(csv[k]), b.p0[k]/2.*dt_done);
            add_cs(&(v0[k]), &(csv[k]), a0[k]*dt_done);
         */   
    }


    //}

    // This is where the particle arrays would be used.
    /*
    for(int i=0; i<N3; i++){
	x[i] = x0[i];
	v[i] = v0[i];
    }
    */

    r->t += dt_done;
    r->dt_last_done = dt_done;
    r->dt = dt;

    encke_hh_acceleration(at,csa0,
			  atr, gravitycsvr,
			  x0,Msun,m,
			  xkep[si],vkep[si],
			  xtwobod2[si],xtwobod3[si],
			  csxhtemp[si],csvhtemp[si],
			  factor1keep[si],
			  N);	    

    // Swap particle buffers
    for(int k=0;k<N;++k) {
        int mk = map[k];
        particles[mk].x = x0[3*k+0];  // Set final position
        particles[mk].y = x0[3*k+1];
        particles[mk].z = x0[3*k+2];

        particles[mk].vx = v0[3*k+0]; // Set final velocity
        particles[mk].vy = v0[3*k+1];
        particles[mk].vz = v0[3*k+2];

        particles[mk].ax = at[3*k+0]; // Set final acceleration
        particles[mk].ay = at[3*k+1];
        particles[mk].az = at[3*k+2];
    }

    // Update the orbits, without reset, for all bodies.
    for(int i=0; i<N; i++){
	getorb(Msun, m[i], &xkep[si][3*i], &vkep[si][3*i], 0, &(orb[i]));	    	    	    	
    }

    // Update the compensated sums for all bodies.    
    for(int i=0; i<N; i++){
	for(int k=0; k<NDIM; k++){
	    orb[i].csxh[k] = csxhtemp[si][3*i+k];
	    orb[i].csvh[k] = csvhtemp[si][3*i+k];	    
	}
    }
    
    copybuffers(e,er,N3);     
    copybuffers(b,br,N3);       
    double ratio = dt/ dt_done;
    predict_next_step(ratio, N3, e, b, e, b);

    return 1; // Success.
}       


// Checks if retification is needed and act accordingly
void update_two_bod(struct reb_simulation* r,
		    int *flg,
		    int *rcount)
{

    int ikeep;    
    static double apertp[NDIM];

    // Should get rid of NMAX and make these dynamic.
    static double csxht[NMAX],csvht[NMAX],xhel[NMAX],vhel[NMAX];
    static double xpert[NMAX],vpert[NMAX];

    orbstruct* restrict const orb = r->ri_encke_hh.orb;

    double* restrict const m = r->ri_encke_hh.m;

    double* restrict const csx = r->ri_encke_hh.csx;
    double* restrict const csv = r->ri_encke_hh.csv;             
    double* restrict const csa0 = r->ri_encke_hh.csa0;

    //double* restrict const gravitycsvr = r->ri_encke_hh.gravitycsvr;
    //double* restrict const atr = r->ri_encke_hh.atr;             

    const struct reb_dpconst7 e  = dpcast(r->ri_encke_hh.e);
    const struct reb_dpconst7 b  = dpcast(r->ri_encke_hh.b);
    const struct reb_dpconst7 er = dpcast(r->ri_encke_hh.er);
    const struct reb_dpconst7 br = dpcast(r->ri_encke_hh.br);

    double Msun = r->ri_encke_hh.Msun;

    int N = r->N;

    *flg = 0;
    int flgrpt = 0;
    *rcount = 0;

    int* map = r->ri_encke_hh.map;
    struct reb_particle* const particles = r->particles;

    for(int k=0;k<N;k++) {
        int mk = map[k];
	xpert[3*k]   = particles[mk].x;
        xpert[3*k+1] = particles[mk].y;
        xpert[3*k+2] = particles[mk].z;
        vpert[3*k]   = particles[mk].vx;
        vpert[3*k+1] = particles[mk].vy;
        vpert[3*k+2] = particles[mk].vz;
    }    

    r->ri_encke_hh.t0 = r->t;

    // Loop over the planets, checking if any
    // satisfy the rectification criterion
    //
    for(int i=0; i<N; i++){
	if(*flg == 0){
	    double rp = 0.;
	    //double vp = 0.;
	    for(int k=0; k<NDIM; k++){      
		rp += xpert[3*i+k] * xpert[3*i+k];
		//vp += vpert[3*i+k] * vpert[3*i+k];
	    }
	    rp = sqrt(rp);
	    //vp = sqrt(vp);
	    double q = 1.0/(orb[i].peri) * (rp);
	    if(q > EPSILONRESET) {
		*rcount = 1;
		if(flgrpt == 1){
		    printf("error resetting, need improve code\n");
		    exit(0);
		}
		*flg = 1;
		ikeep = i;
		flgrpt = 1;
	    }
	}
    }

    for(int i=0; i<N; i++){
	xhel[3*i+0] = orb[i].x;
	xhel[3*i+1] = orb[i].y;
	xhel[3*i+2] = orb[i].z;

	vhel[3*i+0] = orb[i].vx;
	vhel[3*i+1] = orb[i].vy;
	vhel[3*i+2] = orb[i].vz;
    }

    // A reference orbit needs to be reset.
    if(*flg == 1){

	// Compute true positions and velocities
	for(int i=0; i<N; i++){                  
	    for(int k=0; k<NDIM; k++){
		csxht[3*i+k] = orb[i].csxh[k];
		csvht[3*i+k] = orb[i].csvh[k];		
		add_cs(&(xhel[3*i+k]), &(csxht[3*i+k]), xpert[3*i+k]);
		add_cs(&(vhel[3*i+k]), &(csvht[3*i+k]), vpert[3*i+k]);
	    }
	}

	// Reset the reference position, velocity, and compensated sums
	// for particular orbit.
	int i = ikeep;
	for(int k=0; k<NDIM; k++){
	    orb[i].csxh[k] = csxht[3*i+k];
	    orb[i].csvh[k] = csvht[3*i+k];
	}

	// Reset the orbit for body i.
	getorb(Msun, m[i], &xhel[3*i], &vhel[3*i], 1, &(orb[i]));	

	// Recalculate the Encke acceleration.
	encke_hh_acceleration_i(apertp,&csa0[3*i],xhel,m,N,i);	
	
	// Need to reset the particle also.
	struct reb_particle* const particles = r->particles;
	
	particles[i].x = 0.0;
	particles[i].y = 0.0;
	particles[i].z = 0.0;	    	    	    
	particles[i].vx = 0.0;
	particles[i].vy = 0.0;
	particles[i].vz = 0.0;	    	    	    
	particles[i].ax = apertp[0];
	particles[i].ay = apertp[1];
	particles[i].az = apertp[2];

	// Perhaps there is a better way to reset these values.
	for(int k=3*i;k<(3*i+3);k++){
	    csx[k] = 0.; csv[k] = 0.;
	    e.p0[k] = 0.; e.p1[k] = 0.; e.p2[k] = 0.; e.p3[k] = 0.; e.p4[k] = 0.; e.p5[k] = 0.; e.p6[k] = 0.;
	    b.p0[k] = 0.; b.p1[k] = 0.; b.p2[k] = 0.; b.p3[k] = 0.; b.p4[k] = 0.; b.p5[k] = 0.; b.p6[k] = 0.;
	    er.p0[k] = 0.; er.p1[k] = 0.; er.p2[k] = 0.; er.p3[k] = 0.; er.p4[k] = 0.; er.p5[k] = 0.;
	    br.p6[k] = 0.; br.p0[k] = 0.; br.p1[k] = 0.; br.p2[k] = 0.; br.p3[k] = 0.; br.p4[k] = 0.; br.p5[k] = 0.; br.p6[k] = 0.;
	}
    }
}

// Same as getorb but is only for one body.
void getorb(double Msun, double m, double *xhel, double *vhel, int reset, orbstruct *orb)
{
    double v2 = 0;
    double r0 = 0;
    double eta = 0;
    for(int k=0; k<NDIM; k++){
	v2 += vhel[k] * vhel[k];
	r0 += xhel[k] * xhel[k];	    
	eta += xhel[k] * vhel[k];
    }
    r0 = sqrt(r0);

    if(reset){
	double kc,beta,a,ec,es,e;

	kc = GNEWT*(Msun + m);
	beta = 2*kc/r0 - v2;
	a = kc/beta;
	ec = 1-r0/a;
	es = eta/sqrt(kc*a);
	e = sqrt(ec *ec + es*es);

	orb->gm = kc;    
	orb->beta = beta;
	orb->b = sqrt(fabs(beta));
	orb->peri = a*(1-e);
	orb->period = 2*PI*sqrt(a*a*a/kc);
    }

    orb->r0 = r0;
    orb->v2 = v2;
    orb->eta = eta;
    orb->zeta = orb->gm - orb->beta*r0;
    orb->x = xhel[0];
    orb->y = xhel[1];
    orb->z = xhel[2];
    orb->vx = vhel[0];
    orb->vy = vhel[1];
    orb->vz = vhel[2];
}

static void copybuffers(const struct reb_dpconst7 _a, const struct reb_dpconst7 _b, int N3){
    for (int i=0;i<N3;i++){ 
        _b.p0[i] = _a.p0[i];
        _b.p1[i] = _a.p1[i];
        _b.p2[i] = _a.p2[i];
        _b.p3[i] = _a.p3[i];
        _b.p4[i] = _a.p4[i];
        _b.p5[i] = _a.p5[i];
        _b.p6[i] = _a.p6[i];
    }
// The above code seems faster than the code below, probably due to some compiler optimizations. 
//  for (int i=0;i<7;i++){  
//      memcpy(_b[i],_a[i], sizeof(double)*N3);
//  }
}

// This is original from IAS15.
void predict_next_step(double ratio, int N3,  const struct reb_dpconst7 _e, const struct reb_dpconst7 _b, const struct reb_dpconst7 e, const struct reb_dpconst7 b){
    if (ratio>20.){
        // Do not predict if stepsize increase is very large.
        for(int k=0;k<N3;++k) {
            e.p0[k] = 0.; e.p1[k] = 0.; e.p2[k] = 0.; e.p3[k] = 0.; e.p4[k] = 0.; e.p5[k] = 0.; e.p6[k] = 0.;
            b.p0[k] = 0.; b.p1[k] = 0.; b.p2[k] = 0.; b.p3[k] = 0.; b.p4[k] = 0.; b.p5[k] = 0.; b.p6[k] = 0.;
        }
    }else{
        // Predict new B values to use at the start of the next sequence. The predicted
        // values from the last call are saved as E. The correction, BD, between the
        // actual and predicted values of B is applied in advance as a correction.
        //
        double q1 = ratio;   
        double q2 = q1 * q1;
        double q3 = q1 * q2;
        double q4 = q2 * q2;
        double q5 = q2 * q3;
        double q6 = q3 * q3;
        double q7 = q3 * q4;

        for(int k=0;k<N3;++k) {
            double be0 = _b.p0[k] - _e.p0[k];
            double be1 = _b.p1[k] - _e.p1[k];
            double be2 = _b.p2[k] - _e.p2[k];
            double be3 = _b.p3[k] - _e.p3[k];
            double be4 = _b.p4[k] - _e.p4[k];
            double be5 = _b.p5[k] - _e.p5[k];
            double be6 = _b.p6[k] - _e.p6[k];


            e.p0[k] = q1*(_b.p6[k]* 7.0 + _b.p5[k]* 6.0 + _b.p4[k]* 5.0 + _b.p3[k]* 4.0 + _b.p2[k]* 3.0 + _b.p1[k]*2.0 + _b.p0[k]);
            e.p1[k] = q2*(_b.p6[k]*21.0 + _b.p5[k]*15.0 + _b.p4[k]*10.0 + _b.p3[k]* 6.0 + _b.p2[k]* 3.0 + _b.p1[k]);
            e.p2[k] = q3*(_b.p6[k]*35.0 + _b.p5[k]*20.0 + _b.p4[k]*10.0 + _b.p3[k]* 4.0 + _b.p2[k]);
            e.p3[k] = q4*(_b.p6[k]*35.0 + _b.p5[k]*15.0 + _b.p4[k]* 5.0 + _b.p3[k]);
            e.p4[k] = q5*(_b.p6[k]*21.0 + _b.p5[k]* 6.0 + _b.p4[k]);
            e.p5[k] = q6*(_b.p6[k]* 7.0 + _b.p5[k]);
            e.p6[k] = q7* _b.p6[k];

            b.p0[k] = e.p0[k] + be0;
            b.p1[k] = e.p1[k] + be1;
            b.p2[k] = e.p2[k] + be2;
            b.p3[k] = e.p3[k] + be3;
            b.p4[k] = e.p4[k] + be4;
            b.p5[k] = e.p5[k] + be5;
            b.p6[k] = e.p6[k] + be6;
        }
    }
}

void converst(State s, double *x, double *v)
{
  x[0] = s.x;
  x[1] = s.y;
  x[2] = s.z;
  v[0] = s.xd;
  v[1] = s.yd;
  v[2] = s.zd;
}

// Calculate the reference and other quantites at the 7 substeps
void predict_kepler(orbstruct* restrict const orb, 
		    double dt, double del0,
		    double **xkep, double **vkep,
		    double **csxhtemp, double **csvhtemp,
		    double **xtwobod3, double **xtwobod2o,		
		    int n,
		    int si)
{
    double xo[NDIM],vo[NDIM],delx[NDIM],delv[NDIM];

    for(int j=0; j<n; j++){
	// This part of the loop could be a separate function.
	for(int i=1;i<si+1; i++){
	    double delt = dt*h[i] + del0;
	    kepler_stepxv_no_reset_prime(delx,delv,xo,vo,delt,orb[j]);

	    double xtwobod2 = 0.;
	    for(int k=0; k<NDIM; k++){

		csxhtemp[i][3*j+k] = orb[j].csxh[k];
		csvhtemp[i][3*j+k] = orb[j].csvh[k];
		
		add_cs(&(xo[k]), &(csxhtemp[i][3*j+k]), delx[k]);
		add_cs(&(vo[k]), &(csvhtemp[i][3*j+k]), delv[k]);

		xkep[i][3*j+k] = xo[k];
		vkep[i][3*j+k] = vo[k];

		xtwobod2 += xkep[i][3*j+k] * xkep[i][3*j+k];
	    }
	    xtwobod2o[i][j] = xtwobod2;
	    xtwobod3[i][j] = xtwobod2 * sqrt(xtwobod2);
	}
    }
}


// Calculated the perturbed acceleration, with compensated sums for
// all the bodies at substep n0.
// 
void encke_hh_acceleration(double *apert, double *atc, 
			   double *apertr, double *atcr,
			   double *xpert,		     
			   double Msun,
			   double *m,
			   double *xkep, double *vkep,
			   double *xtwobod2in, double *xtwobod3in,
			   double *csxhtemp, double *csvhtemp,		     
			   double *factor1keep,
			   int n)
{

    // These should be allocated.
    double xtrue[NMAX];
    double rij[NDIM];
    double factor1[NMAX];
    double csxh[NMAX];

    int N3 = 3*n;

    // initialization
    for(int i=0; i<N3; i++){
	apert[i] = 0.;
	atc[i] = 0.;
	csxh[i] = csxhtemp[i];
    }

    // Compute true positions and related quantities
    for(int i=0; i<n; i++){    
	double xtrue2 = 0;
	for(int k=0; k<NDIM; k++){	
	    xtrue[3*i+k] = xkep[3*i+k];
	    add_cs(&(xtrue[3*i+k]), &(csxh[3*i+k]), xpert[3*i+k]);	    
	    xtrue2 += xtrue[3*i+k] * xtrue[3*i+k];
	}
	double xtrue3 = xtrue2 * sqrt(xtrue2);

	double fac = GNEWT*m[i]/xtrue3;
	for(int k=0; k<NDIM; k++){		
	    factor1[3*i+k] = xtrue[3*i+k]*fac; // indirect terms
	    factor1keep[3*i+k] = factor1[3*i+k];
	}
    }

    // Double loop over all the pairs of non-sun masses
    // This part could be parallelized
    for(int i=0; i<n; i++){            
	for(int j=i+1;j<n;j++){
	    double r2 = 0.;
	    for(int k=0; k<NDIM; k++){
		rij[k] = xtrue[3*i+k] - xtrue[3*j+k];
		r2 += rij[k] * rij[k];
	    }
	    double r3 = r2 * sqrt(r2);
	    for(int k=0; k<NDIM; k++){
		double factor = GNEWT*rij[k]/r3;
		add_cs(&(apert[3*i+k]), &(atc[3*i+k]), -factor * m[j]);  // direct terms
		add_cs(&(apert[3*i+k]), &(atc[3*i+k]), -factor1[3*j+k]); // indirect terms
		add_cs(&(apert[3*j+k]), &(atc[3*j+k]), factor * m[i]);   // direct terms
		add_cs(&(apert[3*j+k]), &(atc[3*j+k]), -factor1[3*i+k]); // indirect terms
	    }
	}
    }

    // Save the accelerations and compensated sums, in case a reference orbit
    // needs to be updated.
    if(apertr != NULL && atcr != NULL){
	for(int i=0; i<N3; i++){
	    apertr[i] = apert[i];
	    atcr[i] = atc[i];
	}
    }

    // Encke terms.  This part is set to zero when a reference orbit
    // is updated.  
    for(int i=0; i<n; i++){                
	double q1t = 0.;
	double q2t = 0.;
	for(int k=0; k<NDIM; k++){	
	    q1t += (xpert[3*i+k] * xpert[3*i+k]);
	    q2t += (xpert[3*i+k] * xkep[3*i+k]);	    
	}

	double qt = (q1t + 2.0*q2t)/xtwobod2in[i];
	double oneplusq = 1.0 + qt;
	double quant = oneplusq * sqrt(oneplusq);
	double ft = qt*(3.0 + 3.0*qt + qt*qt)/(quant*quant + quant);

	double h = GNEWT*(Msun + m[i])/xtwobod3in[i];
	for(int k=0; k<NDIM; k++){	
	    double temp = h*(-xpert[3*i+k] + xtrue[3*i+k] * ft);
	    add_cs(&(apert[3*i+k]), &(atc[3*i+k]), temp);
	}
    }
	
}

// Calculate the accelation on body ikeep from the non-sun masses, along with
// the associated compensated sums. This is only called within update_twobod.
//
void encke_hh_acceleration_i(double *apert, double *atc, double *xtrue, double *m, int n, int ikeep)
{
    double factor1[NMAX];
    double rij[NDIM];

    int i = ikeep;

    for(int k=0; k<NDIM; k++){
	apert[k] = 0.;
	atc[k] = 0.;
    }

    for(int j=0; j<n; j++){  
	double xtrue2 = 0.;
	for(int k=0; k<NDIM; k++){      
	    xtrue2 += xtrue[3*j+k] * xtrue[3*j+k];
	}
	double xtrue3 = xtrue2 * sqrt(xtrue2);
	//double fac = GNEWT*m[j]/xtrue3;
	for(int k=0; k<NDIM; k++){            
	    factor1[3*j+k] = xtrue[3*j+k]/xtrue3*GNEWT*m[j];
	    //factor1[3*j+k] = xtrue[3*j+k]*fac; // indirect terms	    
	}
    }

    for(int j=0; j<n; j++){    
	if(j == i) continue;
	double r2 = 0.;
	for(int k=0; k<NDIM; k++){          
	    rij[k] = xtrue[3*i+k] - xtrue[3*j+k];
	    r2 += rij[k] * rij[k];
	}
	double r3 = r2*sqrt(r2);
	for(int k=0; k<NDIM; k++){                
	    double factor = GNEWT*rij[k]/r3;
	    add_cs(&(apert[k]), &(atc[k]), -factor*m[j]);    // indirect
	    add_cs(&(apert[k]), &(atc[k]), -factor1[3*j+k]); // direct
	}
    }
}


// Kepler stepper for calculating the reference trajectories
void kepler_stepxv_no_reset_prime(double *delx, double *delv,
				  double *xo, double *vo,
				  double delt, orbstruct orb)
{

    void kepler_step_depth(double kc, double dt, double beta, double b, 
			   State *s0, State *s, int depth, 
			   double r0, double v2, double eta, double zeta);
    
    State s0,s;
    double deltin = delt;

    s0.x = orb.x;
    s0.y = orb.y;
    s0.z = orb.z;
    s0.xd = orb.vx;
    s0.yd = orb.vy;
    s0.zd = orb.vz;
    if(deltin > orb.period){
	deltin = mfmod(deltin, orb.period);
    }

    kepler_step_depth(orb.gm, deltin, orb.beta, orb.b, &s0, &s, 0, orb.r0, orb.v2, orb.eta, orb.zeta);
    converst(s,delx,delv);
    converst(s0,xo,vo);

}

double mfmod(double x,double y) { 
    double a1, a2;
    a1 = x/y;
    a2 = floor(a1); 
    return y*(a1 - a2); 
}

// Calculate time step adaptively
void adjust_step(double **xkep, double **x3,
	    double *apert,
	    double dt_done, double *dt_new,
	    double epsilon, int si, int nprog,
	    double **factor1keep)
{
    
    int N3 = nprog*3;
    static int n3_alloc = 0;

    double a0[NMAX],at2v[NMAX],gravitycsv[NMAX],csa0[NMAX];
    double maxak=0.,maxb6ktmp=0.;

    static struct reb_dp7 g;
    static struct reb_dp7 b;
    static struct reb_dp7 csb;         ///< Compensated summation for b

    if (N3 > n3_alloc){
        realloc_dp7(&g,N3);
        realloc_dp7(&b,N3);
        realloc_dp7(&csb,N3);
	n3_alloc = N3;
    }
    
    for(int k=0;k<N3;++k) {
	b.p0[k] = 0.; b.p1[k] = 0.; b.p2[k] = 0.; 
	b.p3[k] = 0.; b.p4[k] = 0.; b.p5[k] = 0.; b.p6[k] = 0.;
	g.p0[k] = 0.; g.p1[k] = 0.; g.p2[k] = 0.; 
	g.p3[k] = 0.; g.p4[k] = 0.; g.p5[k] = 0.; g.p6[k] = 0.;
	gravitycsv[k] = 0.; csa0[k] = 0.;
	csb.p0[k] = 0.; csb.p1[k] = 0.; csb.p2[k] = 0.; 
	csb.p3[k] = 0.; csb.p4[k] = 0.; csb.p5[k] = 0.; csb.p6[k] = 0.;
    }

    //accelerations at time 0
    for(int i=0; i<N3; i++){
	a0[i] = factor1keep[si][i];
    }

    //intermediate accelerations
    for(int n=1;n<si;n++){
	for(int i=0; i<nprog; i++){	
	    int k0 = 3*i;
	    int k1 = 3*i+1;
	    int k2 = 3*i+2;
	    at2v[k0] = factor1keep[n][k0];
	    at2v[k1] = factor1keep[n][k1];
	    at2v[k2] = factor1keep[n][k2];
	}

	switch (n) {                  // Improve b and g values
	case 1:
	    for(int k=0;k<N3;++k) {
		double gk = at2v[k];
		double gk_cs = gravitycsv[k];
		add_cs(&gk, &gk_cs, -a0[k]);
		add_cs(&gk, &gk_cs, csa0[k]);
		g.p0[k]  = gk/rr[0];
	    } 
	    break;
	case 2:
	    for(int k=0;k<N3;++k) {
		double gk = at2v[k];
		double gk_cs = gravitycsv[k];
		add_cs(&gk, &gk_cs, -a0[k]);
		add_cs(&gk, &gk_cs, csa0[k]);
		g.p1[k] = (gk/rr[1] - g.p0[k])/rr[2];
	    } 
	    break;
	case 3:
	    for(int k=0;k<N3;++k) {
		double gk = at2v[k];
		double gk_cs = gravitycsv[k];
		add_cs(&gk, &gk_cs, -a0[k]);
		add_cs(&gk, &gk_cs, csa0[k]);
		g.p2[k] = ((gk/rr[3] - g.p0[k])/rr[4] - g.p1[k])/rr[5];
	    } break;
	case 4:
	    for(int k=0;k<N3;++k) {
		double gk = at2v[k];
		double gk_cs = gravitycsv[k];
		add_cs(&gk, &gk_cs, -a0[k]);
		add_cs(&gk, &gk_cs, csa0[k]);
		g.p3[k] = (((gk/rr[6] - g.p0[k])/rr[7] - g.p1[k])/rr[8] - g.p2[k])/rr[9];
	    } break;
	case 5:
	    for(int k=0;k<N3;++k) {
		double gk = at2v[k];
		double gk_cs = gravitycsv[k];
		add_cs(&gk, &gk_cs, -a0[k]);
		add_cs(&gk, &gk_cs, csa0[k]);
		g.p4[k] = ((((gk/rr[10] - g.p0[k])/rr[11] - g.p1[k])/rr[12] - g.p2[k])/rr[13] - g.p3[k])/rr[14];
	    } break;
	case 6:
	    for(int k=0;k<N3;++k) {
		double gk = at2v[k];
		double gk_cs = gravitycsv[k];
		add_cs(&gk, &gk_cs, -a0[k]);
		add_cs(&gk, &gk_cs, csa0[k]);
		g.p5[k] = (((((gk/rr[15] - g.p0[k])/rr[16] - g.p1[k])/rr[17] - g.p2[k])/rr[18] - g.p3[k])/rr[19] - g.p4[k])/rr[20];
	    } break;
	case 7:
	    {
		maxak = 0.0;
		maxb6ktmp = 0.0;
		for(int k=0;k<N3;++k) {
		    double tmp = g.p6[k];
		    double gk = at2v[k];
		    double gk_cs = gravitycsv[k];
		    add_cs(&gk, &gk_cs, -a0[k]);
		    add_cs(&gk, &gk_cs, csa0[k]);
		    g.p6[k] = ((((((gk/rr[21] - g.p0[k])/rr[22] - g.p1[k])/rr[23] - g.p2[k])/rr[24] - g.p3[k])/rr[25] - g.p4[k])/rr[26] - g.p5[k])/rr[27];
		    tmp = g.p6[k] - tmp;    
		    add_cs(&(b.p6[k]), &(csb.p6[k]), tmp); 

		    double atot = fabs(at2v[k]); //add these w/ compensated?
		    double b6t = fabs(b.p6[k]);
		    if(atot > maxak) maxak = atot;
		    double temp1 = b6t;
		    if(temp1 > maxb6ktmp) maxb6ktmp = temp1;
		}
		break;
	    }
	}
    }
    double integrator_error = maxb6ktmp/maxak;
    double dttemp = sqrt7(epsilon/integrator_error)*dt_done;
    *dt_new =dttemp;
      
}



// Do nothing here. This is only used in a leapfrog-like DKD integrator. IAS15 performs one complete timestep.
void reb_integrator_encke_hh_part1(struct reb_simulation* r){
    r->gravity_ignore_terms = 0;
}

void reb_integrator_encke_hh_part2(struct reb_simulation* r){
#ifdef GENERATE_CONSTANTS
    integrator_generate_constants();
#endif  // GENERATE_CONSTANTS
    // Try until a step was successful.
    while(!reb_integrator_encke_hh_step(r));
}

void reb_integrator_encke_hh_synchronize(struct reb_simulation* r){
}
