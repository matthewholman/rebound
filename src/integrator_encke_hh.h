#include <math.h>
#include <stdio.h>

#define PI 3.14159265358979323846
/*au, day, solar mass*/
#define GNEWT 0.01720209895*0.01720209895
#define YEAR 1.0
#define NMAX 50
#define NMAX2 2000
#define NDIM 3
#define NSTEP 9

#ifndef _STATE_
#define _STATE_
typedef struct {
  double x, y, z, xd, yd, zd, xdd, ydd, zdd;
} State;
#endif

void update_two_bod(struct reb_simulation* r, int *flg, int *rcount);

#ifdef __cplusplus

// At least GCC and clang support the restrict keyword as an extension.
#if defined(__GNUC__) || defined(__clang__)

#define REBOUND_RESTRICT __restrict__

#else

// For other compilers, we disable it.
#define REBOUND_RESTRICT

#endif

extern "C" {

#else

#define REBOUND_RESTRICT restrict

#endif

/**
 * @brief Struct containing pointers to intermediate values
 */
struct reb_dpconst7 {
    double* const restrict p0;  ///< Temporary values at intermediate step 0 
    double* const restrict p1;  ///< Temporary values at intermediate step 1 
    double* const restrict p2;  ///< Temporary values at intermediate step 2 
    double* const restrict p3;  ///< Temporary values at intermediate step 3 
    double* const restrict p4;  ///< Temporary values at intermediate step 4 
    double* const restrict p5;  ///< Temporary values at intermediate step 5 
    double* const restrict p6;  ///< Temporary values at intermediate step 6 
};

int reb_integrator_encke_hh_step(struct reb_simulation* r);

void predict_next_step(double ratio, int n3,  const struct reb_dpconst7 _e, const struct reb_dpconst7 _b, const struct reb_dpconst7 e, const struct reb_dpconst7 b);

void predict_kepler(orbstruct* restrict const orb, //struct reb_simulation* r,
		double dt, double del0,
		double **xkep,
		double **vkep,
		double **csxhtemp, double **csvhtemp,
		double **xtwobod3, double **xtwobod2o,		
		int n,
		int si);
    
void encke_hh_acceleration(double *apert, double *atc,
		     double *apertr, double *atcr,
		     double *xpert,
		     double Msun,
		     double *m,
		     double *xkeppred, double *vkeppred,
		     double *xtwobod2in, double *xtwobod3in,
		     double *cstempx, double *cstempv,		     
		     double *factor1keep,
		     int n);

void encke_hh_acceleration_i(double *apert, double *atc, double *xtrue, double *m, int n, int ikeep);
    
void kepler_stepxv_no_reset_prime(double *delx, double *delv,
				  double *xo, double *vo,
				  double delt, orbstruct orb);

void adjust_step(double **xkep, double **x3,
	    double *apert,
	    double dt_done, double *dt_new,
	    double epsilon, int si, int nprog,
	    double **factor1keep);

double mfmod(double x,double y);

static inline void add_cs(double* p, double* csp, double inp);
void add_cs2_long(long double* p, long double* csp, long double inp);

void converst(State s, double *x, double *v);

void getorb(double Msun, double m, double *xhel, double *vhel, int reset, orbstruct *orb);    

#ifndef _INTEGRATOR_ENCKE_HH_H
#define _INTEGRATOR_ENCKE_HH_H
void reb_integrator_encke_hh_part1(struct reb_simulation* r);              ///< Internal function used to call a specific integrator
void reb_integrator_encke_hh_part2(struct reb_simulation* r);              ///< Internal function used to call a specific integrator    
void reb_integrator_encke_hh_synchronize(struct reb_simulation* r);        ///< Internal function used to call a specific integrator
void reb_integrator_encke_hh_clear(struct reb_simulation* r);              ///< Internal function used to call a specific integrator
void reb_integrator_encke_hh_alloc(struct reb_simulation* r);              ///< Internal function, alloctes memory for ENCKE
void reb_integrator_encke_hh_reset(struct reb_simulation* r);          ///< Internal function used to call a specific integrator
#endif
    
