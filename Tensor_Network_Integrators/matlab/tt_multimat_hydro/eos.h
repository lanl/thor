#ifndef _EOS_
#define _EOS_

#ifdef __cplusplus
extern "C" {
#endif

#include "globals.h"
#include "minip.h" 

void sspd_solid(int which_solid, double rho, double *cs);
void p_solid(int which_solid, double rho, double ei, double *p);
void e_solid(int which_solid, double rho, double p, double *e);

void rho_polynominal(double p, double *rho); 

#ifdef __cplusplus
}
#endif
#endif
