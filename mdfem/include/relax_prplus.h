// Last Modified : Tue Mar 17 00:54:40 2009

//////////////////////////////////////////////////////
// Stand alone conjugate gradient relaxiation program
//
//    transplanted from LAMMPS min_cg.cpp
//
//  writen by Keonwook Kang at Stanford, Mar 16 2009
//
//////////////////////////////////////////////////////

#ifndef _RELAX_PRPLUS_H
#define _RELAX_PRPLUS_H

/* Special Interface to potential_wrapper in md.h md2d.h */
//extern class MDFrame *__conj_simframe;
//extern class SIHMDFrame *__conj_sihmdframe;
//extern class SimFrame2d *__conj_simframe2d;
//extern class Morse2d *__conj_morse2d;
//extern class SpFrame2d *__conj_spframe2d;
/* End of Special Interface */

#ifdef __cplusplus
#ifdef __gps
void CGRelax_PRplus(void (*func)(int,double*,double *, double*),
             int n,double etol,double ftol,int max_iter,int max_eval,double maxdist,
             double x[],double g[],double *f, double *buffer)
#else
extern "C" void CGRelax_PRplus(void (*func)(int,double*,double *, double*),
                        int n,double etol,double ftol,int max_iter,int max_eval,double maxdist,
            double x[],double g[],double *f, double *buffer);
#endif
#endif

#endif // _RELAX_PRPLUS_H


