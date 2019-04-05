// Last Modified : Sat Aug 13 21:58:13 2005

//////////////////////////////////////////////////////
// Stand alone conjugate gradient relaxiation program
//
//    Based on IMSL's U2CGG.F
//
//  writen by Dongyi Liao at MIT 1999
//
//////////////////////////////////////////////////////

#ifndef _RELAX_ZXCGR_H
#define _RELAX_ZXCGR_H

/* Special Interface to potential_wrapper in md.h md2d.h */
extern class MDFrame *__conj_simframe;
//extern class SIHMDFrame *__conj_sihmdframe;
//extern class SimFrame2d *__conj_simframe2d;
//extern class Morse2d *__conj_morse2d;
//extern class SpFrame2d *__conj_spframe2d;
/* End of Special Interface */

#ifdef __cplusplus
#ifdef __gps
void CGRelax(void (*func)(int,double*,double *, double*),
             int n,double acc,int maxfn,double dfpred,
             double x[],double g[],double *f, double *buffer);
#else
extern "C" void CGRelax(void (*func)(int,double*,double *, double*),
            int n,double acc,int maxfn,double dfpred,
            double x[],double g[],double *f, double *buffer);
#endif
#endif

#endif // _RELAX_ZXCGR_H


