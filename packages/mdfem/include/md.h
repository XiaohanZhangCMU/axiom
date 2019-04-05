/*
  md.h
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Mon Jan  3 18:14:20 2011

  FUNCTION  :  Easy-to-use MD simulation Framefork

  Featuring :  1. Scripting input
               2. X-window display
               3. automatic simulation log
               4. convenient configuration and output file handle
               5. perfect lattice and dislocation config creator
               6. conjugate-gradient relaxation
               7. NVT MD simulation
               8. NPH and NPT simulation (Parrinello-Rahman + Nose-Hoover)
               9. Cell-Verlist combined neighbor list
               
  This is a single CPU code.
  Previous shared-memory option removed.
  May implment MPI parallelization in the future.
*/


#ifndef _MD_H
#define _MD_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <vector>

#ifdef _USEPY
#include "Python.h"
#endif

#ifdef _USE_R
#include <R.h>
#endif

#ifdef _USEOCTAVE
/* The following definition is intended for switching "OCTINTERP_API" to blank 
 in the header files of octave.h, oct.h, and parse.h. */
#define OCTINTERP_API
#include <octave/octave.h>
#include <octave/oct.h>
#include <octave/parse.h>
#endif

#ifdef _USEBOOST
#include "boost/multi_array.hpp"
#endif

#ifdef _USEFFTW
#include "fftw3.h"
#else
#define fftw_plan int
#endif

#ifdef _USEHDF5
#include "hdf5.h"
#endif

#ifdef _GSL
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_coupling.h>
#endif

#include "general.h"
#include "organizer.h"
#include "display.h"
#include "linalg3.h"
#include "relax_zxcgr.h"
#include "relax_prplus.h"

/* Internal units
   Energy:  eV
   Length:  A
   Stress:  eV/A^3

   Stored Internal Variables
   Position:     r               (in A)
   Velocity:     v*TIMESTEP      (in A)
   Accelaration: a*TIMESTEP^2/2  (in A)
*/

/* Physical constants */
#define AVO  6.0221367E+23
#define EV   1.6021772E-19     /* (J)    */
#define BOLZ 1.380658E-23      /* (J/K) BOLZ=KB*EV */
#define KB   0.8617336E-4      /* (eV/K) */
#define PLANCK_H 4.135667e-15  /* eV/s */

/* mass unit conversion from (gram/mol) to (eV/(A/ps)^2) */
#define MASSCONVERT   1e-3/AVO/EV/1e20*1e24

/* physical constants for Ewald summation */
#define P_CLMB_E 1.0 // for the correction of Coulomb potential (dielectric function constant) 
#define P_CLMB (P_E*P_E/P_E*1.0e10)/(4.0*M_PI*P_EPSILON0*P_CLMB_E) // (eV*A) conversion factor for the Coulomb force
#define P_SQRT_CLMB sqrt(P_CLMB) //3.79470135954879
#define M_SQRT_PI sqrt(M_PI)     //1.77245385090552


//#ifndef max
//#define max(a,b) (((a)>(b))?(a):(b))
//#endif
//#ifndef min
//#define min(a,b) (((a)<(b))?(a):(b))
//#endif

#ifndef SQR
#define SQR(A) ((A)*(A))
#endif

#ifndef CUBE
#define CUBE(A) ((A)*(A)*(A))
#endif

#ifndef POW6
#define POW6(A) SQR(CUBE(A))
#endif

#ifndef POW12
#define POW12(A) SQR(POW6(A))
#endif

#define PME_int_mod(x,y) (((x%y)+y)%y) /* Particle Mesh Ewald: mod function */

#define MAXSPECIES 10          /* maximum number of species */
#define MAXCOLORS  10          /* maximum number of colors in display */

#ifndef MAXCONSTRAINATOMS
#define MAXCONSTRAINATOMS 50001 /* maximum number of atoms in constrained minimization */
#endif

#define MAXPAIRCONSTRAINTS 3001 /* maximum number of pair displacement constraints */
#define MAXCMDLEN  1000        /* maximum length of command string */
#define MAXNHCLEN  20          /* maximum length of NH chain */
#define MAXNLSKIP  20000        /* maximum length of nl_skip_pairs array */
#define MAXINPUTLEN 20000
#define MAXOUTPUTDATLEN 2000
#define MAXOUTPUTSTRLEN 10000
#define MAXFILENAMELEN 200

#define NEBSPECSIZE 20
#define ANNEALSPECSIZE 10
#define SURFTENSIONSPECSIZE 10
#define TORQUESPECSIZE 3
#define BENDSPECSIZE 4
#define EXTFORCESIZE 100
#define ENERGYTHRESHOLDSIZE 10

/* Random number generator for platforms not supporting drand48 (e.g. cygwin) */
#ifdef _NODRAND48
#define drand48 (1.0/RAND_MAX)*rand
#endif

//shared memory MP deprecated replaced by MPI
//#define MP_Sync() {}
//#define MP_BeginMasterOnly() {}
//#define MP_EndMasterOnly() {}
//#define MP_Free(s) free(s)
//#define MP_cpu() 0
//#define MP_ncpu() 1

#ifndef _USE_R
#define Realloc(p,t,s) {p=(t *)realloc(p,sizeof(t)*s);}
#else
#define Realloc(p,t,s) {if(p==NULL) p=R_Calloc(s,t); else R_Realloc(p,s,t);}
#endif

#define bindcommand_sync bindcommand



/* Input/output files */
class CNFile : public AUXFile /* configuration(cn) file */
{
public:
    CNFile(int t):AUXFile(t){};
    virtual char * describe();
private:
    virtual int writeblock(void *p);
    virtual int readblock(void *p);
};

class LAMMPSFile : public AUXFile /* LAMMPS file */
{
public:
    LAMMPSFile(int t):AUXFile(t){};
    virtual char * describe();
private:
    virtual int writeblock(void *p);
    virtual int readblock(void *p);
};

class MDCASKconFile : public AUXFile /* MDCASK configuration file */
{
public:
    MDCASKconFile(int t):AUXFile(t){};
    virtual char * describe();
private:
    virtual int writeblock(void *p);
    virtual int readblock(void *p);
};

class MDCASKinputFile : public AUXFile /* MDCASK input file */
{
public:
    MDCASKinputFile(int t):AUXFile(t){};
    virtual char * describe();
private:
    virtual int writeblock(void *p);
    virtual int readblock(void *p);
};

class PropFile : public AUXFile/* Output property file */
{
public:
    PropFile(int t):AUXFile(t){};
    virtual char * describe();
    virtual int writeentry(void *p);
};


//------------------------------------------------------------------- 
// Define EComplex Class (will be used in Classic Ewald summation)
//-------------------------------------------------------------------
struct EComplex
{
    double Re, Im; //Real and imaginary part (only data member)
public:
    double Norm() { return sqrt(Re*Re+Im*Im); }
    double Norm2() { return Re*Re+Im*Im; }
    
    static void Mul(const EComplex &a, const EComplex &b, EComplex &c)//c=a*b
    {
        double tmp=a.Re*b.Re-a.Im*b.Im;
        c.Im=a.Re*b.Im+a.Im*b.Re; 
        c.Re=tmp;
    }
    static void MulConjugate(const EComplex &a, const EComplex &b, EComplex &c)
    {
        c.Re=a.Re*b.Re+a.Im*b.Im; 
        c.Im=-a.Re*b.Im+a.Im*b.Re;
    }
    EComplex &operator *=(double d)
    {
        Re*=d;
        Im*=d;
        return *this;
    }
    EComplex &operator +=(EComplex &c)
    {
        Re+=c.Re;
        Im+=c.Im;
        return *this;
    }
    EComplex &operator -=(EComplex &c)
    {
        Re-=c.Re;
        Im-=c.Im;
        return *this;
    }
    void ExpI(double theta)
    { 
        Re=cos(theta);
        Im=sin(theta);
    }
    EComplex &operator =(double d) 
    { 
        Re=d;
        Im=0.0;
        return *this;
    }
    friend LOStream &operator <<(LOStream &os, EComplex &c)
    {
        return os << "["<<c.Re<<"(+)"<<c.Im<<"(I)]";
    }
};

/*-------------------------------------------------------------------
 * Define Complex Array Type (will be used in PME)
 *------------------------------------------------------------------- */
typedef double ComplexType[2];

/* MD simulation framework */
class MDFrame : public Organizer
{
public:
    int _SAVEMEMORY;    /* controls how much memory is allocated per atom */
    int _NP;            /* total number of atoms */
    int _NP0, _NIMAGES; /* number of original and image atoms */
    int _NPfree,_NPfixed; /* number of free and fixe atoms */
    int _NP_sp[MAXSPECIES]; /* numober of atoms per species */

    class Vector3 *_R, *_SR;   /* real and scaled coordinates of atoms */
    class Vector3 *_SRA; /* running average coordinates of atoms in scaled coordinate */
    class Vector3 *_VR, *_VSR; /* real and scaled velocities of atoms */
    class Vector3 *_F, *_F0;   /* force on atoms */
    class Vector3 *_R0;        /* place to save old coordinates of atoms */
    class Vector3 *_Fext;      /* external force on atoms */
    int _ENABLE_FEXT;

    /* parameter for external flat indentor */ 
    double _FLAT_INDENTOR_POS, _FLAT_INDENTOR_POS0, _FLAT_INDENTOR_V, _FLAT_INDENTOR_K;
    double _FLAT_INDENTOR_F;   /* total force exerted by indentor */
    int    _FLAT_INDENTOR_DIR, _ENABLE_FLAT_INDENTOR;

    int *fixed;                /* 1 if atom is fixed */
    int *image;                /* image atom */
    int *species;              /* species type of atom */
    int *group;                /* group ID */
    int nspecies;              /* number of atom species */
    char element[MAXSPECIES][10]; /* species names */

    class Matrix33 _H, _H0;     /* simulation box vectors and copy */
    class Matrix33 _VIRIAL, _VH;/* Virial stress, and velocity of H */
    class Matrix33 *_VIRIAL_IND;/* Virial stress contribution from each atom */
    double _EPOT, *_EPOT_IND;   /* total potential energy and contribution per atom */
    double *_EPOT_RMV;          /* potential energy change if atom is removed */
    double _EPOT_BOX;           /* potential energy of the box */
    double _EPOT0;              /* potential energy from previous step */
    double _EPOT0_ext;          /* reference potential energy from boundary atoms */
    double *_TOPOL;             /* atomic environment, e.g. central symmetry parameter */
    
    /* Molecular Dynamics parameters */
    double _ESTRAIN;             /* box strain energy (in eV) */
    double _OMEGA;               /* box volume (in A^3) */
    double _KATOM, _KBOX;        /* kinetic energies (in eV) */
    double _ATOMMASS[MAXSPECIES], _WALLMASS; /* in (g/mol) */
    double _ATOMCHARGE[MAXSPECIES]; /* in electron */
    double _TIMESTEP;            /* in picosecond (ps) */
    double _TDES,_T;             /* desired and current temperature */
    int _DOUBLET;                /* 1: if initialize T double desired temperature */
    double _ATOMTCPL,_BOXTCPL;   /* temperature coupling constants */
    double _BOXDAMP;             /* damping coefficient for box PR method to prevent oscillation */
    double zeta,zetav,zetaa,vt2; /* nose hoover thermal stat */
    double _HELM,_HELMP;         /* Helmholtz free energy and a conserved var */
//    int usenosehoover;           /* 1: if use nose hoover thermostat */
//    int usescalevelocity;        /* 1: if scale velocity for temperature control */
    int applyshear;              /* apply constant shear deformation */
    double extforce[EXTFORCESIZE];        /* external forces on different atom groups */
    double forcemul;             /* a multiplication factor for extforce */
    Matrix33 _SHEARRATE;         /* shear deformation rate, if applyshear==1 */
    int runavgposition;
    int _ENABLE_DPD;
    double _DPD_fc;
    double _DPD_RATIO;
    double _DPD_RCUT;
    int _ENABLE_ANDERSON;
    double _ANDERSON_RATIO;
    
    /* Neighbor list */
    double _RLIST,_SKIN;        /* potential cut off radius (in A) */
    int _NIC, _current_NIC;     /* maximum allowable and actual number of atoms per cell */
    int _NNM, _current_NNM;     /* maximum allowable and actual number of neighbors per atom */

#ifndef _USEBOOST
    int *nn;                    /* nn[i]: number of neighbors of atom i */
    int **nindex;               /* nindex[i][j]: the jth neighbor of atom i */
    char *nindex_mem;
    int ****celllist;
    char *cell_mem;
#else
    boost::multi_array<int, 1> nn;
    boost::multi_array<int, 2> nindex;
    boost::multi_array<int, 4> celllist;
#endif

    int *link_head, *link_list; /* link list */
    
    bool firsttime;

    double rc_plot;             /* another neighbor list for visualization */
    int _NNM_plot;
    int *nn_plot;
    int **nindex_plot;
    char *nindex_plot_mem;
    int nl_skip_pairs[MAXNLSKIP];    /* skip certain pairs of atoms when contructing neighbor list */

    
    /* Numerical integrator */
    Vector3 *s2, *s3, *s4, *s5, *s2real; /* Gear 6th predictor-corrector */
    class Matrix33 h2,h3,h4,h5,h2real;
    double zeta2,zeta3,zeta4,zeta5;      /* for Nose-Hoover thermal stat */
    int npold;
    int NHChainLen;                      /* Length Nose-Hoover Chain */
    double zetaNHC[MAXNHCLEN], zetaNHCv[MAXNHCLEN], zetaNHCa[MAXNHCLEN],
           zetaNHC2[MAXNHCLEN], zetaNHC3[MAXNHCLEN], zetaNHC4[MAXNHCLEN], zetaNHC5[MAXNHCLEN];
    double NHMass[MAXNHCLEN];            /* Masses of Nose-Hoover Chain */

                                        /* Nose Hoover chain for box */
//    double zetaBNHC[MAXNHCLEN], zetaBNHCv[MAXNHCLEN], zetaBNHCa[MAXNHCLEN],
//           zetaBNHC2[MAXNHCLEN], zetaBNHC3[MAXNHCLEN], zetaBNHC4[MAXNHCLEN], zetaBNHC5[MAXNHCLEN];    
//    double BNHMass[MAXNHCLEN];         
 
    /* Internal and applied stresses */
    class Matrix33 _PSTRESS;    /* stress due to momenta (in eV/A^3) */
    class Matrix33 _TOTSTRESS;  /* total stress (_PSTRESS+_VIRIAL) */
    class Matrix33 _TOTSTRESSinMPa; /* total stress in MPa */
    class Matrix33 _EXTSTRESS;  /* external applied stress */
    double _TOTPRESSURE;
    double _EXTSTRESSMUL;       /* multiplier of applied stress */
    double _EXTPRESSADD;        /* addition of hydrostatic pressure, not affected by ESTSTRESSMUL */
    double _VACUUMRATIO;        /* stress correction factor when medium is smaller than simulation box */
    class Matrix33 _DVIRIAL_Dexx; /* for computing elastic constants at finite temperature */
    double _DVIRIAL[3][3][3][3];  /* the Born elastic constant C_{ijkl}^B*(Volume) */
    double _Dexx;
    
    /* Conjugate gradient relaxation */   
    int conj_itmax,conj_fevalmax,conj_fixbox,conj_fixshape,
        conj_fixatoms,conj_step,conj_fixdir[3];
    class Matrix33 conj_fixboxvec;
    double conj_dfpred,conj_f,conj_etol,conj_ftol,conj_g2res;
    double *conj_space;
    double _HPRECOND;                     /* preconditioner for H */
    double _XPRECOND,_YPRECOND,_ZPRECOND; /* preconditioner for x,y,z */
    class Matrix33 _SIGMA;                /* stress term in Parilleno-Raman */
    class Matrix33 _GH;                   /* Energy derivative dE_pot / dH */

    /* Conjugate gradient relaxation with constraints */
    double constrainS, constrainS_inst, constrain_dist, constrainF;
    int constrainatoms[MAXCONSTRAINATOMS];
    double pair_displacement_constraints[MAXPAIRCONSTRAINTS];
    char relaxation_algorithm[100]; 

    /* constrained MD*/
    int _constrainedMD;
    
    /* Molecular Dynamics simulation */
    int totalsteps, equilsteps, curstep, step0, continue_curstep;  /* simulation steps */
    char ensemble_type[100], integrator_type[100], initvelocity_type[10], zerorot[4];
    int implementation_type, algorithm_id;
    
    /* Monte Carlo simulation */
    int MC_atom, MC_accept, MC_accept_tot, MC_switchoffatoms[10];
    int MC_dV_freq, MC_dN_freq;
    double MC_accept_ratio, MC_dVmax, MC_dfrac, MC_P_ext, MC_mu_ext, MC_Lambda3, MC_fracatom;
    Vector3 MC_dr;
    double *_EPOT_IND1, *_EPOT_INDUMB;
    Vector3 *_SR_MC_SAV,*_SR_MC_SAVUMB;
    double MC_maxd;


    /* Free energy calculations */
    double _LAMBDA, _LAMBDA0, _LAMBDA1, dEdlambda, dlambdadt, _ECOH, _WTOT, _WAVG;
    int _ENABLE_SWITCH, _SWITCHFREQ, _RANDSEED, _CHAINLENGTH, _SWITCHFUNC, _REFPOTENTIAL;
    int heslen, *hesind; double *hes; /* Hessian matrix */
    double _ECSPRING;                 /* spring constant of Einstein crystal */
    double _I12_RC, _I12_SIGMA, _I12_EPSILON; /* Inverse r12 potential */
    double _GAUSS_RC, _GAUSS_SIGMA, _GAUSS_EPSILON; /* Gaussian potential */

    /* Nudged-Elastic-Band (NEB) method for saddle point search */
    class Vector3 *_SR1, *_SR2, *_SR3;
    class Vector3 dscom12;      /* center of mass difference for constrained atoms between _SR1 and _SR1 */
    class Vector3 _COM;         /* center of mass */
    class Vector3 **_Rc, **_Rctmp, **_Fc; /* atom positions and forces along a chain of states */
    double *_nebinterp;
    int nebspec[NEBSPECSIZE];
    double energythreshold[ENERGYTHRESHOLDSIZE];
    int readrchainspec;
    double annealspec[ANNEALSPECSIZE];
  
    /* Surface Tension Calculation (Seunghwa Ryu), only works for cubic cell */
    int    SURFTEN;
    int    surftensionspec[SURFTENSIONSPECSIZE];
    int    Pndir, Ptdir1, Ptdir2; /* direction normal to slot, and the two other */
    int    Nnorm; /* number of slots */
    double slotvolume; /* volume of slot, determined by Nnorm */
    double slotheight; /* height of slot, determined by above */
    double slot_halfheight_s; /* half height of slot, in scaled coordinate */
    double *_Pn; /* P normal arrays */
    double *_Pt; /* P transverse arrays */
    int    *_Nslot;
    double *_normpoints; /* coordinates for slot center */
    int    _ENABLE_separation_pot;
    int    ST_step,ST_orien,STsepa_flag;
    double ST_K,ST_Kmax,ST_Y0,ST_Y0max,ST_LAMBDA,ST_LMAX;
    double ST_dUdLAM_L,ST_dUdLAM_POT,ST_dUdLAM_TOT;

    /* Forward Flux Sampling parameters (FFS) (Seunghwa Ryu) */ 
    double *_QLM, *_QL, *_QLweighted;         /* Local order parameter for crystaline structure */
    int    *_Lam_array, *_Lam_check; /* lambda arrays for FFS */
    int    YES_FFS;
    int    L_for_QLM;           /* angular momentum number for QLM */
    int    l_for_qlm;           /* angular momentum number for qlm */
    int    N_solid_P;		/* total number of solid-like particle */
    int    N_skin_P;	        /* total number of skin particle */
    int    N_lgst_cluster;      /* number of particles in largest cluster */
    double Rc_for_QLM;          /* cut-off radius for QLM */
    double N_lgst_inDOUBLE;     /* temporary variable for non-integer RC */
    double QLM_cutoff;		/* order param cutoff for solid-like particle*/
    int    *QLM_mark;		/* DFS search algorithm mark 1 if searched, 0 if not */
    int    *QLM_solid;		/* equals 1 if solid, 0 if not */
    int    N_cluster_temp;	/* temporary memory for cluster size counting */
    Vector3 Principal_Inertia;	/* Inertia for lgst_cluster */
    double N_lgst_index;	/* Special TOPOL indexing for N_lgst_cluster */
    double N_lgst_skin;		/* total number of skin atom for lgst_cluster*/
    double N_lgst_skin_index;   /* Special TOPOL indexing for lgst skin */

#ifdef _GSL
    double *_ql,*_wl,*_ql_,*_wl_,*_Ql;
    double *_qlmr,*_qlmi,*_qlm_r,*_qlm_i,*_qlm;
    char   qlm_type[100];
    int    qlm_id;
#endif

    int    wSKIN;               /* 1 for skin inclusion for cluster, 0 for not */  
    int    DLxyz;		/* x=0 y=1 z=1 when using shear stress  */
    int    DLNdiv;		/* DLNdiv number of division when using shear */

    int    saveFFScn;		/* configuration save if saveFFScn=1 */
    int    FFSoption;		/* 0 for q/T, 1 for P(i+1/i) */
    int    FFSpruning;		/* 0 for time-limited-run, 1 for pruning */
    double    FFScn_weight;     /* default = 1, pruning factor */
    double    FFS_Pp;           /* Pruning Probability */
    int    FFS0_check;   	/* check key for FFS0 run */
    int    lambda_A;		/* FFS i th order parameter */
    int    lambda_B;          /* FFS (i+1)th order parameter */
    int    FFScurstep;         /* FFS (i-1)th order parameter */
    int    FFSautoend;	       /* autoend mark */

    int    FFShist;
    int    *FFSfullarray;
    int    *FFSfullhistogram;
    int    Delta_LAM;
    int    MIN_LAM;
    int    MAX_LAM;
    int    LAM_MAX_INDEX;

    int    FFSbackward;
    int    B_lambda_cut;
    int    B_lambda_B;
    int    B_lambda_0;

    int    FFScommitor;

    /* Umbrella Sampling parameters */
    int YES_UMB; /* set to be 0 as a default, 1 if UMB is running */
    double react_coord, react_coord_old;
    int react_coord_type;
    int react_coord_ngrid; /* number of grid points from state A to state B */

    /* To do: move many of the following outside md.cpp, for not being general enough */
    /*int YES_DK; int DKgroup;*/
    int UMB_order_param; /* 0-9: crystalline order, 10: slip, 20: kink pair */
    int UMB_order_group; /* assign group number for local area constriction */
    double KinkDmax; 
    int printUMBorder;
    Vector3 UMB_nvec;
    Vector3 UMB_R0;
    Vector3 UMB_slipvec;
    int UMB_noslipdirection; /* only consider the size of slip, no directino */
    double UMB_thickness;
    int HETERO_DN;
    double UMB_HMIN, UMB_HMAX;
    int UMB_NNBR;

    int **UMB_nindex_ref;
    int *UMB_nn_ref;
    char *UMB_nindex_ref_mem;

    int MCequilstep;
    double UMB_K_in_Kelvin; /* k of the bias function 0.5k(n-n_c)^2 */
    int UMB_equilstep; /* equilibration step before UMB */
    int UMB_curstep; /* curstep check for UMB */
    double UMB_tryrate;
    int UMB_continue;
    int MC_UMB_cal_period;
    int MC_UMB_log_period;
    int MC_UMB_accept_tot;
    double MC_UMB_accept_ratio;
    int MC_UMB_num_of_trials;
    int n_center; /* n_c of bias function */
    int delta_n; /* the width of the sampling */
    int *Narray; /* the array of number including all n in the window */
    int *Occurence; /* row occurence data */
    double *Probability; /* occurence times the weighting factor */
    int acc_UMB; /* 1 if UMB is accepted, 0 if UMB is rejected */ 

    /* Forward flux sampling parameters */
    int MC_RC_cal_period;
 
    int Kinetic; /* set to be 0 as a defualt, 1 if kinetic is running */
    double Kinetic_Time;
    int Kinetic_Swip;
    double *time_array;
    int *n_array;
    int kinetic0;
    int kineticNmax;
    int N_lgst_temp;
    int kinetic_flag;
    int KN_Center;
    Vector3 *_SR4, *_R4, *_VSR4, *_VR4;
    Matrix33 _H4;

    int YES_HMC; /* set to be 0 as a default, 1 if UMB is running */
    int Ns_HMC; /* number of steps before HMC acc evaluation */
    double EPOT_temp; /* potential E of the system of previous step */
    double KATOM_temp; /* kinetic E of the system of previous step */
    int acc_HMC; /* 1 if HMC is accepted, 0 if HMC is rejected */ 
    double T_HMC; /* next target temperature */

    int YES_SEP; /* set to be 0 as a default, 1 if species speration on */
    Vector3 _CLUSTER_CM; /* center of mass for largest cluster */
    Vector3 _SRforCM;
    double SEPA_ORDER; /* larger this, more separation */
    double SEPA_TARGET; /* YES_SEP = 2, this target is used */
    double SEPA_RATIO; /* ratio of separation in reaction coordinate */
    double NORM_REACTION; /* normalized reaction coordinate for each section */
    double NC_TIMES_SEPA; /* N_lgst_cluster x SEPA_order */
    int IMP_INDEX; /* species index for impurity particle */
    double IMP_TOPOL_INDEX; /* visuallization topology index for impurity particle */
    double IMP_R_EXP; /* impurity order-parameter r exponent */

    /* Ewald Summation for Coulomb interaction */
    Vector3 *_F_Ewald, *_F_Ewald_Real, *_F_Ewald_Rec;
    double *_EPOT_IND_Ewald, *_EPOT_IND_Ewald_Real, *_EPOT_IND_Ewald_Rec;
    Matrix33 _VIRIAL_Ewald, _VIRIAL_Ewald_Real, _VIRIAL_Ewald_Rec;
    int Ewald_CE_or_PME;    /*0: CE (classical Ewald), 1: PME (particle mesh Ewald) */
    int Ewald_option_Alpha; /*different options to determine Alpha (spread of Gaussian charge) */
    double Ewald_Rcut;      /* real space cut-off distance for Ewald summation */
    double Ewald_precision; /* cut-off error in both real and reciprocal space */
    double Ewald_time_ratio;/* computational time ratio between real and reciprocal space */
    double _EPOT_Ewald, _EPOT_Ewald_Real, _EPOT_Ewald_Rec, _EPOT_Ewald_Self;
    double Ewald_Alpha, Ewald_qSQRsum;
    double Ewald_H_L1, Ewald_H_L2, Ewald_H_L3, Ewald_cell_V; /* simulation box size */
    
    int CE_MKV, CE_nKV;         /* max and actual number of K-vectors in list */
    int CE_Mka, CE_Mkb, CE_Mkc; /* Largest repeat reciprocal distance in each direction */
    double CE_Kc;               /* cut-off radius in reciprocal space */
    
    struct K_rec { double kx, ky, kz, k2; int a,b,c; } *CE_KV; /* Actual K-vectors */
    EComplex *CE_sc_axby, *CE_sckk;  /* structural factors exp(iK(rx+ry)), qj*exp(iK(rx+ry+rz)) */
    EComplex **CE_scx, **CE_scy, **CE_scz; /* sin-cos table for computing structural factor */
    
    clock_t CE_sREAL, CE_sREC, CE_eREAL, CE_eREC; /* time information for CE calculation */
    double  CE_iREAL, CE_iREC, CE_tREAL, CE_tREC;
    
    /* Particle Mesh Ewald Summation (Hark Lee, 2007/10/18) */
    int PME_bsp_n;                       /* order of B-spline interpolation */

    int PME_K1, PME_K2, PME_K3, PME_K3S; /* number of mesh points */
    int PME_K23, PME_K123, PME_K23S, PME_K123S;
    double PME_K1d, PME_K2d, PME_K3d, PME_K123d;
    
    Vector3 *_UR;                        /* Scaled atom positions */

    double *PME_B;                       /* B array */
    double *PME_C;                       /* C array */
    double *PME_BC;                      /* B*C array */
    double *PME_Q;                       /* Q array */
    double *PME_CONV;                    /* convolution of BC and IQ */
    ComplexType *PME_IQ;                 /* inverse Q array */

    ComplexType PME_R_bsp_Bm;            /* entries of B */
    double *PME_bsp_Mk;                  /* to calculate B entries */
    double *PME_bsp_fac;                 /* factorial array */

    int *PME_m1s, *PME_m2s, *PME_m3s;    /* k-space start indices */
    int *PME_m1, *PME_m2, *PME_m3;       /* k-space indices */
    int *PME_m1mod, *PME_m2mod, *PME_m3mod; /* actual indices after PBC */

    double *PME_x1, *PME_x2, *PME_x3;    /* (u-k) after PBC */
    double *PME_MU1, *PME_MU2, *PME_MU3; /* Mn(u) */
    double *PME_d1MU1, *PME_d1MU2, *PME_d1MU3; /* Mn-1(u) */
    double *PME_d2MU1, *PME_d2MU2, *PME_d2MU3; /* Mn-1(u-1) */

    fftw_plan PME_fft_plan1, PME_fft_plan2; /* fftw plans */
    
    clock_t PME_sREC, PME_eREC;          /* time information for PME calculation */
    double  PME_iREC, PME_tREC;
    
    /* Configuration manipulation */
    class Vector3 *storedr;         /* atomic displacements */
    class Matrix33 dH;              /* change of matrix H (cell size) */
    char  crystalstructure[30];     /* crystal structure for makecrystal() */
    double latticeconst[3];         /* lattice constant  for makecrystal() */
    double latticesize[3][4];       /* lattice size for makecrystal() */
    int _TORSIONSIM;
    double torquespec[TORQUESPECSIZE];           /* torsion specification */
    double _TORQUE, *_TORQUE_IND;
    int _BENDSIM;
    double bendspec[BENDSPECSIZE];             /* bending specification */
    double _BENDMOMENT, *_BENDMOMENT_IND;
    Vector3 _P_COM, _F_COM, _L_COM; /* center of mass momenta */
    double _PTHETA_COM; 

    /* File input and output */
    double input[MAXINPUTLEN];              /* general input variable in script file */
    double output_dat[MAXOUTPUTDATLEN];                /* general output data array */
    char   output_str[MAXOUTPUTSTRLEN];               /* general output string */
    char   output_fmt[MAXOUTPUTSTRLEN];               /* output format */
    class CNFile initcn,intercn,FFScn,finalcn,continuecn;/* configuration file objects */
    class PropFile pf;                      /* property file objects */
    char incnfile[MAXFILENAMELEN], finalcnfile[MAXFILENAMELEN];   /* configuration file name */
    char intercnfile[MAXFILENAMELEN];                  /* configuration file name */
    char continuecnfile[MAXFILENAMELEN];
    char FFScnfile[MAXFILENAMELEN];
    char intercfgfile[MAXFILENAMELEN];                 /* atomeye cfg file name */
    char outpropfile[MAXFILENAMELEN];                  /* property file name */
    char myname[MAXFILENAMELEN];                       /* name of the simulation */
    char potfile[MAXFILENAMELEN];                      /* name of potential file */    
    char command[MAXCMDLEN]; int ncom;      /* shell command */
    int  savecn, savecnfreq;                /* frequency of saving cn files */
    int  savecontinuecnfreq;	/* frequency of saving continuecn files */
    int  savecfg, savecfgfreq, savecfgcount;/* frequency of saving cfg files */
    int  saveprop, savepropfreq;            /* frequency of saving prop files */
    int  printfreq;                         /* frequency of printing on scren */
    int  filecounter;                       /* number of current file */
    int  FFSfilecounter;                    /* number of current FFS file */
    int allocmultiple;                      /* allocate more memory than current number of atoms */
    int writevelocity;                      /* write velocity into cn file */
    int writeall;                           /* write "sx,sy,sz,vx,vy,vz,pot,fixed" into cn file */
    int zipfiles;                           /* zip output files */
    int fixedatomenergypartition;           /* 1: disregard energy contribution from fixed atoms */    

    /* Interface to Fortran code (Relax) */
    char fortranpath[MAXFILENAMELEN], fortranexe[MAXFILENAMELEN];

    /* Interface to AtomEye and other viewers */
    char atomeyepath[MAXFILENAMELEN], atomeyeexe[MAXFILENAMELEN];
    int atomeyerepeat[4];

    /* Interface to LAMMPS */
    class LAMMPSFile initlammps,finallammps;
    
    /* Visualization */
    YWindow *win;                   /* Display window object */

    int win_width, win_height;      /* size of the window */
    double rotateangles[5];         /* rotation angles and scaling parameter */
    int plotfreq;                   /* frequency of re-plotting */
    double atomradius[MAXSPECIES];  /* size of atoms */
    double bondradius, bondlength;  /* thickness and length of bonds */

    char atomcolor[MAXSPECIES][30]; /* color of atom species */
    char bondcolor[30];             /* color of bonds */    
    char fixatomcolor[30];          /* color of fixed atoms, fixed[i]=1 */
    char highlightcolor[30];        /* color of highlighted atoms */
    char backgroundcolor[30];       /* background color of window */
    char colornames[MAXCOLORS][30]; /* names of allocated colors */
    unsigned colors[MAXCOLORS+15];  /* value of allocated colors */

    int plot_highlight_atoms[10000];/* indices of atoms to be highlighed */    
    double plot_limits[10];         /* set x,y,z, limits of plotting regime */
    int    plot_atom_info;          /* print information when atom is clicked */
    int    plot_map_pbc;            /* map all atoms in primary cell when plotting */
    double plot_color_windows[100]; /* only plot atoms whose properties fall into specified windows */
    double plot_color_bar[10];      /* how to color atoms according to their properties */
    int    plot_color_axis;          /* which atomic property specifies color */

    int autowritegiffreq;           /* frequency of outputing .GIF graphics file */
    double *color_ind;              /* property linked to atom color */
    int NCS;                        /* number of neighbors for central symmetry parameter comp. */
    int coloratoms;                 /* determine how to color atoms 
				   0 = color by species, 1 = color by groups */

    /* Constructor (set initial values of member variables) */
    MDFrame(): _SAVEMEMORY(0),_NP(0),_NP0(0),_NIMAGES(0),_NPfree(0),_NPfixed(0),

               _R(0),_SR(0),_SRA(0),_VR(0),_VSR(0),_F(0),_R0(0),_Fext(0),_ENABLE_FEXT(0),

               _FLAT_INDENTOR_POS(0), _FLAT_INDENTOR_POS0(0), _FLAT_INDENTOR_V(0),
               _FLAT_INDENTOR_K(0),   _FLAT_INDENTOR_F(0),    _FLAT_INDENTOR_DIR(0),
               _ENABLE_FLAT_INDENTOR(0),

               fixed(0),species(0),group(0),nspecies(1),
               _EPOT(0),_EPOT_IND(0),_EPOT_RMV(0),_EPOT_BOX(0),_EPOT0(0),_EPOT0_ext(0),_TOPOL(0),

               /* Molecular Dynamics parameters */
               _ESTRAIN(0),_KATOM(0),_KBOX(0),_WALLMASS(1),
               _TIMESTEP(1e-3),_TDES(0),_T(0),_DOUBLET(0),
               _ATOMTCPL(100),_BOXTCPL(100),_BOXDAMP(0),
               zeta(0),zetav(0),zetaa(0),vt2(0.25),_HELM(0),_HELMP(0),
               //usenosehoover(0),usescalevelocity(0),
               applyshear(0),
               forcemul(1), runavgposition(0),
               _ENABLE_DPD(0), _DPD_fc(0), _DPD_RATIO(1), _DPD_RCUT(100),
               _ENABLE_ANDERSON(0), _ANDERSON_RATIO(0),

               /* Neighbor list */
               _RLIST(0),_SKIN(0),_NIC(100),_current_NIC(0),_NNM(60),_current_NNM(0),
#ifndef _USEBOOST
               nn(0), nindex(0), nindex_mem(0), celllist(0), cell_mem(0),
#endif
               link_head(0),link_list(0), firsttime(true),
               rc_plot(0), _NNM_plot(4), nn_plot(0), nindex_plot(0), nindex_plot_mem(0),
               
               /* Numerical integrator */
               s2(0),s3(0),s4(0),s5(0),s2real(0), /* Gear 6th predictor-corrector */
               zeta2(0),zeta3(0),zeta4(0),zeta5(0),
               npold(0),

               /* Nose Hoover chain */
               NHChainLen(2), 

               /* Internal and applied stresses */
               _TOTPRESSURE(0),_EXTSTRESSMUL(1),_EXTPRESSADD(0),_VACUUMRATIO(0),_Dexx(0),

               /* Conjugate gradient relaxation */
               conj_itmax(1000),conj_fevalmax(10000),conj_fixbox(1),
               conj_fixshape(0),conj_fixatoms(0),conj_step(0),
               conj_dfpred(1e-3),conj_f(0),conj_etol(1e-8),conj_ftol(1e-8),conj_g2res(0),conj_space(0),
               _HPRECOND(1),_XPRECOND(1),_YPRECOND(1),_ZPRECOND(1),

               /* Conjugate gradient relxation with constrains */
               constrainS(0),constrain_dist(0),
               
               /* constrained MD */
               _constrainedMD(0),
               
               /* Molecular Dynamics simulation steps */
               totalsteps(100), equilsteps(0), curstep(0), step0(0), continue_curstep(0),
               implementation_type(0), algorithm_id(0),

               /* Monte Carlo simulation */
               MC_atom(0),MC_accept(0),MC_accept_tot(0), 
               MC_dV_freq(1), MC_dN_freq(1), MC_accept_ratio(0), MC_dVmax(0.1), MC_dfrac(0.1),
               MC_P_ext(0), MC_mu_ext(0), MC_Lambda3(0), MC_fracatom(1),
               _EPOT_IND1(0), _EPOT_INDUMB(0), 
               _SR_MC_SAV(0), _SR_MC_SAVUMB(0),
               MC_maxd(0),
               

               /* Free energy calculations */
               _LAMBDA(0),
               _LAMBDA0(0),_LAMBDA1(1),dEdlambda(0),dlambdadt(0),_ECOH(0),_WTOT(0),_WAVG(0),
               _ENABLE_SWITCH(0), _SWITCHFREQ(1),
               _RANDSEED(12345),_CHAINLENGTH(0),_SWITCHFUNC(0),_REFPOTENTIAL(0),
               heslen(0), hesind(0), hes(0),

               /* Switching potentials */
               _ECSPRING(0),
               _I12_RC(0), _I12_SIGMA(0), _I12_EPSILON(0),
               _GAUSS_RC(0), _GAUSS_SIGMA(0), _GAUSS_EPSILON(0),
               
               /* Nudged-Elastic-Band (NEB) method */
               _SR1(0),_SR2(0),_SR3(0),_Rc(0),_Rctmp(0),_Fc(0),_nebinterp(0), //_Fnebinterp(0),
               readrchainspec(-1),
 
               /* Surface tension calculation */
	       SURFTEN(0), Pndir(0), Ptdir1(1), Ptdir2(2), Nnorm(1), 
               slotvolume(1), slotheight(1),
               slot_halfheight_s(1), _Pn(0), _Pt(0), _Nslot(0), _normpoints(0),
               _ENABLE_separation_pot(0), ST_step(0),ST_orien(0),
               STsepa_flag(0), ST_K(0), ST_Kmax(0),
               ST_Y0(0), ST_Y0max(0), ST_LAMBDA(0), ST_LMAX(0),
               ST_dUdLAM_L(0), ST_dUdLAM_POT(0), ST_dUdLAM_TOT(0),


               /* Forward Flux Sampling (FFS) parameters (Seunghwa Ryu) */ 
               _QLM(0),_QL(0), _QLweighted(0), _Lam_array(0),_Lam_check(0), YES_FFS(0), L_for_QLM(0), l_for_qlm(0),
               N_solid_P(0),N_skin_P(0),N_lgst_cluster(0),Rc_for_QLM(0),
               N_lgst_inDOUBLE(0),
               QLM_cutoff(0),QLM_mark(0), QLM_solid(0), N_cluster_temp(0), 
               N_lgst_index(0), N_lgst_skin(0), N_lgst_skin_index(0), 
#ifdef _GSL
    _ql(0),_wl(0),_ql_(0),_wl_(0),_Ql(0),
    _qlmr(0),_qlmi(0),_qlm_r(0),_qlm_i(0),_qlm(0),qlm_id(1),
#endif
	       wSKIN(0), DLxyz(0), DLNdiv(0),
               saveFFScn(0), FFSoption(0), FFSpruning(0), FFScn_weight(1), 
               FFS_Pp(0), FFS0_check(1), lambda_A(0), lambda_B(0), 
               FFScurstep(0), FFSautoend(0), FFShist(0), FFSfullarray(0), 
               FFSfullhistogram(0), Delta_LAM(0), MIN_LAM(0), MAX_LAM(0), LAM_MAX_INDEX(0),
               FFSbackward(0), B_lambda_cut(0), B_lambda_B(0), B_lambda_0(0),
               FFScommitor(0),
   
               /* Umbrella Sampling parameters */
               YES_UMB(0), /* YES_DK(0), DKgroup(0), */
               react_coord(0), react_coord_old(0), react_coord_type(0), react_coord_ngrid(0),
               UMB_order_param(0), UMB_order_group(0), 
               KinkDmax(1), printUMBorder(0), UMB_noslipdirection(0),
               UMB_thickness(0), HETERO_DN(0),
               UMB_HMIN(0), UMB_HMAX(0), UMB_NNBR(0), 
               UMB_nindex_ref(0), UMB_nn_ref(0), UMB_nindex_ref_mem(0), 
               MCequilstep(0), UMB_K_in_Kelvin(0), UMB_equilstep(0), UMB_curstep(0), 
	       UMB_tryrate(1), UMB_continue(0),
               MC_UMB_cal_period(1000), MC_UMB_log_period(1000),
               MC_UMB_accept_tot(0), MC_UMB_accept_ratio(0),
               MC_UMB_num_of_trials(0), n_center(0),
               delta_n(0), Narray(0), Occurence(0), Probability(0), acc_UMB(0),
               MC_RC_cal_period(1000),
               Kinetic(0), Kinetic_Time(10),
               Kinetic_Swip(0), time_array(0), n_array(0), kinetic0(0), kineticNmax(0), N_lgst_temp(0), 
	       kinetic_flag(0), KN_Center(0), _SR4(0), _R4(0), _VSR4(0), _VR4(0), YES_HMC(0), Ns_HMC(10), 
	       EPOT_temp(0), KATOM_temp(0), acc_HMC(0), T_HMC(0),
               YES_SEP(0),SEPA_ORDER(0),SEPA_TARGET(0), SEPA_RATIO(0), NORM_REACTION(0),NC_TIMES_SEPA(0), 
               IMP_INDEX(0), IMP_TOPOL_INDEX(0), IMP_R_EXP(2.5),
               
               /* Ewald Summation for Coulomb interaction */
               _F_Ewald(0), _F_Ewald_Real(0), _F_Ewald_Rec(0),
               _EPOT_IND_Ewald(0), _EPOT_IND_Ewald_Real(0), _EPOT_IND_Ewald_Rec(0),
               Ewald_CE_or_PME(0),
               Ewald_option_Alpha(0),
               Ewald_Rcut(9),
               Ewald_precision(4.2),
               Ewald_time_ratio(3.6),
               _EPOT_Ewald(0), _EPOT_Ewald_Real(0),
               _EPOT_Ewald_Rec(0), _EPOT_Ewald_Self(0),               
               Ewald_Alpha(0), Ewald_qSQRsum(0),
               Ewald_H_L1(0), Ewald_H_L2(0), Ewald_H_L3(0), Ewald_cell_V(0), 
               
               CE_MKV(0), CE_nKV(0),
               CE_Mka(0), CE_Mkb(0), CE_Mkc(0),
               CE_Kc(0),
               CE_KV(0),
               CE_sc_axby(0), CE_sckk(0),
               CE_scx(0), CE_scy(0), CE_scz(0),
               CE_sREAL(0), CE_sREC(0), CE_eREAL(0), CE_eREC(0),
               CE_iREAL(0), CE_iREC(0), CE_tREAL(0), CE_tREC(0),    
               
               /* Particle Mesh Ewald Summation (Hark Lee, 2007/10/18) */ 
               PME_bsp_n(8),
               PME_K1(64), PME_K2(64), PME_K3(64), PME_K3S(0),
               PME_K23(0), PME_K123(0), PME_K23S(0), PME_K123S(0),
               PME_K1d(0), PME_K2d(0), PME_K3d(0), PME_K123d(0),
               _UR(0),

               PME_B(0), PME_C(0), PME_BC(0), PME_Q(0), PME_CONV(0), PME_IQ(0),
               PME_bsp_Mk(0), PME_bsp_fac(0),
               PME_m1s(0), PME_m2s(0), PME_m3s(0),
               PME_m1(0), PME_m2(0), PME_m3(0),
               PME_m1mod(0), PME_m2mod(0), PME_m3mod(0),

               PME_x1(0), PME_x2(0), PME_x3(0),
               PME_MU1(0), PME_MU2(0), PME_MU3(0),
               PME_d1MU1(0), PME_d1MU2(0), PME_d1MU3(0),
               PME_d2MU1(0), PME_d2MU2(0), PME_d2MU3(0),

               PME_sREC(0), PME_eREC(0), PME_iREC(0), PME_tREC(0),
               
               /* Configuration manipulation */
               storedr(0),

               /* Torsion simulation */
               _TORSIONSIM(0),_TORQUE(0),

               /* Bending simulation */
               _BENDSIM(0), _BENDMOMENT(0),
        
               /* Center of mass momenta */
               _PTHETA_COM(0),               

               /* File input and output */
               initcn(AUXFile::BLOCK),intercn(AUXFile::SERIES),
               FFScn(AUXFile::SERIES),
               finalcn(AUXFile::BLOCK),
	       continuecn(AUXFile::BLOCK),
               pf(AUXFile::ENTRIES),
               ncom(0),savecn(0),savecnfreq(100),
               savecontinuecnfreq(0), savecfg(0),
               savecfgfreq(100),savecfgcount(1),
               saveprop(0),savepropfreq(100),
               printfreq(100),filecounter(1),FFSfilecounter(1),
               allocmultiple(1),writevelocity(0),writeall(0),
               zipfiles(0),fixedatomenergypartition(0),

               /* Interface to LAMMPS */
               initlammps(AUXFile::BLOCK),finallammps(AUXFile::BLOCK),

               /* Visualization */
               win(0),win_width(350),win_height(350),
               plotfreq(1),bondradius(0.1),bondlength(0),plot_atom_info(1),plot_map_pbc(0),
               plot_color_axis(0),autowritegiffreq(0),color_ind(0),NCS(12),coloratoms(0)
               
    {
        /* Input/output control */
        input[0]=0;
        output_str[0]=0;
        output_fmt[0]=0;
        _ATOMMASS[0]=1;
        
        /* Configuration manipulation */
        latticeconst[0]=latticeconst[1]=latticeconst[2]=1.0;
        sprintf(crystalstructure,"simple-cubic");

        /* Workspace for Monte Carlo simulation */
        MC_switchoffatoms[0]=0;

        /* Plot settings */
        plot_highlight_atoms[0]=0;
        plot_limits[0]=0;
        plot_color_windows[0]=0;
        plot_color_bar[0]=0;

        
    };

    virtual ~MDFrame() { delete(win); }
            
    /* Parser */    
    void initvars();
    void clear_input();       /* set all entries in input array to zero */
    
    /* Coordinate transformation */
    void SHtoR();
    void RHtoS();
    void VSHtoVR();
    void VRHtoVS();
    void RtoR0();
    void R0toR();
    void clearR0(); /* clear R0 so that neighbor list will be reconstructed */
    bool Bond(int I, int J) const;
    
    virtual void Alloc();

    /* Potentials */    
    virtual void call_potential(); /* by default calls potential but in mdparallel calls potential_parallel */
    void call_ANDERSON();
    virtual void potential();
    virtual void potential_energyonly();
    virtual double potential_energyonly(int iatom);
    void ECpotential();            /* Harmonic potential of Einstein crystal */
    void ECpotential_energyonly();
    double ECpotential_energyonly(int iatom);
    void HApotential();
    void HApotential_energyonly(); /* Harmonic potential based on Hessian matrix */
    void I12potential_energyonly(); /* Inverse r-12 potential as liquid reference */
    void I12potential();
    void GAUSSpotential();          /* Gaussian potential as liquid reference */
    void VASPpotential();           /* call VASP to compute energy and force */
    void SWITCHpotential_energyonly(double lambda); /* for adiabatic switching */
    void SWITCHpotential(double lambda); 
    double SWITCH_Find_Lambda();
    void eval_insertparticle();      /* compute potential energy of test particle */
    void eval_removeparticle();      /* compute potential energy of test particle */

    /* Neighbor list */
    virtual void NbrList_reconstruct(int iatom=-1);  
    void NbrList_reconstruct_use_cell_list(int iatom=-1); /* use cell list with Verlist */
    void NbrList_reconstruct_use_link_list(int iatom=-1); /* use link list with Verlist */
    
    void NbrList_print();        /* for debug use */
    void NbrList_refresh();
    bool NbrList_needrefresh();
    void NbrList_init(int, int);
    void NbrList_free();

    void NbrList_initlinklist(int);
    void NbrList_freelinklist();
    
    void NbrList_initcell(int,int,int,int);
    void NbrList_freecell();

    void NbrListPlot_reconstruct();
    void NbrListPlot_init(int,int);
#define refreshneighborlist NbrList_refresh
    
    /* Numerical integrator */
    void NVE_Gear6();
    void NVT_Gear6();
    void NVTC_Gear6();
    void NPH_Gear6();
    void NPT_Gear6();
    void NPTC_Gear6();
    void Gear6_init();
    void VVerlet_Get_s2();
    void VVerlet_Get_h2();
    void VVerlet_Get_Tinst();
    void VVerlet_NHCINT();

    void NVE_VVerlet(int);
    void NVT_VVerlet_Implicit(int);
    void NVT_VVerlet_Explicit_1(int);
    void NVT_VVerlet_Explicit_2(int);
    void NVTC_VVerlet_Explicit(int);
    void NPH_VVerlet(int);
    void NPH_VVerlet_Explicit_1(int);
    /* void NPH_VVerlet_Explicit_2(int); */ /* not working yet */
    void NPT_VVerlet_Implicit(int);
    void NPT_VVerlet_Explicit_1(int);
    void NPT_VVerlet_Explicit_2(int);

    /* Conjugate gradient relaxation */
    int calconstrainS_inst(); /* calculate constrainS_inst */
    int applyconstraint();
    int relax();           /* CG relax */
    int constrainedrelax();/* CG relax with constraints */
    void relax_by_group(double *);
    int rigid_relaxation();
    static void potential_wrapper(int *n,double x[],double *f,double df[]);
    static void potential_wrapper_fixbox(int *n,double x[],double *f,double df[]);
    static void potential_wrapper_c(int n,double x[],double *f,double df[]);
    static void potential_wrapper_rigid_translate(int n,double x[],double *f,double df[]);
    static void cstr_potential_wrapper(int *n,double x[],double *f,double df[]);
    static void cstr_potential_wrapper_c(int n,double x[],double *f,double df[]);
    void strelax();         /* steepest descent relaxation */
    
    /* Molecular Dynamics simulation */
    virtual void integrator();       /* MD integrator */

    virtual void run();              /* run the simulation */
    virtual void calcprop();         /* calculate physical properties */
    virtual void calcoutput();       /* calculate output string */
    void printResult();
    void calDVIRIALDexx();           /* calculate strain derivative of stress */
    void calcentralsymmetry();       /* calculate central symmetry deviation */
    void caldisregistry();           /* for slip (dislocation) visualization */
    void calHessian();               /* numerical calculation of Hessian matrix */
    void readHessian();              /* read Hessian matrix from file */
    void calModeHessian();           /* calculate Hessian matrix for a set of deformation modes */
    void calphonondisp();            /* calculate phonon dispersion relation */
    void calmisfit();                /* calculate misfit potential */
    void calmisfit2();               /* calculate misfit potential */
    void calmisfit2_rigidrlx();      /* calculate misfit potential. Relax normal to the surface */
    void testpotential();            /* test consistency between force and energy */
    void findcore();                 /* find center of mass of defect core atoms */
    void initvelocity();
    void randomposition();
    //void velscaler(); /* disable velocity scaler */
    void perturbevelocity(); 	     /* gaussian perturbation */
    void MCperturbevelocity();       /* shooting for MC TPS for NVE, NVT */
    void zero_angmom();
    void zero_totmom();
    void multiplyvelocity();
    void step();
    void eval();
    void multieval();                /* for testing/timing purposes */

    virtual double potential_energyonly_before(int iatom);
    virtual double potential_energyonly_after(int iatom);
    int MCstep();
    int MCstep_moveatom();
    int MCstep_dV();

#ifdef _MC_CHANGE_NUM_ATOMS
    int MCstep_dN();
    int MCstep_fracatom();
    int MCstep_moveatom_fracatom();
    int MCstep_dV_fracatom();
    int MCstep_dN_fracatom();    
#endif

    int MCSWITCHstep(double lambda);
    int runMDSWITCH();
    
    /* Nudged-Elastic-Band relaxation */
    void AllocChain();
    void allocChainVars(int n, double **pEc, double **pFm, double **pFt, double **pFs, double **pdR, double **dLc, double **pTangMag2, Vector3 ***pdTang, void **pc);
    void calChainForce(int n, int relax_surround, double *Ec);
    void calChainTang(int n, double *Ec, double *dR, double *TangMag2, Vector3 **dTang);
    void orthoChainForce(int moveleftend, int moverightend, int yesclimbimage, int EmaxDomain, int n, double *Fm, double *Ft, double *TangMag2, Vector3 **dTang);
    void reparamChain(int moveleftend, int moverightend, int yesclimbimage, int EmaxDomain, int n, double *plavg, double *plavg0, double *dR);
    void caldscom12();
    void initRchain();
    int readRchain();
    int readRchain_parallel_toCN(int j);
    int writeRchain();
    int writeRchain_parallel_fromCN(int j);
    void stringrelax();
#ifdef _STRINGMETHOD
    void runstringmethod();
#endif
    void nebrelax();   
    void annealpath();       /* simulated annealing of paths */
    void cutpath();
    int statedistance();     /* compute distance between two states _SR1 and _SR */
    int interpCN(double s); /* 2007/3/6 Wei Cai */
    int interpCN_Quadratic(double);
    void copyRchaintoCN(int j);
    void copyCNtoRchain(int j);
    void moveRchain(int, int);
 
    //void Fold_into_Unitcell(); /* the same as maptoprimarycell() */
    void STseparation();

#ifdef _GSL
    void calqlwl();
    void allocqlm();
#endif
    void cal_react_coord();
    void cal_react_coord_from_MEP();
    void caldislocationorder();
    void allocQLM();
    void allocKinetic();
    void printhist();
    void readhist();
    void kineticprint();
    void initKinetic();
    int Kineticprocess();
    void setconfig4();
    void copyS4toS();
    void copyStoS4();

    
    /* Configuration manipulation */
    void makecrystal();         /* create perfect crystal structure */
    void makecut();             /* create cut surfaces in lattice */
    void makedipole();          /* create dislocation dipole structure */
    void makedislocation();     /* create an arbitrarily oriented dislocation line */
    void makedislpolygon();     /* create a polygonal dislocation loop */
    void makedislellipse();     /* create an elliptical dislocation loop */
    void makecylinder();        /* cut a cylinder out of the crystal */
    void setfixbufferatoms();   /* fix atoms around a cylinder */
    void makedislcylinder();    /* create a single dislocation in a cylinder */
    void commit_storedr();      /* apply the displacement field in storedr to R */
    void commit_removeatoms();  /* remove atoms marked with fixed = -1 */
    
    /*void makedisloop();*/         /* create dislocation loop */
    void calloopdisp            /* compute the displacement field */
         (Vector3 *,double,double,Vector3 *,Vector3 *,Vector3 *);
    void makegrainboundary();   /* create grain boundary */
    void makewave();            /* initiate wave form */
    void cutbonds();            /* remove bonds across specified cut plane */
    void cutbonds_by_ellipse(); /* remove bonds across specified cut plane */
    void maketorquehandle();    /* prepare for torque application */
    void removetorquehandle();  /* restore usual PBC */
    void makebendhandle();      /* prepare for bend application */
    void copytorqueatoms();     /* copy end atoms for torsion simulation */
    void copybendatoms();       /* copy end atoms for bending simulation */
    void writeimagefile(const char *fname);/* write image information into file */
    void readimagefile(const char *fname); /* read image information from file */
    
    void scaleH();              /* hydrostatic scale of simulation box */
    void setH();
    void saveH();
    void restoreH();            /* restore H from H0 */
    void scaleVel();            /* scale velocity */
    
    void shiftbox();            /* shift PBC box */
    void redefinepbc();
    void applystrain();         /* apply homogeneous strain */
    void extendbox();           /* make multiple copies of the box*/
    void switchindex();         /* switch two box directions, e.g. 1(x), 2(y) */ 
    void maptoprimarycell();    /* map all atoms to primary cell */

    void cutslice();            /* cut a portion from box, inverse of extendbox */
    void splicecn();            /* put two configurations together */
    void cutpastecn();          /* cut and paste a region of the config */

    void setconfig1();          /* copy current _SR to _SR1 */
    void setconfig2();          /* copy current _SR to _SR2 */
    void setconfig3();          /* copy current _SR to _SR3 */
    void SR1toSR();
    void switchconfig();        /* R - R1 + R2 */
    void copytoconfig1(int);
    void copytoconfig2(int);
    void replacefreeatom();     /* replace free (not fixed) atoms from another config */
    void relabelatom();
    
    void moveatom();
    void movegroup();
    void rotategroup();
    void setgroupcomvel();      /* set the center of mass velocity of specified groups */
    void addtorque();
    void addbend();
    void printatoms_in_sphere();/* print the indices of all atoms in a sphere */
    void pbcshiftatom();        /* shift atoms in PBC box */
    
    void fixatoms_by_ID();      /* fix a set of atoms whose ID are specified */
    void fixatoms_by_position();/* fix all atoms whose position falls within a specified regime */
    void fixatoms_by_group(double *);
    void RtoR0_by_group(double *);
    void R0toR_by_group(double *);
    void SR1toSR_by_group(double *);
    void SR2toSR_by_group(double *);
    void fixatoms_by_species(double *);
    void fixatoms_by_pos_topol();/* fix atoms whose position and TOPOL value are in a specified regime */
    void fixallatoms();
    void freeallatoms();
    void reversefixedatoms();
    void constrain_fixedatoms();
    void fix_constrainedatoms();
    void fix_imageatoms();

    void setfixedatomsspecies();
    void setfixedatomsgroup();
    void reversespecies();      /* species: 0->1 , 1->0 */
    void movefixedatoms();      /* move atoms that are marked as fixed */
    void removefixedatoms();    /* remove atoms that are marked as fixed */
    void markremovefixedatoms();/* mark fixed atoms for removal (fix=-1) */
    void makeellipsoid();       /* mark or remove all atoms within an ellipsoid */    
    void removeellipsoid();     /* remove all atoms within an ellipsoid */    
    void removerectbox();       /* remove all atoms within a parallelpepid */
    void removeoverlapatoms();  /* remove atoms that are too close to each other*/
    void findcenterofmass(class Vector3 *);
    void translate_com();       /* apply rigid-body translation so that COM matches that of config1 */
    void rotate_com();          /* apply rigid-body rotation so that COM matches that of config1 */

    void clearFext();           /* reset Fext array */
    void addFext_to_group();    /* add Fext to atoms belonging to specific group */

#ifdef _GSL
    void compute_XRD_intensity(); /* compute spherically averaged XRD intensity I(2theta) */
#endif

    /* File input and output */
    virtual void saveintercn(int);   /* save intermediate cn files */
    virtual void saveintercfg(int);   /* save intermediate cfg files */
    int readcn();
    int readcontinuecn();
    int readXYZ();
    int readLAMMPS();
    int convertXDATCAR();
    int setfilecounter();
    int setFFSfilecounter();
    int openintercnfile();
    int openFFScnfile(); 
    int openpropfile();
    int closepropfile();
        
    int writefinalcnfile(int zip=1,bool bg=true);
    void savecontinuecn(int step);
    int writecontinuecnfile(int zip=1,bool bg=true);
    int writeavgcnfile(int zip=1,bool bg=true);
    int writeLAMMPS();
    void writePINYMD(const char *fname);
    int writeintercnfile(int zip=1,bool bg=false);
    int writeFFScnfile(int zip=1,bool bg=false);

    /* Interface with Fortran code (Relax) */
    void writefortraninifile(const char *fname);

    /* Interface with AtomEye and other viewers */
    virtual void writeatomeyecfgfile(const char *fname);
    void convertCNtoCFG();
    void atomeye();
    void writepovray(const char *fname);
    void writeMDCASKXYZ(const char *fname);
    void writeRASMOLXYZ(const char *fname);
    void readRASMOLXYZ(const char *fname);
    void writeatomtv(const char *fname);
   
#ifdef _USEHDF5
    /* Interface with HDF5 (HDFView) */
    int writecn_to_HDF5();
#endif
     
    /* Print out energies of current configuration */
    void writeENERGY(const char *fname);
    void writeFORCE(const char *fname);
    void writePOSITION(const char *fname);
    void GnuPlotHistogram();

    /* Visualization */
    virtual void winplot();        
    virtual void winplot(int);
    virtual void openwin();
    virtual void plot();
    int     openwindow(int w,int h,const char *n);
    void    closewindow();
    void    wintogglepause();
    void    alloccolors();
    void    alloccolorsX();
    void    rotate();
    void    saverot();
   
    /* set variables from arg list */
    void set_dirname(const char *);
    void set_randseed(const char *); 

    //bondcolor
    inline std::string get_bondcolor() { std::string str(bondcolor); return str; }
    inline void set_bondcolor(std::string s){ strcpy(bondcolor, s.c_str()); };

    //backgroundcolor
    inline std::string get_backgroundcolor() { std::string str(backgroundcolor); return str; }
    inline void set_backgroundcolor(std::string s){ strcpy(backgroundcolor, s.c_str()); };

    //fixatomcolor
    inline std::string get_fixatomcolor() {std::string str(fixatomcolor); return str; }
    inline void set_fixatomcolor(std::string s){ strcpy(fixatomcolor, s.c_str()); };

    //hightlightcolor
    inline std::string get_highlightcolor() {std::string str(highlightcolor); return str; }
    inline void set_highlightcolor(std::string s){strcpy(highlightcolor,s.c_str());};

    //incnfile 
    inline std::string get_incnfile() { std::string str(incnfile); return str; }
    inline void set_incnfile(std::string s){ strcpy(incnfile, s.c_str()); };

    //finalcnfile 
    inline std::string get_finalcnfile(){ std::string str(finalcnfile); return str; };
    inline void set_finalcnfile(std::string s){ strcpy(finalcnfile, s.c_str()); };

    //crystalstructure
    inline std::string get_crystalstructure(){ std::string str(crystalstructure); return str; };
    inline void set_crystalstructure(std::string s){ strcpy(crystalstructure, s.c_str()); };

    //latticeconst
    inline std::vector<double> get_latticeconst() {
        return std::vector<double>(latticeconst, latticeconst+3); }
    inline void set_latticeconst(std::vector<double> v){ 
        latticeconst[0] = v[0]; latticeconst[1] = v[1]; latticeconst[2] = v[2];
    }

    //latticesize
    inline std::vector<double> get_latticesize() { 
        std::vector<double> v;
        for(int i = 0;i<3;i++) for(int j =0;j<4;j++) v.push_back(latticesize[i][j]); 
        return v;
    }
    inline void set_latticesize(std::vector<double> v){ 
        int cnt = 0;
        for(int i = 0;i<3;i++) for(int j =0;j<4;j++) latticesize[i][j] = v[cnt++];
    }

    //rotateangles
    inline std::vector<double> get_rotateangles() { 
        std::vector<double> v;
        for(int i = 0;i<5;i++)  v.push_back(rotateangles[i]); 
        return v;
    }
    inline void set_rotateangles(std::vector<double> v){ 
        for(int i = 0;i<5;i++) rotateangles[i] = v[i];
    }

    //plot_color_windows
    inline std::vector<double> get_plot_color_windows() { 
        std::vector<double> v;
        for(unsigned int i = 0;i<100;i++)  v.push_back(plot_color_windows[i]); 
        return v;
    }
    inline void set_plot_color_windows(std::vector<double> v){ 
        for(unsigned int i = 0;i<v.size();i++) plot_color_windows[i] = v[i];
    }

    //plot_limits
    inline std::vector<double> get_plot_limits() { 
        std::vector<double> v;
        for(unsigned int i = 0;i<100;i++)  v.push_back(plot_color_windows[i]); 
        return v;
    }
    inline void set_plot_limits(std::vector<double> v){ 
        for(unsigned int i = 0;i<v.size();i++) plot_limits[i] = v[i];
    }

    //atomcolor
    inline std::vector<std::string, std::allocator<std::string>> get_atomcolor() { 
        std::vector<std::string, std::allocator<std::string>> v;
    //    for(int i = 0;i<MAXSPECIES;i++) { 
    //        std::string str(atomcolor[i]);
    //        //v.push_back(std::string str(atomcolor[i])); 
    //    }
        return v;
    }
    inline void set_atomcolor(std::vector<std::string, std::allocator<std::string> > v) { 
        for(int i = 0;i<nspecies;i++) strcpy(atomcolor[i], v[i].c_str());
    }

    //atomradius
    inline std::vector<double> get_atomradius() { 
        std::vector<double> v;
        for(int i = 0;i<MAXSPECIES;i++) v.push_back(atomradius[i]); 
        return v;
    }
    inline void set_atomradius(std::vector<double> v) { 
        for(int i = 0;i<MAXSPECIES;i++) atomradius[i] = v[i];
    }

};

#endif //_MD_H
