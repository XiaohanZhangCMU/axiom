/*
  md.cpp
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Mon Jan  3 18:28:35 2011

  FUNCTION  :  Easy-to-use MD simulation Framefork

  Featuring :  1. Scripting input
               2. X-window display
               3. automatic simulation log
               4. convenient configuration and output file handle
               5. perfect lattice and dislocation config creator
               6. conjugate-gradient relaxation
               7. NVT MD simulation (Nose-Hoover thermostat)
               8. NPH and NPT simulation (Parrinello-Rahman + Nose-Hoover)
               9. Cell-Verlist combined neighbor list

  This is a single CPU code.
  Previous shared-memory option removed.
  May implment MPI parallelization in the future.

  TO DO:
  1. implement new potential function (call SW, FS, EAM, etc)
  2. implement new neighbor list function (call Verlist)
  3. implement new integrator (call Gear6 or VVerlet)

  4. re-organize example script files
  5. cannot run in tcl mode if script file too long

  new link list for neighbor list construction
  benchmark Si NW 47580 atoms meam-lammps_mpicc
        link list (new)      cell list (old)
  ------------------------------------------
  1 cpu   103.15 s             99.5 s
  2 cpu    62.70 s             59.5 s

History:
2007.07.10 (CW & WC) add refreshneighborlist(); in EAMFrame::potential() in eam.cpp
2007.04.17 (WC) Add option savecfg, savecfgfreq
2007.04.13 (WF) new definition of meam_force_ in meam-lammps.h
                meam-lammps.cpp
                meam_force.f
2007.03.20 (KW) added Nose-Hoover chain for velocity Verlet,
                NVTC_VVerlet_Explicit()
2007.03.20 (WC) added Nose Hoover chain, ensemble type "NVTC", "NPTC"
                introduced MASSCONVERT from g/mol to eV ps^2/A^2
                new thermal NHMass, which if nonzero, replaces vt2
2007.03.15 (WC) added continue_curstep
----
2007.02.28 (WC) comment out repetition in cstr_potential_wrapper()
2007.02.19 (KW) Implement explicit NPH integrator,
                NPH_VVerlet_Explicit_1()
2007.02.14 (KW) Add _VH*=factor in multiplyvelocity()
2007.02.12 (KW) Correct Stormer Verlet integrator. integrator.cpp
                NVT_VVerlet_Explicit_1(), added VVerlet_Get_h2
2007.02.05 (KW)
           1. calcoutput() and pf.write(this) interchanged in run()
           2. small change in PropFile::describe()
           3. calcoutput() is added in the last part of run().
2007.01.23 (KW)
           Implement Liouville-formulation of explicit integrator for NVT.
           NVT_VVerlet_Explicit_2() now working
2007.01.22 (KW)
           1. New variable "step0" in run()
           2. Add zeta and zetav in CNFile:writeblock(void *p)
              and CNFile:readblock(void *p)
           3. New func. VVerlet_Get_s2() in integrators.cpp
2007.01.19 (KW) 1. small change in run(), put savecn() after step()
                2. NbrList_needrefresh() changed.
2007.01.16 (KW) modification in src/Makefile
2006.08.14 (KW) call_potential(), calcprop() added in run()
*/

#include "md.h"

/********************************************************************/
/* Begining of class MDFrame definition */
/********************************************************************/


/********************************************************************/
/* Initialize script file parser */
/********************************************************************/
void MDFrame::clear_input()
{
    int i;
    for(i = 0;i < MAXINPUTLEN; i++)
       input[i] = 0;
}

void MDFrame::initvars()
{
    int i;
    
    INFO("MDFrame::initvars()");
    
    strcpy(command,"echo Hello World");
    
    strcpy(incnfile,"con.cn");
    strcpy(intercnfile,"inter.cn");
    strcpy(FFScnfile,"FFS.cn");
    strcpy(intercfgfile,"auto.cfg");
    strcpy(finalcnfile,"final.cn");
    strcpy(continuecnfile,"continue.cn");
    strcpy(potfile,"../Mo.pot");

    strcpy(ensemble_type,"NVE");
    strcpy(integrator_type,"Gear6");
    strcpy(initvelocity_type,"Uniform");
    strcpy(zerorot,"none");

    strcpy(outpropfile,"prop.out");
#ifdef _CALPHONONSPECTRUM
    strcpy(phonondispfile,"dispersion.m");
#endif

    strcpy(bondcolor,"red");
    strcpy(highlightcolor,"purple");

    for(i=0;i<MAXSPECIES;i++)
    {
        atomradius[i]=0.1;
        strcpy(atomcolor[i],"orange");
        strcpy(element[i],"Mo");
    }
    for(i=0;i<MAXCOLORS;i++)
    {
        strcpy(colornames[i],"gray50");
    }

    conj_fixboxvec.clear();

    strcpy(element[1],"Mo");
    strcpy(element[2],"Si");
    strcpy(element[3],"Cu");

    strcpy(relaxation_algorithm, "ZXCGR");

    atomeyerepeat[0]=0;

    extforce[0] = 0;
    _SIGMA.clear();
    _H.clear();
    _H0.clear();
    
    clear_input();

    for(i=0;i<MAXNHCLEN;i++) NHMass[i] = 0;

    for(i=0;i<NEBSPECSIZE;i++) nebspec[i] = 0;
    for(i=0;i<ANNEALSPECSIZE;i++) annealspec[i] = 0;
    for(i=0;i<SURFTENSIONSPECSIZE;i++) surftensionspec[i] = 0;
    for(i=0;i<ENERGYTHRESHOLDSIZE;i++) energythreshold[i] = 0;

    for(i=0;i<TORQUESPECSIZE;i++) torquespec[i] = 0;
    for(i=0;i<BENDSPECSIZE;i++) bendspec[i] = 0;
}





/********************************************************************/
/* Coordinate transformation */
/********************************************************************/

void MDFrame::SHtoR()
{ for(int i=0;i<_NP;i++) _R[i]=_H*_SR[i]; }
void MDFrame::RHtoS()
{ Matrix33 hinv=_H.inv();
 for(int i=0;i<_NP;i++) _SR[i]=hinv*_R[i]; }
void MDFrame::VSHtoVR()
{ for(int i=0;i<_NP;i++) _VR[i]=_H*_VSR[i]; }
void MDFrame::VRHtoVS()
{ Matrix33 hinv=_H.inv();
 for(int i=0;i<_NP;i++) _VSR[i]=hinv*_VR[i]; }
void MDFrame::RtoR0()
{ int i; SHtoR(); for(i=0;i<_NP;i++) _R0[i]=_R[i]; }
void MDFrame::R0toR()
{ int i; for(i=0;i<_NP;i++) _R[i]=_R0[i]; RHtoS(); }
void MDFrame::clearR0()
{ int i; for(i=0;i<_NP;i++) _R0[i].clear(); }

bool MDFrame::Bond(int I, int J) const
{/* Returns true if bond(i,j) is stored in I */
    return (I<=J)? ((I^J)&1)  : !((I^J)&1);
}

/* Memory Allocation */
void MDFrame::Alloc()
{
    int i, size;
    size=_NP*allocmultiple;

    /* _SAVEMEMORY == 9 : can allocate structure and write to cn file */ 
    Realloc(_SR,Vector3,size);
    Realloc(fixed,int,size);

    memset(fixed,0,sizeof(int)*size);
    bindvar("SR",_SR,DOUBLE);
    bindvar("fixed",fixed,INT);

    if(_SAVEMEMORY>=9) return;
    /* The following is only allocated if _SAVEMEMORY < 9 */

    /* _SAVEMEMORY == 8 : can allocate structure and write real space coordinates to file */ 
    Realloc(_R,Vector3,size);
    Realloc(_R0,Vector3,size);

    bindvar("R",_R,DOUBLE);
    bindvar("R0",_R0,DOUBLE);

    if(_SAVEMEMORY>=8) return;
    /* The following is only allocated if _SAVEMEMORY < 8 */

    /* _SAVEMEMORY == 7 : can perform energy minimization */ 
    Realloc(_F,Vector3,size);
    Realloc(_EPOT_IND,double,size);
    Realloc(species,int,size);
    
    memset(_EPOT_IND,0,sizeof(double)*size);
    memset(species,0,sizeof(int)*size);

    bindvar("F",_F,DOUBLE);
    bindvar("EPOT_IND",_EPOT_IND,DOUBLE);
    bindvar("species",species,INT);

    if(_SAVEMEMORY>=7) return;
    /* The following is only allocated if _SAVEMEMORY < 7 */

    /* _SAVEMEMORY == 6 : can perform single-component MD simulation */ 
    Realloc(_VR,Vector3,size);
    Realloc(_VSR,Vector3,size);

    /*memset(_VSR,0,sizeof(Vector3)*size);*/    // See KW Kang 2007.06.01

    bindvar("VR",_VR,DOUBLE);
    bindvar("VSR",_VSR,DOUBLE);

    if(_SAVEMEMORY>=6) return;
    /* The following is only allocated if _SAVEMEMORY < 6 */

    /* _SAVEMEMORY == 5 : can perform multi-species simulation and compute CSD parameter */ 
    Realloc(image,int,size);
    Realloc(_TOPOL,double,size);
    Realloc(group,int,size);

    //memset(image,0,sizeof(int)*size);
    for(i=0;i<size;i++) image[i]=-1;

    memset(_TOPOL,0,sizeof(double)*size);
    memset(group,0,sizeof(int)*size);

    bindvar("image",image,INT);
    bindvar("TOPOL",_TOPOL,DOUBLE);
    bindvar("group",group,INT);

    if(_SAVEMEMORY>=5) return;
    /* The following is only allocated if _SAVEMEMORY < 5 */

    /* _SAVEMEMORY == 4 : can perform free energy calculation and apply external force */ 
    Realloc(_F0,Vector3,size);
    Realloc(_Fext,Vector3,size);

    memset(_Fext,0,sizeof(Vector3)*size);
    bindvar("F0",_F0,DOUBLE);
    bindvar("Fext",_Fext,DOUBLE);

    if(_SAVEMEMORY>=1) return;
    /* The following is only allocated if _SAVEMEMORY < 1 */
    
    /* _SAVEMEMORY == 1 : (default) can perform torsion and bending simulation */ 
    Realloc(_EPOT_RMV,double,size);
    Realloc(_VIRIAL_IND,Matrix33,size);
    Realloc(_TORQUE_IND,double,size);
    Realloc(_BENDMOMENT_IND,double,size);

    memset(_EPOT_RMV,0,sizeof(double)*size);
    memset(_VIRIAL_IND,0,sizeof(Matrix33)*size);
    memset(_TORQUE_IND,0,sizeof(double)*size);
    memset(_BENDMOMENT_IND,0,sizeof(double)*size);

    bindvar("EPOT_RMV",_EPOT_RMV,DOUBLE);
    bindvar("VIRIAL_IND",_VIRIAL_IND,DOUBLE);
    bindvar("TORQUE_IND",_TORQUE_IND,DOUBLE);
    bindvar("BENDMOMENT_IND",_BENDMOMENT_IND,DOUBLE);
    
}


/********************************************************************/                         
/* Conjugate gradient relaxation */
/********************************************************************/                         
void MDFrame::potential_wrapper_fixbox(int *n,double x[],double *f,double df[])
{
    int i, np0, np1;
    Matrix33 h, hinv;
    Vector3 gr;

    __conj_simframe->winplot(__conj_simframe->conj_step);
    __conj_simframe->saveintercn(__conj_simframe->conj_step);
    __conj_simframe->saveintercfg(__conj_simframe->conj_step);
    
    __conj_simframe->conj_step++;

    h=__conj_simframe->_H;
    hinv = h.inv();
    np0=0;
    np1=__conj_simframe->_NP;
    /* assign _SR, _R */
    for(i=np0;i<np1;i++)
    {
        __conj_simframe->_R[i].set(x+i*3);
        __conj_simframe->_SR[i]=hinv*__conj_simframe->_R[i];
    }

    /* multi process functions */
    __conj_simframe->call_potential();
    __conj_simframe->calcprop();
    
#ifndef PAIR_DISPLACEMENT_CONSTRAINTS
#define PAIR_DISPLACEMENT_CONSTRAINTS
#endif
#ifdef  PAIR_DISPLACEMENT_CONSTRAINTS
    int nconstraints, ipt1, ipt2;
    double dFdotb, bx, by, bz;
    Vector3 dF, burgV;

    /* apply pair displacement constraints:  each constraint affects two atoms */
    nconstraints = (int) __conj_simframe->pair_displacement_constraints[0];
    //INFO_Printf(" nconstraints = %d\n", nconstraints);
    
    /* enforce constraint on atom positions (to implement later) */

    /* orthogonalize force (in real space) */
    for(i=0;i<nconstraints;i++)
    { /* applying each constraint separately */
      ipt1 = (int)__conj_simframe->pair_displacement_constraints[i*5+1];
      ipt2 = (int)__conj_simframe->pair_displacement_constraints[i*5+2];
      bx   =      __conj_simframe->pair_displacement_constraints[i*5+3];
      by   =      __conj_simframe->pair_displacement_constraints[i*5+4];
      bz   =      __conj_simframe->pair_displacement_constraints[i*5+5];
      burgV.set(bx,by,bz);

      dF = __conj_simframe->_F[ipt2] - __conj_simframe->_F[ipt1];
      dFdotb = dot( dF, burgV ) ;
      dF = burgV; dF *= (dFdotb / burgV.norm2() / 2.0);

      __conj_simframe->_F[ipt2] -= dF;
      __conj_simframe->_F[ipt1] += dF;
    }

#endif
    
    for(i=np0;i<np1;i++)
    {
        if((__conj_simframe->fixed[i])
           ||(__conj_simframe->conj_fixatoms))
            gr.clear();
        else
            gr=__conj_simframe->_F[i] * (-1.0);

        gr.copytoarray(df+i*3);   //forces to _R are obtained
    }
    *f=__conj_simframe->_EPOT;
    
    /* computing residual gradient */
    __conj_simframe->conj_g2res=0;
    for(i=0;i<*n;i++)
        __conj_simframe->conj_g2res+=df[i]*df[i];

}

void MDFrame::potential_wrapper(int *n,double x[],double *f,double df[])
{/* ensure Sync in the end */
    int i, j, np0, np1;
    Matrix33 h,htran,htraninv,s,hinv,hinvtran,hsigma,sigmahtranh;
    double pressure, pres,vb;
    Vector3 gs;

    __conj_simframe->winplot(__conj_simframe->conj_step);
    __conj_simframe->saveintercn(__conj_simframe->conj_step);
    __conj_simframe->saveintercfg(__conj_simframe->conj_step);
    
    __conj_simframe->conj_step++;
    __conj_simframe->_H.set(x);

    /* x = H / HPRECOND */
    __conj_simframe->_H*=__conj_simframe->_HPRECOND;
    
    h=__conj_simframe->_H;
    htran=h.tran();
    htraninv=htran.inv();

    np0=0;
    np1=__conj_simframe->_NP;
    /* assign _SR, _R */
    for(i=np0;i<np1;i++)
    {
        __conj_simframe->_SR[i].set(x+i*3+9);
        __conj_simframe->_SR[i].x*=__conj_simframe->_XPRECOND;
        __conj_simframe->_SR[i].y*=__conj_simframe->_YPRECOND;
        __conj_simframe->_SR[i].z*=__conj_simframe->_ZPRECOND;
        __conj_simframe->_R[i]=h*__conj_simframe->_SR[i];
    }

    /* multi process functions */
    __conj_simframe->call_potential();
    __conj_simframe->calcprop();
    vb=0;
    
    if(__conj_simframe->conj_fixbox) __conj_simframe->_GH.clear();
    else
    {
        /* fixshape relaxation still need to be implemented */
        if(__conj_simframe->conj_fixshape)
        {
            pressure=__conj_simframe->_VIRIAL.trace()/3;
            __conj_simframe->_GH.clear();
        }
        else
        {
            /* PR in finite stress _EXTSTRESS */
            pres=__conj_simframe->_EXTSTRESS.trace()/3;
            pres*=__conj_simframe->_EXTSTRESSMUL;
            pres+=__conj_simframe->_EXTPRESSADD;
            pres/=160.2e3;/*convert from MPa to eV/A^3 */
            s=__conj_simframe->_VIRIAL;
            s[0][0]-=__conj_simframe->_OMEGA*(1-__conj_simframe->_VACUUMRATIO)*pres;
            s[1][1]-=__conj_simframe->_OMEGA*(1-__conj_simframe->_VACUUMRATIO)*pres;
            s[2][2]-=__conj_simframe->_OMEGA*(1-__conj_simframe->_VACUUMRATIO)*pres;
            
            hinv=h.inv(); hinvtran=hinv.tran();
            s=s*hinvtran;
            hsigma=h*__conj_simframe->_SIGMA;
            __conj_simframe->_GH=s-hsigma;
            __conj_simframe->_GH*=-1;
            for(i=0;i<3;i++)
                for(j=0;j<3;j++)
                    if(__conj_simframe->conj_fixboxvec[i][j]==1)
                        __conj_simframe->_GH[i][j]=0;
            sigmahtranh=(__conj_simframe->_SIGMA*htran)*h;
            vb=pres*__conj_simframe->_OMEGA*(1-__conj_simframe->_VACUUMRATIO)+sigmahtranh.trace()*0.5;
        }
    }

    /* dE/dx = _GH * HPRECOND */
    __conj_simframe->_GH*=__conj_simframe->_HPRECOND;
    
    __conj_simframe->_GH.copytoarray(df); // forces to _H are obtained

#ifndef PAIR_DISPLACEMENT_CONSTRAINTS
#define PAIR_DISPLACEMENT_CONSTRAINTS
#endif
#ifdef  PAIR_DISPLACEMENT_CONSTRAINTS
    int nconstraints, ipt1, ipt2;
    double dFdotb, bx, by, bz;
    Vector3 dF, burgV;

    /* apply pair displacement constraints:  each constraint affects two atoms */
    nconstraints = (int) __conj_simframe->pair_displacement_constraints[0];
    //INFO_Printf(" nconstraints = %d\n", nconstraints);
    
    /* enforce constraint on atom positions (to implement later) */

    /* orthogonalize force (in real space) */
    for(i=0;i<nconstraints;i++)
    { /* applying each constraint separately */
      ipt1 = (int)__conj_simframe->pair_displacement_constraints[i*5+1];
      ipt2 = (int)__conj_simframe->pair_displacement_constraints[i*5+2];
      bx   =      __conj_simframe->pair_displacement_constraints[i*5+3];
      by   =      __conj_simframe->pair_displacement_constraints[i*5+4];
      bz   =      __conj_simframe->pair_displacement_constraints[i*5+5];
      burgV.set(bx,by,bz);

      dF = __conj_simframe->_F[ipt2] - __conj_simframe->_F[ipt1];
      dFdotb = dot( dF, burgV ) ;
      dF = burgV; dF *= (dFdotb / burgV.norm2() / 2.0);

      __conj_simframe->_F[ipt2] -= dF;
      __conj_simframe->_F[ipt1] += dF;
    }

#if 0 /* debug */
    for(i=0;i<nconstraints;i++)
    { /* applying each constraint separately */
      ipt1 = (int)__conj_simframe->pair_displacement_constraints[i*5+1];
      ipt2 = (int)__conj_simframe->pair_displacement_constraints[i*5+2];
      bx   =      __conj_simframe->pair_displacement_constraints[i*5+3];
      by   =      __conj_simframe->pair_displacement_constraints[i*5+4];
      bz   =      __conj_simframe->pair_displacement_constraints[i*5+5];
      burgV.set(bx,by,bz);

      dF = __conj_simframe->_F[ipt2] - __conj_simframe->_F[ipt1];
      dFdotb = dot( dF, burgV ) ;
      INFO_Printf(" constraint[%d] = (%d %d %e %e %e), dFdotb = %e (should be zero)\n", 
                    i, ipt1, ipt2, bx, by, bz, dFdotb);
    }
#endif
#endif
    
    for(i=np0;i<np1;i++)
    {
        if((__conj_simframe->fixed[i])
           ||(__conj_simframe->conj_fixatoms))
            gs.clear();
        else
            gs=-htran*__conj_simframe->_F[i];

        gs.x*=__conj_simframe->_XPRECOND;
        gs.y*=__conj_simframe->_YPRECOND;
        gs.z*=__conj_simframe->_ZPRECOND;
        gs.copytoarray(df+i*3+9);   //forces to _SR are obtained
    }
    *f=__conj_simframe->_EPOT + vb;
    
    /* computing residual gradient */
    __conj_simframe->conj_g2res=0;
    for(i=0;i<*n;i++)
        __conj_simframe->conj_g2res+=df[i]*df[i];
}

void MDFrame::potential_wrapper_c(int n,double x[],double *f,double df[])
{/* ensure Sync in the end */

  if (__conj_simframe->conj_fixbox == 1) {
    potential_wrapper_fixbox(&n,x,f,df);
  }
  else {
    potential_wrapper(&n,x,f,df);
  }
}

void MDFrame::potential_wrapper_rigid_translate(int n,double x[],double *f,double df[])
{
    int i, j, np0, np1, ngrp, grp_id, match;
    int ngrpatomA, ngrpatomB;
    Matrix33 h, hinv;
    Vector3 gr, F;

    __conj_simframe->winplot(__conj_simframe->conj_step);
    __conj_simframe->saveintercn(__conj_simframe->conj_step);
    __conj_simframe->saveintercfg(__conj_simframe->conj_step);
    
    __conj_simframe->conj_step++;

    ngrp = (int)__conj_simframe->input[0];
    //INFO("ngrp= "<<ngrp);

    h=__conj_simframe->_H;
    hinv = h.inv();
    np0=0;
    np1=__conj_simframe->_NP;
    /* assign _SR, _R */
    for(i=np0;i<np1;i++)
    {
        __conj_simframe->_R[i].set(x+i*3);
        __conj_simframe->_SR[i]=hinv*__conj_simframe->_R[i];
    }
    __conj_simframe->call_potential();

    ngrpatomA = 0; ngrpatomB = 0;
    F.set(0.0,0.0,0.0); 
    for(i=np0;i<np1;i++)
    {
        match = 0;
        for(j=0;j<ngrp;j++)
        {
            grp_id = (int)__conj_simframe->input[j+1];
            if ( __conj_simframe->group[i] == grp_id)
            {
                match = 1; break;
            }
        }
      
        if (match == 1)
        {
            F+=__conj_simframe->_F[i];
            ngrpatomA++;
        } else 
        {
            ngrpatomB++;
        }
    }

    /* zero out the force in desired direction */
    if (__conj_simframe->conj_fixdir[0]==1) F.x=0;
    if (__conj_simframe->conj_fixdir[1]==1) F.y=0;
    if (__conj_simframe->conj_fixdir[2]==1) F.z=0;

    for(i=np0;i<np1;i++)
    {
        match = 0;
        for(j=0;j<ngrp;j++)
        {
            grp_id = (int)__conj_simframe->input[j+1];
            if ( __conj_simframe->group[i] == grp_id)
            {
                match = 1; break;
            }
        }
      
        if (match == 1)
            gr=F*(-1.0/ngrpatomA);
        else
            gr=F*(1.0/ngrpatomB);

        gr.copytoarray(df+i*3);
    }
    *f=__conj_simframe->_EPOT;
    
    /* computing residual gradient */
    __conj_simframe->conj_g2res=0;
    for(i=0;i<n;i++)
        __conj_simframe->conj_g2res+=df[i]*df[i];
}

int MDFrame::relax()
{
    /* This function will be executed by multiprocesses synchronously */
    
  int i,n,_NP0,_NP1;
  double *conj_space,*conj_x,*conj_g,*conj_w, pres;
  Matrix33 prp, hinv, hinvtran;
    
  INFO("Relax");

  if (conj_fixbox == 1) {

    _NP0=0;
    _NP1=_NP;
    n=_NP*3;
    conj_step=0;
    
    conj_space=NULL;
    Realloc(conj_space,double,n*8);
    
    conj_x=conj_space;    /* Array length n (=3*NP+9) for H and _SR */
    conj_g=conj_space+n;  /* Array length n for corrsponding forces to H and _SR */
    conj_w=conj_space+n*2;/* Array length 6n that will be used as a buffer */

    SHtoR();
    for(i=_NP0;i<_NP1;i++)
        _R[i].copytoarray(conj_x+i*3);

    __conj_simframe=this;

    /* All processes enter CGRelax */
    if ( strcmp(relaxation_algorithm, "ZXCGR") == 0 )
    {
        CGRelax(potential_wrapper_c,n,
                conj_ftol,conj_fevalmax,conj_dfpred,
                conj_x,conj_g,&conj_f,conj_w);
    } else if ( strcmp(relaxation_algorithm, "PRPLUS") == 0 )
    {
        /* CG relax by Polak-Ribiere plus method */
        CGRelax_PRplus(potential_wrapper_c,n,conj_etol,conj_ftol,
                       conj_itmax,conj_fevalmax,conj_dfpred,
                       conj_x,conj_g,&conj_f,conj_w);
    } else
    {
        ERROR("invalid relaxation_algorithm. Valid choices are ZXCGR or PRPLUS");
    }

    for(i=_NP0;i<_NP1;i++)
    {
        _R[i].set(conj_x+i*3);
    }
    RHtoS();
    INFO_Printf("\tconj_f = %24.18g\n",conj_f);

    free(conj_space);

  } 
  else {

    _NP0=0;
    _NP1=_NP;
    n=_NP*3+9;
    conj_step=0;
    
    conj_space=NULL;
    Realloc(conj_space,double,n*8);
    
    conj_x=conj_space;    /* Array length n (=3*NP+9) for H and _SR */
    conj_g=conj_space+n;  /* Array length n for corrsponding forces to H and _SR */
    conj_w=conj_space+n*2;/* Array length 6n that will be used as a buffer */

    _H.copytoarray(conj_x);

    /* preconditioning */
    for(i=0;i<9;i++) conj_x[i]/=_HPRECOND;

    for(i=_NP0;i<_NP1;i++)
        _SR[i].copytoarray(conj_x+i*3+9);

    for(i=_NP0;i<_NP1;i++)
    {
        conj_x[9+i*3]  /=_XPRECOND;
        conj_x[9+i*3+1]/=_YPRECOND;
        conj_x[9+i*3+2]/=_ZPRECOND;
    }
    __conj_simframe=this;

    /* prepare _SIGMA for PR relaxation */
    pres=_EXTSTRESS.trace()/3;
    prp=_EXTSTRESS;
    prp[0][0]-=pres;
    prp[1][1]-=pres;
    prp[2][2]-=pres;
    prp*=_EXTSTRESSMUL;
    pres*=_EXTSTRESSMUL;
    pres+=_EXTPRESSADD;

    /* Use H0 as initial H if H0 is specified */
    if(fabs(_H0.det())<1e-10)
    {
        INFO("H0 has not been specified.  Use current H as reference.");
        hinv=_H.inv();
        hinvtran=hinv.tran();
        _OMEGA=_H.det();
    }
    else
    {
        INFO("H0 has been specified as reference.");
        hinv=_H0.inv();
        hinvtran=hinv.tran();
        _OMEGA=_H0.det();
    }
    _SIGMA=((hinv*prp)*hinvtran)*_OMEGA*(1-_VACUUMRATIO);
    _SIGMA/=160.2e3;/* convert MPa to eV/A^3 */

    /* All processes enter CGRelax */
    if ( strcmp(relaxation_algorithm, "ZXCGR") == 0 )
    {
        CGRelax(potential_wrapper_c,n,
                conj_ftol,conj_fevalmax,conj_dfpred,
                conj_x,conj_g,&conj_f,conj_w);
    } else if ( strcmp(relaxation_algorithm, "PRPLUS") == 0 )
    {
        /* CG relax by Polak-Ribiere plus method */
        CGRelax_PRplus(potential_wrapper_c,n,conj_etol,conj_ftol,
                       conj_itmax,conj_fevalmax,conj_dfpred,
                       conj_x,conj_g,&conj_f,conj_w);
    } else
    {
        ERROR("invalid relaxation_algorithm. Valid choices are ZXCGR or PRPLUS");
    }


    _H.set(conj_x);

    /* preconditioing */
    _H*=_HPRECOND;

    for(i=_NP0;i<_NP1;i++)
    {
        _SR[i].set(conj_x+i*3+9);
        _SR[i].x*=_XPRECOND;
        _SR[i].y*=_YPRECOND;
        _SR[i].z*=_ZPRECOND;
    }
    INFO_Printf("\tconj_f = %24.18g\n",conj_f);

    free(conj_space);
  }

  return 0;
}

void MDFrame::strelax()
{
    int i, ignoreE;
    double Ediff, Eold, Edrop;
    
    INFO("Run Steepest Descent Relax");

    Edrop = input[0];    
    ignoreE = (int) input[1];
    Ediff = 0.0;
    INFO("Tolerance = "<<Edrop<<" (eV)");
    INFO("Time step = "<<_TIMESTEP);

    call_potential();
    //INFO("Epot = "<<_EPOT);
    SHtoR();
    /* store old configuration and energy */
    if(_SRA==NULL) Realloc(_SRA,Vector3,_NP);
    for(i=0;i<_NP;i++) _SRA[i]=_R[i];
    Eold = _EPOT;
    //INFO("Eold = "<<Eold);

    if(continue_curstep)
        step0 = curstep;
    else
        step0 = 0;
    INFO("curstep="<<curstep<<", Epot="<<_EPOT);
    for(curstep=step0;curstep<totalsteps;curstep++)
    {
        for(i=0;i<_NP;i++)
            _R[i] += _F[i]*_TIMESTEP;

        RHtoS();
        call_potential();
        Ediff = _EPOT - Eold;

        if ((curstep%plotfreq) == 0)
            winplot();

        if (fabs(Ediff) < Edrop)
        {
            INFO("Delta(E) = "<<Ediff<<" (eV)");
            INFO("Energy tolerance = "<<Edrop<<" (eV) achieved ");
            break;        
        }
        else if (Ediff>0)
        {
            INFO("curstep = "<<curstep<<" Energy is increasing, Delta(E)="<<Ediff<<" (eV)");
            if (!ignoreE)
            {
                _TIMESTEP/=10.0;
                /* restore old configuration */
                for(i=0;i<_NP;i++) _R[i]=_SRA[i];
                RHtoS();
                call_potential();
                INFO("Decrease timestep as "<<_TIMESTEP);
                INFO("Restore old configuration.");
            }
        }
        else 
        {
            for(i=0;i<_NP;i++) _SRA[i]=_R[i];
        }
        

        if (((curstep%printfreq)==0)&&(curstep>step0))
            INFO("curstep="<<curstep<<", Epot="<<_EPOT);

        /* store old configuration and energy */
        RtoR0();         
        Eold = _EPOT;
    }
    INFO("curstep="<<curstep<<", Epot="<<_EPOT);
    
    if ((curstep>=totalsteps) && (fabs(Ediff)>=Edrop))
        INFO("Maximum iterations reached!!");
}

void MDFrame::relax_by_group(double *myinput)   // 2007/03/09 Keonwook Kang
{
    /* relax atoms that belong to specified groups
     *
     * input = [NgrpID, grpID_1 grpID_2 ...]
     *
     * NgrpID : number of groups to be relaxed
     * grpID_1, grpID_2, ... : group ID's to be relaxed
     */
    
    int i, j, NgrpID, flag;
    NgrpID = (int) myinput[0];         

    if (NgrpID <= 0 )
    {
        ERROR("Number of groups should be greater than zero.");
        return;
    }
    
    for (i=0; i<_NP; i++)
    {
        if (fixed[i]==1) /* If already fixed, skip */
            continue;
        else
        {
            flag = 1;
            for (j=0; j<NgrpID; j++)
            {// If group[i] is one of groups to be relaxed,
                if (group[i]==(int) myinput[j+1]) 
                {// flag becomes zero.
                    flag = 0; break;
                }
            }// Atoms with flag=1 should be fixed.
            if (flag) fixed[i]=1;
        }
    }
    INFO("test");
    relax(); //call relax
}

int MDFrame::rigid_relaxation()
{
  int i,n;
  double *conj_space,*conj_x,*conj_g,*conj_w;
  Vector3  disp;
  Matrix33 prp, hinv, hinvtran;
    
  INFO("Rigid Relaxation");

  n=_NP*3;
  conj_step=0;

  conj_space=NULL;
  Realloc(conj_space,double,n*8);
    
  conj_x=conj_space;    /* Array length n for rigid displacement dx, dy, and dz */
  conj_g=conj_space+n;  /* Array length n for corrsponding forces to dx, dy, and dz */
  conj_w=conj_space+n*2;/* Array length 6n that will be used as a buffer */

  SHtoR();
  //memset(conj_x,0.0,sizeof(double)*n);
  for(i=0;i<_NP;i++)
    _R[i].copytoarray(conj_x+i*3);

  __conj_simframe=this;

  /* All processes enter CGRelax */
  if ( strcmp(relaxation_algorithm, "ZXCGR") == 0 )
  {
      CGRelax(potential_wrapper_rigid_translate,n,
              conj_ftol,conj_fevalmax,conj_dfpred,
              conj_x,conj_g,&conj_f,conj_w);
  } else if ( strcmp(relaxation_algorithm, "PRPLUS") == 0 )
  {
      /* CG relax by Polak-Ribiere plus method */
      CGRelax_PRplus(potential_wrapper_rigid_translate,n,conj_etol,conj_ftol,
                     conj_itmax,conj_fevalmax,conj_dfpred,
                     conj_x,conj_g,&conj_f,conj_w);
  } else
  {
      ERROR("invalid relaxation_algorithm. Valid choices are ZXCGR or PRPLUS");
  }

  //INFO("conj_x = ["<<conj_x[0]<<", "<<conj_x[1]<<", "<<conj_x[2]<<"]");
  RHtoS();
  INFO_Printf("\tconj_f = %24.18g\n",conj_f);

  free(conj_space);
  return 0;
}

/* Constrained Conjugate gradient relaxation */
/* cannot run with multi-processors */
/* not sure about pre-conditioning */
/* no fixed atoms */
void MDFrame::cstr_potential_wrapper(int *n,double x[],double *f,double df[])
{/* ensure Sync in the end */
    int i, j, np0, np1, ipt;
    Matrix33 h,htran,htraninv,s,hinv,hinvtran,hsigma,sigmahtranh, G;
    double pressure, pres, vb, gn, dist;
    Vector3 gs, ds0, dF, dsstarstar, dr0;

    __conj_simframe->winplot(__conj_simframe->conj_step);
    __conj_simframe->saveintercn(__conj_simframe->conj_step);
    __conj_simframe->saveintercfg(__conj_simframe->conj_step);
    
    __conj_simframe->conj_step++;
    __conj_simframe->_H.set(x);

    /* x = H / HPRECOND */
    __conj_simframe->_H*=__conj_simframe->_HPRECOND;
    
    h=__conj_simframe->_H;
    htran=h.tran();
    htraninv=htran.inv();

    np0=0;
    np1=__conj_simframe->_NP;
    /* assign _SR, _R in Shared memory */
    for(i=np0;i<np1;i++)
    {
        __conj_simframe->_SR[i].set(x+i*3+9);
        __conj_simframe->_SR[i].x*=__conj_simframe->_XPRECOND;
        __conj_simframe->_SR[i].y*=__conj_simframe->_YPRECOND;
        __conj_simframe->_SR[i].z*=__conj_simframe->_ZPRECOND;
        __conj_simframe->_R[i]=h*__conj_simframe->_SR[i];
    }

    /* multi process functions */
    __conj_simframe->call_potential();
    __conj_simframe->calcprop();
    vb=0;
    
    if(__conj_simframe->conj_fixbox) __conj_simframe->_GH.clear();
    else
    {
        /* fixshape relaxation still need to be implemented */
        if(__conj_simframe->conj_fixshape)
        {
            pressure=__conj_simframe->_VIRIAL.trace()/3;
            __conj_simframe->_GH.clear();
        }
        else
        {
            /* PR in finite stress _EXTSTRESS */
            pres=__conj_simframe->_EXTSTRESS.trace()/3;
            pres*=__conj_simframe->_EXTSTRESSMUL;
            pres+=__conj_simframe->_EXTPRESSADD;
            pres/=160.2e3;/*convert from MPa to eV/A^3 */
            s=__conj_simframe->_VIRIAL;
            s[0][0]-=__conj_simframe->_OMEGA*(1-__conj_simframe->_VACUUMRATIO)*pres;
            s[1][1]-=__conj_simframe->_OMEGA*(1-__conj_simframe->_VACUUMRATIO)*pres;
            s[2][2]-=__conj_simframe->_OMEGA*(1-__conj_simframe->_VACUUMRATIO)*pres;
            
            hinv=h.inv(); hinvtran=hinv.tran();
            s=s*hinvtran;
            hsigma=h*__conj_simframe->_SIGMA;
            __conj_simframe->_GH=s-hsigma;
            __conj_simframe->_GH*=-1;
            for(i=0;i<3;i++)
                for(j=0;j<3;j++)
                    if(__conj_simframe->conj_fixboxvec[i][j]==1)
                        __conj_simframe->_GH[i][j]=0;
            sigmahtranh=(__conj_simframe->_SIGMA*htran)*h;
            vb=pres*__conj_simframe->_OMEGA*(1-__conj_simframe->_VACUUMRATIO)+sigmahtranh.trace()*0.5;
        }
    }

    /* dE/dx = _GH * HPRECOND */
    __conj_simframe->_GH*=__conj_simframe->_HPRECOND;
    
    __conj_simframe->_GH.copytoarray(df);

    /* apply constraint: a single constraint affecting many atoms */
    /* compute out-of-plane gradient */
    gn=0; G = htran*h; dist = 0;
    for(i=1;i<=__conj_simframe->constrainatoms[0];i++)
    {
        ipt=__conj_simframe->constrainatoms[i];
        ds0=(__conj_simframe->_SR2[ipt]-__conj_simframe->dscom12)-__conj_simframe->_SR1[ipt];
        dsstarstar = G*ds0;
        gs=htran*__conj_simframe->_F[ipt];
        gn+=dot(gs,dsstarstar);
        dist+=dot(dsstarstar,dsstarstar);
    }
    /* gn : out-of-constraint-plain component */
    gn/=dist;

    /* orthogonalize force */
    for(i=1;i<=__conj_simframe->constrainatoms[0];i++)
    {
        ipt=__conj_simframe->constrainatoms[i];
        ds0=(__conj_simframe->_SR2[ipt]-__conj_simframe->dscom12)-__conj_simframe->_SR1[ipt];
        dsstarstar = h*ds0;
        dF=dsstarstar*gn;
        /* find in-plain force by subtracting out-of-plain force*/
        __conj_simframe->_F[ipt]-=dF;
    }
    
    for(i=np0;i<np1;i++)
    {
        if((__conj_simframe->fixed[i])
           ||(__conj_simframe->conj_fixatoms))
            gs.clear();
        else
            gs=-htran*__conj_simframe->_F[i];
        gs.x*=__conj_simframe->_XPRECOND;
        gs.y*=__conj_simframe->_YPRECOND;
        gs.z*=__conj_simframe->_ZPRECOND;
        gs.copytoarray(df+i*3+9);
    }
    *f=__conj_simframe->_EPOT + vb;
}

void MDFrame::cstr_potential_wrapper_c(int n,double x[],double *f,double df[])
{/* ensure Sync in the end */
//    static int ncalls=0;
    cstr_potential_wrapper(&n,x,f,df);
//    INFO_Printf("%10d : %22.16e\n",ncalls++,*f);
}

int MDFrame::calconstrainS_inst()
{
    int i,ipt;
    double sumds02, ds_star2;
    Vector3 ds, ds0, dr0, ds_star;
    Matrix33 hinv, htraninv;

    if(constrainatoms[0]<=0)
    {
        ERROR("constrainedatoms not set!");
        return -1;
    }
    if(constrainatoms[0]>MAXCONSTRAINATOMS)
    {
        ERROR("Number of constrained atoms is too big!");
        return -1;
    }
    if(_SR1 == NULL || _SR2 == NULL)
    {
        if(_SR1 == NULL)
            ERROR("config1 is not set!");
        else
            ERROR("config2 is not set!");
        return -1;
    }
    
    /* center of mass : sum(m_i*r_i)/sum(m_i)
                         = sum(r_i)/N         if all atoms have same mass.*/
    /* shift of c.o.m between two states A and B:
             sum(r_i^B)/N - sum(r_i^A)/N = sum(r_i^B - r_i^A)/N           */
    dscom12.clear();
    for(i=1;i<=constrainatoms[0];i++)
    {
        ipt=constrainatoms[i];
        dscom12+=_SR2[ipt];
        dscom12-=_SR1[ipt];
    }
    dscom12 /= constrainatoms[0];
    //INFO_Printf("center-of-mass shift: [%f %f %f]\n",dscom12.x,dscom12.y,dscom12.z);

    /* do not shift center of mass of configB */
#if 0    
    for(i=1;i<=constrainatoms[0];i++)
    {
        ipt=constrainatoms[i];
        _SR2[ipt]-=dscom;
    }
#endif
    hinv=_H.inv(); htraninv=hinv.tran(); 
    constrainS_inst=0; sumds02=0; ds_star2=0;
    
    /* compute current constrain value */
    constrain_dist=0; 
    for(i=1;i<=constrainatoms[0];i++)
    {
        ipt=constrainatoms[i];
        ds0=(_SR2[ipt]-dscom12)-_SR1[ipt];
        ds_star=htraninv*ds0;
        ds =_SR[ipt]-_SR1[ipt];
        constrainS_inst+=dot(ds,ds0);
        constrain_dist+=ds0.norm2();
        sumds02+=ds0.norm2();
        ds_star2+=ds_star.norm2();
    }
    constrainS_inst/=constrain_dist;
    //INFO_Printf("(R_B - R_A)*S^star/|S^star| = %25.18e\n",sumds02/sqrt(ds_star2));
    return 0;
}

int MDFrame::applyconstraint() /* Keonwook Kang, May 03, 2007 */
{
    int i,ipt;
    Vector3 ds, ds0, dr0;
    Matrix33 hinv, htraninv;

    /* compute current constrain value */
    calconstrainS_inst();
    INFO_Printf("Before displacement, s=%f constrain_dist=%e\n",constrainS_inst,constrain_dist);

    /* shift current config to desired S */
    for(i=1;i<=constrainatoms[0];i++)
    {
    /* If you enforce interpCN externally, thils loop does no work,
       because (constrainS-constrainS_inst) = 0. */
        ipt=constrainatoms[i];
        ds0=(_SR2[ipt]-dscom12)-_SR1[ipt];
        ds0*=(constrainS-constrainS_inst);
        
        _SR[ipt]+=ds0;
        /* Don't we need to add dscom12? (2007.04.20 Keonwook Kang)
           e.g. _SR[ipt]+=ds0+dscom12;      */
    }
    
    /* compute new constrain value */
    calconstrainS_inst();

    INFO_Printf("after displacement s=%f constrain_dist=%e\n",constrainS_inst,constrain_dist);

    return 0;
}

int MDFrame::constrainedrelax()
{
    /* This function will be executed by multiprocesses synchronously */
    
    int i,n,ipt;
    double *conj_space,*conj_x,*conj_g,*conj_w, pres;
    Matrix33 prp, htran, hinv, hinvtran, G;
    Vector3 ds, ds0, dsstarstar;
    FILE *fp;
    char filename[100];
    
    INFO("Constrained Relax");

    if(constrainatoms[0]<=0)
    {
        ERROR("constrainedatoms not set!");
        return -1;
    }
    if(constrainatoms[0]>MAXCONSTRAINATOMS)
    {
        ERROR("Number of constrained atoms is too big!");
        return -1;
    }
    if(_SR1 == NULL || _SR2 == NULL)
    {
        if(_SR1 == NULL)
            ERROR("config1 is not set!");
        else
            ERROR("config2 is not set!");
        return -1;
    }

    caldscom12();

    htran=_H.tran(); G=htran*_H;
    /* compute current constrain value */
    constrain_dist=0; constrainS_inst=0;
    for(i=1;i<=constrainatoms[0];i++)
    {
        ipt=constrainatoms[i];
        ds0=(_SR2[ipt]-dscom12)-_SR1[ipt];
        dsstarstar = G*ds0;
        ds =_SR[ipt]-_SR1[ipt];
        constrainS_inst+=dot(ds,dsstarstar);
        constrain_dist+=dot(ds0,dsstarstar);
    }
    constrainS_inst/=constrain_dist;
    INFO_Printf("before displacement s=%f constrain_dist=%e\n",constrainS_inst,constrain_dist);

    /* shift current config to desired S */
    for(i=1;i<=constrainatoms[0];i++)
    {
        ipt=constrainatoms[i];
        ds0=(_SR2[ipt]-dscom12)-_SR1[ipt];
        ds0*=(constrainS-constrainS_inst);
        _SR[ipt]+=ds0;
        /* Don't we need to add dscom12? (2007.04.20 Keonwook Kang)
           e.g. _SR[ipt]+=ds0+dscom12;      */
    }

    /* compute new constrain value */
    constrainS_inst=0;
    for(i=1;i<=constrainatoms[0];i++)
    {
        ipt=constrainatoms[i];
        ds0=(_SR2[ipt]-dscom12)-_SR1[ipt];
        dsstarstar=G*ds0;
        ds =_SR[ipt]-_SR1[ipt];
        constrainS_inst+=dot(ds,dsstarstar);
    }
    constrainS_inst/=constrain_dist;
    INFO_Printf("after  displacement s=%f constrain_dist=%e\n",constrainS_inst,constrain_dist);
    
    n=_NP*3+9;
    conj_step=0;
    
    conj_space=NULL;
    Realloc(conj_space,double,n*8);
    
    conj_x=conj_space;    /* Array length n  */
    conj_g=conj_space+n;  /* Array length n  */
    conj_w=conj_space+n*2;/* Array length 6n */

    _H.copytoarray(conj_x);

    /* preconditioning */
    for(i=0;i<9;i++) conj_x[i]/=_HPRECOND;

    //for(i=_NP;i<_NP;i++)  //!Bug
    for(i=0;i<_NP;i++)
        _SR[i].copytoarray(conj_x+i*3+9);

    for(i=0;i<_NP;i++)
    {
        conj_x[9+i*3]  /=_XPRECOND;
        conj_x[9+i*3+1]/=_YPRECOND;
        conj_x[9+i*3+2]/=_ZPRECOND;
    }
    __conj_simframe=this;

    /* prepare _SIGMA for PR relaxation */
    pres=_EXTSTRESS.trace()/3;
    prp=_EXTSTRESS;
    prp[0][0]-=pres;
    prp[1][1]-=pres;
    prp[2][2]-=pres;
    prp*=_EXTSTRESSMUL;
    pres*=_EXTSTRESSMUL;
    pres+=_EXTPRESSADD;

    hinv=_H.inv();
    hinvtran=hinv.tran();
    _SIGMA=((hinv*prp)*hinvtran)*__conj_simframe->_OMEGA*(1-__conj_simframe->_VACUUMRATIO);
    _SIGMA/=160.2e3;/* convert MPa to eV/A^3 */
    
    /* All processes enter CGRelax */
    if ( strcmp(relaxation_algorithm, "ZXCGR") == 0 )
    {
        CGRelax(cstr_potential_wrapper_c,n,
                conj_ftol,conj_fevalmax,conj_dfpred,
                conj_x,conj_g,&conj_f,conj_w);
    } else if ( strcmp(relaxation_algorithm, "PRPLUS") == 0 )
    {
        /* CG relax by Polak-Ribiere plus method */
        CGRelax_PRplus(cstr_potential_wrapper_c,n,conj_etol,conj_ftol,
                       conj_itmax,conj_fevalmax,conj_dfpred,
                       conj_x,conj_g,&conj_f,conj_w);
    } else
    {
        ERROR("invalid relaxation_algorithm. Valid choices are ZXCGR or PRPLUS");
    }

    _H.set(conj_x);

    /* preconditioing */
    _H*=_HPRECOND;

    for(i=0;i<_NP;i++)
    {
        _SR[i].set(conj_x+i*3+9);
        _SR[i].x*=_XPRECOND;
        _SR[i].y*=_YPRECOND;
        _SR[i].z*=_ZPRECOND;
    }
    INFO_Printf("\tconj_f = %24.18g\n",conj_f);

    /* compute final constrain value */
    constrainS_inst=0;
    for(i=1;i<=constrainatoms[0];i++)
    {
        ipt=constrainatoms[i];
        ds0=(_SR2[ipt]-dscom12)-_SR1[ipt];
        dsstarstar = G*ds0;
        ds =_SR[ipt]-_SR1[ipt];
        constrainS_inst+=dot(ds,dsstarstar);
    }
    constrainS_inst/=constrain_dist;
    INFO_Printf("after relaxation   s=%f E=%25.15e \n",constrainS_inst,conj_f);
    
    /* write configuration in NEB format */
    sprintf(filename,"chain%5.3f.out",constrainS_inst);
    fp=fopen(filename,"w");
    for(i=1;i<=constrainatoms[0];i++)
    {
        ipt=constrainatoms[i];
        fprintf(fp,"%25.15e %25.15e %25.15e\n",_R[ipt].x,_R[ipt].y,_R[ipt].z);
    }
    fclose(fp);
        
    free(conj_space);
    return 0;
}



/********************************************************************/
/* Molecular Dynamics simulation */
/********************************************************************/

void MDFrame::call_potential()
{
    double lambda, vn, gn, ds_star2, sigma, theta, fd, fr, rijnorm, y, ewall, fwall;
    int i, j, ipt, jpt;
    Vector3 ds, dr, dF, ds0, dv, gs, ds_star, sij, rij, rij_unit, vij, fij, vrand, si, ri;
    Matrix33 htran, htraninv, hinv;
#ifdef _TORSION_OR_BENDING
    if (_NIMAGES>0)
    {
       if (_TORSIONSIM) copytorqueatoms();
       if (_BENDSIM)    copybendatoms();
    }
#endif
    if (_Dexx>0)
    {
        calDVIRIALDexx();
    }

    if (_ENABLE_SWITCH)
    {
    }
    else {
         potential();
    }
#ifdef _TORSION_OR_BENDING
    if (_NIMAGES>0)
    {
        for(i=0;i<_NIMAGES;i++)
        {
             _EPOT -= _EPOT_IND[_NP0+i];
             _VIRIAL = _VIRIAL - _VIRIAL_IND[_NP0+i];
             _TORQUE -= _TORQUE_IND[_NP0+i];
             _BENDMOMENT -= _BENDMOMENT_IND[_NP0+i];
        }             
    }
#endif
    if (_constrainedMD)
    {
        htran = _H.tran();
        htraninv = htran.inv();
        /* compute out-of-constraint-plane component*/
        gn=0; vn=0; ds_star2=0;
        for(i=1;i<=constrainatoms[0];i++)
        {
            ipt=constrainatoms[i];
            ds0=(_SR2[ipt]-dscom12)-_SR1[ipt];
            ds_star=htraninv*ds0;
            gn+=dot(_F[ipt],ds_star);
            vn+=dot(_VSR[ipt],ds0);
            ds_star2+=ds_star.norm2();
        }
        /* gn, vn : out-of-constraint-plain component */
        gn/=ds_star2;
        vn/=constrain_dist;

        constrainF = gn*sqrt(ds_star2); 
        
        /* orthogonalize force and velocity*/
        /* find in-plain force by subtracting out-of-plain force*/
        for(i=1;i<=constrainatoms[0];i++)
        {
            ipt=constrainatoms[i];
            ds0=(_SR2[ipt]-dscom12)-_SR1[ipt];
            ds_star=htraninv*ds0;
            dF = ds_star*gn;
            dv = ds0*vn;
            _F[ipt]-=dF;
            _VSR[ipt]-=dv;
        }        
    }

    /* add external forces to surface atoms, 
     * for self-consistency, a linear potential is also added to potential energy
     */
    if(extforce[0] > 0)
    {
        for(i=0;i<_NP;i++)
        {
            if( (group[i]>0) && (group[i]<=extforce[0]) )
            {
                Vector3 f;
                f.set(extforce+(group[i]-1)*3+1); f*=forcemul;
                _F[i]+=f;
                _R[i]=_H*_SR[i];
                _EPOT-=dot(f,_R[i]);
            }
        }
    }

    if (_ENABLE_FEXT)
    {
        SHtoR();
        for(i=0;i<_NP;i++)
            _F[i] += _Fext[i];

        if(_SR1==NULL)
        {
            for(i=0;i<_NP;i++)
                _EPOT -= dot(_Fext[i],_R[i]);
        }
        else /* if there is a reference structure */
        {
            for(i=0;i<_NP;i++)
            {
                ds = _SR[i] - _SR1[i];
                dr = _H*ds;
                _EPOT -= dot(_Fext[i],dr);
            }
        }
        _EPOT += _EPOT0_ext;
    }

/* surface tension calculation by separation */
    if (_ENABLE_separation_pot)
    {
        ST_dUdLAM_POT = 0;
	for (i=0;i<_NP;i++)
	{
	    si=_SR[i]; si.subint();
	    ri=_H*si;
	    y=ri.y;
//            if (i==0)
//		INFO_Printf("test y=%f ST_Y0=%f ST_LAMBDA=%f\n",y,ST_Y0,ST_LAMBDA);
	    if ( fabs(y) < ST_Y0 )
	    {
//            INFO_Printf("YESYES!!! id=%d test y=%f ST_Y0=%f ST_LAMBDA=%f\n",i,y,ST_Y0,ST_LAMBDA);
		ewall = ST_K*pow( (pow(ST_Y0,2)-pow(y,2)),2)/pow(ST_Y0,2);
		fwall = ST_K*4*( pow(ST_Y0,2)-pow(y,2) )*y/pow(ST_Y0,2);
		_F[i].y += fwall;
		_EPOT_IND[i] += ewall;
		_EPOT += ewall;
		_VIRIAL[1][1] += fwall*y;
		if (ST_LAMBDA >0)
                {
		    if (ST_step == 1 ) ST_dUdLAM_POT += fwall*pow(ST_Y0,2)/(y*ST_LAMBDA)-2*ewall/ST_LAMBDA;
		    else if (ST_step == 2 ) ST_dUdLAM_POT += ewall/ST_LAMBDA;
		    else ST_dUdLAM_POT = 0;
                }
	    }
        }
    }

/* Disspative Particle Dynamics, Stochastic Thermostat */
    if (_ENABLE_DPD && drand48() < _DPD_RATIO )
    {
	SHtoR();
        VSHtoVR();
        sigma=sqrt(2.0*KB*_TDES*_DPD_fc);
        for(ipt=0;ipt<_NP;ipt++)
        {
	    for(j=0;j<nn[ipt];j++)
	    {
		jpt=nindex[ipt][j];
		if (ipt>jpt) continue;
                sij=_SR[ipt]-_SR[jpt]; sij.subint();
                vij=_VR[ipt]-_VR[jpt];
                rij=_H*sij;
                rijnorm=rij.norm();
                if (rijnorm>_DPD_RCUT) continue;
                rij_unit=rij/rijnorm;
                vij=vij/_TIMESTEP;

                theta=randnorm48();
                fd=-1.0*_DPD_fc*dot(rij_unit,vij); 
                fr=1.0*theta*sigma/sqrt(_TIMESTEP); 
                fij=rij_unit*(fd+fr);
//if (ipt<1) INFO_Printf("i=%d j=%d fdx = %e frx = %e\n",ipt,jpt,fd.x,fr.x);
                theta=_F[ipt].norm();
                _F[ipt]+=fij;
                _F[jpt]-=fij;
//if (ipt<1) INFO_Printf("curstep=%d i=%d j=%d forig = %e df =%e fd=%e fr=%e \n",curstep,ipt,jpt,theta,fij.norm(),fd,fr);
                _VIRIAL.addnvv(1.,fij,rij);
	    }
	}
    }

    if (_ENABLE_FLAT_INDENTOR)
    {
        SHtoR();
        _FLAT_INDENTOR_F= 0; 
        _FLAT_INDENTOR_POS = _FLAT_INDENTOR_POS0 + _FLAT_INDENTOR_V * (curstep*_TIMESTEP);
        int ind;
        ind = _FLAT_INDENTOR_DIR - 1;
        for(i=0;i<_NP;i++)
        {
            if(_R[i][ind] > _FLAT_INDENTOR_POS) 
            {
               _FLAT_INDENTOR_F += -_FLAT_INDENTOR_K*(_R[i][ind] - _FLAT_INDENTOR_POS);
               _F[i][ind]       += -_FLAT_INDENTOR_K*(_R[i][ind] - _FLAT_INDENTOR_POS);
               _EPOT            += 0.5*_FLAT_INDENTOR_K*SQR(_R[i][ind] - _FLAT_INDENTOR_POS);
            }
        }
    }
}

void MDFrame::call_ANDERSON()
{
    int ipt;
    double sigma,theta;
    Vector3 vrand;
    Matrix33 hinv;

    if ( _ENABLE_ANDERSON )
    {
	sigma = sqrt((KB*_TDES)*SQR(_TIMESTEP)
		/(_ATOMMASS[0]*MASSCONVERT));
        hinv=_H.inv();
	for (ipt=0;ipt<_NP;ipt++)
	{
	    theta=1.0;
	    if (drand48()<_ANDERSON_RATIO)
            {
		if (species[ipt]!=0) 
		    theta=sqrt(_ATOMMASS[0]/_ATOMMASS[species[ipt]]);
		vrand.x=theta*sigma*randnorm48();
		vrand.y=theta*sigma*randnorm48();
		vrand.z=theta*sigma*randnorm48();
                _VSR[ipt]=hinv*vrand;
                _VR[ipt]=vrand;
            }
	}
    }
    
}

void MDFrame::potential()
{
    INFO("MDFrame::potential() not implemented yet.");
}

void MDFrame::potential_energyonly()
{    
    INFO("MDFrame::potential_energyonly() not implemented yet.");
    /* if not implemented then simply call potential */
}

double MDFrame::potential_energyonly(int iatom)
{
    INFO("MDFrame::potential_energyonly(int) not implemented yet.");
    /* if not implemented then simply call potential_energyonly()
     */
    return _EPOT;
}

void MDFrame::integrator()
{
    switch(algorithm_id) {
    case(10000): NVE_Gear6();                     break;
    case(20000): NVT_Gear6();                     break;
    case(21000): NVTC_Gear6();                    break;
    case(30000): NPH_Gear6();                     break;
    case(40000): NPT_Gear6();                     break;
    case(41000): NPTC_Gear6();                    break;
    case(10100): NVE_VVerlet(curstep);            break;
    case(20100): NVT_VVerlet_Implicit(curstep);   break;
    case(20101): NVT_VVerlet_Explicit_1(curstep); break;
    case(20102): NVT_VVerlet_Explicit_2(curstep); break;
    case(21100): NVTC_VVerlet_Explicit(curstep);  break; //temporary fix. Keonwook Kang
    case(21101): NVTC_VVerlet_Explicit(curstep);  break;
/*    case(30100): NPH_VVerlet(curstep);            break; */
    case(30100): NPH_VVerlet_Explicit_1(curstep); break; //temporary fix. Keonwook Kang Mar. 03 2008
    case(30101): NPH_VVerlet_Explicit_1(curstep); break;
 /* case(30102): NPH_VVerlet_Explicit_2(curstep); break; */
    case(40100): NPT_VVerlet_Implicit(curstep);   break;
    case(40101): NPT_VVerlet_Explicit_1(curstep); break;
    case(40102): NPT_VVerlet_Explicit_2(curstep); break;
    default: FATAL("unknown algorithm_id ("<<algorithm_id<<")"); break;
    }
}

/********************************************************************/
/* Numerical integrators now defined in integrators.cpp             */
/********************************************************************/
#include "integrators.in"

int MDFrame::runMDSWITCH()
{
    int n;
    n = _ENABLE_SWITCH;
    _ENABLE_SWITCH = 1;

    _WTOT = _WAVG = 0;

    run();

    _ENABLE_SWITCH = n;

    return 0;
}

//void MDFrame::runABF()
//{
//    int n;
//    n = _ENABLE_ABF;
//    _ENABLE_ABF = 1;
//
//    run();
//
//    _ENABLE_ABF = n;
//
//    return 0;
//}

void MDFrame::run()
{
    int itmp, i, key;
    double pres, omega0;
    class Matrix33 prp, h0inv, h0invtran;

    /* decide which algorithm to run */
    if (strcmp(ensemble_type,"NVE")==0)
        algorithm_id = 10000;
    else if (strcmp(ensemble_type,"NVT")==0)
        algorithm_id = 20000;
    else if (strcmp(ensemble_type,"NVTC")==0)
        algorithm_id = 21000;
    else if (strcmp(ensemble_type,"NPH")==0)
        algorithm_id = 30000;
    else if (strcmp(ensemble_type,"NPT")==0)
        algorithm_id = 40000;
    else if (strcmp(ensemble_type,"NPTC")==0)
        algorithm_id = 41000;
    else
    {
        algorithm_id = 100;
        ERROR("unknown ensemble_type("<<ensemble_type<<" use NVE\n"
            "choices: NVE, NVT, NPH, NPT");        
    }
    if (strcmp(integrator_type,"Gear6")==0)
        algorithm_id += 0;
    else if (strcmp(integrator_type,"VVerlet")==0)
        algorithm_id += 100;
    else
    {
        ERROR("unknown ensemble_type("<<ensemble_type<<" use Gear6\n"
              "choices: Gear6 VVerlet");
    }

    algorithm_id += implementation_type;
    INFO("algorithm_id = "<<algorithm_id);
    
    itmp = algorithm_id / 10000; /* 1: NVE, 2: NVT, 3: NPH, 4: NPT */
        
    SHtoR();
    /* prepare _SIGMA for Parrinello-Rahman */
    
    if((itmp==3)||(itmp==4)) /* NPH or NPT */
    {
        //call_potential(); /* do not need to call potential here */
        prp=_EXTSTRESS;
        prp*=_EXTSTRESSMUL;
//        pres=(-1)*prp.trace()/3;
        pres=prp.trace()/3;
        pres+=_EXTPRESSADD;
        prp[0][0]-=pres;
        prp[1][1]-=pres;
        prp[2][2]-=pres;

        if(fabs(_H0.det())<1e-10)
            FATAL("you need to specify H0 by calling saveH "
                  "before using the Parrinello-Rahman method");
        
        h0inv=_H0.inv();
        h0invtran=h0inv.tran();
        omega0=_H0.det()*(1-_VACUUMRATIO);
        _SIGMA=((h0inv*prp)*h0invtran)*omega0;
        _SIGMA/=160.2e3;/* convert MPa to eV/A^3 */
    }
    
    if(runavgposition) // 2007.06.04 Keonwook Kang
    {
        if(_SRA==NULL)
        {
            INFO("Running average position, _SRA allocated and initialized");
            Realloc(_SRA,Vector3,_NP);
        }
        else
            INFO("Running average position, _SRA reinitialized");
        memcpy(_SRA,_SR,sizeof(Vector3)*_NP);
    }
            
    if(continue_curstep)
       step0 = curstep;
    else
       step0 = 0;

    /* Surface tension calculation initialization */ 
    // specify surftensionspec(1st Pndir, 2nd Nnorm), SURFTEN, slotheight

    for(curstep=step0;curstep<(step0 + totalsteps);curstep++)
    //for(curstep=0;curstep<totalsteps;curstep++)
    {
        key=0;
        while(win!=NULL)
        {
            if(win->IsPaused()) sleep(1);
            else break;
        }

        if(curstep%savepropfreq==0) INFO("curstep="<<curstep);
        winplot();
        
        if(curstep>=equilsteps)
        {
            if(saveprop)
            {
                if((curstep%savepropfreq)==0)
                {
		    if ( (algorithm_id%1000)==0 || curstep==step0 || curstep==0 ) /* added by Seunghwa Ryu May/20/2009 */
		    {
                    	call_potential();/* added by Keonwook Kang, Aug 14, 2006 */
                    	calcprop();      /* added by Keonwook Kang, Aug 14, 2006 */
		    }
                    calcoutput();
                    pf.write(this);
                }
            }
           
        }
        /* Integrator */
        step();

        if(runavgposition) /* calculate running average position */
        {/* May 15 2007 Keonwook Kang */
            for(i=0;i<_NP;i++)
            {
                if(fixed[i]) continue;
                _SRA[i]=(_SRA[i]*(curstep-step0+1)+_SR[i])/(curstep-step0+2);
            }
        }

	/* save configuration for continuation run */        
	savecontinuecn(curstep);

        /* put savecn after step() to save latest config, KW Jan 19 2007 */
        if(curstep>=equilsteps)
        {
            /* save intermediate cfg (and cn files when writeall == 1) */
            if(   (savecfg &&  ((curstep%savecfgfreq)==0) )
               || (savecn  &&  ((curstep%savecnfreq )==0) && (writeall==1) ) )
            {
                if((curstep%plotfreq) != 0) // To avoid redundant calculation
                {
                    if(plot_color_axis==2)
                        calcentralsymmetry();
#ifdef _GSL
                    if(plot_color_axis==7)
                        calqlwl();
#endif
                }
            }
            //calcrystalorder(); /*only for test*/
            saveintercn(curstep);
            saveintercfg(curstep);
        }

        if(win!=NULL)
            if(win->IsAlive())
            if(autowritegiffreq>0)
                if((curstep%autowritegiffreq)==0)
                {
                    win->LockWritegif(); //lock and wait for writegif
                }
    }
/*************************/
    if(saveprop)
   {
        INFO("curstep="<<curstep); /*/ Keonwook Kang, May 03, 2007 */
        call_potential();/* added by Keonwook Kang, Aug 14, 2006 */
        calcprop();      /* added by Keonwook Kang, Aug 14, 2006 */
        calcoutput();    /* added by Keonwook Kang, Feb 05, 2007 */
        pf.write(this);
    }
}

void MDFrame::calcprop()  /*calculate physical properties*/
{
    int i, j, itmp;
    Matrix33 dH, s, hinv, hinvtran, htran, hsigma, sigmahtranh;
    Vector3 r0, r, dr, f, v;
    double tmp, pres;
#ifdef _TORSION_OR_BENDING        
    double tmp2, tmp4, curvature, Radius;
#endif
    
    tmp=0.5*(_ATOMMASS[0]*MASSCONVERT)/(_TIMESTEP*_TIMESTEP);
    _KATOM=0; _PSTRESS.clear();
    _EPOT_BOX=0; _KBOX=0;
    _OMEGA=_H.det();
    
#if 0 /* disabled, also see FLAT_INDENTOR */
    /* add indenter forces (spherical indentor)*/
    if(indenterspec[0]!=0)
    {
        /* use real coordinates to calculate indentor force, no PBC */
        x0=indenterspec[1];
        y0=indenterspec[2];
        z0=indenterspec[3];
        R0=indenterspec[4];
        K =indenterspec[5];
        r0.set(x0,y0,z0);
        for(i=0;i<_NP;i++)
        {
            r=_H*_SR[i];
            dr=r-r0;
            dR=dr.norm();
            if(dR<R0)
            {
                phi=0.5*K*(R0-dR)*(R0-dR);
                f=dr*(K*(R0-dR)/R0);
                _F[i]+=f;
                _EPOT_IND[i]+=phi;
                _EPOT+=phi;
            }
        }
    }
#endif

    /* every process repeat the same procedure */
    _NPfree = 0;
  if(_VR!=NULL)
  {
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        _NPfree ++;
        _VR[i]=_H*_VSR[i];
        if(species[i]==0)
        {
            _KATOM+=tmp*_VR[i].norm2();
            _PSTRESS.addnvv(2.0*tmp,_VR[i],_VR[i]);
        }
        else
        {
            if(_ATOMMASS[species[i]]<=0)
            {
                FATAL("_ATOMMASS (for species "<<species[i]<<") = "<<_ATOMMASS[species[i]]<<" is not valid");
            }
            _KATOM+=tmp*_VR[i].norm2()*_ATOMMASS[species[i]]/_ATOMMASS[0];
            _PSTRESS.addnvv(2.0*tmp*_ATOMMASS[species[i]]/_ATOMMASS[0],
                            _VR[i],_VR[i]);
        }
    }
  }
    /* _T=_KATOM*EV/(1.5*BOLZ*(_NPfree-1)); */
    _T=_KATOM/(1.5*KB*_NPfree);
    _TOTSTRESS=_PSTRESS+_VIRIAL;
    _TOTSTRESS/=_OMEGA*(1-_VACUUMRATIO); /* in eV/A^3 */
    _TOTSTRESSinMPa=_TOTSTRESS*160.2e3;  /* in MPa */
    _TOTPRESSURE = _TOTSTRESS.trace()/3;

    if (_ENABLE_separation_pot==1)
    {
	if (ST_step == 1) ST_dUdLAM_L = -ST_LMAX*_TOTSTRESS[1][1]*_H[0][0]*_H[2][2];
        else if (ST_step == 2) ST_dUdLAM_L = 0;
        else ST_dUdLAM_L = 0;
	ST_dUdLAM_TOT = ST_dUdLAM_L + ST_dUdLAM_POT;
    }
    
    /* compute linear and angular momenta */
    if(_VR!=NULL)
    {
      _P_COM.clear(); _L_COM.clear();
        for(i=0;i<_NP;i++)
        {
            if(fixed[i]) continue;
            r = _H*(_SR[i] - _COM);
            if(species[i]==0)
	    {
	        _P_COM += _VR[i]; _L_COM += cross(r,_VR[i]);
	    }
            else
	    {
	        v = _VR[i] * _ATOMMASS[species[i]] / _ATOMMASS[0];
                _P_COM += v; _L_COM += cross(_R[i],v);
	    }
        }
        _P_COM *= (_ATOMMASS[0]*MASSCONVERT/_TIMESTEP);
        _L_COM *= (_ATOMMASS[0]*MASSCONVERT/_TIMESTEP);
    }

#ifdef _TORSION_OR_BENDING    
    if (_TORSIONSIM)
    {
        /* Compute _PTHETA_COM */
        _PTHETA_COM = 0;
        tmp2 = tmp4 = 0;
        for(i=0;i<_NP0;i++)
        {
            r = _H*_SR[i];
            v = _H*_VSR[i];
            if(species[i]!=0)
                v *= _ATOMMASS[species[i]] / _ATOMMASS[0];
            tmp2 += (v.y*r.x - v.x*r.y);
            tmp4 += (r.x*r.x + r.y*r.y);
        }
        _PTHETA_COM = tmp2 / tmp4;
        
        /* Recompute _P_COM */
        _P_COM.clear(); _F_COM.clear();
        for(i=0;i<_NP0;i++)
        {
            if(species[i]==0)
                _P_COM += _VR[i];
            else
                _P_COM += _VR[i] * _ATOMMASS[species[i]] / _ATOMMASS[0];
            _F_COM += _F[i];
        }        
    }

    if (_BENDSIM)
    {
        /* Compute _PTHETA_COM */
        _PTHETA_COM = 0;
        tmp2 = tmp4 = 0;
        curvature = bendspec[3];
        Radius = 1.0/curvature;
        for(i=0;i<_NP0;i++)
        {
            r = _H*_SR[i];
            r.y += Radius;
            v = _H*_VSR[i];
            if(species[i]!=0)
                v *= _ATOMMASS[species[i]] / _ATOMMASS[0];
            tmp2 += (v.y*r.x - v.x*r.y);
            tmp4 += (r.x*r.x + r.y*r.y);
        }
        _PTHETA_COM = tmp2 / tmp4;

        /* Recompute _P_COM */
        _P_COM.clear(); _F_COM.clear();
        for(i=0;i<_NP0;i++)
        {
            if(species[i]==0)
                _P_COM += _VR[i];
            else
                _P_COM += _VR[i] * _ATOMMASS[species[i]] / _ATOMMASS[0];
            _F_COM += _F[i];            
        }        
    }
#endif
    
    itmp = algorithm_id / 10000; /* 1: NVE, 2: NVT, 3: NPH, 4: NPT */
    
    if((itmp==3)||(itmp==4)) 
    {
        /* PR in finite stress _EXTSTRESS */
        pres=_EXTSTRESS.trace()/3;
        pres*=_EXTSTRESSMUL;
        pres+=_EXTPRESSADD;
        pres/=160.2e3;/*convert from MPa to eV/A^3 */
        s=_TOTSTRESS;
        s*=_OMEGA*(1-_VACUUMRATIO);
        s[0][0]-=_OMEGA*(1-_VACUUMRATIO)*pres;
        s[1][1]-=_OMEGA*(1-_VACUUMRATIO)*pres;
        s[2][2]-=_OMEGA*(1-_VACUUMRATIO)*pres;
        
        hinv=_H.inv(); hinvtran=hinv.tran(); htran=_H.tran();
        s=s*hinvtran;
        hsigma=_H*_SIGMA;
        _GH=s-hsigma;
        //_GH*=-1;  /* bug fixed. Keonwook Kang Feb. 23 2008*/

        for(i=0;i<3;i++)
            for(j=0;j<3;j++)
                if(conj_fixboxvec[i][j]==1)
                    _GH[i][j]=0;        
        sigmahtranh=(_SIGMA*htran)*_H;
        
        /* potential energy of box */
        _EPOT_BOX=pres*_OMEGA*(1-_VACUUMRATIO)+sigmahtranh.trace()*0.5;

        /* kinetic energy of box */
        for(i=0;i<3;i++)
            for(j=0;j<3;j++)
                if(conj_fixboxvec[i][j]==0)
                    _KBOX+=tmp*_VH[i][j]*_VH[i][j]*_WALLMASS/_ATOMMASS[0];
    }

    /* only master keep track of Nose-Hoover zeta, HELM */

    itmp = algorithm_id / 1000; /* 1: NVE, 2: NVT, 3: NPH, 4: NPT */

    /* _HELMP is the conserved quantity */
    if(itmp==10) /* NVE */
    {
        _HELM  = _EPOT + _KATOM;
        _HELMP = _HELM;
    }
    else if(itmp==20) /* NVT */
    {
        if((NHMass[0]==0)&&(vt2!=0))
        {
            _HELM  = _EPOT + _KATOM +
                 3*_NPfree*_TDES*KB*zetav*zetav/_TIMESTEP/_TIMESTEP/1e-24/2/vt2;
        }
        else
        {
            _HELM  = _EPOT + _KATOM + 0.5*NHMass[0]*SQR(zetav/_TIMESTEP);
        }
        _HELMP = _HELM + 3*_NPfree*_TDES*KB*zeta;
    }
    else if(itmp==21) /* NVTC */
    {
        _HELM  = _EPOT + _KATOM;
        for(i=0;i<NHChainLen;i++)
            _HELM += 0.5*NHMass[i]*SQR(zetaNHCv[i]/_TIMESTEP);
        _HELMP = _HELM + 3*_NPfree*_TDES*KB*zetaNHC[0];
        for(i=1;i<NHChainLen;i++)
            _HELMP += KB*_TDES*zetaNHC[i];
    }
    else if(itmp==30) /* NPH */
    {
        _HELM  = _EPOT + _KATOM + _KBOX;
        _HELMP = _HELM + _EPOT_BOX;
    }
    else if(itmp==40) /* NPT */
    {
        if((NHMass[0]==0)&&(vt2!=0))
        {
            _HELM  = _EPOT + _KATOM + _KBOX +
                3*_NPfree*_TDES*KB*zetav*zetav/_TIMESTEP/_TIMESTEP/1e-24/2/vt2;            
        }
        else
        {
            _HELM  = _EPOT + _KATOM + _KBOX + 0.5*NHMass[0]*SQR(zetav/_TIMESTEP);
        }
        _HELMP = _HELM + _EPOT_BOX + 3*_NPfree*_TDES*KB*zeta;
    }
    else if(itmp==41)  /* NPTC */
    {
        _HELM  = _EPOT + _KATOM + _KBOX;
        for(i=0;i<NHChainLen;i++)
            _HELM += 0.5*NHMass[i]*SQR(zetaNHCv[i]/_TIMESTEP);
        _HELMP = _HELM + _EPOT_BOX + 3*_NPfree*_TDES*KB*zetaNHC[0];
        for(i=1;i<NHChainLen;i++)
            _HELMP += KB*_TDES*zetaNHC[i];
    }

    /* apply constant shear rate on box H */
    if(applyshear)
    {
        dH=_H*_SHEARRATE; dH*=_TIMESTEP;
        _H+=dH;
    }

    if(plot_color_axis==2)
    {
        if(((curstep%plotfreq)==0)&&(conj_step%plotfreq)==0)
            calcentralsymmetry();
    }

    if(plot_color_axis==4)
    {
        if(((curstep%plotfreq)==0)&&(conj_step%plotfreq)==0)
            caldisregistry();
    }

#ifdef _GSL
    if(plot_color_axis==7)
    {
        if (strcmp(integrator_type,"VVerlet")==0)
        {
            if((( (curstep+1)%plotfreq)==0)&&((conj_step)%plotfreq)==0)
            calqlwl();
            if (curstep == step0 || curstep == 0 ) calqlwl();
        }
        else
        {
            if((( (curstep)%plotfreq)==0)&&((conj_step)%plotfreq)==0)
            calqlwl();
        }
    }
#endif  

    if(plot_color_axis==8)
    {
	}
	else
	{
    }
 
    if(plot_color_axis==9)
    {
	}
	else
	{
    }
    
}

void MDFrame::calcoutput()
{
    char names[10000], *p, *q;
    /* separate output_fmt into short strings between white spaces */

    if(strlen(output_fmt)==0)
    {
        output_str[0]=0;
        return;
    }
    
    strcpy(names,output_fmt);
    p=names;

    output_str[0]=0;
    while(strlen(p)>0)
    {
        switch(*p)
        {
        case ' ':
        case '\t': /* white spaces */
            p++; break;
        default:
            q=strchr(p,' ');
            if(q!=NULL) *q=0;
            identify(p);
            if(curn>=0)
                switch(vartype[curn])                                    
                {                                                            
                case(INT): sprintf(output_str+strlen(output_str),"%10d ",*((int *)varptr[curn]+shift)); break; 
                case(LONG): sprintf(output_str+strlen(output_str),"%10ld ",*((long *)varptr[curn]+shift)); break;
                case(DOUBLE): sprintf(output_str+strlen(output_str),"%23.16e ",*((double *)varptr[curn]+shift)); break;
                case(STRING): sprintf(output_str+strlen(output_str),"%s ",(char *)varptr[curn]+shift); break;
                default: FATAL("unknown vartype ("<<vartype[curn]<<")");
                }
            if(q!=NULL) p=q+1;
            else goto it;
            break;
        }
    }
 it:
    sprintf(output_str+strlen(output_str),"%s","\n");
}

void MDFrame::calDVIRIALDexx()
{/* Specify a positive value for Dexx (e.g. Dexx = 1e-4)       
    before starting MD simulation will enable the calculation of
    d(_VIRIAL)/d(_exx) during MD simulation by numerical differentiation
    (centered difference)
  */
    Matrix33 virial1, virial2, htmp;
    
    htmp = _H;
    
    _H[0][0]+=_H[0][0]*_Dexx;
    _H[0][1]+=_H[0][1]*_Dexx;
    _H[0][2]+=_H[0][2]*_Dexx;

    potential();

    virial1 = _VIRIAL;

    _H = htmp;
    
    _H[0][0]-=_H[0][0]*_Dexx;
    _H[0][1]-=_H[0][1]*_Dexx;
    _H[0][2]-=_H[0][2]*_Dexx;

    potential();

    virial2 = _VIRIAL;

    _DVIRIAL_Dexx = virial1 - virial2;
    _DVIRIAL_Dexx /= (2*_Dexx);

    _DVIRIAL[0][0][0][0] = _DVIRIAL_Dexx[0][0]; // C^B_1111
    _DVIRIAL[0][0][1][1] = _DVIRIAL_Dexx[1][1]; // C^B_1122
    _DVIRIAL[0][0][2][2] = _DVIRIAL_Dexx[2][2]; // C^B_1133
    _DVIRIAL[0][0][0][1] = _DVIRIAL_Dexx[0][1]; // C^B_1112
    _DVIRIAL[0][0][1][0] = _DVIRIAL_Dexx[1][0]; // C^B_1121
    _DVIRIAL[0][0][1][2] = _DVIRIAL_Dexx[1][2]; // C^B_1123
    _DVIRIAL[0][0][2][1] = _DVIRIAL_Dexx[2][1]; // C^B_1132
    _DVIRIAL[0][0][2][0] = _DVIRIAL_Dexx[2][0]; // C^B_1131
    _DVIRIAL[0][0][0][2] = _DVIRIAL_Dexx[0][2]; // C^B_1113

    _H = htmp;

    _H[1][0]+=_H[1][0]*_Dexx;
    _H[1][1]+=_H[1][1]*_Dexx;
    _H[1][2]+=_H[1][2]*_Dexx;

    potential();

    virial1 = _VIRIAL;

    _H = htmp;
    
    _H[1][0]-=_H[1][0]*_Dexx;
    _H[1][1]-=_H[1][1]*_Dexx;
    _H[1][2]-=_H[1][2]*_Dexx;

    potential();

    virial2 = _VIRIAL;

    _DVIRIAL_Dexx = virial1 - virial2;
    _DVIRIAL_Dexx /= (2*_Dexx);

    _DVIRIAL[1][1][0][0] = _DVIRIAL_Dexx[0][0]; // C^B_2211
    _DVIRIAL[1][1][1][1] = _DVIRIAL_Dexx[1][1]; // C^B_2222
    _DVIRIAL[1][1][2][2] = _DVIRIAL_Dexx[2][2]; // C^B_2233
    _DVIRIAL[1][1][0][1] = _DVIRIAL_Dexx[0][1]; // C^B_2212
    _DVIRIAL[1][1][1][0] = _DVIRIAL_Dexx[1][0]; // C^B_2221
    _DVIRIAL[1][1][1][2] = _DVIRIAL_Dexx[1][2]; // C^B_2223
    _DVIRIAL[1][1][2][1] = _DVIRIAL_Dexx[2][1]; // C^B_2232
    _DVIRIAL[1][1][2][0] = _DVIRIAL_Dexx[2][0]; // C^B_2231
    _DVIRIAL[1][1][0][2] = _DVIRIAL_Dexx[0][2]; // C^B_2213

    _H = htmp;

    _H[2][0]+=_H[2][0]*_Dexx;
    _H[2][1]+=_H[2][1]*_Dexx;
    _H[2][2]+=_H[2][2]*_Dexx;

    potential();

    virial1 = _VIRIAL;

    _H = htmp;
    
    _H[2][0]-=_H[2][0]*_Dexx;
    _H[2][1]-=_H[2][1]*_Dexx;
    _H[2][2]-=_H[2][2]*_Dexx;

    potential();

    virial2 = _VIRIAL;

    _DVIRIAL_Dexx = virial1 - virial2;
    _DVIRIAL_Dexx /= (2*_Dexx);

    _DVIRIAL[2][2][0][0] = _DVIRIAL_Dexx[0][0]; // C^B_3311
    _DVIRIAL[2][2][1][1] = _DVIRIAL_Dexx[1][1]; // C^B_3322
    _DVIRIAL[2][2][2][2] = _DVIRIAL_Dexx[2][2]; // C^B_3333
    _DVIRIAL[2][2][0][1] = _DVIRIAL_Dexx[0][1]; // C^B_3312
    _DVIRIAL[2][2][1][0] = _DVIRIAL_Dexx[1][0]; // C^B_3321
    _DVIRIAL[2][2][1][2] = _DVIRIAL_Dexx[1][2]; // C^B_3323
    _DVIRIAL[2][2][2][1] = _DVIRIAL_Dexx[2][1]; // C^B_3332
    _DVIRIAL[2][2][2][0] = _DVIRIAL_Dexx[2][0]; // C^B_3331
    _DVIRIAL[2][2][0][2] = _DVIRIAL_Dexx[0][2]; // C^B_3313

    _H = htmp;

    _H[1][0]+=_H[2][0]*_Dexx;
    _H[1][1]+=_H[2][1]*_Dexx;
    _H[1][2]+=_H[2][2]*_Dexx;

    potential();

    virial1 = _VIRIAL;

    _H = htmp;
    
    _H[1][0]-=_H[2][0]*_Dexx;
    _H[1][1]-=_H[2][1]*_Dexx;
    _H[1][2]-=_H[2][2]*_Dexx;

    potential();

    virial2 = _VIRIAL;

    _DVIRIAL_Dexx = virial1 - virial2;
    _DVIRIAL_Dexx /= (2*_Dexx);

    _DVIRIAL[1][2][0][0] = _DVIRIAL_Dexx[0][0]; // C^B_2311
    _DVIRIAL[1][2][1][1] = _DVIRIAL_Dexx[1][1]; // C^B_2322
    _DVIRIAL[1][2][2][2] = _DVIRIAL_Dexx[2][2]; // C^B_2333
    _DVIRIAL[1][2][0][1] = _DVIRIAL_Dexx[0][1]; // C^B_2312
    _DVIRIAL[1][2][1][0] = _DVIRIAL_Dexx[1][0]; // C^B_2321
    _DVIRIAL[1][2][1][2] = _DVIRIAL_Dexx[1][2]; // C^B_2323
    _DVIRIAL[1][2][2][1] = _DVIRIAL_Dexx[2][1]; // C^B_2332
    _DVIRIAL[1][2][2][0] = _DVIRIAL_Dexx[2][0]; // C^B_2331
    _DVIRIAL[1][2][0][2] = _DVIRIAL_Dexx[0][2]; // C^B_2313

    _H = htmp;

    _H[2][0]+=_H[1][0]*_Dexx;
    _H[2][1]+=_H[1][1]*_Dexx;
    _H[2][2]+=_H[1][2]*_Dexx;

    potential();

    virial1 = _VIRIAL;

    _H = htmp;
    
    _H[2][0]-=_H[1][0]*_Dexx;
    _H[2][1]-=_H[1][1]*_Dexx;
    _H[2][2]-=_H[1][2]*_Dexx;

    potential();

    virial2 = _VIRIAL;

    _DVIRIAL_Dexx = virial1 - virial2;
    _DVIRIAL_Dexx /= (2*_Dexx);

    _DVIRIAL[2][1][0][0] = _DVIRIAL_Dexx[0][0]; // C^B_3211
    _DVIRIAL[2][1][1][1] = _DVIRIAL_Dexx[1][1]; // C^B_3222
    _DVIRIAL[2][1][2][2] = _DVIRIAL_Dexx[2][2]; // C^B_3233
    _DVIRIAL[2][1][0][1] = _DVIRIAL_Dexx[0][1]; // C^B_3212
    _DVIRIAL[2][1][1][0] = _DVIRIAL_Dexx[1][0]; // C^B_3221
    _DVIRIAL[2][1][1][2] = _DVIRIAL_Dexx[1][2]; // C^B_3223
    _DVIRIAL[2][1][2][1] = _DVIRIAL_Dexx[2][1]; // C^B_3232
    _DVIRIAL[2][1][2][0] = _DVIRIAL_Dexx[2][0]; // C^B_3231
    _DVIRIAL[2][1][0][2] = _DVIRIAL_Dexx[0][2]; // C^B_3213

    _H = htmp;

    _H[0][0]+=_H[2][0]*_Dexx;
    _H[0][1]+=_H[2][1]*_Dexx;
    _H[0][2]+=_H[2][2]*_Dexx;

    potential();

    virial1 = _VIRIAL;

    _H = htmp;
    
    _H[0][0]-=_H[2][0]*_Dexx;
    _H[0][1]-=_H[2][1]*_Dexx;
    _H[0][2]-=_H[2][2]*_Dexx;

    potential();

    virial2 = _VIRIAL;

    _DVIRIAL_Dexx = virial1 - virial2;
    _DVIRIAL_Dexx /= (2*_Dexx);

    _DVIRIAL[0][2][0][0] = _DVIRIAL_Dexx[0][0]; // C^B_1311
    _DVIRIAL[0][2][1][1] = _DVIRIAL_Dexx[1][1]; // C^B_1322
    _DVIRIAL[0][2][2][2] = _DVIRIAL_Dexx[2][2]; // C^B_1333
    _DVIRIAL[0][2][0][1] = _DVIRIAL_Dexx[0][1]; // C^B_1312
    _DVIRIAL[0][2][1][0] = _DVIRIAL_Dexx[1][0]; // C^B_1321
    _DVIRIAL[0][2][1][2] = _DVIRIAL_Dexx[1][2]; // C^B_1323
    _DVIRIAL[0][2][2][1] = _DVIRIAL_Dexx[2][1]; // C^B_1332
    _DVIRIAL[0][2][2][0] = _DVIRIAL_Dexx[2][0]; // C^B_1331
    _DVIRIAL[0][2][0][2] = _DVIRIAL_Dexx[0][2]; // C^B_1313

    _H = htmp;

    _H[2][0]+=_H[0][0]*_Dexx;
    _H[2][1]+=_H[0][1]*_Dexx;
    _H[2][2]+=_H[0][2]*_Dexx;

    potential();

    virial1 = _VIRIAL;

    _H = htmp;
    
    _H[2][0]-=_H[0][0]*_Dexx;
    _H[2][1]-=_H[0][1]*_Dexx;
    _H[2][2]-=_H[0][2]*_Dexx;

    potential();

    virial2 = _VIRIAL;

    _DVIRIAL_Dexx = virial1 - virial2;
    _DVIRIAL_Dexx /= (2*_Dexx);

    _DVIRIAL[2][0][0][0] = _DVIRIAL_Dexx[0][0]; // C^B_3111
    _DVIRIAL[2][0][1][1] = _DVIRIAL_Dexx[1][1]; // C^B_3122
    _DVIRIAL[2][0][2][2] = _DVIRIAL_Dexx[2][2]; // C^B_3133
    _DVIRIAL[2][0][0][1] = _DVIRIAL_Dexx[0][1]; // C^B_3112
    _DVIRIAL[2][0][1][0] = _DVIRIAL_Dexx[1][0]; // C^B_3121
    _DVIRIAL[2][0][1][2] = _DVIRIAL_Dexx[1][2]; // C^B_3123
    _DVIRIAL[2][0][2][1] = _DVIRIAL_Dexx[2][1]; // C^B_3132
    _DVIRIAL[2][0][2][0] = _DVIRIAL_Dexx[2][0]; // C^B_3131
    _DVIRIAL[2][0][0][2] = _DVIRIAL_Dexx[0][2]; // C^B_3113

    _H = htmp;

    _H[0][0]+=_H[1][0]*_Dexx;
    _H[0][1]+=_H[1][1]*_Dexx;
    _H[0][2]+=_H[1][2]*_Dexx;

    potential();

    virial1 = _VIRIAL;

    _H = htmp;
    
    _H[0][0]-=_H[1][0]*_Dexx;
    _H[0][1]-=_H[1][1]*_Dexx;
    _H[0][2]-=_H[1][2]*_Dexx;

    potential();

    virial2 = _VIRIAL;

    _DVIRIAL_Dexx = virial1 - virial2;
    _DVIRIAL_Dexx /= (2*_Dexx);

    _DVIRIAL[0][1][0][0] = _DVIRIAL_Dexx[0][0]; // C^B_1211
    _DVIRIAL[0][1][1][1] = _DVIRIAL_Dexx[1][1]; // C^B_1222
    _DVIRIAL[0][1][2][2] = _DVIRIAL_Dexx[2][2]; // C^B_1233
    _DVIRIAL[0][1][0][1] = _DVIRIAL_Dexx[0][1]; // C^B_1212
    _DVIRIAL[0][1][1][0] = _DVIRIAL_Dexx[1][0]; // C^B_1221
    _DVIRIAL[0][1][1][2] = _DVIRIAL_Dexx[1][2]; // C^B_1223
    _DVIRIAL[0][1][2][1] = _DVIRIAL_Dexx[2][1]; // C^B_1232
    _DVIRIAL[0][1][2][0] = _DVIRIAL_Dexx[2][0]; // C^B_1231
    _DVIRIAL[0][1][0][2] = _DVIRIAL_Dexx[0][2]; // C^B_1213

    _H = htmp;

    _H[1][0]+=_H[0][0]*_Dexx;
    _H[1][1]+=_H[0][1]*_Dexx;
    _H[1][2]+=_H[0][2]*_Dexx;

    potential();

    virial1 = _VIRIAL;

    _H = htmp;
    
    _H[1][0]-=_H[0][0]*_Dexx;
    _H[1][1]-=_H[0][1]*_Dexx;
    _H[1][2]-=_H[0][2]*_Dexx;

    potential();

    virial2 = _VIRIAL;

    _DVIRIAL_Dexx = virial1 - virial2;
    _DVIRIAL_Dexx /= (2*_Dexx);

    _DVIRIAL[1][0][0][0] = _DVIRIAL_Dexx[0][0]; // C^B_2111
    _DVIRIAL[1][0][1][1] = _DVIRIAL_Dexx[1][1]; // C^B_2122
    _DVIRIAL[1][0][2][2] = _DVIRIAL_Dexx[2][2]; // C^B_2133
    _DVIRIAL[1][0][0][1] = _DVIRIAL_Dexx[0][1]; // C^B_2112
    _DVIRIAL[1][0][1][0] = _DVIRIAL_Dexx[1][0]; // C^B_2121
    _DVIRIAL[1][0][1][2] = _DVIRIAL_Dexx[1][2]; // C^B_2123
    _DVIRIAL[1][0][2][1] = _DVIRIAL_Dexx[2][1]; // C^B_2132
    _DVIRIAL[1][0][2][0] = _DVIRIAL_Dexx[2][0]; // C^B_2131
    _DVIRIAL[1][0][0][2] = _DVIRIAL_Dexx[0][2]; // C^B_2113

    _H = htmp;

}

void MDFrame::calcentralsymmetry()
{
    int i, j, k, kpt, *csflags, mk1, mk2;
    double mind1, mind2, d1, d2; Vector3 sik1, sik, rik;
    
    csflags = (int *)malloc(sizeof(int)*_NNM);
    
    for(i=0;i<_NP;i++)
    {
        _TOPOL[i]=0;
    }
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==-1) continue;
        for(j=0;j<_NNM;j++)
            csflags[j]=0;

        /* handling the exceptional case of not enough nearst neighbors */
        if(nn[i]<NCS) {
            _TOPOL[i]=NCS+1.0;
            continue;
        }
        
        for(j=0;j<NCS/2;j++)
        {
            mind1=1e10; mk1=-1; sik1.set(-1,-1,-1);
            /* finding the nearest neighor */
            for(k=0;k<nn[i];k++)
            {
                if(csflags[k]!=0) continue;
                kpt=nindex[i][k];
                if(kpt==i) continue;
                if(fixed[kpt]==-1) continue;
                sik=_SR[kpt]-_SR[i]; sik.subint();
                rik=_H*sik; d1=rik.norm2();
                if(mind1>d1)
                {
                    mind1=d1; mk1=k; sik1=sik;
                }
            }
            csflags[mk1]=1;
            mind2=1e10; mk2=-1; sik.set(-1,-1,-1);
            /* finding the opposite atom */
            for(k=0;k<nn[i];k++)
            {
                if(csflags[k]!=0) continue;
                kpt=nindex[i][k];
                if(kpt==i) continue;
                if(fixed[kpt]==-1) continue;
                sik=_SR[kpt]-_SR[i]; sik+=sik1; sik.subint();
                rik=_H*sik; d2=rik.norm2();
                if(mind2>d2)
                {
                    mind2=d2; mk2=k;
                }
            }
            csflags[mk2]=1;
            _TOPOL[i]+=mind2;
            
        }
    }

    free(csflags);
}

void MDFrame::caldisregistry()
{
    int i, j, jpt;
    Vector3 r, s;
    double d;

    for(i=0;i<_NP;i++)
        _TOPOL[i]=0;

    for(i=0;i<_NP;i++)
    {
        d=0;
        for(j=0;j<nn_plot[i];j++)
        {
            jpt=nindex_plot[i][j];
            s=_SR[i]-_SR[jpt];
            s.subint();
            r=_H*s;
            d+= (r.norm() - latticeconst[0]*sqrt(3)/4.0);
        }
        _TOPOL[i]=d;
    }
}

void MDFrame::findcore()
{ /* find the center of mass of defect atoms
   *  assuming that there is only one defect
   *  and that the defect can be characterized by CSD parameter
   */
    int i, j, n, calshape;
    double tmin, tmax;
    Vector3 s0, s1, ds, savg, ravg, ri0;
    Matrix33 Inertia, Unit, Temp;
    double m,q,p;
    double rri0;
    double phi, I1, I2, I3;

    tmin = input[0];
    tmax = input[1];
    calshape = (int) input[2];

    // Matrix Initialize
    Inertia.clear();
    Unit.clear();
    Unit[0][0]=1.0;Unit[1][1]=1.0;Unit[2][2]=1.0;

    n = 0;  savg.clear();
    for(i=0;i<_NP;i++)
    {
      if(plot_limits[0]==1)
        if((_SR[i].x<plot_limits[1])||(_SR[i].x>plot_limits[2])
         ||(_SR[i].y<plot_limits[3])||(_SR[i].y>plot_limits[4])
         ||(_SR[i].z<plot_limits[5])||(_SR[i].z>plot_limits[6]))
        {
            continue;
        }
        if( (_TOPOL[i]>=tmin) && (_TOPOL[i]<tmax) )
        {
            if(n==0)
            {
                s0 = _SR[i];
                savg = s0;
            }
            else
            {
                s1 = _SR[i];
                ds = s1 - s0;  ds.subint();
                s1 = s0 + ds;
                savg+=s1;

 		// Inertia Tensor
                if (calshape == 1)
                {
                    ri0 = _H*ds; rri0=ri0.norm();
                    Inertia[0][0]+=rri0*rri0-ri0[0]*ri0[0];
                    Inertia[0][1]+=         -ri0[0]*ri0[1];
                    Inertia[0][2]+=         -ri0[0]*ri0[2];
                    Inertia[1][0]+=         -ri0[1]*ri0[0];
                    Inertia[1][1]+=rri0*rri0-ri0[1]*ri0[1];
                    Inertia[1][2]+=         -ri0[1]*ri0[2];
                    Inertia[2][0]+=         -ri0[2]*ri0[0];
                    Inertia[2][1]+=	        -ri0[2]*ri0[1];
                    Inertia[2][2]+=rri0*rri0-ri0[2]*ri0[2];
                }

            }
            n++;
        }
    }
    if (n>0)
       savg /= n;
    else
       ERROR("no defect atoms found!");

    ravg = _H*savg;
    INFO_Printf("findcore: avg pos =  %20.10e %20.10e %20.10e  %20.10e %20.10e %20.10e\n",
                savg.x,savg.y,savg.z, ravg.x,ravg.y,ravg.z);

    if (calshape == 1)
    { /* compute the principle axes of the nucleus */
        m = Inertia.trace()/3;
        Temp = Inertia - (Unit*m);
        q = Temp.det()/2;   
        p = 0;
        for (i=0;i<3;i++)
            for (j=0;j<3;j++)
            {
                p += Temp[i][j]*Temp[i][j];
            }
        p = p/6; 
        phi = atan( sqrt(p*p*p-q*q)/q )/3;
        if (phi<0) phi += M_PI; 
        I1 = m + 2*sqrt(p)*cos(phi);
        I2 = m - sqrt(p)*( cos(phi) + sqrt(3) * sin(phi));
        I3 = m - sqrt(p)*( cos(phi) - sqrt(3) * sin(phi));
        Principal_Inertia[0]=I1;
        Principal_Inertia[1]=I2;
        Principal_Inertia[2]=I3;
        INFO_Printf("Inertia Eigenvalues = %15.7e %15.7e %15.7e\n",I1,I2,I3);
    }
}

void MDFrame::calHessian()
{
    FILE *fp;
    int i, ind, j, ipt;
    Matrix33 hinv;
    Vector3 df, ds, dr;

    /* prints out df/dx which is the negative of the Hessian matrix */
    fp=fopen("hessian.out","w");

    INFO("open hessian.out");
    SHtoR();
    call_potential();
    memcpy(_VR,_F,_NP*sizeof(Vector3));
    memcpy(_VSR,_SR,_NP*sizeof(Vector3));
    hinv=_H.inv();
    if(input[0]==0)
    {
        for(ipt=0;ipt<_NP;ipt++)
        {
            INFO("atom "<<ipt);
            for(ind=0;ind<3;ind++)
            {
                _R[ipt]=_H*_VSR[ipt];
                _R[ipt][ind]-=_TIMESTEP;
                _SR[ipt]=hinv*_R[ipt];
                call_potential();
                memcpy(_VR,_F,_NP*sizeof(Vector3));
                
                _R[ipt]=_H*_VSR[ipt];
                _R[ipt][ind]+=_TIMESTEP;
                _SR[ipt]=hinv*_R[ipt];
                call_potential();

                _R[ipt]=_H*_VSR[ipt];
                _SR[ipt]=_VSR[ipt];

                for(j=0;j<_NP;j++)
                {
                    df=_F[j]-_VR[j]; df/=(2*_TIMESTEP);
                    fprintf(fp,"%20.13e %20.13e %20.13e\n",df.x,df.y,df.z);
                }
            }
        }
    }
    else
    {
        for(i=0;i<(int)input[0];i++)
        {
            ipt=(int)input[i+1];
            //INFO("atom "<<ipt<<" has "<<nn[ipt]<<" neighbors:");
            //for(j=0;j<nn[ipt];j++)
            //    INFO_Printf(" %d ",nindex[ipt][j]);
            //INFO_Printf("\n");
            for(ind=0;ind<3;ind++)
            {
                /* use centered difference */
                _R[ipt]=_H*_VSR[ipt];
                _R[ipt][ind]-=_TIMESTEP;
                _SR[ipt]=hinv*_R[ipt];
                call_potential();
                memcpy(_VR,_F,_NP*sizeof(Vector3));

                //if(ipt==4) INFO(_F[ipt]<<"  "<<_F[167]);
                
                _R[ipt]=_H*_VSR[ipt];
                _R[ipt][ind]+=_TIMESTEP;
                _SR[ipt]=hinv*_R[ipt];
                call_potential();

                //if(ipt==4) INFO(_F[ipt]<<"  "<<_F[167]);
                
                _R[ipt]=_H*_VSR[ipt];
                _SR[ipt]=_VSR[ipt];
                
                for(j=0;j<_NP;j++)
                {
                    df=_F[j]-_VR[j]; df/=(2*_TIMESTEP);
                    ds=_SR[j]-_SR[ipt]; ds.subint(); dr=_H*ds;
                    if(df.norm()>1e-10)
                        fprintf(fp,"%6d %3d %6d %20.13e %20.13e %20.13e   %20.13e %20.13e %20.13e\n",
                                ipt,ind+1,j,df.x,df.y,df.z,dr.x,dr.y,dr.z);
                }
            }
        }        
    }
    fclose(fp);
}

void MDFrame::readHessian()
{
    FILE *fp;
    int i, id, j;
    Vector3 df;

    if(_NP<=0) FATAL("readHessian: NP="<<_NP);
    
    /* prints out df/dx which is the negative of the Hessian matrix */
    fp=fopen(incnfile,"r");

    Realloc(hes,double,(_NP*3)*(_NNM*3));
    Realloc(hesind,int,(_NP*3)*(_NNM*3)*4);

    heslen = 0;
    for(i=0;i<_NP;i++)
    {
        for(id=0;id<3;id++)
        {
            for(j=0;j<_NP;j++)
            {
                fscanf(fp,"%le %le %le",&df.x,&df.y,&df.z);
                if(df.norm2()>1e-16)
                {
                    hes[heslen] = -df.x;
                    hesind[heslen*4  ] = i;
                    hesind[heslen*4+1] = id;
                    hesind[heslen*4+2] = j;
                    hesind[heslen*4+3] = 0;                    
                    heslen ++;

                    hes[heslen] = -df.y;
                    hesind[heslen*4  ] = i;
                    hesind[heslen*4+1] = id;
                    hesind[heslen*4+2] = j;
                    hesind[heslen*4+3] = 1;                    
                    heslen ++;

                    hes[heslen] = -df.z;
                    hesind[heslen*4  ] = i;
                    hesind[heslen*4+1] = id;
                    hesind[heslen*4+2] = j;
                    hesind[heslen*4+3] = 2;                    
                    heslen ++;

                    if(heslen>(_NP*3)*(_NNM*3))
                    {
                        INFO_Printf("heslen = %d _NP = %d _NNM = %d\n",heslen,_NP,_NNM);
                        FATAL("heslen too large, increase NNM");
                    }
                }
            }
        }
    }
    INFO_Printf("readHessian: heslen = %d\n",heslen);
    fclose(fp);
}

void MDFrame::calModeHessian()
{
    //FILE *fp;
    int i, ind, j, ipt, k, p;
    Matrix33 hinv;
    Vector3 df, ds, dr;

    int *mode_atomID, natoms, groupID, n, m, nmodes, pmode;
    Vector3 *mode_value;
    double s, sum2, tmp, mag;
    
    /* atoms with group == groupID will be included for mode */
    groupID = (int) input[0];
    n = (int) input[1];
    m = (int) input[2];
    pmode = (int) input[3];
    mag = input[4];
    
    nmodes = n*m*3;

    natoms = 0;
    for(i=0;i<_NP;i++)
        if(group[i]==groupID)
            natoms++;
    mode_atomID = (int *) malloc(sizeof(int)*natoms);
    mode_value = (Vector3 *) malloc(sizeof(Vector3)*natoms*nmodes);

    k=0;
    for(i=0;i<_NP;i++)
        if(group[i]==groupID)
        {
            mode_atomID[k]=i;
            k++;
        }    

    /* initialize modes */
    SHtoR();
    for(i=0;i<n;i++)
        for(j=0;j<m;j++)
        {
            for(p=0;p<natoms;p++)
            {
                ipt=mode_atomID[p];
                s = pow(_R[ipt].x,(double)i)*pow(_R[ipt].y,(double)j);

                for(k=0;k<3;k++)
                {
                    ind = i*(m*3)+j*3+k;
                    mode_value[ind*natoms+p][k] = s;
                }
            }
        }

    /* orthonormalize modes */
    for(i=0;i<nmodes;i++)
    {
        for(j=0;j<i;j++)
        { /* orthongalize mode i wrt mode j */
            sum2 = 0;
            for(p=0;p<natoms;p++)
                sum2 += dot(mode_value[i*natoms+p],mode_value[j*natoms+p]);
            for(p=0;p<natoms;p++)
                mode_value[i*natoms+p] -= mode_value[j*natoms+p]*sum2;            
        }

        /* normalize mode i */
        sum2 = 0;
        for(p=0;p<natoms;p++)
            sum2 += mode_value[i*natoms+p].norm2();
        tmp = sqrt(sum2);
        for(p=0;p<natoms;p++)
            mode_value[i*natoms+p]/=tmp; 
    }

    if(mag>0)
    { /* visualize mode */
        for(p=0;p<natoms;p++)
            _R[mode_atomID[p]] += mode_value[pmode*natoms+p]*mag;
        RHtoS();
        winplot();        
    }
    
    free(mode_atomID);
    free(mode_value);    
}

void MDFrame::calphonondisp()
{
    INFO("MDFrame::calphonondisp not implemented");
}

void MDFrame::calmisfit()
{
    int nx, ny, nz, i, j, k, n, ind, inda, indb;
    double x0, y0, z0, x1, y1, z1, dx, dy, dz, area;
    FILE *fp;
    class Vector3 c1, c2, c3, ca, cb, caxcb, dc1, dc2, dc3; /* *S0 */
    class Matrix33 H_origin; //, StreSS;

    /* direction of the surface normal */
    ind = (int)input[0];
    if(ind < 1 || ind > 3)
    {
        ERROR("Surface normal should be one of 1, 2, and 3.");
        return;
    }
    ind--;

    /* the other two directions orthogonal to ind */
    inda = (ind+1)%3; indb = (inda+1)%3;
    
    x0 = input[1]; dx = input[2]; x1 = input[3];
    y0 = input[4]; dy = input[5]; y1 = input[6];
    z0 = input[7]; dz = input[8]; z1 = input[9];

    INFO_Printf("%%        x0 = %e, dx = %e, x1 = %e\n",x0,dx,x1);
    INFO_Printf("%%        y0 = %e, dy = %e, y1 = %e\n",y0,dy,y1);
    INFO_Printf("%%        z0 = %e, dz = %e, z1 = %e\n",z0,dz,z1);
    
    if(dx!=0) nx = (int) round((x1-x0)/dx); else nx = 0;
    if(dy!=0) ny = (int) round((y1-y0)/dy); else ny = 0;
    if(dz!=0) nz = (int) round((z1-z0)/dz); else nz = 0;

    INFO_Printf("%%        nx = %d, ny = %d, nz = %d\n",nx,ny,nz);

    fp = fopen("Emisfit.out","w");
    if(fp==NULL) FATAL("calmisfit: open file failure");
        
    fprintf(fp,"%% Data file for misfit energy E(i,j,k) (eV/A^2)\n");
    fprintf(fp,"%% Format: i, j, k, E(i,j,k)\n");
    fprintf(fp,"%%         i=1:nx, j=1:ny, k=1:nz\n");
    fprintf(fp,"%%         nx = %d, ny = %d, nz = %d\n",nx,ny,nz);
    fprintf(fp,"\n");    
    fprintf(fp,"%%         x0 = %e, dx = %e, x1 = %e\n",x0,dx,x1);
    fprintf(fp,"%%         y0 = %e, dy = %e, y1 = %e\n",y0,dy,y1);
    fprintf(fp,"%%         z0 = %e, dz = %e, z1 = %e\n",z0,dz,z1);
    fprintf(fp,"\n");

    c1.set(_H[0][0],_H[1][0],_H[2][0]);
    c2.set(_H[0][1],_H[1][1],_H[2][1]);
    c3.set(_H[0][2],_H[1][2],_H[2][2]);
    c1 /= latticesize[0][3]; // one lattice vector in each direction
    c2 /= latticesize[1][3];
    c3 /= latticesize[2][3];
    
    ca.set(_H[0][inda],_H[1][inda],_H[2][inda]); 
    cb.set(_H[0][indb],_H[1][indb],_H[2][indb]);
    caxcb = cross(ca, cb);
    area = caxcb.norm();

    fprintf(fp,"%% area = %20.14e\n", area);
    fprintf(fp,"\n");

    H_origin = _H;
    for(n=0;n<_NP;n++) _VSR[n]=_R[n]; // Store original position

    call_potential();
    
    fprintf(fp,"%% Epot_0/area = %20.14e (eV/A^2)\n",_EPOT/area);
    fprintf(fp,"\n");    
    fprintf(fp,"%% i\t j\t k\t Epot/area (eV/A^2)\n");
    fprintf(fp,"%%-------------------------------------------\n");

    SHtoR();
    
    for(k=0;k<=nz;k++)
        for(j=0;j<=ny;j++)
            for(i=0;i<=nx;i++)
            {
                dc1 = c1*(x0+i*dx);
                dc2 = c2*(y0+j*dy);
                dc3 = c3*(z0+k*dz);

                _H[0][ind]= H_origin[0][ind]+dc1.x + dc2.x + dc3.x;
                _H[1][ind]= H_origin[1][ind]+dc1.y + dc2.y + dc3.y;
                _H[2][ind]= H_origin[2][ind]+dc1.z + dc2.z + dc3.z;

                for(n=0;n<_NP;n++) _R[n]=_VSR[n];
                RHtoS();

                /* forces neighborlist reconstruction every time (important!!) */
                /* clear Verlet neighbor list => enforce neighbor list updated at every time */
                for(n=0;n<_NP;n++) _R0[n].clear();
                call_potential();

                fprintf(fp, "%d\t %d\t %d\t %20.14e\n",i,j,k,_EPOT/area);

                //StreSS = _VIRIAL/(_H.det()*(1-_VACUUMRATIO))*160.2e3;
                //fprintf(fp,"%d %d %d %18.14e %lf %lf %lf %lf %lf %lf\n",p,q,r,_EPOT/area,StreSS[0][0],StreSS[0][1],StreSS[0][2],StreSS[1][1],StreSS[1][2],StreSS[2][2]);
            }
    fclose(fp);

}

void MDFrame::calmisfit2()
{
    int nx, ny, nz, i, j, k, n, ind, inda, indb;
    double x0, y0, z0, x1, y1, z1, dx, dy, dz, xmin, xmax, ymin, ymax, zmin, zmax, mx, my, mz, Lx, Ly, Lz, area;
    FILE *fp;
    class Vector3 c1, c2, c3, ca, cb, caxcb; 

    WARNING("**************************************************************");
    WARNING(" calmisfit2() is implemented only for rectangular supercells!");
    WARNING(" Use calmisfit() for tilted supercells.");
    WARNING("**************************************************************");

    /* direction of the surface normal */
    ind = (int)input[0]; 
    if(ind < 0 || ind > 3)
    {
        ERROR("Surface normal should be one of 1, 2, and 3.");
        return;
    }
    ind--;

    /* the other two directions orthogonal to ind */
    inda = (ind+1)%3; indb = (inda+1)%3;

    x0 = input[1]; dx = input[2]; x1 = input[3];
    y0 = input[4]; dy = input[5]; y1 = input[6];
    z0 = input[7]; dz = input[8]; z1 = input[9];
    xmin = input[10]; xmax = input[11];
    ymin = input[12]; ymax = input[13];
    zmin = input[14]; zmax = input[15]; 

    INFO_Printf("%%        x0 = %e, dx = %e, x1 = %e\n",x0,dx,x1);
    INFO_Printf("%%        y0 = %e, dy = %e, y1 = %e\n",y0,dy,y1);
    INFO_Printf("%%        z0 = %e, dz = %e, z1 = %e\n",z0,dz,z1);
    
    if(dx!=0) nx = (int) round((x1-x0)/dx); else nx = 0;
    if(dy!=0) ny = (int) round((y1-y0)/dy); else ny = 0;
    if(dz!=0) nz = (int) round((z1-z0)/dz); else nz = 0;

    INFO_Printf("%%        nx = %d, ny = %d, nz = %d\n",nx,ny,nz);

    fp = fopen("Emisfit.out","w");
    if(fp==NULL) FATAL("calmisfit2: open file failure");
        
    fprintf(fp,"%% Data file for misfit energy E(i,j,k) (eV/A^2)\n");
    fprintf(fp,"%% Format: i, j, k, E(i,j,k)\n");
    fprintf(fp,"%%         i=1:nx, j=1:ny, k=1:nz\n");
    fprintf(fp,"%%         nx = %d, ny = %d, nz = %d\n",nx,ny,nz);
    fprintf(fp,"\n");    
    fprintf(fp,"%%         x0 = %e, dx = %e, x1 = %e\n",x0,dx,x1);
    fprintf(fp,"%%         y0 = %e, dy = %e, y1 = %e\n",y0,dy,y1);
    fprintf(fp,"%%         z0 = %e, dz = %e, z1 = %e\n",z0,dz,z1);
    fprintf(fp,"\n");

    c1.set(_H[0][0],_H[1][0],_H[2][0]); 
    c2.set(_H[0][1],_H[1][1],_H[2][1]);
    c3.set(_H[0][2],_H[1][2],_H[2][2]);
    c1 /= latticesize[0][3]; // one lattice vector in each direction
    c2 /= latticesize[1][3];
    c3 /= latticesize[2][3];
    
    ca.set(_H[0][inda],_H[1][inda],_H[2][inda]); 
    cb.set(_H[0][indb],_H[1][indb],_H[2][indb]);
    caxcb = cross(ca, cb);
    area = caxcb.norm();

    fprintf(fp,"%% area = %20.14e\n", area);
    fprintf(fp,"\n");

    call_potential();
    
    fprintf(fp,"%% Epot_0/area = %20.14e (eV/A^2)\n",_EPOT/area);
    fprintf(fp,"\n");    
    fprintf(fp,"%% i\t j\t k\t Epot/area (eV/A^2)\n");
    fprintf(fp,"%%-------------------------------------------\n");

    Lx = c1.norm(); if(nx!=0) mx = Lx*(x1-x0)/nx; else mx=0;
    Ly = c2.norm(); if(ny!=0) my = Ly*(y1-y0)/ny; else my=0;
    Lz = c3.norm(); if(nz!=0) mz = Lz*(z1-z0)/nz; else mz=0;
    
    SHtoR();
    for(n=0;n<_NP;n++)
    {
        if((_SR[n].x>=xmin)&&(_SR[n].x<xmax)
           &&(_SR[n].y>=ymin)&&(_SR[n].y<ymax)
           &&(_SR[n].z>=zmin)&&(_SR[n].z<zmax))
        { // store the initial position
            _VSR[n].x = _R[n].x + Lx*x0;  
            _VSR[n].y = _R[n].y + Ly*y0; 
            _VSR[n].z = _R[n].z + Lz*z0;
        }
    }
    RHtoS();

    for(k=0;k<=nz;k++)
        for(j=0;j<=ny;j++)
            for(i=0;i<=nx;i++)
            {
                SHtoR();
                for(n=0;n<_NP;n++)
                {
                    if((_SR[n].x>=xmin)&&(_SR[n].x<xmax)
                       &&(_SR[n].y>=ymin)&&(_SR[n].y<ymax)
                       &&(_SR[n].z>=zmin)&&(_SR[n].z<zmax))
                    {
                        _R[n].x = _VSR[n].x + i*mx; 
                        _R[n].y = _VSR[n].y + j*my; 
                        _R[n].z = _VSR[n].z + k*mz;
                    }
                }
                //for(n=0;n<_NP;n++) _R[n]=_VSR[n];
                RHtoS();
                
                /* forces neighborlist reconstruction every time (important!!) */
                /* clear Verlet neighbor list => enforce neighbor list updated at every time */
                for(n=0;n<_NP;n++) _R0[n].clear();
                call_potential();

                fprintf(fp, "%d\t %d\t %d\t %20.14e\n",i,j,k,_EPOT/area);
                INFO_Printf("%d\t %d\t %d\t %20.14e\n",i,j,k,_EPOT/area);
            }
    fclose(fp);
}

void MDFrame::calmisfit2_rigidrlx()
{
    int nx, ny, nz, i, j, k, n, ind, inda, indb, ngrp, *grp_id;
    double x0, y0, z0, x1, y1, z1, dx, dy, dz, xmin, xmax, ymin, ymax, zmin, zmax, mx, my, mz, Lx, Ly, Lz, area;
    FILE *fp;
    class Vector3 c1, c2, c3, ca, cb, caxcb; 

    WARNING("**************************************************************");
    WARNING(" calmisfit2_rigidrlx() is implemented only for rectangular supercells!");
    WARNING(" Use calmisfit() for tilted supercells.");
    WARNING("**************************************************************");

    /* direction of the surface normal */
    ind = (int)input[0]; 
    if(ind < 0 || ind > 3)
    {
        ERROR("Surface normal should be one of 1, 2, and 3.");
        return;
    }
    ind--;

    /* the other two directions orthogonal to ind */
    inda = (ind+1)%3; indb = (inda+1)%3;

    x0 = input[1]; dx = input[2]; x1 = input[3];
    y0 = input[4]; dy = input[5]; y1 = input[6];
    z0 = input[7]; dz = input[8]; z1 = input[9];
    xmin = input[10]; xmax = input[11];
    ymin = input[12]; ymax = input[13];
    zmin = input[14]; zmax = input[15];
    ngrp = input[16]; 

    grp_id = new int [ngrp];
    memset(grp_id,0,sizeof(int)*ngrp);
    for(i=0;i<ngrp;i++) grp_id[i] = input[17+i];

    INFO_Printf("%%        x0 = %e, dx = %e, x1 = %e\n",x0,dx,x1);
    INFO_Printf("%%        y0 = %e, dy = %e, y1 = %e\n",y0,dy,y1);
    INFO_Printf("%%        z0 = %e, dz = %e, z1 = %e\n",z0,dz,z1);

    if (ind==0)
    {
        conj_fixdir[0]=0; conj_fixdir[1]=1; conj_fixdir[2]=1;
    }
    else if (ind==1)
    {
        conj_fixdir[0]=1; conj_fixdir[1]=0; conj_fixdir[2]=1;
    }
    else if (ind==2)
    {
        conj_fixdir[0]=1; conj_fixdir[1]=1; conj_fixdir[2]=0;
    }

    /* Toss ngrp and grp_id for rigid_relaxation */
    input[0]=ngrp;
    for(i=0;i<ngrp;i++) input[i+1]=grp_id[i];

    if(dx!=0) nx = (int) round((x1-x0)/dx); else nx = 0;
    if(dy!=0) ny = (int) round((y1-y0)/dy); else ny = 0;
    if(dz!=0) nz = (int) round((z1-z0)/dz); else nz = 0;

    INFO_Printf("%%        nx = %d, ny = %d, nz = %d\n",nx,ny,nz);

    fp = fopen("Emisfit.out","w");
    if(fp==NULL) FATAL("calmisfit2: open file failure");
        
    fprintf(fp,"%% Data file for misfit energy E(i,j,k) (eV/A^2)\n");
    fprintf(fp,"%% Format: i, j, k, E(i,j,k)\n");
    fprintf(fp,"%%         i=1:nx, j=1:ny, k=1:nz\n");
    fprintf(fp,"%%         nx = %d, ny = %d, nz = %d\n",nx,ny,nz);
    fprintf(fp,"\n");    
    fprintf(fp,"%%         x0 = %e, dx = %e, x1 = %e\n",x0,dx,x1);
    fprintf(fp,"%%         y0 = %e, dy = %e, y1 = %e\n",y0,dy,y1);
    fprintf(fp,"%%         z0 = %e, dz = %e, z1 = %e\n",z0,dz,z1);
    fprintf(fp,"\n");

    c1.set(_H[0][0],_H[1][0],_H[2][0]); 
    c2.set(_H[0][1],_H[1][1],_H[2][1]);
    c3.set(_H[0][2],_H[1][2],_H[2][2]);
    c1 /= latticesize[0][3]; // one lattice vector in each direction
    c2 /= latticesize[1][3];
    c3 /= latticesize[2][3];
    
    ca.set(_H[0][inda],_H[1][inda],_H[2][inda]); 
    cb.set(_H[0][indb],_H[1][indb],_H[2][indb]);
    caxcb = cross(ca, cb);
    area = caxcb.norm();

    fprintf(fp,"%% area = %20.14e\n", area);
    fprintf(fp,"\n");

    call_potential();
    
    fprintf(fp,"%% Epot_0/area = %20.14e (eV/A^2)\n",_EPOT/area);
    fprintf(fp,"\n");    
    fprintf(fp,"%% i\t j\t k\t Epot/area (eV/A^2)\n");
    fprintf(fp,"%%-------------------------------------------\n");

    Lx = c1.norm(); if(nx!=0) mx = Lx*(x1-x0)/nx; else mx=0;
    Ly = c2.norm(); if(ny!=0) my = Ly*(y1-y0)/ny; else my=0;
    Lz = c3.norm(); if(nz!=0) mz = Lz*(z1-z0)/nz; else mz=0;
    
    SHtoR();
    for(n=0;n<_NP;n++)
    {
        if((_SR[n].x>=xmin)&&(_SR[n].x<xmax)
           &&(_SR[n].y>=ymin)&&(_SR[n].y<ymax)
           &&(_SR[n].z>=zmin)&&(_SR[n].z<zmax))
        { // store the initial position
            _VSR[n].x = _R[n].x + Lx*x0;  
            _VSR[n].y = _R[n].y + Ly*y0; 
            _VSR[n].z = _R[n].z + Lz*z0;
            group[n] = grp_id[0];
        }
    }
    RHtoS();

    for(k=0;k<=nz;k++)
        for(j=0;j<=ny;j++)
            for(i=0;i<=nx;i++)
            {
                SHtoR();
                for(n=0;n<_NP;n++)
                {
                    if(group[n]==grp_id[0])
                    {
                        _R[n].x = _VSR[n].x + i*mx; 
                        _R[n].y = _VSR[n].y + j*my; 
                        _R[n].z = _VSR[n].z + k*mz;
                    }
                }
                //for(n=0;n<_NP;n++) _R[n]=_VSR[n];
                RHtoS();
                
                /* forces neighborlist reconstruction every time (important!!) */
                /* clear Verlet neighbor list => enforce neighbor list updated at every time */
                for(n=0;n<_NP;n++) _R0[n].clear();
                rigid_relaxation();
                call_potential();

                INFO_Printf("%d\t %d\t %d\t %20.14e\n",i,j,k,_EPOT/area);
                fprintf(fp, "%d\t %d\t %d\t %20.14e\n",i,j,k,_EPOT/area);
                fflush(fp);
            }
    fclose(fp);
    delete[] grp_id;
}

void MDFrame::testpotential()
{
    int i, n, ipt;
    double E0, fx0, fy0, fz0, fx, fy, fz, dx, f2, fmax;
    int imax;
    Vector3 r0;

    call_potential();
    E0=_EPOT;
    f2=_F[0].norm2(); imax=0; fmax=f2;
    for(i=0;i<_NP;i++)
    {
        f2=_F[i].norm2();
        if(f2>1e4) INFO_Printf("atom [%d] f=(%e %e %e)\n",i,_F[i].x,_F[i].y,_F[i].z);
        if(fmax<f2)
        {
            imax=i;
            fmax=f2;
        }
    }
    INFO_Printf("max force\n");
    INFO_Printf("atom [%d] f=(%e %e %e)\n",imax,_F[imax].x,_F[imax].y,_F[imax].z);
    n=(int)input[0];
    dx=input[1];
    for(i=1;i<=n;i++)
    {
        ipt=(int)input[i+1];
        SHtoR();
        call_potential();
        E0=_EPOT;
        fx0=_F[ipt].x;
        fy0=_F[ipt].y;
        fz0=_F[ipt].z;
        
        r0=_R[ipt];
        
        _R[ipt].x+=dx; RHtoS();
        call_potential();
        fx=(E0-_EPOT)/dx;
        _R[ipt]=r0;    RHtoS();
        
        _R[ipt].y+=dx; RHtoS();
        call_potential();
        fy=(E0-_EPOT)/dx;
        _R[ipt]=r0;    RHtoS();
        
        _R[ipt].z+=dx; RHtoS();
        call_potential();
        fz=(E0-_EPOT)/dx;
        _R[ipt]=r0;    RHtoS();

        INFO_Printf("atom [%d] \tf0  =(%e %e %e)\n"
                    "          \tfnum=(%e %e %e)\n",
                    ipt, fx0, fy0, fz0, fx, fy, fz);
    }
}

void MDFrame::initvelocity()
{
    int i, n, n_c=0;
    Vector3 vavg,v2avg,r,v,omega,ds;
    Matrix33 hinv, I, Inertia, dInertia, Iinv, rotimesr;
    double v2desired, v2, totalmass, tmp;
#ifdef _TORSION_OR_BENDING    
    double tmp2, tmp4, curvature, Radius;
#endif
    
    v2desired = 2*(KB*_TDES)/(_ATOMMASS[0]*MASSCONVERT)
                *SQR(_TIMESTEP);

    if(_DOUBLET!=1) v2desired*=.5;
    
    vavg.clear(); v2avg.clear();
    
    n=0; totalmass=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i])
        {
            _VR[i].clear();
            continue;
        }

        if ( strcmp(initvelocity_type, "Uniform") == 0 )
        {// uniform distribution
            _VR[i].x=drand48();
            _VR[i].y=drand48();
            _VR[i].z=drand48();
        } 
        else if ( strcmp(initvelocity_type,"Gaussian") == 0)
        {// Maxwell-Boltzmann distribution
            _VR[i].x=randnorm48();
            _VR[i].y=randnorm48();
            _VR[i].z=randnorm48();
        } 
        else 
        {
            ERROR("unknown initvelocity_type "<<initvelocity_type);
        } 

        if(species[i]!=0)
        {
            tmp = sqrt(_ATOMMASS[0]/_ATOMMASS[species[i]]);
            _VR[i] *= tmp;
        }
        vavg+=_VR[i]*_ATOMMASS[species[i]];
        totalmass+=_ATOMMASS[species[i]];
        n++;
    }

    /* zero out COM translation velocity */
    vavg/=totalmass;
    for(i=0;i<_NP;i++)
    {//_VR has no unit here, until multiplied with sqrt(v2desired/v2avg)
        if(fixed[i]) continue;
        _VR[i]-=vavg;
    }
    
    /* zero out COM angular velocity */
    if ((strcmp(zerorot,"x")==0) || (strcmp(zerorot,"y")==0) || (strcmp(zerorot,"z")==0) || (strcmp(zerorot,"all")==0))
    {
        INFO_Printf("Zero out %s component(s) of the angular momentum\n",zerorot);
        VRHtoVS();
        zero_angmom();

        // number of additional fixed degree of freedom
        if (strcmp(zerorot,"all")==0) 
	    n_c=3;
        else
	    n_c=1;
    }

#ifdef _TORSION_OR_BENDING
    /* zero out COM rotation velocity */
    if (_TORSIONSIM)
    {
        /* Compute _PTHETA_COM */
        _PTHETA_COM = 0;
        tmp2 = tmp4 = 0;
	for(i=0;i<_NP0;i++)
        {
            r = _H*_SR[i];
            v = _H*_VSR[i];
            if(species[i]!=0)
                v *= _ATOMMASS[species[i]] / _ATOMMASS[0];
            tmp2 += (v.y*r.x - v.x*r.y);
            tmp4 += (r.x*r.x + r.y*r.y);
        }
        _PTHETA_COM = tmp2 / tmp4;

        for(i=0;i<_NP0;i++)
        {
            _VR[i].x += r.y*_PTHETA_COM;
            _VR[i].y -= r.x*_PTHETA_COM;
        }
        _PTHETA_COM = 0;
    }

    if (_BENDSIM)
    {
        /* Compute _PTHETA_COM */
        _PTHETA_COM = 0;
        tmp2 = tmp4 = 0;
        curvature = bendspec[3];
        Radius = 1.0/curvature;
        for(i=0;i<_NP0;i++)
        {
            r = _H*_SR[i];
            r.y += Radius;
            v = _H*_VSR[i];
            if(species[i]!=0)
                v *= _ATOMMASS[species[i]] / _ATOMMASS[0];
            tmp2 += (v.y*r.x - v.x*r.y);
            tmp4 += (r.x*r.x + r.y*r.y);
        }
        _PTHETA_COM = tmp2 / tmp4;
        for(i=0;i<_NP0;i++)
        {
            _VR[i].x += r.y*_PTHETA_COM;
            _VR[i].y -= r.x*_PTHETA_COM;
        }
        _PTHETA_COM = 0;
    }
#endif
    
    if(_NIMAGES==0) _NP0=_NP;
    v2 = 0.0;
    for(i=0;i<_NP0;i++)
    {
        if(species[i]==0)
	{
            v2avg+=_VR[i].sq(); v2+=_VR[i].norm2();
        } 
        else
	{
            v2avg+=_VR[i].sq() * _ATOMMASS[species[i]]/_ATOMMASS[0];
            v2+=_VR[i].norm2() * _ATOMMASS[species[i]]/_ATOMMASS[0];
        }
    }
    v2avg/=(n-1)*v2desired;v2avg=v2avg.sqroot();
    v2/=(3*n-3-n_c)*v2desired; v2=sqrt(v2);

    _KATOM=0;
    for(i=0;i<_NP0;i++)
    {/* Now _VR has unit of length */
        if (strcmp(zerorot,"none")==0)
        {
            if(v2avg.x>0)_VR[i].x/=v2avg.x;
            if(v2avg.y>0)_VR[i].y/=v2avg.y;
            if(v2avg.z>0)_VR[i].z/=v2avg.z;
        }
        else 
            _VR[i]/=v2;
	
	tmp=0.5*(_ATOMMASS[0]*MASSCONVERT)/(_TIMESTEP*_TIMESTEP);
        if (species[i]==0)
	    _KATOM+=tmp*_VR[i].norm2();
	else
	    _KATOM+=tmp*_VR[i].norm2()*_ATOMMASS[species[i]]/_ATOMMASS[0];
    }
    _T=_KATOM/(1.5*KB*n);

    /* compute linear and angular momenta */
    findcenterofmass(_R);    // center of mass _COM in real coordinate
    _P_COM.clear(); _L_COM.clear();
    for (i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        r = _R[i] - _COM;
        if(species[i]==0)
        {
            _P_COM += _VR[i]; _L_COM += cross(r,_VR[i]);
        }
        else
        {
	    v = _VR[i] * _ATOMMASS[species[i]] / _ATOMMASS[0];
            _P_COM += v; _L_COM += cross(_R[i],v);
	}
    }
    _P_COM *= (_ATOMMASS[0]*MASSCONVERT/_TIMESTEP);
    _L_COM *= (_ATOMMASS[0]*MASSCONVERT/_TIMESTEP);
    INFO("P_COM = [ "<<_P_COM<<" ] in eV/(A/ps)");
    INFO("L_COM = [ "<<_L_COM<<" ] in eV/(1/ps)");


    hinv=_H.inv();
    for(i=0;i<_NP0;i++) /* _VSR has no unit. */
        _VSR[i]=hinv*_VR[i];

    _KATOM=3.0*BOLZ*_TDES/EV*(n-1);
    if(_DOUBLET!=1) _KATOM*=.5;

    _VH.clear();

}

void MDFrame::perturbevelocity()
{
    int i, n;
    double fraction ; 
    Vector3 vavg, r, v, w, L, Rcm, scm,v2avg;
    Matrix33 hinv,Inertia,iinv;
    double v2desired, totalmass, totalmass1, tmp, tmp1;
    double rnorm;
    
    fraction = input[0];
    if (YES_HMC!=1)
        v2desired = 2* (KB*_TDES)/(_ATOMMASS[0]*MASSCONVERT)*SQR(_TIMESTEP);
    else
        v2desired = 2* (KB*T_HMC)/(_ATOMMASS[0]*MASSCONVERT)*SQR(_TIMESTEP);
    if(_DOUBLET!=1) v2desired*=.5;
    
    //
//  INFO_Printf("vdesire %f frac %f \n",v2desired, fraction);

    vavg.clear(); scm.clear(); Rcm.clear(); 
    n=0; totalmass=0; totalmass1=0;
    for(i=0;i<_NP;i++)
    {
        scm+=_SR[i]*_ATOMMASS[species[i]];
        totalmass1+=_ATOMMASS[species[i]];

	if(fixed[i])
	{
	    continue;
	}
        _VR[i]=_H*_VSR[i];
//        INFO_Printf("i %d vsr %f vr %f rn %f\n",i,_VSR[i].x,_VR[i].x,randnorm());
        tmp = v2desired*sqrt(_ATOMMASS[0]/_ATOMMASS[species[i]]);
	_VR[i].x = tmp * fraction*(randnorm48());
        _VR[i].y = tmp * fraction*(randnorm48());
        _VR[i].z = tmp * fraction*(randnorm48());
//        INFO_Printf("tmp %f vx %f\n",tmp,_VR[i].x);
	vavg+=_VR[i]*_ATOMMASS[species[i]];
	totalmass+=_ATOMMASS[species[i]];
        n++;
    }
    scm /= totalmass1;
    Rcm = _H*scm;
    // zero out COM translational velocity 
//INFO_Printf("before perturbe vavgx=%e, vavgy=%e, vavgz=%e\n",vavg.x,vavg.y,vavg.z);
    vavg/=totalmass;
    for(i=0;i<_NP;i++)
    {
	if(fixed[i]) continue;
	_VR[i]-=vavg;
    }
//    INFO_Printf("after perturbe\nvavgx %10.8e vavgy %10.8e vavgz %10.8e\n",vavg[0],vavg[1],vavg[2]); 
//    INFO_Printf("atommass %10.8e totmass %10.8e totmass1 %10.8e\n",_ATOMMASS[0],totalmass,totalmass1);

//    INFO_Printf("Rcmx = %10.8e Rcmy = %10.8e Rcmz = %10.8e \n\n", Rcm[0],Rcm[1],Rcm[2]);
    hinv=_H.inv();
    // zero out COM rotational velocity
    Inertia.clear();
    L.clear();w.clear();vavg.clear();

    for (i=0;i<_NP;i++)
    {
        r = _H*_SR[i] - Rcm;
        rnorm = r.norm(); 
        v = _VR[i];
        vavg+=_VR[i]*_ATOMMASS[species[i]];
	tmp = (_ATOMMASS[species[i]]);
        Inertia[0][0]+= (rnorm*rnorm - r[0]*r[0])*tmp;
        Inertia[0][1]+= (            - r[0]*r[1])*tmp;
        Inertia[0][2]+= (            - r[0]*r[2])*tmp;
        Inertia[1][0]+= (            - r[1]*r[0])*tmp;
        Inertia[1][1]+= (rnorm*rnorm - r[1]*r[1])*tmp;
        Inertia[1][2]+= (            - r[1]*r[2])*tmp;
        Inertia[2][0]+= (            - r[2]*r[0])*tmp;
        Inertia[2][1]+= (            - r[2]*r[1])*tmp;
        Inertia[2][2]+= (rnorm*rnorm - r[2]*r[2])*tmp;
        L[0]         += (r[1]*v[2]  -  r[2]*v[1])*tmp;   
        L[1]         += (r[2]*v[0]  -  r[0]*v[2])*tmp;
        L[2]         += (r[0]*v[1]  -  r[1]*v[0])*tmp;
    }
//        INFO_Printf("before Lx = %10.8e, Ly = %10.8e, Lz = %10.8e\n",L[0],L[1],L[2]);
//        INFO_Printf("i %d rn %f vx %f rx %f In00 %f L0 %f\n", i, rnorm, v.x, r.x, Inertia[0][0], L[0]); 
    
    vavg/=totalmass;
    iinv=Inertia.inv();
    w=iinv*L;  
//    INFO_Printf("wx %10.8e wy %10.8e wz %10.8e\n", w[0],w[1],w[2]);
    L.clear(); vavg.clear(); v2avg.clear();
    for (i=0;i<_NP;i++)
    {
        r = _H*_SR[i] - Rcm;
        _VR[i][0]    -= (w[1]*r[2]  -  w[2]*r[1]);
        _VR[i][1]    -= (w[2]*r[0]  -  w[0]*r[2]);
        _VR[i][2]    -= (w[0]*r[1]  -  w[1]*r[0]);
        vavg+=_VR[i]*_ATOMMASS[species[i]];
       
	if (species[i]==0)
  	    v2avg+=_VR[i].sq();
	else
	    v2avg+=_VR[i].sq()*_ATOMMASS[species[i]]/_ATOMMASS[0];

        _VSR[i] = hinv * _VR[i];
//        INFO_Printf("i %d vr %f vsr %f\n",i,_VR[i].x,_VSR[i].x); 
        v = _VR[i];
        L[0]         += (r[1]*v[2]  -  r[2]*v[1])*_ATOMMASS[species[i]];
        L[1]         += (r[2]*v[0]  -  r[0]*v[2])*_ATOMMASS[species[i]];
        L[2]         += (r[0]*v[1]  -  r[1]*v[0])*_ATOMMASS[species[i]];
    }
    vavg/=totalmass;
//    v2avg/=(n-1)*v2desired;v2avg=v2avg.sqroot();
    if (YES_HMC!=1)
        tmp=(v2avg.x+v2avg.y+v2avg.z)/3/( (n-1)*v2desired );
    else 
	tmp=(v2avg.x+v2avg.y+v2avg.z)/3/( n*v2desired );
    tmp=sqrt(tmp);
//INFO_Printf("tmp=%e _VR0x=%e _Rx=%e wx=%e",tmp,_VR[0].x,_R[0].x,w[0]);

    _KATOM=0;
    for (i=0;i<_NP;i++)
    {
	if(tmp>0) {
	_VR[i].x/=tmp;
	_VR[i].y/=tmp;
	_VR[i].z/=tmp;
	}
 
        tmp1=0.5*(_ATOMMASS[0]*MASSCONVERT)/(_TIMESTEP*_TIMESTEP);
	if (species[i]==0)
	    _KATOM+=tmp1*_VR[i].norm2();
	else
	    _KATOM+=tmp1*_VR[i].norm2()*_ATOMMASS[species[i]]/_ATOMMASS[0];
    }
//    _T=_KATOM/(1.5*KB*n);

    for (i=0;i<_NP;i++) _VSR[i]=hinv*_VR[i];
//    INFO_Printf("after Lx %10.8e Ly %10.8e Lz %10.8e\n", L[0],L[1],L[2]);
//    INFO_Printf("vavgx %10.8e vavgy %10.8e vavgz %10.8e\n",vavg[0],vavg[1],vavg[2]);
//    L.clear();
//    for (i=0;i<_NP;i++)
//   {
//	r=_H*_SR[i]-Rcm;
//        v=_VR[i];
//        L[0]+=(r[1]*v[2]-r[2]*v[1])*_ATOMMASS[species[i]];
//        L[1]+=(r[2]*v[0]-r[0]*v[2])*_ATOMMASS[species[i]];
//        L[2]+=(r[0]*v[1]-r[1]*v[0])*_ATOMMASS[species[i]];
//   } 
//    INFO_Printf("Lx =%e Ly= %e Lz= %e",L[0],L[1],L[2]);
}  

void MDFrame::MCperturbevelocity()
{
    int i;
    double tmp, Ktmp;
    double Ko,Kn,v2;
    int MC_accept, MC_trial;
    Matrix33 hinv;
    
    Ktmp=0.5*(_ATOMMASS[0]*MASSCONVERT)/(_TIMESTEP*_TIMESTEP);
    hinv=_H.inv();
    Ko=0;
    MC_trial=0;

    for (i=0;i<_NP;i++)
    {
        _R[i]=_H*_VSR[i];
        v2=_R[i].norm2(); 
        Ko+=Ktmp*v2*_ATOMMASS[species[i]]/_ATOMMASS[0];
    }

    while(1) {
    MC_trial++;
    Kn=0;
    perturbevelocity();
    for (i=0;i<_NP;i++)
    {
        v2=_VR[i].norm2();
        Kn+=Ktmp*v2*_ATOMMASS[species[i]]/_ATOMMASS[0];
    }

    if (strcmp(ensemble_type,"NVE")==0) 
    {   
        tmp=sqrt(Ko/Kn);
        Kn=0;
        for (i=0;i<_NP;i++)
        {
	    _R[i]=_H*_SR[i];
            _VR[i]*=tmp;
            _VSR[i]*=tmp;
            v2=_VR[i].norm2();
            Kn+=Ktmp*v2*_ATOMMASS[species[i]]/_ATOMMASS[0];
        }
//        INFO_Printf("NVE Kfinal = %10.8e (%f K)\n",Kn,Kn/(1.5*KB*(_NP-1)));
        break;
    }
    else
    {
        if(Kn<Ko)
            MC_accept=1;
        else
        {
            tmp=(Ko-Kn)/(KB*_TDES);
//    INFO_Printf("KBT = %10.8e tmp = %10.8e\n",KB*_TDES,tmp);
            if (tmp<-40) MC_accept=0;
            else
            {
                if(drand48()<exp(tmp))
		    MC_accept=1;
	        else
		    MC_accept=0;
	    }
        } 
//    INFO_Printf("MC_accept = %d\n",MC_accept);

        if (MC_accept==1) 
        {
            for(i=0;i<_NP;i++)
	    {
	        _R[i]=_H*_SR[i];
            }
            INFO_Printf("MC_trial = %d\n",MC_trial);
            break;
        }
        else
        {
            for(i=0;i<_NP;i++)
            {
                _VR[i]=_R[i];
                _VSR[i]=hinv*_VR[i];
  	    }
        }
    }

    }
}

void MDFrame::zero_angmom()
{ /* Zero out angular momentum.  July 01 2009, Keonwook Kang */
    int i;
    double rnorm2, massratio;
    Vector3 r, v, va, omega, vavg;
    Matrix33 hinv, I, Inertia, dInertia, Iinv, rotimesr;

    // 3x3 identity matrix I
    I.eye(); 
    
    SHtoR(); VSHtoVR();
    findcenterofmass(_R);    // center of mass in scaled coordinate
    
    Inertia.clear(); _L_COM.clear(); rotimesr.clear();
    for (i=0;i<_NP;i++)
    {// calculate intertia matrix Inertia and angular momentum _L_COM

        if(fixed[i]) continue;
    
        r = _R[i] - _COM; v = _VR[i];
        rnorm2 = r.norm2(); rotimesr.dyad(r,r);
        if(species[i]==0)
	{
            Inertia += (I*rnorm2 - rotimesr); _L_COM += cross(r,v);
        }
        else
	{
	    massratio = _ATOMMASS[species[i]] / _ATOMMASS[0];
            v *= massratio;
            dInertia = (I*rnorm2 - rotimesr) * massratio;
            Inertia += dInertia;
            _L_COM  += cross(r,v);
	}
    }
    Iinv=Inertia.inv(); omega=Iinv*_L_COM;
    // Don't need to include the following lines before omega, 
    // cause mass will be cancelled each other.
    //Inertia *= (_ATOMMASS[0]*MASSCONVERT); 
    //_L_COM  *= (_ATOMMASS[0]*MASSCONVERT/_TIMESTEP);
    //INFO("Before zero_angmom,");
    //INFO("    L_COM = [ "<<_L_COM<<" ] in eV/(1/ps)");
    //INFO("  Inertia = [\n"<<Inertia<<" ] in eV*ps^2");
    //INFO("    omega = [ "<<omega/_TIMESTEP<<" ] in 1/ps");

    //_L_COM.clear();
    for (i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        r = _R[i] - _COM;

        va = cross(omega,r);
        if (strcmp(zerorot,"x")==0) {_VR[i].y -= va.y; _VR[i].z -= va.z;}
        if (strcmp(zerorot,"y")==0) {_VR[i].z -= va.z; _VR[i].x -= va.x;}
        if (strcmp(zerorot,"z")==0) {_VR[i].x -= va.x; _VR[i].y -= va.y;}
        if (strcmp(zerorot,"all")==0)  _VR[i] -= va;        

    //	if(species[i]==0)
    //    {
    //        _P_COM += _VR[i]; _L_COM += cross(r,_VR[i]);
    //    }
    //    else
    //	{
    //        v = _VR[i] * _ATOMMASS[species[i]] / _ATOMMASS[0];
    //        _P_COM += v; _L_COM += cross(_R[i],v);
    //	}
    }
    //_P_COM*=(_ATOMMASS[0]*MASSCONVERT/_TIMESTEP); // in eV/(A/ps)
    //_L_COM*=(_ATOMMASS[0]*MASSCONVERT/_TIMESTEP); // in eV/(1/ps)
    //INFO("After zero_angmom,");
    //INFO("    P_COM = [ "<<_P_COM<<" ] in eV/(A/ps)");
    //INFO("    L_COM = [ "<<_L_COM<<" ] in eV/(1/ps)");

    VRHtoVS();

}

void MDFrame::zero_totmom()
{

    int i, n;
    Vector3 vavg;
    double totalmass;

    vavg.clear();
    n=0; totalmass=0;

    for(i=0;i<_NP;i++)
    {

        if(fixed[i]) continue;
        vavg+=_VR[i]*_ATOMMASS[species[i]];
        totalmass+=_ATOMMASS[species[i]];
        n++;
    }

    /* zero out COM translation velocity */
    vavg/=totalmass;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        _VR[i]-=vavg;
    }
    VRHtoVS();
}


void MDFrame::randomposition()   
{
    int i;
    for(i=0;i<_NP;i++)
    {
        _SR[i].x=drand48()-0.5;
        _SR[i].y=drand48()-0.5;
        _SR[i].z=drand48()-0.5;
    }
}

#if 0 /* disable velocity scaler */
void MDFrame::velscaler()
{
//    INFO("use velscaler");
    /* Box temperature rescaling not implemented */
    double katomdesired,eratio,vratio; int i,n;

    _NPfree=0;
    for(i=0;i<_NP;i++)
    {
        if(!fixed[i]) _NPfree++;
    }
    katomdesired=1.5*KB*_TDES*(_NPfree-1);
    eratio=katomdesired/_KATOM;
    eratio=min(eratio,1000);
    vratio=(eratio-1.)/_ATOMTCPL+1.;

    for(i=0;i<_NP;i++)
    {
        _VSR[i]*=vratio;
    }
}
#endif

void MDFrame::scaleVel()
{
    double vratio; int i;

    vratio = input[0];
    for(i=0;i<_NP;i++)
    {
        if(!fixed[i])
            _VSR[i]*=vratio;
    }
}


void MDFrame::multiplyvelocity()
{
    double factor; int i, itmp;

    factor = input[0];
    for(i=0;i<_NP;i++)
    {
        _VSR[i]*=factor;
    }
    itmp = algorithm_id / 1000; /* 10: NVE, 20: NVT, 30: NPH, 40: NPT */
    if((itmp==20)||(itmp==40))
    {
        zetav*=factor;
    }
    if(itmp==21) /* NVTC integrator, Keonwook Kang Mar 03 2008 */
        for(i=0;i<NHChainLen;i++)
            zetaNHCv[i]*=factor;

    if((itmp==30)||(itmp==40))   /* Feb 14 2007 Keonwook Kang, for reversibility test */
    {
        _VH*=factor;
    }
}



void MDFrame::step()
{
    /* call external function */
    integrator();
    call_ANDERSON();

    /* velscaler is disabled, use NVT, or NPT ensemble instead (Nose Hoover) */
//    if(usescalevelocity&&(!usenosehoover))
//        velscaler();

    /* accumulate reversible work during switching */
    if (_ENABLE_SWITCH)
    {
       _WTOT += dEdlambda * dlambdadt * (_LAMBDA1-_LAMBDA0);
       //_WAVG = _WTOT/(curstep+1);
       _WAVG = _WTOT / totalsteps;
    }
}
    
void MDFrame::eval()
{
    /* multi process functions */
    call_potential();
    calcprop();
    printResult();
}

void MDFrame::eval_insertparticle()
{
    int particle_id;

    particle_id = _NP; 
    _SR[particle_id].x = drand48()-0.5;
    _SR[particle_id].y = drand48()-0.5;
    _SR[particle_id].z = drand48()-0.5;
    fixed[particle_id]=0;
    
    /* make sure allocmultiple > 1 before calling this function */
    if(allocmultiple<=1)
        FATAL("eval_insertparticle cannot run with allocmultiple = "<<allocmultiple);
    NbrList_reconstruct_use_link_list();
    NbrList_reconstruct_use_link_list(particle_id); /* append list */
    potential_energyonly(particle_id);
    INFO_Printf("EPOT_RMV[%d] = %20.12e\n",particle_id,_EPOT_RMV[particle_id]);

#if 0 /* debug, compare with brute force subtraction */
    double E0, E1, dE;

    INFO("------");
    _NP ++;
    NbrList_reconstruct_use_link_list();
    potential_energyonly();

    E1 = _EPOT;
    INFO("------");
    
    _NP --;
    NbrList_reconstruct_use_link_list();
    potential_energyonly();
    E0 = _EPOT;

    INFO("------");
    NbrList_reconstruct_use_link_list();
    NbrList_reconstruct_use_link_list(particle_id); /* append list */
    potential_energyonly(particle_id);
    dE = _EPOT_RMV[particle_id];
    
    INFO_Printf("E1=%20.12e E0=%20.12e E1-E0=%20.12e dE=%20.12e",E1,E0,E1-E0,dE);        
#endif
}

void MDFrame::eval_removeparticle()
{
    int particle_id;
    
    particle_id = (int) input[0]; 
    potential_energyonly();
    INFO_Printf("EPOT_RMV[%d] = %20.12e\n",particle_id,_EPOT_RMV[particle_id]);

#if 0 /* debug, compare with brute force subtraction */    
    double E0, E1, dE;
    dE = _EPOT_RMV[particle_id];
    E1 = _EPOT;

    fixed[particle_id]=-1;
    potential_energyonly();
    E0 = _EPOT;

    fixed[particle_id]=0;
    INFO_Printf("E1=%20.12e E0=%20.12e E1-E0=%20.12e dE=%20.12e",E1,E0,E1-E0,dE);
#endif
}

void MDFrame::printResult()
{
    calcoutput();
    if(strlen(output_str)==0)
    {
        INFO_Printf("EPOT=%18.14e\n",_EPOT);
        INFO_Printf("KATOM=%18.14e\n",_KATOM);
        INFO_Printf("Tinst=%18.14e\n",_T);
        INFO_Printf("PRESSURE=%18.14e in eV/Angstrom^3\n",_TOTSTRESS.trace()/3.);
        INFO("\nVIRIAL=["<<_VIRIAL<<"]\n"
             "\nH=["<<_H<<"]\n");
        INFO("\nStress (in MPa, not including momentum term)=\n["<<_VIRIAL/(_H.det()*(1-_VACUUMRATIO))*160.2e3<<"]\n");
        INFO("\nTotal Stress (in MPa, including momentum term)=\n["<<_TOTSTRESSinMPa<<"]\n");
    }
    else
    {
        INFO_Printf("Format: %s\n",output_fmt);
        INFO_Printf("Result: %s\n",output_str);
    }
}
    
void MDFrame::multieval()
{
    int i;
    for(i=0;i<totalsteps;i++)
        call_potential();
}


/********************************************************************/
/* Functions for Umbrella Sampling in MC simulations                */
/********************************************************************/


/********************************************************************/
/* Functions for Forward Flux Sampling in MC simulations            */
/********************************************************************/


double MDFrame::potential_energyonly_before(int iatom)
{
    return potential_energyonly(iatom);
}

double MDFrame::potential_energyonly_after(int iatom)
{
    return potential_energyonly(iatom);
}







/**********************************************************/
/* Function for MC simulations                            */
/**********************************************************/


/****************************************************************/
/* Function for Adiabatic Switching simulations for Free Energy */
/****************************************************************/

double MDFrame::SWITCH_Find_Lambda()
{
     double lambda, t, s, t5, t4, t3, t2;
     
     t = ((double) curstep)/ (totalsteps-1);
     if(_SWITCHFUNC == 0)
     {
         t2 = t*t;
         t3 = t2*t;
         t4 = t3*t;
         t5 = t4*t;
         s = (t5*(70.*t4-315.*t3+540.*t2-420.*t+126.));
         dlambdadt = (t4*(630.*t4-2520.*t3+3780.*t2-2520.*t+630.));
     }
     else if(_SWITCHFUNC == 1)
     { /* a linear switching function */
         s = t;
         dlambdadt = 1;
     }
     else
     {
         s = t;
         dlambdadt = 1;
         ERROR("unknown switchfunc!");
     }
     lambda = _LAMBDA0 +  s * (_LAMBDA1-_LAMBDA0) ;
     return lambda;
}




/********************************************************************/
/* Nudged-Elastic-Band relaxation */
/********************************************************************/

/* Memory Allocation */
void MDFrame::AllocChain()
{/* need to fix this allocation later */
    int n,size,i;
    void *c;
    
    n=constrainatoms[0];
    size=sizeof(Vector3 *)*(_CHAINLENGTH+1)+sizeof(Vector3)*(_CHAINLENGTH+1)*n + 10;

    if(n<=0)
    {
            ERROR("Constrained atoms not set! (atoms not selected)");
            return;
    }
    if(n>MAXCONSTRAINATOMS)
    {
        ERROR("Number of constrained atoms is too big!");
        return;
    }
   
    INFO_Printf("AllocChain: malloc size = %d\n",size);
    c=malloc(size); memset(c,0,size);
    _Rc=(Vector3 **) c;   _Rc[0]=(Vector3 *)(_Rc+(_CHAINLENGTH+1));
    for(i=1;i<=_CHAINLENGTH;i++) _Rc[i]=_Rc[i-1]+n;

    c=malloc(size); memset(c,0,size);
    _Rctmp=(Vector3 **) c;   _Rctmp[0]=(Vector3 *)(_Rctmp+(_CHAINLENGTH+1));
    for(i=1;i<=_CHAINLENGTH;i++) _Rctmp[i]=_Rctmp[i-1]+n;

    c=malloc(size); memset(c,0,size);
    _Fc=(Vector3 **) c;   _Fc[0]=(Vector3 *)(_Fc+(_CHAINLENGTH+1));
    for(i=1;i<=_CHAINLENGTH;i++) _Fc[i]=_Fc[i-1]+n;
    
    _nebinterp=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1));
    memset(_nebinterp,0,sizeof(double)*(_CHAINLENGTH+1));
    
    INFO("Memory for _Rc and _Fc are allocated.");
}

void MDFrame::caldscom12()
{
    /* center of mass : sum(m_i*r_i)/sum(m_i)
                         = sum(r_i)/N         if all atoms have same mass.*/
    /* shift of c.o.m between two states A and B:
             sum(r_i^B)/N - sum(r_i^A)/N = sum(r_i^B - r_i^A)/N           */
    int i, n, ipt;

    n=constrainatoms[0];

    /* calculate dscom12 */
    dscom12.clear();
    for(i=1;i<=n;i++)
    {
        ipt=constrainatoms[i];
        dscom12+=_SR2[ipt];
        dscom12-=_SR1[ipt];
    }
    dscom12/=n;

#if 0
    /* rid center of mass shift in configB */
    for(i=0;i<_NP;i++)
    {
        _SR2[i]-=dscom12;
        _SR2[i].subint();
    }
    dscom12.clear();
#endif    
    INFO_Printf("center-of-mass shift: [%f %f %f]\n",dscom12.x,dscom12.y,dscom12.z);
}
    
int MDFrame::readRchain()
{
    int i, j, nchain;
    char *buffer; char *pp, *q; 
    char fname[300];
 
    LFile::SubHomeDir(incnfile,fname);
    INFO("filename="<<fname);
    LFile::LoadToString(fname,buffer,0);

    pp=buffer;
    
    q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
    sscanf(q, "%d", &nchain);

    /* Default value of readrchainspec = -1.
       If you want to read only first n chains, set readrchainspec=n. */
    if(readrchainspec>=0 && readrchainspec<nchain)
        _CHAINLENGTH = readrchainspec;
    else  
        _CHAINLENGTH = nchain;
    INFO("CHAINLEGNTH="<<_CHAINLENGTH);

    q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
    sscanf(q, "%d", constrainatoms);

    INFO("constrainatoms[0]="<<constrainatoms[0]);

    /* Bug fixed by Keonwook Kang Apr/20/2010
       : Read constrainatoms before caldscom12() and initRchain() */
    for(i=0;i<constrainatoms[0];i++)
    {
            q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
            sscanf(q, "%d", constrainatoms+i+1);
    }

    if (_Rc==NULL) // If _Rc is not allocated
    {
        INFO("call AllocChain()");
        AllocChain();     // Mar 08 2007 Keonwook Kang         
        caldscom12();     // Apr 02 2010 Keonwook Kang
        initRchain();     // Apr 02 2010 Keonwook Kang
    } else 
        INFO("Warning: _Rc will be overwritten!!");

    for(j=0;j<=_CHAINLENGTH;j++)
    {
        for(i=0;i<constrainatoms[0];i++)
        {
            q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
            sscanf(q, "%lf %lf %lf",&(_Rc[j][i].x),&(_Rc[j][i].y),&(_Rc[j][i].z));
        }
    }

    Free(buffer);
    INFO("readRchain finished");
    return 0;
}

int MDFrame::readRchain_parallel_toCN(int j)
{
    int i;
    char fname[200], fullname[200];
    char *buffer; char *pp, *q; 
    
    sprintf(fname,"%s.cpu%02d",incnfile,j);
    LFile::SubHomeDir(fname, fullname);
    INFO("filename="<<fullname);
    LFile::LoadToString(fullname,buffer,0);
    
    pp=buffer;

    q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
    sscanf(q, "%d", &_CHAINLENGTH);

    q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
    sscanf(q, "%d", constrainatoms);

    INFO("constrainatoms[0]="<<constrainatoms[0]);

    /* Bug fixed by Keonwook Kang Apr/20/2010
       : Read constrainatoms before initRchain_parallel() */
    for(i=0;i<constrainatoms[0];i++)
    {
            q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
            sscanf(q, "%d", constrainatoms+i+1);
    }

    for(i=0;i<constrainatoms[0];i++)
    {
        q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
        sscanf(q, "%lf %lf %lf",&(_R[constrainatoms[i+1]].x),&(_R[constrainatoms[i+1]].y),&(_R[constrainatoms[i+1]].z));
    }

    Free(buffer);
    INFO_Printf("readRchain_parallel_toCN(%d) finished\n",j);
    return 0;
}


int MDFrame::writeRchain()
{
    int i, j;
    FILE *fp;

    fp=fopen(finalcnfile,"w");

    fprintf(fp,"%d\n",_CHAINLENGTH);
    fprintf(fp,"%d\n",constrainatoms[0]);

    for(i=0;i<constrainatoms[0];i++)
        fprintf(fp,"%d\n",constrainatoms[i+1]);

    for(j=0;j<=_CHAINLENGTH;j++)
        for(i=0;i<constrainatoms[0];i++)
            fprintf(fp,"%25.17e %25.17e %25.17e\n",_Rc[j][i].x,_Rc[j][i].y,_Rc[j][i].z);

    fclose(fp);  /* 2007/03/22 Keonwook Kang */
    DUMP("writeRchain finished");
    return 0;
}

int MDFrame::writeRchain_parallel_fromCN(int j)
{
    int i;
    char fname[200];
    FILE *fp;

    sprintf(fname,"%s.cpu%02d",finalcnfile,j);
    fp=fopen(fname,"w");

    SHtoR();

    fprintf(fp,"%d\n",_CHAINLENGTH);
    fprintf(fp,"%d\n",constrainatoms[0]);

    for(i=0;i<constrainatoms[0];i++)
        fprintf(fp,"%d\n",constrainatoms[i+1]);

    for(i=0;i<constrainatoms[0];i++)
        fprintf(fp,"%25.17e %25.17e %25.17e\n",_R[constrainatoms[i+1]].x,_R[constrainatoms[i+1]].y,_R[constrainatoms[i+1]].z);

    fclose(fp);
    INFO_Printf("writeRchain_parallel_fromCN(%d) finished\n",j);
    return 0;
}

#if 0
void MDFrame::initRchain()
{
    int ipt,i,j,n; double s; Vector3 ds;
    
    n=constrainatoms[0];    
    /* initialize path */
    for(j=0;j<=_CHAINLENGTH;j++)
    {
        s=(1.0*j)/_CHAINLENGTH;
        _nebinterp[j] = s;
        for(i=0;i<n;i++)
        {
             ipt=constrainatoms[i+1];
             ds=_SR2[ipt]-dscom12; ds-=_SR1[ipt]; ds.subint(); /* 2007/3/6 Wei Cai */
	     ds*=s; ds+=_SR1[ipt];
            _Rc[j][i]=_H*ds;
        }
    }
    INFO("Chain _Rc is initialized");
}
#endif

/* incnfile for middle point needs to be specified */
/* question about calChainForce. Why s=_nebinterp[j] and interpCN(s) is called? */
void MDFrame::initRchain()
{
    int j; double s; Vector3 ds;

    for(j=0;j<=_CHAINLENGTH;j++)
    {
        s=(1.0*j)/_CHAINLENGTH;
        _nebinterp[j] = s;
        if (_SR3 == NULL)
        {
           interpCN(s);
        } else 
        {
           interpCN_Quadratic(s);
        }
        copyCNtoRchain(j);
    }
    INFO("Chain _Rc is initialized");
}


int MDFrame::interpCN(double s)
{
  int i; 
  Vector3 ds;

  // 2007.04.20 Keonwook Kang
  if(_SR1 == NULL || _SR2 == NULL)
  {
      if(_SR1 == NULL)
          ERROR("config1 is not set!");
      else
          ERROR("config2 is not set!");
      return -1;
  }
  
#if 0 /* 2007/3/7 Wei Cai, for debug */
  if(s==0)
  { 
    INFO_Printf("interpCN s=%g special case\n",s);
    for(i=0;i<_NP;i++)
      _SR[i]=_SR1[i];
  }
  else if(s==1)
  {
    INFO_Printf("interpCN s=%g special case\n",s);
    for(i=0;i<_NP;i++)
      _SR[i]=_SR2[i];
  }
  else 
#else
  if(1)
#endif
  {
    for(i=0;i<_NP;i++)
      {
        //ds=_SR2[i]; ds-=_SR1[i]; ds.subint(); /* 2007/3/6 Wei Cai */
        ds=_SR2[i]-dscom12; ds-=_SR1[i]; ds.subint(); /* 2010/4/6 Keonwook */
        ds*=s; ds+=_SR1[i];
        _SR[i]=ds;
      }
  }
  return 0;
}

int MDFrame::interpCN_Quadratic(double x)
{
  int i; 
  Vector3 ds, ds2, ds3;

  // 2007.04.20 Keonwook Kang
  if(_SR1 == NULL || _SR2 == NULL)
  {
      if(_SR1 == NULL)
          ERROR("config1 is not set!");
      else 
          ERROR("config2 is not set!");
      return -1;
  }
  if (_SR3 == NULL) { 
      ERROR("config3 is not set!");
      return -1;
  }

#if 1
  x = 1-x;
  for(i=0;i<_NP;i++)
  {
      ds = _SR1[i]; ds *= (0.5-x)*(1-x)/(0.5);
      ds2 = _SR2[i]-dscom12; ds2 *= (-x)*(0.5-x)/(1)/(0.5);
      ds3 = _SR3[i]; ds3 *= (-x)*(1-x)/(-0.5)/(0.5);
      ds += ds2 + ds3;
      //ds.subint(); 
      _SR[i] = ds;
  }
#else //Straight path. For Debug only
  for(i=0;i<_NP;i++)
  {
      ds = _SR1[i]; ds *= (1-x);
      ds2 = _SR2[i]-dscom12; ds2 *= x;
      ds += ds2;
      _SR[i] = ds;
  }
#endif

  return 0;
}


void MDFrame::copyRchaintoCN(int j)
{
  int i, n, ipt;
  Matrix33 hinv;
  
  hinv = _H.inv();
  n=constrainatoms[0];
  for(i=0;i<n;i++)
  {
    ipt=constrainatoms[i+1];
    _SR[ipt]=hinv*_Rc[j][i];
    //_SR[ipt].subint();                    /* 2009/10/09 Keonwook Kang */
  }
}

void MDFrame::copyCNtoRchain(int j)
{
  int i, n, ipt;
  
  n=constrainatoms[0];
  for(i=0;i<n;i++)
  {
    ipt=constrainatoms[i+1];
    _Rc[j][i]= _H*_SR[ipt];
  }
}

void MDFrame::moveRchain(int jdes, int jsrc)
{
  int i, n;
  
  n=constrainatoms[0];
  for(i=0;i<n;i++)
  {
    _Rc[jdes][i]= _Rc[jsrc][i];
  }
}


/********************************************************************/
/* Weinan E's String method                                         */
/* J. of Chem. Phys. Vol.126 164103, 2007                           */
/********************************************************************/
#ifdef _STRINGMETHOD
#include "stringmethod.cpp"
#endif

void MDFrame::allocChainVars(int n, double **pEc, double **pFm, double **pFt, double **pFs, double **pdR, double **pLc, double **pTangMag2, Vector3 ***pdTang, void **pc)
{
    /* Allocate workspace variables for nebrelax 
     * Ec:       potential energy of each replica
     * Fm:       magnitude of the orthogonalized force (perpendicular to local tangent) 
     * Ft:       dot product between atomic forces and the local tangent vector
     * Fs:       magnitude of spring force (along the local tangent)
     * dR:       segment length (along path) in Angstrom dR_j = |R_j - R_{j-1}| 
     * dTang:    2-dimensional array storing the local tangent vector at each replica
     * TangMag2: square of the norm of dTang vector
     */
    int size, i; 

    *pEc=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1));
    *pFm=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1));
    *pFt=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1));
    *pFs=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1));
    *pdR=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1));
    *pLc=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1));
    *pTangMag2=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1));
    memset(*pEc,0,sizeof(double)*(_CHAINLENGTH+1));
    memset(*pFm,0,sizeof(double)*(_CHAINLENGTH+1));
    memset(*pFt,0,sizeof(double)*(_CHAINLENGTH+1));
    memset(*pFs,0,sizeof(double)*(_CHAINLENGTH+1));
    memset(*pdR,0,sizeof(double)*(_CHAINLENGTH+1));
    memset(*pLc,0,sizeof(double)*(_CHAINLENGTH+1));
    memset(*pTangMag2,0,sizeof(double)*(_CHAINLENGTH+1));

    size=sizeof(Vector3 *)*(_CHAINLENGTH+1)+sizeof(Vector3)*(_CHAINLENGTH+1)*n + 10;
    *pc=malloc(size); memset(*pc,0,size);
    *pdTang = (Vector3 **) (*pc); (*pdTang)[0]=(Vector3 *)((*pdTang)+(_CHAINLENGTH+1));
    for(i=1;i<_CHAINLENGTH+1;i++) (*pdTang)[i]=(*pdTang)[i-1]+n;

}

void MDFrame::calChainForce(int n, int relax_surround, double *Ec)
{
    int i, j, ipt;
    double s;
    /* Get atomic forces to _Fc array along the entire chain */
    for(j=0;j<=_CHAINLENGTH;j++)
    {
        /* Interpolate all atoms, replace constrained atoms from Rchain */
        s=_nebinterp[j];
        if (_SR3 == NULL)
        {
            interpCN(s);       /* interpolate all atoms for a given s */
        } else
        {
            interpCN_Quadratic(s);       /* interpolate all atoms for a given s */
        }
        copyRchaintoCN(j); /* constrained atoms at chain value    */

        if(relax_surround==1) { INFO("Relax surrounding atoms in replica "<<j<<"."); relax(); }

        clearR0(); /* enforces neighborlist reconstruction every step (important!!) */

        /* compute potential energy and forces for this replica */
        call_potential();
        Ec[j]=_EPOT;
 
        /* Put atomic forces to _Fc array */
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            _Fc[j][i]=_F[ipt]; //Force acting on the i-th constrained atom in replica j
        }
    } 
}

void MDFrame::calChainTang(int n, double *Ec, double *dR, double *TangMag2, Vector3 **dTang)
{
    int i, j;
    double dEmax, dEmin, r2;
    Vector3 dt;
    /* calculate tangential vector, dTang[j][i] */
    for(j=1;j<_CHAINLENGTH;j++) // j; index of chains
    {
        if (((Ec[j+1]-Ec[j]>0) && (Ec[j-1]-Ec[j]>0)) ||
            ((Ec[j+1]-Ec[j]<0) && (Ec[j-1]-Ec[j]<0))) // convex or concave energy profile
        {
            dEmax = std::max(fabs(Ec[j+1]-Ec[j]),fabs(Ec[j-1]-Ec[j]));
            dEmin = std::min(fabs(Ec[j+1]-Ec[j]),fabs(Ec[j-1]-Ec[j]));
            
            if (Ec[j+1]>Ec[j-1])                
                for(i=0;i<n;i++) // loop through all constrained atoms
                {
                    dt=(_Rc[j+1][i]-_Rc[j][i])*dEmax + (_Rc[j][i]-_Rc[j-1][i])*dEmin;
                    dTang[j][i]=dt;
                }
            else
                for(i=0;i<n;i++) // loop through all constrained atoms
                {
                    dt=(_Rc[j+1][i]-_Rc[j][i])*dEmin + (_Rc[j][i]-_Rc[j-1][i])*dEmax;
                    dTang[j][i]=dt;
                }
        }
        else
        {
            if (Ec[j+1]>Ec[j-1]) // uphill
                for(i=0;i<n;i++) // loop through all constrained atoms
                {
                    dt=_Rc[j+1][i]-_Rc[j][i];
                    dTang[j][i]=dt;
                }
            else // downhill
                for(i=0;i<n;i++) // n: number of constrained atoms
                {
                    dt=_Rc[j][i]-_Rc[j-1][i];
                    dTang[j][i]=dt;
                }
        }
        
        r2=0;        
        for(i=0;i<n;i++) // compute magnitude of tangential vector
            r2+=dTang[j][i].norm2();
        TangMag2[j]=r2;
    }
    /* calculate tangentail vector at the left end */
    j=0; r2=0;
    for(i=0;i<n;i++)
    {
        dt= _Rc[j+1][i]-_Rc[j][i];
        dTang[j][i]=dt;
        r2+=dTang[j][i].norm2();
    }
    TangMag2[j]=r2;
    /* calculate tangentail vector at the right end */
    j=_CHAINLENGTH; r2=0;
    for(i=0;i<n;i++)
    {
        dt= _Rc[j][i]-_Rc[j-1][i];
        dTang[j][i]= dt;
        r2+=dTang[j][i].norm2();
    }
    TangMag2[j]=r2;

    /* calculate chain segment length dR[j] = | Rc[j] - Rc[j-1] | */
    memset(dR,0,sizeof(double)*(_CHAINLENGTH+1));
    for(j=1;j<=_CHAINLENGTH;j++)
    {
        r2=0;
        for(i=0;i<n;i++)
        {
            dt=_Rc[j][i]-_Rc[j-1][i];
            r2+=dt.norm2();
        }
        dR[j] = sqrt(r2);
    }
    
}

void MDFrame::orthoChainForce(int moveleftend, int moverightend, int yesclimbimage, int EmaxDomain, int n, double *Fm, double *Ft, double *TangMag2, Vector3 **dTang)
{
    int i, j;
    double fr, fm2;
    Vector3 dt;

    /* orthogonalize forces _Fc against local tangent vector dTang */
    for(j=0;j<=_CHAINLENGTH;j++)
    {
        fr=0;            
        for(i=0;i<n;i++)
            fr+=dot(_Fc[j][i],dTang[j][i]);

        Ft[j] = fr/sqrt(TangMag2[j]);
        fr/=TangMag2[j]; // normalizing
        
        fm2=0;
        for(i=0;i<n;i++)
        {
            dt=dTang[j][i]*fr;   // tangential force vector

            if ((j==0 && moveleftend==3) || (j==_CHAINLENGTH && moverightend==3))
            { 
               /* ends move in free fall, do not orthogonalize forces */
            } else if ((j==0 && moveleftend==2) || (j==_CHAINLENGTH && moverightend==2))
            {
               /* set force along local tangent direction */
               _Fc[j][i] = dt; 
            } else {
               /* orthogonalize forces */
               _Fc[j][i] -= dt; 
            }
            fm2 += _Fc[j][i].norm2();

            if (yesclimbimage && (curstep>=equilsteps) && j==EmaxDomain)
    	    _Fc[j][i]-=dt;   // revert the force along tangent direction for climbing image
        }
        Fm[j]=sqrt(fm2); // store orthgonal force magnitude
    }

}

void MDFrame::reparamChain(int moveleftend, int moverightend, int yesclimbimage, int EmaxDomain, int n, double *plavg, double *plavg0, double *dR)
{
    int i, j, leftendDomain, rightendDomain, left_index, right_index;
    double r2, Ltot, new_position, alpha;
    Vector3 dt;

    leftendDomain = 0; rightendDomain = _CHAINLENGTH;

    if (moveleftend) 
    {
        ERROR("reparamChain: moveleftend not implemented yet!");
        /* update leftendDomain */
        /* perhaps we can allow the user to specify leftendDomain */
    }
    if (moverightend) 
    {
        ERROR("reparamChain: moverightend not implemented yet!");
        /* update rightendDomain */
        /* perhaps we can allow the user to specify leftendDomain */
    }
    if (yesclimbimage)
    {
        ERROR("reparamChain: not compatible with yesclimbimage!");
    }

    /* calculate chain segment length dR[j] = | Rc[j] - Rc[j-1] | */
    memset(dR,0,sizeof(double)*(_CHAINLENGTH+1));
    for(j=1;j<=_CHAINLENGTH;j++)
    {
        r2=0;
        for(i=0;i<n;i++)
        {
            dt=_Rc[j][i]-_Rc[j-1][i];
            r2+=dt.norm2();
        }
        dR[j] = sqrt(r2);
    }

    /* calculate the total length from leftendDomain to rightendDomain */
    Ltot = 0.0;
    for(j=leftendDomain+1;j<=rightendDomain;j++) Ltot += dR[j];

    /* find left_index, right_index, alpha, for each replica j */
    for(j=0;j<=_CHAINLENGTH;j++)
    {
        new_position = (1.0*j*Ltot)/_CHAINLENGTH;
        left_index = leftendDomain;
        while(new_position > dR[left_index+1] && left_index < _CHAINLENGTH) {
            left_index++;
            new_position -= dR[left_index];
        }
        if (left_index == _CHAINLENGTH) {
       	    left_index = _CHAINLENGTH-1; 
            right_index = _CHAINLENGTH;
            alpha = 1;
        } else { 
            right_index = left_index + 1; 
            alpha = new_position/dR[right_index];
        }
        for(i=0;i<n;i++) { 
	    _Rctmp[j][i] = _Rc[left_index][i]*(1.0-alpha)+_Rc[right_index][i]*alpha;
        }
    }
    for(j=0;j<=_CHAINLENGTH;j++)
    {
        for(i=0;i<n;i++) { 
            _Rc[j][i] = _Rctmp[j][i];
        }
    }
}

void MDFrame::nebrelax()
{/*
  * New Algorithm described in
  *
  * Improved tangent estimate in the nudged elastic band method for finding 
  * minimum energy paths and saddle points, Henkelman and Johnsson, 
  * Journal of Chemical Physics Vol.113 9978-9985, 2000
  *
  * relax _Rc by nudged elastic band method (steepest descent)
  *
  * Notes:  works for fix box only
  * See:    http://micro.stanford.edu/MDpp/entries?utf8=%E2%9C%93&search=nebrelax
  * 
  * Settings specified in nebspec
  * See:    http://micro.stanford.edu/MDpp/entries?utf8=%E2%9C%93&search=nebspec
  */
    int n, i, ipt, j;
    int relax_surround, moveleftend, moverightend, yesclimbimage, EmaxDomain;
    Vector3 ds0;
    double minFm, maxFm, Fspring, spring_const, springK, E0, Emax;
    double *Ec, *Fm, *Ft, *Fs, *dR, *Lc, *TangMag2; Vector3 **dTang; void *c;
    FILE *fp;
    
    INFO("NEB Relax");

    relax_surround= nebspec[0];
    spring_const  = (double) nebspec[1]; /* suggested spring constant, if <= 0 then automatically adjusted */ 
    moveleftend   = nebspec[2]; 
    moverightend  = nebspec[3];
    yesclimbimage = nebspec[4];

    /* if moveleftend = 0, fix the left end replica
                      = 1, move the left end replica by the spring force in the normal to the real force direction
                           (PNAS vol.104 3031-3036, 2007)
                      = 2, move the left end replica along the chain tangential direction (constrained free fall)
                      = 3, move the left end replica along the real force direction (free fall) */
    if (yesclimbimage)
        INFO_Printf("NEB: The climbing image method will be applied after %d steps.\n",equilsteps);

    /* Makes sure _SR1 and _SR2 (for configs A and B) are allocated */
    if(_SR1 == NULL || _SR2 == NULL)
    {
        if(_SR1 == NULL) ERROR("config1 is not set!");
        if(_SR2 == NULL) ERROR("config2 is not set!");
        return;
    }

    /* Makes sure _Rc (positions of constrained atoms) and _Fc (forces on therm)
     *  are allocated (as 2-dimensional arrays) */
    if(_Rc==NULL || _Fc==NULL)
    {
        AllocChain();
        caldscom12();
        initRchain();
    }
   
    n=constrainatoms[0]; 

    fp=fopen("nebeng.out","w");

    /* Allocate workspace variables for nebrelax */
    Ec = NULL; Fm = NULL; Ft = NULL; Fs = NULL; Fs = NULL; dR = NULL; TangMag2 = NULL; dTang = NULL; c = NULL;
    allocChainVars(n, &Ec, &Fm, &Ft, &Fs, &dR, &Lc, &TangMag2, &dTang, &c);
   
    /* calculate center-of-mass shift from _SR1 (state A) to _SR2 (state B) */ 
    caldscom12();
    /* compute current constrain value */
    constrain_dist=0; 
    for(i=0;i<n;i++)
    {
        ipt=constrainatoms[i+1];
        ds0=(_SR2[ipt]-dscom12)-_SR1[ipt];
        constrain_dist+=ds0.norm2();
    }
    INFO_Printf("NEB: constrain_dist=%e\n",constrain_dist);

    if(relax_surround==1)
    {/* fix constrained atoms for relaxing surrounding atoms */
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            fixed[ipt]=1;
        }
    }

    /* Store the energy of the left end for future reference. 
       This was implemented for moving left end, but is also 
       applicable when the left end is fixed. Oct/04/2010 KW */
    for(i=0;i<_NP;i++) _SR[i]=_SR1[i];          
    clearR0(); call_potential(); E0=_EPOT; /* E0 is the energy of left end (state A) */

    /* ?? Do we need to update E0 if moveleftend = 1 ?? */
               
    step0 = curstep;
    for(curstep=step0;curstep<=(step0 + totalsteps);curstep++)    
    {
        /* Get atomic forces to _Fc array along the entire chain */
        calChainForce(n, relax_surround, Ec);

        /* Find max-energy replica (i.e. the saddle point) */
        EmaxDomain = 0; Emax=E0; 
        for(j=1;j<=_CHAINLENGTH;j++)
        {
            if(Ec[j]>Emax) { Emax=Ec[j]; EmaxDomain=j; }
        }

        /* calculate tangential vector, dTang[j][i] */
        calChainTang(n, Ec, dR, TangMag2, dTang);

        /* orthogonalize forces _Fc against local tangent vector dTang */
        orthoChainForce(moveleftend, moverightend, yesclimbimage, EmaxDomain, n, Fm, Ft, TangMag2, dTang);

        /* adjust spring constant if the suggested value is 0 or negative*/
        if (spring_const>0.0)
        {
            springK = spring_const; 
            if (curstep==step0) INFO_Printf("use manually set spring constant k = %20.12e\n", springK);
        }
        else
        {
            minFm = maxFm = Fm[0]; // Min. and Max. magnitude of orthogonal force over j
            for(j=1;j<=_CHAINLENGTH;j++)
            {
                if (Fm[j]<minFm) minFm = Fm[j];
                if (Fm[j]>maxFm) maxFm = Fm[j];
            }
            springK = maxFm*_CHAINLENGTH;
            INFO_Printf("automatically chosen spring constant k = %20.12e\n", springK);
        }

        /* add tangential component of spring force (not including left and right ends) */
        for(j=1;j<_CHAINLENGTH;j++)
        {
            Fspring = springK*(dR[j+1] - dR[j]);

            Fs[j] = Fspring;
            
            /* If climbing image, skip */
            if (yesclimbimage && (curstep>=equilsteps) && j==EmaxDomain) continue;

            for(i=0;i<n;i++)
                _Fc[j][i]+=dTang[j][i]*(Fspring/sqrt(TangMag2[j]));
        }

        if(curstep%printfreq==0)  // Mar. 08 2007 Keonwook Kang 
        {
            INFO_Printf("curstep = %d\n",curstep);
            for(j=0;j<=_CHAINLENGTH;j++)
	        INFO_Printf("%20.12e %25.15e %25.15e %25.15e %25.15e %25.15e\n",
                            _nebinterp[j],Ec[j]-E0,Fm[j],Ft[j],Fs[j],dR[j]);
            fprintf(fp,"%d ",curstep);
            for(j=0;j<=_CHAINLENGTH;j++)
                fprintf(fp,"%20.12e %25.15e %25.15e %25.15e %25.15e %25.15e ",
                            _nebinterp[j],Ec[j]-E0,Fm[j],Ft[j],Fs[j],dR[j]);
            fprintf(fp,"\n");
	    fflush(fp); 
        }
        
        if(curstep<(step0 + totalsteps))
        {
            /* move chains along force */
            for(j=0;j<=_CHAINLENGTH;j++)
            {
                if ( (j==0) && (moveleftend<=0) )             continue;
                if ( (j==_CHAINLENGTH) && (moverightend<=0) ) continue;
                for(i=0;i<n;i++)
                    _Rc[j][i]+=_Fc[j][i]*_TIMESTEP;
            }
        }    
    }

    free(Ec); free(Fm); free(Ft); free(Fs);
    free(dR); free(dTang); free(TangMag2);
    fclose(fp);
}

void MDFrame::stringrelax()
{/*
  * relax _Rc by string method (steepest descent)
  *
  * Note: fix box only
  */
    int n, i, ipt, j;
    int relax_surround, redistri_freq, moveleftend, moverightend, yesclimbimage, EmaxDomain;
    Vector3 ds0; 
    double s, lavg0, lavg, E0, Emax;
    double *Ec, *Fm, *Ft, *Fs, *dR, *Lc, *TangMag2; Vector3 **dTang; void *c;
    int *left_index, *right_index;
    int rightendEnergyTooHigh, myEnergyTooHigh, manualcut, manualleftend, manualrightend, cgstepsmid,cgstepsrightend;
    double energySpike, energySlope;
    FILE *fp;
    
    INFO("String Relax");

    relax_surround= nebspec[0];
    redistri_freq = nebspec[1];
    moveleftend   = nebspec[2]; 
    moverightend  = nebspec[3];
    yesclimbimage = nebspec[4];
    /* islengthconstant=nebspec[5]; deprecated */

    if (moveleftend) { ERROR("string: moveleftend not yet implemented"); }

    if (yesclimbimage)
        INFO_Printf("string: The climbing image method will be applied after %d steps.\n",equilsteps);

    /* Makes sure _SR1 and _SR2 (for configs A and B) are allocated */
    if(_SR1 == NULL || _SR2 == NULL)
    {
        if(_SR1 == NULL) ERROR("config1 is not set!");
        if(_SR2 == NULL) ERROR("config2 is not set!");
        return;
    }

    /* Makes sure _Rc (positions of constrained atoms) and _Fc (forces on therm)
     *  are allocated (as 2-dimensional arrays) */
    if(_Rc==NULL || _Fc==NULL)
    {
        AllocChain();
        caldscom12();
        initRchain();
    }

    n = constrainatoms[0];

    fp=fopen("stringeng.out","w");
    
    /* Allocate workspace variables for stringrelax */
    Ec = NULL; Fm = NULL; Ft = NULL; Fs = NULL; Fs = NULL; dR = NULL; TangMag2 = NULL; dTang = NULL; c = NULL;
    allocChainVars(n, &Ec, &Fm, &Ft, &Fs, &dR, &Lc, &TangMag2, &dTang, &c);

    /* calculate center-of-mass shift from _SR1 (state A) to _SR2 (state B) */ 
    caldscom12();
    /* compute current constrain value */
    constrain_dist=0; s=0;
    for(i=0;i<n;i++)
    {
        ipt=constrainatoms[i+1];
        ds0=(_SR2[ipt]-dscom12)-_SR1[ipt];
        constrain_dist+=ds0.norm2();
    }
    INFO_Printf("string: constrain_dist=%e\n",constrain_dist);

    if(relax_surround==1)
    {/* fix constrained atoms for relaxing surrounding atoms */
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            fixed[ipt]=1;
        }
    }

    /* Store the energy of the left end for future reference. 
       This was implemented for moving left end, but is also 
       applicable when the left end is fixed. Oct/04/2010 KW */
    for(i=0;i<_NP;i++) _SR[i]=_SR1[i];          
    clearR0(); call_potential(); E0=_EPOT;
               
    step0 = curstep; lavg = 0;
    for(curstep=step0;curstep<=(step0 + totalsteps);curstep++)    
    {
        /* Get atomic forces to _Fc array along the entire chain */
        calChainForce(n, relax_surround, Ec);
           
        /* Find max-energy replica for climbing image method */
        EmaxDomain = 0; Emax=E0; 
        for(j=1;j<=_CHAINLENGTH;j++)
        {
            if(Ec[j]>Emax) { Emax=Ec[j]; EmaxDomain=j; }
        }

        /* calculate tangential vector, dTang[j][i] */
        calChainTang(n, Ec, dR, TangMag2, dTang);

        /* orthogonalize forces _Fc against local tangent vector dTang */
        orthoChainForce(moveleftend, moverightend, yesclimbimage, EmaxDomain, n, Fm, Ft, TangMag2, dTang);
       
        /* calculate the initial averaged chain length, lavg0 */
        if(curstep==step0)
        {
            lavg0 = 0;
            for(j=1;j<=_CHAINLENGTH;j++) lavg0 += dR[j];
            lavg0 /= _CHAINLENGTH;
        }

        if(curstep%printfreq==0)
        {
           INFO_Printf("curstep = %d\n",curstep);
           for(j=0;j<=_CHAINLENGTH;j++)
                   INFO_Printf("%8d %25.15e %25.15e %25.15e %25.15e\n",j,Ec[j]-E0,Fm[j],Ft[j],dR[j]);
           fprintf(fp,"%d ",curstep);
           for(j=0;j<=_CHAINLENGTH;j++)
                   fprintf(fp,"%8d %25.15e %25.15e %25.15e %25.15e ",j,Ec[j]-E0,Fm[j],Ft[j],dR[j]);
           fprintf(fp,"\n");
           fflush(fp);
        }

        if(curstep<(step0 + totalsteps))
        {
            /* move chains along force */
            for(j=0;j<=_CHAINLENGTH;j++)
            {
                if ( (j==0) && (moveleftend<=0) )             continue;
                if ( (j==_CHAINLENGTH) && (moverightend<=0) ) continue;
                for(i=0;i<n;i++)
                    _Rc[j][i]+=_Fc[j][i]*_TIMESTEP;
            }

            /* reparameterize path (redistribution) */
            if(curstep%redistri_freq==0) 
                 reparamChain(moveleftend, moverightend, yesclimbimage, EmaxDomain, n, &lavg, &lavg0, dR);

        }    
    }
      
    free(Ec); free(Fm); free(Ft); free(Fs);
    free(dR); free(dTang); free(TangMag2);
    fclose(fp);
}

void MDFrame::annealpath()
{/* simulated annealing of transition path by Monte Carlo */
 /* fix box only */
    int i,n,j,k, ipt, naccept, alg;
    Vector3 ds, dr, *Rold; Matrix33 hinv;
    double s, rmax, T, T0, lambda, r0, r1, r2, r02, r12, r22, rt;
    double Emax, Enewmax, Eold, Epen, Enewpen, Etot, Enewtot;
    double xi, crit, accratio, *Ec, *Lc;
    FILE *fp;
    
    INFO("annealpath");

    if(constrainatoms[0]<=0)
    {
        ERROR("constrainedatoms not set! (atoms not selected)");
        return;
    }

    rmax  = annealspec[0];
    T0    = annealspec[1];
    lambda= annealspec[2];
    alg   = (int)annealspec[3];
    
    if((rmax<=0)||(T0<=0)||(lambda<=0))
    {
        ERROR("Invalide annealspec !");
        return;
    }
    
    fp=fopen("nebeng.out","w");
    
    Ec=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1));
    Lc=(double *)malloc(sizeof(double)*(_CHAINLENGTH+1));
    Rold=(Vector3 *)malloc(sizeof(Vector3)*constrainatoms[0]);
                           
    hinv=_H.inv();
    
    /* rid center of mass shift in configB */
    dscom12.clear();
    constrain_dist=0; s=0;
    n=constrainatoms[0];
    for(i=0;i<n;i++)
    {
        ipt=constrainatoms[i+1];
        dscom12+=_SR2[ipt];
        dscom12-=_SR1[ipt];
    }
    dscom12/=n;

#if 0
    for(i=0;i<n;i++)
    {
        ipt=constrainatoms[i+1];
        _SR2[ipt]-=dscom12;
    }
#endif

    INFO("find energy along the entire chain");
    /* get potential energy along the chain */
    for(j=0;j<=_CHAINLENGTH;j++)
    {
        s=(1.0*j)/_CHAINLENGTH;
        /* surrounding atoms at interpolated value */
        for(i=0;i<_NP;i++)
        {
            ds=_SR2[i]-dscom12; ds-=_SR1[i]; ds*=s; ds+=_SR1[i];
            _SR[i]=ds;
        }
        /* constained atoms at chain value */
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            _SR[ipt]=hinv*_Rc[j][i];
        }
        potential_energyonly();
        Ec[j]=_EPOT;        
    }
    
    Emax=Ec[0]; Epen=0; Etot=0;
    for(j=1;j<=_CHAINLENGTH;j++)
    {
        /* find maximum energy along the path*/
        if(Emax<Ec[j])    Emax=Ec[j];
        /* find penalty energy */
        if(Ec[j]>Ec[j-1]) Epen+=Ec[j]-Ec[j-1];
        Etot+=Ec[j]-Ec[0];
    }
    
    INFO("begin annealing simulation");
    /* start annealing loop */
    naccept = 0;

    if(continue_curstep)
       step0 = curstep;
    else
       step0 = 0;
    for(curstep=step0;curstep<(step0 + totalsteps);curstep++)    
    //for(curstep=0;curstep<totalsteps;curstep++)
    {
        T = T0*exp(-lambda*curstep/totalsteps);
        
        /* randomly pick a state in the chain */
        j = (int) ceil( drand48()*(_CHAINLENGTH-1) );

        /* save old coordinate and energy */
        for(i=0;i<n;i++)
            Rold[i]=_Rc[j][i];
        Eold=Ec[j];
        
        /* randomly perturb the atoms in state j */
        /* replace constraint by 2-norm, rejection free algorithm */
        r02=0;
        for(i=0;i<n;i++)
        {
            dr=_Rc[j+1][i]-_Rc[j-1][i];
            r02+=dr.norm2();
        }        
        r0=sqrt(r02);
        if(r0>1e-6)
        {
            xi=(2*drand48()-1)*(rmax-r0*0.5);
            r12=rmax*rmax-(r0*0.5+fabs(xi))*(r0*0.5+fabs(xi));
            if(r12<0)
            {
                WARNING("r12<0");
                r12=0;
            }                
            r1=sqrt(r12);
            /* randomize Rc[j] */
            for(i=0;i<n;i++)
                for(k=0;k<3;k++)
                    _Rc[j][i][k]=2*drand48()-1;
            /* orthogonize with t[j]=Rc[j+1]-Rc[j-1] */
            rt=0;
            for(i=0;i<n;i++)
            {
                dr=_Rc[j+1][i]-_Rc[j-1][i];
                rt+=dot(_Rc[j][i],dr);
            }
            rt/=r02;
            for(i=0;i<n;i++)
            {
                dr=_Rc[j+1][i]-_Rc[j-1][i];
                dr*=rt;
                _Rc[j][i]-=dr;
            }
            /* normalize it to r1 */
            r22=0;
            for(i=0;i<n;i++)
                r22+=_Rc[j][i].norm2();
            r2=sqrt(r22);
            for(i=0;i<n;i++)
                _Rc[j][i]*=r1/r2;
            /* shift to center + xi along t[j] */
            for(i=0;i<n;i++)
            {
                dr=_Rc[j+1][i]-_Rc[j-1][i];
                dr/=r0;
                dr*=(r0*0.5+xi);
                dr+=_Rc[j-1][i];
                _Rc[j][i]+=dr;
            }
        }
        else
        {/* r0 = 0 */
            /* randomize Rc[j] */
            for(i=0;i<n;i++)
                for(k=0;k<3;k++)
                    _Rc[j][i][k]=2*drand48()-1;
            /* normalize it to rmax */
            r22=0;
            for(i=0;i<n;i++)
                r22+=_Rc[j][i].norm2();
            r2=sqrt(r22);
            for(i=0;i<n;i++)
                _Rc[j][i]*=(rmax-2e-6)/r2;
            /* shift to center */
            for(i=0;i<n;i++)
            {
                _Rc[j][i]+=(_Rc[j-1][i]+_Rc[j+1][i])*0.5;
            }
        }
        /* old version: using infinite-norm */
        /*
        for(i=0;i<n;i++)
            for(k=0;k<3;k++)
            {
                rold=_Rc[j][i][k];
                rnew=rold+(drand48()-0.5)*_TIMESTEP;
                if((fabs(rnew-_Rc[j-1][i][k])<rmax)&&(fabs(rnew-_Rc[j+1][i][k])<rmax))
                    _Rc[j][i][k]=rnew;
            }
        */
        
        s=(1.0*j)/_CHAINLENGTH;
        /* surrounding atoms at interpolated value */
        for(i=0;i<_NP;i++)
        {
            ds=_SR2[i]-dscom12; ds-=_SR1[i]; ds*=s; ds+=_SR1[i];
            _SR[i]=ds;
        }
        /* constained atoms at chain value */
        for(i=0;i<n;i++)
        {
            ipt=constrainatoms[i+1];
            _SR[ipt]=hinv*_Rc[j][i];
        }
        potential_energyonly();
        Ec[j]=_EPOT;

        Enewmax=Ec[0]; Enewpen=0; Enewtot=0;
        for(k=1;k<=_CHAINLENGTH;k++)
        {
            /* find maximum energy along the path*/
            if(Enewmax<Ec[k]) Enewmax=Ec[k];
            /* find penalty energy */
            if(Ec[k]>Ec[k-1]) Enewpen+=Ec[k]-Ec[k-1];
            Enewtot+=Ec[k]-Ec[0];            
        }

        switch(alg) {
        case(0): crit=exp((Emax-Enewmax)/T); break;
        case(1): crit=exp((Epen-Enewpen)/T); break;
        case(2): crit=exp((Etot-Enewtot)/T); break;
        default: crit=0; ERROR("annealpath: MC alg not set!"); break;
        }
        xi = drand48();
        if(xi<crit)
        {/* accept the move */
            Emax=Enewmax;
            Epen=Enewpen;
            Etot=Enewtot;
            naccept ++;
        }
        else
        {/* reject the move */
            for(i=0;i<n;i++)
                _Rc[j][i]=Rold[i];
            Ec[j]=Eold;
        }
        accratio = 1.0*naccept/(curstep+1);
            
        if(curstep%printfreq==0)
        {
           INFO_Printf("curstep = %d\n",curstep);
           for(k=0;k<=_CHAINLENGTH;k++)
                   INFO_Printf("%e  %25.15e\n",1.0*k/_CHAINLENGTH,Ec[k]-Ec[0]);
           INFO_Printf("alg=%d Emax=%e Epen=%e Etot=%e T=%e j=%d acceptance ratio = %d/%d = %e\n",
                       alg,Emax-Ec[0],Epen,Etot,T,j,naccept,curstep+1,accratio);
           fprintf(fp,"%d ",curstep);
           for(k=0;k<=_CHAINLENGTH;k++)
                   fprintf(fp," %25.15e ",Ec[k]-Ec[0]);
           fprintf(fp,"\n");
        }
    }

    free(Ec);
    free(Lc);
    free(Rold);
    fclose(fp);
}

void MDFrame::cutpath()
{/* relax _Rc by nudged elastic band method (steepest descent) */
 /* fix box only */
    int i,n,j,n0,n1,step,jnew;
    
    INFO("cutpath");

    if(constrainatoms[0]<=0)
    {
        ERROR("constrainedatoms not set! (atoms not selected)");
        return;
    }

    step  = (int) annealspec[0];
    n0   = (int) annealspec[1];
    n1   = (int) annealspec[2];

    if(n0<=0) n0=1;
    if(n1>=_CHAINLENGTH) n1=_CHAINLENGTH-1;
    n=constrainatoms[0];
    jnew=1;
    for(j=n0;j<=n1;j+=step)
    {
        /* constained atoms at chain value */
        for(i=0;i<n;i++)
            _Rc[jnew][i]=_Rc[j][i];
        jnew++;
    }
    for(i=0;i<n;i++)
        _Rc[jnew][i]=_Rc[_CHAINLENGTH][i];
    _CHAINLENGTH=jnew;
}

int MDFrame::statedistance()
{/* compute the distance between two states _SR1 and _SR (current) */
    int i, n, j, k, ipt;
    Vector3 ds, dr;
    double r, r2, rinf, dx;

    n=constrainatoms[0];
    if(n<=0)
    {/* distance involving all atoms */
        INFO_Printf("distance involving all %d atoms\n",_NP);
        r2=0; rinf=0;
        for(i=0;i<_NP;i++)
        {
            ds=_SR[i]-_SR1[i];
            dr=_H*ds;
            r2+=dr.norm2();
            for(k=0;k<3;k++)
            {
                dx=fabs(dr[k]);
                if(dx>rinf) rinf=dx;
            }                
        }
        r=sqrt(r2);
    }
    else
    {/* distance only involving constrainatoms */
        INFO_Printf("distance involving %d constrain atoms\n",n);
        r2=0; rinf=0;
        for(j=1;j<=n;j++)
        {
            ipt=constrainatoms[j];
            ds=_SR[ipt]-_SR1[ipt];
            dr=_H*ds;
            r2+=dr.norm2();
            for(k=0;k<3;k++)
            {
                dx=fabs(dr[k]);
                if(dx>rinf) rinf=dx;
            }            
        }
        r=sqrt(r2);
    }
    INFO_Printf("2-norm r = %20.12e  infinite-norm r = %20.12e\n",r,rinf);
    return 0;
}



/********************************************************************/
/* Configuration manipulations */
/********************************************************************/

/* basis data for different crystal structures */
#include "lattice.h"

void MDFrame::makecrystal()
{   /* create perfect crystal structure */
  UnitCell my_unitcell,my_unitcell2;
    Matrix33 myorient, lattice;
    Vector3 tmpv;
    int i,j,k,m,n;

    INFO("makecrystal");
    myorient.set(latticesize[0][0],latticesize[0][1],latticesize[0][2],
                 latticesize[1][0],latticesize[1][1],latticesize[1][2],
                 latticesize[2][0],latticesize[2][1],latticesize[2][2]);

    /* http://cst-www.nrl.navy.mil/lattice/prototype.html */
    /*                          Common Name          Strukturbericht Designation*/
    if(strcmp(crystalstructure,"simple-cubic")==0 || strcmp(crystalstructure,"Ah")==0)
    {
        latticeconst[1]=latticeconst[2]=latticeconst[0];
        my_unitcell.set(1,sc_basis);   /* simple-cubic */
        nspecies=1;
    }
    else if(strcmp(crystalstructure,"body-centered-cubic")==0 || strcmp(crystalstructure,"A2")==0)
    {
        latticeconst[1]=latticeconst[2]=latticeconst[0];
        my_unitcell.set(2,bcc_basis); /* body-centered-cubic */
        nspecies=1;
    }
    else if(strcmp(crystalstructure,"face-centered-cubic")==0 || strcmp(crystalstructure,"A1")==0)
    {
        latticeconst[1]=latticeconst[2]=latticeconst[0];
        my_unitcell.set(4,fcc_basis); /* face-centered-cubic */
        nspecies=1;
    }
    else if(strcmp(crystalstructure,"NaCl")==0 || strcmp(crystalstructure,"B1")==0) 
    {
        latticeconst[1]=latticeconst[2]=latticeconst[0];
        my_unitcell.set(8,nacl_basis,nacl_species); /* NaCl (B1) */
        nspecies=2;
    }
    else if(strcmp(crystalstructure,"CsCl")==0 || strcmp(crystalstructure,"B2")==0)
    {
        latticeconst[1]=latticeconst[2]=latticeconst[0];
        my_unitcell.set(2,bcc_basis,cscl_species); /* NaCl (B1) */
        nspecies=2;
    }
    else if(strcmp(crystalstructure,"AuCu3")==0 || strcmp(crystalstructure,"L1_2")==0)
    {
        latticeconst[1]=latticeconst[2]=latticeconst[0];
        my_unitcell.set(4,fcc_basis,l12_species);   /* L1_2 */
        nspecies=2;
    }
    else if(strcmp(crystalstructure,"AuCu")==0 || strcmp(crystalstructure,"L1_0")==0)
    {
        latticeconst[1]=latticeconst[2]=latticeconst[0];
        my_unitcell.set(4,fcc_basis,l10_species);   /* L1_0 */
        nspecies=2;
    }
    else if(strcmp(crystalstructure,"diamond-cubic")==0 || strcmp(crystalstructure,"A4")==0)
    {
        latticeconst[1]=latticeconst[2]=latticeconst[0];
        my_unitcell.set(8,dc_basis);   /* diamond-cubic */
        nspecies=1;
    }
    else if(strcmp(crystalstructure,"beta-tin")==0 || strcmp(crystalstructure,"A5")==0)
    {
        if(latticeconst[1]<=0) latticeconst[1]=latticeconst[0];
        if(latticeconst[2]<=0) latticeconst[2]=latticeconst[0];
        //latticeconst[2]=latticeconst[0]*sqrt(2); /* when c=a*sqrt(2) is equivalent to diamond-cubic */
        //latticeconst[2]=latticeconst[0]*0.6;     /* when c=a*0.6 for beta-tin phase of Silicon */        
        my_unitcell.set(4,betatin_basis);   /* beta-tin */
        nspecies=1;
    }
    else if(strcmp(crystalstructure,"beta-tin-2")==0)
    {
        if(latticeconst[1]<=0) latticeconst[1]=latticeconst[0];
        if(latticeconst[2]<=0) latticeconst[2]=latticeconst[0];
        my_unitcell.set(4,betatin_basis,betatin_2_species);   /* beta-tin-2 */
        nspecies=2;
    }
    else if(strcmp(crystalstructure,"graphite")==0 || strcmp(crystalstructure,"A9")==0)
    {
        latticeconst[1]=latticeconst[0]*M_SQRT3; 
        latticeconst[2]=latticeconst[0]*2.7257;  /* a = 2.461 A, c = 6.708 A */
        my_unitcell.set(8,graphite_basis);       /* graphite */
        nspecies=1;
    }
    else if(strcmp(crystalstructure,"graphene")==0)
    {
        latticeconst[1]=latticeconst[0]*M_SQRT3; 
        latticeconst[2]=latticeconst[0]*2.7257;  /* a = 2.461 A, c = 6.708 A */
        my_unitcell.set(4,graphene_basis);       /* graphene */
        nspecies=1;
    }
    else if(strcmp(crystalstructure,"zinc-blende")==0 || strcmp(crystalstructure,"B3")==0)
    {
        latticeconst[1]=latticeconst[2]=latticeconst[0];
        my_unitcell.set(8,dc_basis,zb_species);   /* zinc-blende or Spharelite*/
        nspecies=2;
    }
    else if(strcmp(crystalstructure,"wurtzite")==0 || strcmp(crystalstructure,"B4")==0)
    {/* ideal wurtzite structure u = 3/8 
        (see cst-www.nrl.navy.mil/lattice/struk/b4.html) Keonwook Kang Feb 12 2008 */
        latticeconst[1]=latticeconst[0]*M_SQRT3;
        my_unitcell.set(8,wurtzite_basis,wurtzite_species);   /* Wurtzite */
        nspecies=2;
    }
    else if(strcmp(crystalstructure,"CaF2")==0 || strcmp(crystalstructure,"C1")==0)
    {
        latticeconst[1]=latticeconst[2]=latticeconst[0];
        my_unitcell.set(12,caf2_basis,caf2_species);   /* Perovskite */
        nspecies=2;
    }
    else if(strcmp(crystalstructure,"hexagonal-ortho")==0)
    {
        latticeconst[1]=latticeconst[0]*M_SQRT3;
        my_unitcell.set(4,hex_ortho_basis);   /* hexagonal-ortho (hcp) */
        nspecies=1;
    }
    else if(strcmp(crystalstructure,"rhombohedral-Si-Ge")==0)
    {
        latticeconst[1]=latticeconst[0]*M_SQRT3;
        latticeconst[2]=latticeconst[0]*4.89897948556636;
        my_unitcell.set(24,rhsige_basis,rhsige_species);   /* rhomohedral Si_Ge */
        nspecies=2;
    }
    else if(strcmp(crystalstructure,"bc8")==0)
    {
        latticeconst[1]=latticeconst[2]=latticeconst[0];
        my_unitcell.set(16,bc8_basis);  /* Si BC-8 */
        nspecies=1;
    }
    else if(strcmp(crystalstructure,"alpha-quartz")==0)
    {
      latticeconst[1]=latticeconst[0]*M_SQRT3;
      my_unitcell.set(18,alpha_quartz_basis,alpha_quartz_species); /* alpha_quartz orthogonal unit cell */
      nspecies = 2;
    }
    else if(strcmp(crystalstructure,"beta-quartz")==0)
    {
      latticeconst[1]=latticeconst[0]*M_SQRT3;
      my_unitcell.set(18,beta_quartz_basis,beta_quartz_species);  /* beta_quartz orthogonal unit cell */
      nspecies = 2;
    }
    else if(strcmp(crystalstructure,"tridymite")==0)
    {
      latticeconst[1]=latticeconst[0]*M_SQRT3;
      my_unitcell.set(24,tridymite_basis,tridymite_species);  /* Tridymite orthogonal unit cell */
      nspecies = 2;
    }
    else if(strcmp(crystalstructure,"cristobalite")==0)
    {
      latticeconst[1]=latticeconst[0];
      my_unitcell.set(12,cristobalite_basis,cristobalite_species);  /* Tridymite orthogonal unit cell */
      nspecies = 2;
    }
    else if(strcmp(crystalstructure,"ice-Ih-O")==0)
    {
        latticeconst[1]=latticeconst[0]*M_SQRT3;
        my_unitcell.set(8,ice_Ih_O_basis);   /* ice-Ih-Oxygen */
        nspecies=1;
    }
    else if(strcmp(crystalstructure,"ice-Ih")==0)
    {
        latticeconst[1]=latticeconst[0]*M_SQRT3;
        my_unitcell.set(24,ice_Ih_basis,ice_Ih_species);   /* ice-Ih */
        nspecies=2;
    }
    else if(strcmp(crystalstructure,"alpha-Ga")==0 || strcmp(crystalstructure,"A11")==0)
    {
        if(latticeconst[1]==0) latticeconst[1]=latticeconst[0];
        if(latticeconst[2]==0) latticeconst[2]=latticeconst[0];
        my_unitcell.set(8,A11_basis);  /* A11 (orthorhombic) */
	                               /* Ga: a=4.52A b=7.66A c=4.53A) */
        nspecies=1;
    }
    else
    {
        ERROR("unknown lattice structure :"<<crystalstructure
              <<" for makecn");
    }
    my_unitcell=my_unitcell*myorient.tran();

    INFO("my_unitcell\n"<<my_unitcell);
    n=0;
    _NP=(int)(latticesize[0][3]*latticesize[1][3]
              *latticesize[2][3]*my_unitcell.n);

    //printf("%f %f %f %d\n",latticesize[0][3],latticesize[1][3],
    //	   latticesize[2][3],my_unitcell.n);

    INFO("makecn: _NP="<<_NP);
    
    /* Multi process function */
    Alloc();
    
    for(i=0;i<(int)latticesize[0][3];i++)
        for(j=0;j<(int)latticesize[1][3];j++)
            for(k=0;k<(int)latticesize[2][3];k++)
            {
                for(m=0;m<my_unitcell.n;m++)
                {
                    tmpv.set(i,j,k);
                    _SR[n+m]=my_unitcell.basis[m]+tmpv;
                    /* Wei Cai, 8/9/2008, 8/20/2008 */
                    if(_SAVEMEMORY<8)
                            species[n+m]=my_unitcell.species[m];
                }
                n+=my_unitcell.n;
            }
    for(i=0;i<_NP;i++)
    {
        _SR[i].x/=latticesize[0][3];
        _SR[i].y/=latticesize[1][3];
        _SR[i].z/=latticesize[2][3];
    }
    for(i=0;i<_NP;i++)
    {
        _SR[i]-=0.5;
        fixed[i]=0;
        
        /* Wei Cai, 8/9/2008 */
        if(_SAVEMEMORY<7)
                _VSR[i].clear();
    }

    lattice.set(latticeconst[0],0.,0., 0.,latticeconst[1],0., 0.,0.,latticeconst[2]);
    _H=lattice*myorient.tran();
    lattice.set(latticesize[0][3],0.,0., 0.,latticesize[1][3],0., 0.,0.,latticesize[2][3]);
    _H=_H*lattice;
    
    _H=_H.reorient();

    /* Wei Cai, 8/9/2008 */
    if(_SAVEMEMORY<9)
    {
        for(i=0;i<_NP;i++)
            _R[i] = _H*_SR[i];
    }
}


void MDFrame::makecut()
{ /* input = [ ind nx ny nz x0 y0 z0 ]
   * fix all atoms that are on one side of the plane
   */
    int i, ind, nx, ny, nz;
    double x0, y0, z0, rdotn;
    Vector3 a, b, c, r, n;
    Matrix33 M;

    ind = (int) input[0];
    nx  = (int) input[1];
    ny  = (int) input[2];
    nz  = (int) input[3];
    x0 = input[4];
    y0 = input[5];
    z0 = input[6];

    a.set(latticesize[0][0],latticesize[0][1],latticesize[0][2]);
    b.set(latticesize[1][0],latticesize[1][1],latticesize[1][2]);
    c.set(latticesize[2][0],latticesize[2][1],latticesize[2][2]);

    n.set(nx,ny,nz);
    
    a/=a.norm();
    b/=b.norm();
    c/=c.norm();

    M.setcol(a,b,c);

    SHtoR();
    for(i=0;i<_NP;i++)
    {
        r = M*_R[i]; /* atom position in lattice frame */
        r.x -= x0; r.y -= y0; r.z -= z0;
        rdotn = dot(r,n);

        if(rdotn*ind>0) fixed[i]=1; /* mark atoms as fixed */
    }    
}

void MDFrame::cutbonds()
{ /* identify all pairs of atoms that are on two sides of a cut-plane
     register them in nl_skip_pairs array
     input parameters similar to that of makedipole
     input = [ enable zind  yind  x0 y0 y1 z0 z1 ]
  */
    int n, enable, xind, yind, zind, ipt, jpt, j, k, A_inside, B_inside;
    Vector3 sr0, sr1, sri, srj, dsij;
    
    /* read in parameters */
    enable=(int)input[0];
    zind=(int)input[1]; /* 1x/2y/3z */
    yind=(int)input[2]; /* direction contained in the cut-plane */
    xind=zind+yind;if(xind==5)xind=1;if(xind==4)xind=2;
    xind--;yind--;zind--;

    if (enable == 1) 
    {
        sr0[xind]=input[3];sr0[yind]=input[4];sr0[zind]=input[6];
        sr1[xind]=input[3];sr1[yind]=input[5];sr1[zind]=input[7];

        for(ipt=0;ipt<_NP;ipt++)
        {
            for(j=0;j<nn[ipt];j++)
            {
                jpt=nindex[ipt][j];
                if(ipt>=jpt) continue;

                /* handle bonds across PBC */
                sri = _SR[ipt]; srj = _SR[jpt]; dsij = srj - sri; dsij.subint();
                /* ipt, jpt on two sides of cut plane */
                if( (_SR[ipt][xind]-sr0[xind])*(_SR[ipt][xind]+dsij[xind]-sr0[xind]) < 0 )
                {
                    if ((sr1[zind]-sr0[zind])==0)
                    { /* ipt, jpt within the strip area inside the cut plane */
                        if ( (_SR[ipt][yind]>=sr0[yind]) && (_SR[ipt][yind]<sr1[yind])
                             && (_SR[jpt][yind]>=sr0[yind]) && (_SR[jpt][yind]<sr1[yind]) )
                        {
                            if (nl_skip_pairs[0]*2+2>=MAXNLSKIP)
                            {
                                WARNING("MAXNLSKIP ("<<MAXNLSKIP<<") exceeded.");
                                break;
                            }
                            nl_skip_pairs[2*nl_skip_pairs[0]+1] = ipt;
                            nl_skip_pairs[2*nl_skip_pairs[0]+2] = jpt;
                            nl_skip_pairs[0] ++;
                        }
                    }
                    else
                    {
                        if ( (_SR[ipt][yind]>=sr0[yind]) && (_SR[ipt][yind]<sr1[yind])
                             && (_SR[jpt][yind]>=sr0[yind]) && (_SR[jpt][yind]<sr1[yind]) 
                             && (_SR[ipt][zind]>=sr0[zind]) && (_SR[ipt][zind]<sr1[zind]) 
                             && (_SR[jpt][zind]>=sr0[zind]) && (_SR[jpt][zind]<sr1[zind]) )
                        {
                            if (nl_skip_pairs[0]*2+2>=MAXNLSKIP)
                            {
                                WARNING("MAXNLSKIP ("<<MAXNLSKIP<<") exceeded.");
                                break;
                            }
                            nl_skip_pairs[2*nl_skip_pairs[0]+1] = ipt;
                            nl_skip_pairs[2*nl_skip_pairs[0]+2] = jpt;
                            nl_skip_pairs[0] ++;
                        }

                    }
                }
            }
        }
    }
    else if (enable > 1)
    {// multiple cut-windows on the same plane.
        n = enable; // number of rectangular cut-windows (or cut-blocks)
        class Vector3 *sr0_cw, *sr1_cw;

        sr0_cw = 0; sr1_cw = 0;
        Realloc(sr0_cw, Vector3, n);
        Realloc(sr1_cw, Vector3, n);        
        for(k=0;k<n;k++) { sr0_cw[k].clear(); sr1_cw[k].clear(); }
        for(k=0;k<n;k++)
        {
            sr0_cw[k][xind]=input[3]; sr0_cw[k][yind]=input[4+4*k]; sr0_cw[k][zind]=input[6+4*k];
            sr1_cw[k][xind]=input[3]; sr1_cw[k][yind]=input[5+4*k]; sr1_cw[k][zind]=input[7+4*k];
            //INFO("sr0["<<k<<"] = "<<sr0_cw[k]); INFO("sr1["<<k<<"] = "<<sr1_cw[k]);
        }
     
        for(ipt=0;ipt<_NP;ipt++)
        {
            for(j=0;j<nn[ipt];j++)
            {
                jpt=nindex[ipt][j];
                if(ipt>=jpt) continue;

                A_inside = 0; B_inside = 0;
                /* handle bonds across PBC */
                sri = _SR[ipt]; srj = _SR[jpt]; dsij = srj - sri; dsij.subint();
                /* ipt, jpt on two sides of cut plane */
                if( (_SR[ipt][xind]-sr0[xind])*(_SR[ipt][xind]+dsij[xind]-sr0[xind]) < 0 )
                {
                    for(k=0;k<n;k++)
                    {
                        if ( (_SR[ipt][yind]>=sr0_cw[k][yind]) && (_SR[ipt][yind]<sr1_cw[k][yind])
                             && (_SR[ipt][zind]>=sr0_cw[k][zind]) && (_SR[ipt][zind]<sr1_cw[k][zind]))
                            A_inside = 1;
                            
                        if ( (_SR[jpt][yind]>=sr0_cw[k][yind]) && (_SR[jpt][yind]<sr1_cw[k][yind])  
                             && (_SR[jpt][zind]>=sr0_cw[k][zind]) && (_SR[jpt][zind]<sr1_cw[k][zind]))
                            B_inside = 1;

                    }

                    if (A_inside && B_inside)
                    {
                        if (nl_skip_pairs[0]*2+2>=MAXNLSKIP)
                        {
                            WARNING("MAXNLSKIP ("<<MAXNLSKIP<<") exceeded.");
                            break;
                        }
                        nl_skip_pairs[2*nl_skip_pairs[0]+1] = ipt;
                        nl_skip_pairs[2*nl_skip_pairs[0]+2] = jpt;
                        nl_skip_pairs[0] ++;
                    }
                }
            }
        }

        free(sr0_cw); free(sr1_cw);
    }
    /* print out resulting nl_skip_pairs array */
    INFO_Printf("nl_skip_pairs = [ %d", nl_skip_pairs[0]);
    for(j=0;j<nl_skip_pairs[0];j++)
        INFO_Printf(" %d %d   ", nl_skip_pairs[2*j+1], nl_skip_pairs[2*j+2]);
    INFO_Printf(" ]\n");
}

void MDFrame::cutbonds_by_ellipse()
{ /* identify all pairs of atoms that are on two sides of an elliptical cut-plane
     register them in nl_skip_pairs array
     input parameters similar to that of makedipole
     input = [ enable zind  yind  x0 y0 y1 z0 z1 ]
  */
    int enable, xind, yind, zind, ipt, jpt, j;
    Vector3 sr0, sr1, sri, srj, dsij;
    double yc, zc, a, b;
    
    /* read in parameters */
    enable=(int)input[0];
    zind=(int)input[1]; /* 1x/2y/3z */
    yind=(int)input[2]; /* direction contained in the cut-plane */
    xind=zind+yind;if(xind==5)xind=1;if(xind==4)xind=2;
    xind--;yind--;zind--;

    sr0[xind]=input[3];sr0[yind]=input[4];sr0[zind]=input[6];
    sr1[xind]=input[3];sr1[yind]=input[5];sr1[zind]=input[7];

    if (enable == 1) 
    { /* check whether bonds are within the ellipitical area in the cut plane */
        for(ipt=0;ipt<_NP;ipt++)
        {
            for(j=0;j<nn[ipt];j++)
            {
                jpt=nindex[ipt][j];
                if(ipt>=jpt) continue;

                /* handle bonds across PBC */
                sri = _SR[ipt]; srj = _SR[jpt]; dsij = srj - sri; dsij.subint();
                /* ipt, jpt on two sides of cut plane */
                if( (_SR[ipt][xind]-sr0[xind])*(_SR[ipt][xind]+dsij[xind]-sr0[xind]) < 0 )
                {
                    yc = (sr1[yind] + sr0[yind]) / 2;  a = (sr1[yind] - sr0[yind]) / 2;
                    zc = (sr1[zind] + sr0[zind]) / 2;  b = (sr1[zind] - sr0[zind]) / 2;
                    if ( ( SQR(_SR[ipt][yind]-yc)/SQR(a) + SQR(_SR[ipt][zind]-zc)/SQR(b) < 1 ) &&
                         ( SQR(_SR[jpt][yind]-yc)/SQR(a) + SQR(_SR[jpt][zind]-zc)/SQR(b) < 1 ) )
                    {
                        if (nl_skip_pairs[0]*2+2>=MAXNLSKIP)
                        {
                            WARNING("MAXNLSKIP ("<<MAXNLSKIP<<") exceeded.");
                            break;
                        }
                        nl_skip_pairs[2*nl_skip_pairs[0]+1] = ipt;
                        nl_skip_pairs[2*nl_skip_pairs[0]+2] = jpt;
                        nl_skip_pairs[0] ++;
                    }

                }
            }
        }
    }
    /* print out resulting nl_skip_pairs array */
    INFO_Printf("nl_skip_pairs = [ %d", nl_skip_pairs[0]);
    for(j=0;j<nl_skip_pairs[0];j++)
        INFO_Printf(" %d %d   ", nl_skip_pairs[2*j+1], nl_skip_pairs[2*j+2]);
    INFO_Printf(" ]\n");
}

void MDFrame::makedipole()
{
    /* Create Dislocation dipole by using approximated
       isotropic linear elasticity solution under
       Periodic Boundary Condition with image summations */
    int i,j,k,minx,maxx,miny,maxy;
    int xind,yind,zind;
    double theta[50], pre, frac, x, y, r21, r22, bx, by, nu, nua, nub;
    class Vector3 sr0,sr1,r0,r1,ds,dr,sb,b,du,st,due;
    class Vector3 sra,srb,src,ra,rb,rc,dua,dub,duc,ux,uy;
    class Vector3 px, py, pz, pa, pb, pc;
    double pxpa, pxpb, pypa, pypb;
    int removenum,insertnum;
    double dx, dy, dxc, dyc;
    int store, gcut;  
    
    if(input[0]==0)
    {
        ERROR("makedipole(): no dipole geometry set");
        return;
    }
    store = (int) input[17];
    gcut = (int) input[18];
    
    /* read in parameters */
    zind=(int)input[0]; /* 1x/2y/3z */
    yind=(int)input[1];
    xind=zind+yind;if(xind==5)xind=1;if(xind==4)xind=2;
    xind--;yind--;zind--;
    sb.set(input[2],input[3],input[4]);
    b=_H*sb;
    if(fabs(b.x)<1e-8) b.x=0;
    if(fabs(b.y)<1e-8) b.y=0;
    if(fabs(b.z)<1e-8) b.z=0;
    INFO("sb="<<sb<<"\nb="<<b);
    if(!gcut)
    {
        sr0[xind]=input[5];sr0[yind]=input[6];sr0[zind]=0;
        sr1[xind]=input[5];sr1[yind]=input[7];sr1[zind]=0;
    }
    else
    {
        sr0[xind]=input[5];sr0[yind]=input[6];sr0[zind]=0;
        sr1[xind]=input[7];sr1[yind]=input[8];sr1[zind]=0;
        INFO("sr0="<<sr0);
        INFO("sr1="<<sr1);
    }
    r0=_H*sr0;r1=_H*sr1;
    pz.set(_H[0][zind],_H[1][zind],_H[2][zind]);
    pz/=pz.norm();
    py=r1-r0;
    py.orth(pz); py/=py.norm();
    px=cross(py,pz);
    INFO("px="<<px);
    INFO("py="<<py);
    INFO("pz="<<pz);
    pa.set(_H[0][xind],_H[1][xind],_H[2][xind]); //pa/=pa.norm();
    pb.set(_H[0][yind],_H[1][yind],_H[2][yind]); //pb/=pb.norm();
    pc.set(_H[0][zind],_H[1][zind],_H[2][zind]); //pc/=pc.norm();
    pxpa=px*pa; pxpb=px*pb; pypb=py*pb; pypa=py*pa;
    bx=b*px; by=b*py;
    INFO("bx="<<bx<<"  by="<<by);
    if(fabs(pxpa)<1e-8) pxpa=0;
    if(fabs(pxpb)<1e-8) pxpb=0;
    if(fabs(pypb)<1e-8) pypb=0;
    if(fabs(pypa)<1e-8) pypa=0;
    if(pxpa<0) { pxpa=-pxpa; pxpb=-pxpb; px*=-1; INFO("redefine x");}
    INFO_Printf("pxpa=%e pxpb=%e pypb=%e pypa=%e\n",pxpa,pxpb,pypb,pypa);
    
    /* reference (ghost) points */
    sra[xind]=-0.5;sra[yind]=-0.5;sra[zind]=0;
    srb[xind]= 0.5;srb[yind]=-0.5;srb[zind]=0;
    src[xind]=-0.5;src[yind]= 0.5;src[zind]=0;
    ra=_H*sra;rb=_H*srb;rc=_H*src;
    nu=input[8]; nua=1.0/4/(1-nu)/2/M_PI; nub=(1.0-2*nu)*nua;
    minx=(int)input[9]; maxx=(int)input[10];
    miny=(int)input[11];maxy=(int)input[12];
#define getdu_sav(SR) protect(\
                ds[xind]=j; ds[yind]=k; ds[zind]=0;\
                ds+=(SR); ds-=sr0; \
                x=ds[xind]*pxpa+ds[yind]*pxpb;\
                y=ds[xind]*pypa+ds[yind]*pypb; r21=x*x+y*y;\
                theta[0]=atan2(-x,y);\
                ds[xind]=j; ds[yind]=k; ds[zind]=0;\
                ds+=(SR); ds-=sr1; \
                x=ds[xind]*pxpa+ds[yind]*pxpb;\
                y=ds[xind]*pypa+ds[yind]*pypb; r22=x*x+y*y;\
                theta[1]=atan2(-x,y);\
                pre=0.5/M_PI*(theta[1]-theta[0]);\
                du=b*pre; if((fabs(bx)>1e-8)||(fabs(by)>1e-8)){due.clear();\
                due-=px*(-bx*sin(theta[0]*2)*nua+by*(-log(r21)*nub+cos(theta[0]*2)*nua));\
                due-=py*( by*sin(theta[0]*2)*nua-bx*(-log(r21)*nub-cos(theta[0]*2)*nua));\
                due+=px*(-bx*sin(theta[1]*2)*nua+by*(-log(r22)*nub+cos(theta[1]*2)*nua));\
                due+=py*( by*sin(theta[1]*2)*nua-bx*(-log(r22)*nub-cos(theta[1]*2)*nua));\
                /*INFO("due="<<due);*/ du+=due;} \
                )
#define getdu_s(SR) protect(\
                ds[xind]=j; ds[yind]=k; ds[zind]=0;\
                ds+=SR; ds-=sr0; \
                theta[0]=atan2(-ds[xind],ds[yind]);\
                ds[xind]=j; ds[yind]=k; ds[zind]=0;\
                ds+=SR; ds-=sr1; \
                theta[1]=atan2(-ds[xind],ds[yind]);\
                pre=0.5/M_PI*(theta[1]-theta[0]);\
                du=b*pre;)
#define getdu getdu_sav
    
    if(storedr==NULL)
    {
        INFO("allocate storedr");
        Realloc(storedr,Vector3,_NP*allocmultiple);
        bindvar("storedr",storedr,DOUBLE);
        /*storedr = (Vector3 *)malloc(sizeof(Vector3)*_NP*allocmultiple);*/
        SHtoR();
        for(i=0;i<_NP;i++) storedr[i].clear();
        /*for(i=0;i<_NP;i++) storedr[i]=_R[i];*/
    }
    
    /* remove/insert atoms */
    removenum=insertnum=0;
    if(fabs(b*px)>1e-8)
    {
        if((b*px)>0)
        {
            INFO("Need to remove atoms");
            removenum=0;
            dr=r1-r0; dyc=dr*py; dxc=b*px;
            for(i=0;i<_NP;i++)
            {
                /* limit dislocation length (cut depth) */
                if((int)input[14])
                    if((_SR[i][zind]<input[15])||(_SR[i][zind]>=input[16]))
                        continue;
                ds=_SR[i]; ds-=sr0; ds.subint();
                dr=_H*ds; dx=dr*px; dy=dr*py;
                if((dx>=-dxc*0.5)&&(dx<dxc*0.5)&&(dy>=0)&&(dy<dyc))
                {
                    INFO_Printf("need remove atom %d\n",i);
                    fixed[i] = -1; /* flag for removal */
                    removenum++;
                }
            }
            INFO("need removenum="<<removenum);
        }
        else
        {
            INFO("Need to insert atoms");
            insertnum=0;
            dr=r1-r0; dyc=dr*py; dxc=b*px; dxc=dxc*(-1.0);
            for(i=0;i<_NP;i++)
            {
                /* limit dislocation length (cut depth) */
                if((int)input[14])
                    if((_SR[i][zind]<input[15])||(_SR[i][zind]>=input[16]))
                        continue;
                ds=_SR[i]; ds-=sr0; ds.subint();
                dr=_H*ds; dx=dr*px; dy=dr*py;
                if((dx>=-dxc*0.5)&&(dx<dxc*0.5)&&(dy>=0)&&(dy<dyc))
                {
                    INFO_Printf("need replicate atom %d\n",i);
                    _SR[insertnum+_NP] = _SR[i];
                    _R[insertnum+_NP] = _R[i];
                    storedr[insertnum+_NP]=storedr[i];
                    insertnum++;
                }
            }
            INFO("need insertnum="<<insertnum);
        }
    }

    dua.clear(); dub.clear(); duc.clear();

    INFO_Printf("minx = %d  maxx = %d  miny = %d maxy = %d\n",minx,maxx,miny,maxy);
    if((minx!=maxx)||(miny!=maxy))
    {/* pbc correction only if pbc is applied */
        for(j=minx;j<=maxx;j++)
            for(k=miny;k<=maxy;k++)
            {
                getdu(sra); dua+=du;
                getdu(srb); dub+=du;
                getdu(src); duc+=du;
                //INFO("dua="<<dua<<" dub="<<dub<<" duc="<<duc);
            }
    }
    ux=dub-dua;
    uy=duc-dua;
    INFO("dua="<<dua<<" dub="<<dub<<" duc="<<duc);
    INFO("overall tilt ux="<<ux<<" uy="<<uy);

    if ((ux.norm()==0) && (uy.norm()==0))
    {
        FATAL("ux and uy = 0.  This should not happen!  Try compile with build=D.");
    }

    SHtoR();
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        /* limit dislocation length (cut depth) */
        if((int)input[14])
            if((_SR[i][zind]<input[15])||(_SR[i][zind]>=input[16]))
                continue;
        for(j=minx;j<=maxx;j++)
            for(k=miny;k<=maxy;k++)
            {
                st=_SR[i]; st[zind]=0; getdu(st);
                storedr[i]+=du;
            }
        storedr[i]-=ux*_SR[i][xind];
        storedr[i]-=uy*_SR[i][yind];
        if(i==0) INFO("storedr[0]="<<storedr[0]);
    }
    INFO("originally "<<_NP<<" atoms");
    INFO("insert "<<insertnum<<" atoms");
    _NP+=insertnum;
    INFO("now "<<_NP<<" atoms");


    /* shift box */
    if(input[13])
    {
        INFO("shiftbox to accommodate plastic strain");
        frac=-(input[7]-input[6]);
        dH[0][xind]+=(_H[0][0]*input[2]+_H[0][1]*input[3]
                      +_H[0][2]*input[4])*frac;
        dH[1][xind]+=(_H[1][0]*input[2]+_H[1][1]*input[3]
                      +_H[1][2]*input[4])*frac;
        dH[2][xind]+=(_H[2][0]*input[2]+_H[2][1]*input[3]
                      +_H[2][2]*input[4])*frac;
    }

    INFO_Printf("makedipole: store=%d\n",store);
    /* if store == 0, then commit dr, dH, remove atoms
       if store == 1, then commit nothing
       if store == 2, then commit remove atoms, but not dr, dH 
          (for creating reference structure to compute displacement field) */
    if(store==0)
    {
        /* commit displacement field */
        commit_storedr();

        INFO("free storedr");
        free(storedr);
        storedr = NULL;
     }

    if(store==0 || store==2)
    {
        /* commit removing atoms */
        commit_removeatoms();
    }

}

void MDFrame::commit_storedr()
{
    int i;
    /* commit displacement field */
    SHtoR();
    for(i=0;i<_NP;i++)
    {
        _R[i]+=storedr[i];
    }
    free(storedr); storedr=NULL; /* Wei 2014/8/29 */

    RHtoS();

    /* commit shift box */
    _H+=dH;
    dH.clear();
    SHtoR();
}

void MDFrame::commit_removeatoms()
{
    int i, n, removenum;

        n=0; SHtoR();
        removenum=0;
        for(i=0;i<_NP;i++)
        {
            if(fixed[i]!=-1)
            {
                _R0[n]=_R[i];
                species[n]=species[i]; //0;
                group[n]=group[i]; //0;
                fixed[n]=fixed[i]; //0;
                if (_SR1!=NULL) _SR1[n]=_SR1[i];
                if (_SR2!=NULL) _SR2[n]=_SR2[i];
                n++;
            }
            else
            {
                removenum++;
            }
        }
        INFO("originally "<<_NP<<" atoms");
        INFO("remove "<<removenum<<" atoms");
        _NP-=removenum;
        INFO("now "<<_NP<<" atoms");
        R0toR();  RHtoS();
        INFO("NP="<<_NP<<"  n="<<n);
        /* end of remove atoms */
}

#include "solidangle.c"

void MDFrame::makedislpolygon()
{/*
  * Make polygonal dislocation loop
  * input = [ activate store nu
  *           a        bx by bz
  *           Nnode    x1 y1 z1
  *                    x2 y2 z2
  *                    ...
  *                    xn yn zn ]
  *
  * Note: 1. displacement field only proportional to solid angle
  *       2. all cut planes pass through first node (x1,y1,z1)
  *       3. cannot remove atoms yet, so only for glide dislocation loops
  *
  * Ref.1 D. M. Barnett, "The displacement field of a triangular dislocation loop"
  *       Phil. Mag. A. (1985) vol.51 383--387
  * Ref.2 D. M. Barnett and R. W. Balluffi, "The displacement field of a triangular
  *       dislocation loop -- a correction with commentary", Phil. Mag. Lett.
  *       (2007) vol.87 943--944
  * Ref.3 A. Van Oosterom and J. Strackee, "The Solid angle of a plane trangle",
  *       IEEE transactions on biomedical engr. (1983) vol.bme-30 125--126
  */
    int store, Nnode, i;
    Vector3 burg, f, g, du;
    Vector3 sra, srb, src, srd, ra, rb, rc, rd, dua, dub, duc, dud;
    double nu, a, Omega;

    if(input[0]==0)
    {
        INFO_Printf("input = [%f %f %f ]\n",
                    input[0],input[1],input[2]);
        ERROR("makedislpolygon(): no dislocation geometry set");
        return;
    }

    store = (int) input[1];
    nu = input[2]; a = input[3];
    burg.set(input[4]*a,input[5]*a,input[6]*a);
    Nnode = (int) input[7];
        
    if(storedr==NULL)
    {
        INFO("allocate storedr");
        storedr = (Vector3 *)malloc(sizeof(Vector3)*_NP*allocmultiple);
        bindvar("storedr",storedr,DOUBLE);
        for(i=0;i<_NP;i++) storedr[i].clear();
    }

#define getdu_polyloop(R) protect(\
        Omega = SolidAnglePolygon(Nnode, input+8, R.x, R.y, R.z);\
        f = get_f_Polygon(Nnode, input+8, R.x, R.y, R.z, &burg);\
        g = get_g_Polygon(Nnode, input+8, R.x, R.y, R.z, &burg);\
        f *= ((1 - 2*nu)/(8*M_PI*(1-nu)));\
        g /= (8*M_PI*(1-nu));\
        du = burg* (Omega/4.0/M_PI);\
        du = du - f; du = du + g;)
    
    sra.set(-0.5,-0.5,-0.5);
    srb.set( 0.5,-0.5,-0.5);
    src.set(-0.5, 0.5,-0.5);
    srd.set(-0.5,-0.5, 0.5);
    
    ra=_H*sra;rb=_H*srb;rc=_H*src;rd=_H*srd;
    getdu_polyloop(ra); dua=du; //INFO("dua = "<<dua);
    getdu_polyloop(rb); dub=du; //INFO("dub = "<<dub);
    getdu_polyloop(rc); duc=du; //INFO("duc = "<<duc);
    getdu_polyloop(rd); dud=du; //INFO("dud = "<<dud);
    dub -= dua; duc -= dua; dud -= dua;
    output_dat[0] = dub.x; output_dat[1] = dub.y; output_dat[2] = dub.z;
    output_dat[3] = duc.x; output_dat[4] = duc.y; output_dat[5] = duc.z;
    output_dat[6] = dud.x; output_dat[7] = dud.y; output_dat[8] = dud.z;

    SHtoR();
    for(i=0;i<_NP;i++)
    {
        getdu_polyloop(_R[i]);
        storedr[i] += du;
    }
    
    if(!store)
    {
        /* commit displacement field */
        commit_storedr();

        free(storedr);
        storedr = NULL;
    }
}

void MDFrame::makedislellipse()
{   /* create an arbitrarily oriented elliptical dislocation loop
     * Input format:
     *            input = [ 0/1      : activation
     *                      Ra Rb    : two axes of the ellipse (in Angstrom)
     *                      a        : length scale (in Angstrom)
     *                      bx by bz : Burgers vector (in unit of a)
     *                      lx ly lz : line direction vector for the first ellipse axis
     *                      nx ny nz : normal vector of cut plane
     *                      x0 y0 z0 : center of the ellipse
     *                      0/1      : 0 commit, 1 store
     *                    ]
     * Notice:
     *       1. displacement field only proportional to solid angle
     *       2. all cut planes pass through first node (x1,y1,z1)
     *       3. cannot remove atoms yet, so only for glide dislocation loops
     *       4. l and n are vectors which will be automatically normalized
     *       5. still need to cancel mismatch at pbc boundary
     */
    double Ra, Rb, a, x, y, z, Omega;
    Vector3 b, l, n, r0, e1, e2, e3, du;
    Matrix33 M;
    int store, i;
    
    if(input[0]==0)
    {
        INFO_Printf("input = [%f %f %f ]\n",
                    input[0],input[1],input[2]);
        ERROR("makedislellipse(): no dislocation geometry set");
        return;
    }
    
    Ra = input[1]; Rb = input[2]; a = input[3];
    b.set(input[4]*a,input[5]*a,input[6]*a);
    l.set(input[7],input[8],input[9]);
    n.set(input[10],input[11],input[12]);
    r0.set(input[13],input[14],input[15]);
    store = (int) input[16];
    
    /* find real coordinates of b */
    if(fabs(b.x)<1e-8) b.x=0;
    if(fabs(b.y)<1e-8) b.y=0;
    if(fabs(b.z)<1e-8) b.z=0;
    
    /* set up coordinate system */
    e3=n;
    e3/=e3.norm();
    e1=l;
    e1/=e1.norm();
    e2=cross(e3,e1);
    M.setcol(e1,e2,e3);

    INFO("M = "<<M);

    if(storedr==NULL)
    {
        INFO("allocate storedr");
        Realloc(storedr,Vector3,_NP*allocmultiple);
        /*storedr = (Vector3 *)malloc(sizeof(Vector3)*_NP*allocmultiple);*/
        SHtoR();
        for(i=0;i<_NP;i++) storedr[i].clear();
    }

    /* compute displacement field of all atoms */
    SHtoR();
    for(i=0;i<_NP;i++)
    {
        x = (_R[i]-r0)*e1; y = (_R[i]-r0)*e2; z = (_R[i]-r0)*e3;
        Omega = SolidAngleEllipse(Ra, Rb, x, y, z);
        du = b * (Omega/4.0/M_PI);
        storedr[i] += du;
        if(i==0) INFO("storedr[0]="<<storedr[0]);
    }

    if(!store)
    {
        /* commit displacement field */
        commit_storedr();

        free(storedr);
        storedr = NULL;
    }

    return;
}

void MDFrame::makedislocation()
{
    /* create an arbitrarily oriented dislocation line
     * Input format:
     *            input = [ 0/1      : activation
     *                      a        : lattice constant (in Angstrom)
     *                      bx by bz : Burgers vector (in unit of a)
     *                      lx ly lz : line direction vector
     *                      nx ny nz : normal vector of cut plane
     *                      x0 y0 z0 : a point on dislocation (in Angstrom)
     *                      nu       : Poisson's ratio (for nonscrew dislocation)
     *                      0/1      : 0 commit, 1 store
     *                    ]
     * Notice:
     *   b, l, n are now specified in real coordinates (not reduced coordinates s)
     *     (b, l, n are all specified in Lab coordinates (not crystal coordinates)
     *     (we can use crystal coordinates only if the ex=[1 0 0], ey=[0 1 0], ez=[0 0 1]
     *   b is given in unit of a
     *   l and n are vectors which will be automatically normalized
     */
    double a, nu, nua, nub, bx, by, bz, x, y, r2, theta, pre;
    Vector3 b, l, n, r0, e1, e2, e3, b0, du, due;
    Vector3 sra, srb, src, srd, ra, rb, rc, rd, dua, dub, duc, dud;
    Vector3 ds, dr;
    Matrix33 M;
    int store, cut, i, k, removenum, insertnum;
    
    Vector3 *Rtmp;

    Rtmp = NULL;
    Realloc(Rtmp,Vector3,_NP);

    if(input[0]==0)
    {
        INFO_Printf("input = [%f %f %f ]\n",
                    input[0],input[1],input[2]);
        ERROR("makedislocation(): no dislocation geometry set");
        return;
    }
    a = input[1];
    b.set(input[2],input[3],input[4]);
    l.set(input[5],input[6],input[7]);
    n.set(input[8],input[9],input[10]);
    r0.set(input[11],input[12],input[13]);
    nu = input[14];
    store = (int) input[15];
    cut = (int) input[16];
    
    /* find real coordinates of b */
    b*=a;
    if(fabs(b.x)<1e-8) b.x=0;
    if(fabs(b.y)<1e-8) b.y=0;
    if(fabs(b.z)<1e-8) b.z=0;
    

    /* set up coordinate system */    
    e3=l;
    e3/=e3.norm();
    e2=n;
    e2/=e2.norm();
    e1=cross(e2,e3);
    M.setcol(e1,e2,e3);

    bx=b*e1; by=b*e2; bz=b*e3;
    b0.set(bx,by,bz);
    nua=1.0/4/(1-nu)/2/M_PI; nub=(1.0-2*nu)*nua;
    
    if(storedr==NULL)
    {
        INFO("allocate storedr");
        Realloc(storedr,Vector3,_NP*allocmultiple);
        /*storedr = (Vector3 *)malloc(sizeof(Vector3)*_NP*allocmultiple);*/

        SHtoR();
        for(i=0;i<_NP;i++) storedr[i].clear();

        /* Note: the following will not be compatible with makedislpolygon */
        /* for(i=0;i<_NP;i++)  storedr[i]=_R[i]; */
    }

    /* remove/insert atoms */
    removenum=insertnum=0;
    if(fabs(by)>1e-8)
    {
        INFO("Remove/insert atoms not implemented yet");
        if((by)<0)
        {
            INFO("Need to remove atoms");
            removenum=0;

            for(i=0;i<_NP;i++)
            {
                /* limit dislocation length (cut depth) */
                ds=_SR[i]; dr=_H*ds; dr-=r0; x=dr*e1; y=dr*e2;
                if((y>=-fabs(by)*0.5)&&(y<fabs(by)*0.5)&&(x<=0))
                {
                    INFO_Printf("need remove atom %d\n",i);
                    fixed[i] = -1; /* flag for removal */
                    removenum++;
                }
            }
            INFO("need removenum="<<removenum);
        }
        else
        {
            INFO("Need to insert atoms, not implemented yet");
            insertnum=0;
            INFO("need insertnum="<<insertnum);
        }
    }

#define getdu_r(R) protect(\
        x = (R-r0)*e1;\
        y = (R-r0)*e2;\
        if(cut)\
        {\
            if(y>0)\
                du=b0*0.5;\
            else\
                du=b0*(-0.5);\
        }\
        else\
        {\
            r2 = x*x + y*y;\
            theta = atan2(y, x);\
            pre = 0.5/M_PI*theta;\
            du = b0*pre;\
            if ((fabs(bx)>1e-8)||(fabs(by)>1e-8))\
            {\
                due.clear();\
                due.x-=(-bx*sin(theta*2)*nua+by*(-log(r2)*nub+cos(theta*2)*nua));\
                due.y-=( by*sin(theta*2)*nua-bx*(-log(r2)*nub-cos(theta*2)*nua));\
                du += due;\
            }\
        }\
        du = M*du;)

    sra.set(-0.5,-0.5,-0.5);
    srb.set( 0.5,-0.5,-0.5);
    src.set(-0.5, 0.5,-0.5);
    srd.set(-0.5,-0.5, 0.5);
    
    ra=_H*sra;rb=_H*srb;rc=_H*src;rd=_H*srd;
    getdu_r(ra); dua=du; //INFO("dua = "<<dua);
    getdu_r(rb); dub=du; //INFO("dub = "<<dub);
    getdu_r(rc); duc=du; //INFO("duc = "<<duc);
    getdu_r(rd); dud=du; //INFO("dud = "<<dud);
    dub -= dua; duc -= dua; dud -= dua;
    output_dat[0] = dub.x; output_dat[1] = dub.y; output_dat[2] = dub.z;
    output_dat[3] = duc.x; output_dat[4] = duc.y; output_dat[5] = duc.z;
    output_dat[6] = dud.x; output_dat[7] = dud.y; output_dat[8] = dud.z;
    
    SHtoR();
    for(i=0;i<_NP;i++)
    {
        getdu_r(_R[i]);
        storedr[i]+=du;
        if(i==0) INFO("storedr[0]="<<storedr[0]);
    }

    if(!store)
    {
        /* commit displacement field */
        commit_storedr();

        free(storedr);
        storedr = NULL;

        /* commit removing atoms */
        k=0; SHtoR();
        removenum=0;
        for(i=0;i<_NP;i++)
        {
            if(fixed[i]!=-1)
            {
                Rtmp[k]=_R[i];
                fixed[k]=fixed[i];
                k++;
            }
            else
            {
                removenum++;
            }
        }
        INFO("originally "<<_NP<<" atoms");
        INFO("remove "<<removenum<<" atoms");
        _NP-=removenum;
        INFO("now "<<_NP<<" atoms");
        for(i=0;i<_NP;i++) _R[i]=Rtmp[i];
        RHtoS();
        INFO("NP="<<_NP<<"  n="<<n);
    }

    /* by default fixedatomenergypartition == 0 */
    if (fixedatomenergypartition) setfixbufferatoms();

    free(Rtmp);
}


void MDFrame::makecylinder()
{
    /* Cut a cylinder out of the current configuration
     *
     * input = [ zind yind x0 y0 rrem flag szmin szmax ]
     *
     *  zind: direction of the cylinder axis
     *  yind: direction of another axis (perpendicular to zind)
     *  x0, y0: center of the cylinder axis
     *  rrem: radius of the cylinder
     *  flag: if flag==0, then all atoms outside the cylinder is removed
     *        if flag==1, then all atoms outside the cylinder is fixed
     */
    
    int i, xind, yind, zind, flag;
    double r, rrem;
    double x, y, bx, by;
    class Vector3 sr0,r0,ds,dr,sb,b,du,st,due;
    class Vector3 ux,uy;
    class Vector3 px, py, pz, pa, pb, pc;
    double pxpa, pxpb, pypa, pypb, szmin, szmax;
    int removenum,n;
    Vector3 *Rtmp;

    Rtmp = NULL;
    Realloc(Rtmp,Vector3,_NP);

    if(input[0]==0)
    {
        INFO_Printf("input = [%f %f %f ]\n",
                    input[0],input[1],input[2]);
        ERROR("makecylinder(): no dislocation geometry set");
        return;
    }

    /* read in parameters */
    zind=(int)input[0]; /* 1x/2y/3z */
    yind=(int)input[1];
    xind=zind+yind;if(xind==5)xind=1;if(xind==4)xind=2;
    xind--;yind--;zind--;
    
    /* set up coordinate system */
    sr0[xind]=input[2];sr0[yind]=input[3];sr0[zind]=0;
    
    r0=_H*sr0;
    pz.set(_H[0][zind],_H[1][zind],_H[2][zind]);
    pz/=pz.norm();
    py.set(_H[0][yind],_H[1][yind],_H[2][yind]);
    py.orth(pz); py/=py.norm();
    px=cross(py,pz);
    INFO("px="<<px);
    INFO("py="<<py);
    INFO("pz="<<pz);
    pa.set(_H[0][xind],_H[1][xind],_H[2][xind]); //pa/=pa.norm();
    pb.set(_H[0][yind],_H[1][yind],_H[2][yind]); //pb/=pb.norm();
    pc.set(_H[0][zind],_H[1][zind],_H[2][zind]); //pc/=pc.norm();
    pxpa=px*pa; pxpb=px*pb; pypb=py*pb; pypa=py*pa;
    bx=b*px; by=b*py;
    INFO("bx="<<bx<<"  by="<<by);
    if(fabs(pxpa)<1e-8) pxpa=0;
    if(fabs(pxpb)<1e-8) pxpb=0;
    if(fabs(pypb)<1e-8) pypb=0;
    if(fabs(pypa)<1e-8) pypa=0;
    if(pxpa<0) { pxpa=-pxpa; pxpb=-pxpb; px*=-1; INFO("redefine x");}
    INFO_Printf("pxpa=%e pxpb=%e pypb=%e pypa=%e\n",pxpa,pxpb,pypb,pypa);

    rrem=input[4];
    flag=(int) input[5];
    szmin = input[6];
    szmax = input[7];

    if(szmax>szmin)
    {
        INFO_Printf("only atoms with r[%d] between %e and %e will be cut\n",
                    zind,szmin,szmax);
    }
    /* remove all atoms outside rrem */
#define getdr(SR) protect(\
                ds=SR; ds-=sr0; ds[zind]=0; \
                x=ds[xind]*pxpa+ds[yind]*pxpb;\
                y=ds[xind]*pypa+ds[yind]*pypb; r=sqrt(x*x+y*y);)
    
    SHtoR();
    for(i=0;i<_NP;i++)
    {
        getdr(_SR[i]);

        if(szmax<=szmin)
        {
            if(r>rrem)
            {
                if(flag==1)
                    fixed[i]=1;
                else
                    fixed[i]=-1;
            }
        }
        else
        {                       
            if((_SR[i][zind]>=szmin)&&(_SR[i][zind]<szmax))
            {
                if(r>rrem)
                {
                    if(flag==1)
                        fixed[i]=1;
                    else
                        fixed[i]=-1;
                }
            }
        }
    }

    /* commit removing atoms */
    n=0; SHtoR();
    removenum=0;  _NPfixed = 0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]!=-1)
        {
            Rtmp[n]=_R[i];
            fixed[n]=fixed[i];
            
            /* Wei Cai, 8/9/2008 */
            if(_SAVEMEMORY<8) species[n]=species[i];
            n++;
        }
        else
        {
            removenum++;
        }
        if(fixed[i]==1) _NPfixed++;
    }
    INFO("originally "<<_NP<<" atoms");
    INFO("remove "<<removenum<<" atoms");
    _NP-=removenum;
    INFO("now "<<_NP<<" atoms");
    INFO("now "<<_NPfixed<<" atoms fixed");

    for(i=0;i<_NP;i++) _R[i]=Rtmp[i];
    RHtoS();
    INFO("NP="<<_NP<<"  n="<<n);

    /* by default fixedatomenergypartition == 0 */
    if (fixedatomenergypartition) setfixbufferatoms();
    free(Rtmp);
}

void MDFrame::setfixbufferatoms()
{
/* find all fixed atoms that are interacting with free atoms
 * set them fixed[i]=2
 */
    int i, j, jpt, ninter, nfree, nfixed;
    Vector3 sij, rij;
    double r;
    
    INFO("setfixbufferatoms");
    refreshneighborlist();
    
    for(i=0;i<_NP;i++)
    {
        if(!fixed[i]) continue;
        fixed[i]=1;
        for(j=0;j<nn[i];j++)
        {
            jpt=nindex[i][j];
            if(fixed[jpt]) continue;
            sij=_SR[jpt]-_SR[i];
            sij.subint();
            rij=_H*sij;
            r=rij.norm();
            if(r<_RLIST-_SKIN) fixed[i]=2;
        }                
    }
    ninter=0; nfree=0; nfixed=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]!=1) ninter++;
        if(fixed[i]==0) nfree++;
        if(fixed[i]) nfixed++;
    }
    INFO("number of interacting atoms: ninter="<<ninter);
    INFO("number of free        atoms: nfree ="<<nfree);
    INFO("number of fixed       atoms: nfixed="<<nfixed);
}

void MDFrame::makedislcylinder()
{
    /* Create single dislocation in a cylindrical fixed boundary */
    
    int i, xind, yind, zind;
    double r, rfix, rgfbc, rrem;
    
    /*double theta[2], pre, x, y, r21, r22, bx, by, nu, nua, nub;*/
    double theta[2], pre, x, y, r21, bx, by, nu, nua, nub;
    class Vector3 sr0,r0,ds,dr,sb,b,du,st,due;
    class Vector3 ux,uy;
    class Vector3 px, py, pz, pa, pb, pc;
    double pxpa, pxpb, pypa, pypb;
    int removenum,n,insertnum;
    double dx, dxc;
    static class Vector3 *storedr;
    int store;

    INFO("outdated function, use makecylinder and makedislocation instead");
    INFO("see Computer Simulations of Dislocations, Section 3.1, ta-screw.script");
    
    if(input[0]==0)
    {
        INFO_Printf("input = [%f %f %f ]\n",
                    input[0],input[1],input[2]);
        ERROR("makedislocation(): no dislocation geometry set");
        return;
    }
    store = (int) input[11];
    
    /* read in parameters */
    zind=(int)input[0]; /* 1x/2y/3z */
    yind=(int)input[1];
    xind=zind+yind;if(xind==5)xind=1;if(xind==4)xind=2;
    xind--;yind--;zind--;
    
    sb.set(input[2],input[3],input[4]);
    b=_H*sb;
    if(fabs(b.x)<1e-8) b.x=0;
    if(fabs(b.y)<1e-8) b.y=0;
    if(fabs(b.z)<1e-8) b.z=0;

    /* set up coordinate system */
    sr0[xind]=input[5];sr0[yind]=input[6];sr0[zind]=0;
    
    r0=_H*sr0;
    pz.set(_H[0][zind],_H[1][zind],_H[2][zind]);
    pz/=pz.norm();
    py.set(_H[0][yind],_H[1][yind],_H[2][yind]);
    py.orth(pz); py/=py.norm();
    px=cross(py,pz);
    INFO("px="<<px);
    INFO("py="<<py);
    INFO("pz="<<pz);
    pa.set(_H[0][xind],_H[1][xind],_H[2][xind]); //pa/=pa.norm();
    pb.set(_H[0][yind],_H[1][yind],_H[2][yind]); //pb/=pb.norm();
    pc.set(_H[0][zind],_H[1][zind],_H[2][zind]); //pc/=pc.norm();
    pxpa=px*pa; pxpb=px*pb; pypb=py*pb; pypa=py*pa;
    bx=b*px; by=b*py;
    INFO("bx="<<bx<<"  by="<<by);
    if(fabs(pxpa)<1e-8) pxpa=0;
    if(fabs(pxpb)<1e-8) pxpb=0;
    if(fabs(pypb)<1e-8) pypb=0;
    if(fabs(pypa)<1e-8) pypa=0;
    if(pxpa<0) { pxpa=-pxpa; pxpb=-pxpb; px*=-1; INFO("redefine x");}
    INFO_Printf("pxpa=%e pxpb=%e pypb=%e pypa=%e\n",pxpa,pxpb,pypb,pypa);

    nu=input[7]; nua=1.0/4/(1-nu)/2/M_PI; nub=(1.0-2*nu)*nua;
    rfix=input[8]; rgfbc=input[9]; rrem=input[10];
    
    if(storedr==NULL)
    {
        INFO("allocate storedr");
        storedr = (Vector3 *)malloc(sizeof(Vector3)*_NP*allocmultiple);
        bindvar("storedr",storedr,DOUBLE);
        for(i=0;i<_NP;i++) storedr[i].clear();
    }

    /* remove/insert atoms */
    removenum=insertnum=0;
    if(fabs(b*px)>1e-8)
    {
        if((b*px)>0)
        {
            INFO("Need to remove atoms");
            removenum=0;
            dxc=b*px;
            for(i=0;i<_NP;i++)
            {
                ds=_SR[i]; ds-=sr0; ds.subint();
                dr=_H*ds; dx=dr*px; 
                if((dx>=-dxc*0.5)&&(dx<dxc*0.5))
                {
                    INFO_Printf("need remove atom %d\n",i);
                    fixed[i] = -1; /* flag for removal */
                    removenum++;
                }
            }
            INFO("need removenum="<<removenum);
        }
        else
        {
            INFO("Need to insert atoms");
            insertnum=0;
            dxc=b*px; dxc=dxc*(-1.0);
            for(i=0;i<_NP;i++)
            {
                ds=_SR[i]; ds-=sr0; ds.subint();
                dr=_H*ds; dx=dr*px; 
                if((dx>=-dxc*0.5)&&(dx<dxc*0.5))
                {
                    INFO_Printf("need replicate atom %d\n",i);
                    _SR[insertnum+_NP] = _SR[i];
                    _R[insertnum+_NP] = _R[i];
                    storedr[insertnum+_NP]=storedr[i];
                    insertnum++;
                }
            }
            INFO("need insertnum="<<insertnum);
        }
    }

#define getdu_single(SR) protect(\
                ds=SR; ds-=sr0; ds[zind]=0; \
                x=ds[xind]*pxpa+ds[yind]*pxpb;\
                y=ds[xind]*pypa+ds[yind]*pypb; r21=x*x+y*y;\
                theta[0]=atan2(-x,y);\
                pre=0.5/M_PI*(theta[0]);\
                du=b*pre; if((fabs(bx)>1e-8)||(fabs(by)>1e-8)){due.clear();\
                due-=px*(-bx*sin(theta[0]*2)*nua+by*(-log(r21)*nub+cos(theta[0]*2)*nua));\
                due-=py*( by*sin(theta[0]*2)*nua-bx*(-log(r21)*nub-cos(theta[0]*2)*nua));\
                /*due+=px*(-bx*sin(theta[1]*2)*nua+by*(-log(r22)*nub+cos(theta[1]*2)*nua));\*/\
                /*due+=py*( by*sin(theta[1]*2)*nua-bx*(-log(r22)*nub-cos(theta[1]*2)*nua));\*/\
                /*INFO("due="<<due);*/ du+=due;} \
                )
    
    SHtoR();
    for(i=0;i<_NP;i++)
    {
        getdr(_SR[i]);
        if(r>rrem) fixed[i]=-1;
        else if(r>rfix) fixed[i]=1;
        if(fixed[i]==-1) continue;
        st=_SR[i]; st[zind]=0; getdu_single(st);
        storedr[i]+=du;
        
        storedr[i]-=ux*_SR[i][xind];
        storedr[i]-=uy*_SR[i][yind];
        if(i==0) INFO("storedr[0]="<<storedr[0]);
    }
    INFO("originally "<<_NP<<" atoms");
    INFO("insert "<<insertnum<<" atoms");
    _NP+=insertnum;
    INFO("now "<<_NP<<" atoms");


    if(!store)
    {
        /* commit displacement field */
        commit_storedr();

        free(storedr);
        storedr = NULL;

        /* commit removing atoms */
        n=0; SHtoR();
        removenum=0;
        for(i=0;i<_NP;i++)
        {
            if(fixed[i]!=-1)
            {
                _R0[n]=_R[i];
                fixed[n]=fixed[i];
                n++;
            }
            else
            {
                removenum++;
            }
        }
        INFO("originally "<<_NP<<" atoms");
        INFO("remove "<<removenum<<" atoms");
        _NP-=removenum;
        INFO("now "<<_NP<<" atoms");
        R0toR();  RHtoS();
        INFO("NP="<<_NP<<"  n="<<n);

    }

    /* by default fixedatomenergypartition == 0 */
    if (fixedatomenergypartition) setfixbufferatoms();
}

#if 0 /* disable this function (replaced by makedislpolygon() or makedislellipse() */
void MDFrame::makedisloop()
{
    /* Create Dislocation loop by using approximated
       isotropic linear elasticity solution under
       Periodic Boundary Condition with image summations */
    /*
      !!! Not fully implemented
          Not ready to use !!!
    */
    class Vector3 normal, sb, b, du, sr, rt; 
    double radius, height, V;
    int minx, maxx, miny, maxy, minz, maxz, i, j, k, m, n;
    int removenum, store;
    Matrix33 epsilon;
    static class Vector3 *storedr;
    static Matrix33 dH;
        
    if(input[0]==0)
    {
        ERROR("makedisloop(): no dipole geometry set");
        return;
    }
    store = (int) input[16];
    
    /* read in parameters */
    normal.set(input[0],input[1],input[2]);
    normal/=normal.norm();
    INFO("normal="<<normal);
    radius=input[3];
    height=input[4];
    sb.set(input[5],input[6],input[7]);
    b=_H*sb;
    if(fabs(b.x)<1e-8) b.x=0;
    if(fabs(b.y)<1e-8) b.y=0;
    if(fabs(b.z)<1e-8) b.z=0;
    INFO("sb="<<sb<<"\nb="<<b);
    
    /* nu=input[8]; */
    minx=(int)input[9]; maxx=(int)input[10];
    miny=(int)input[11];maxy=(int)input[12];
    minz=(int)input[13];maxz=(int)input[14];
    
    
    if(storedr==NULL)
    {
        INFO("allocate storedr");
        storedr = (Vector3 *)malloc(sizeof(Vector3)*_NP*allocmultiple);
        bindvar("storedr",storedr,DOUBLE);
        for(i=0;i<_NP;i++) storedr[i].clear();
    }
    
    /* remove/insert atoms */
    removenum=0;
    if(fabs(b*normal)>1e-4)
    {
        INFO("b*normal="<<b*normal);
        if((b*normal)>0)
        {
            INFO("Need to remove atoms: not implemented!");
            return;
        }
        else
        {
            INFO("Need to insert atoms: not implemented!");
            return;
        }
    }   

    /* ignore c.c. image contributions for now
       (need to implement !!)
     */
    
    SHtoR();
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        for(j=minx;j<=maxx;j++)
            for(k=miny;k<=maxy;k++)
                for(m=minz;m<=maxz;m++)
                {
                    /* calculate the displacement field */
                    sr=_SR[i]; sr.x+=j; sr.y+=k; sr.z+=m;
                    rt=_H*sr;
                    calloopdisp(&normal,radius,height,&b,&rt,&du);
//                    INFO("du="<<du);
                    storedr[i]+=du;
                }
        if(i==0) INFO("storedr[0]="<<storedr[0]);
    }

    /* ignore shift box for now */
    if(input[15])
    {
        INFO("shiftbox");
        epsilon.clear();
        V=_H.det();
        INFO("V="<<V);
        epsilon.addnvv
            (M_PI*radius*radius/(sqrt(3.0)*_H[0][0]*_H[0][0]),
             b,normal);
        INFO("epsilon="<<epsilon);
//        dH=epsilon*_H;
        dH=epsilon;
    }

    if(!store)
    {
        /* commit displacement field */
        commit_storedr();

        free(storedr);
        storedr = NULL;

        /* commit removing atoms */
        n=0; SHtoR();
        removenum=0;
        for(i=0;i<_NP;i++)
        {
            if(fixed[i]!=-1)
            {
                _R0[n]=_R[i];
                fixed[n]=0;
                n++;
            }
            else
            {
                removenum++;
            }
        }
        INFO("originally "<<_NP<<" atoms");
        INFO("remove "<<removenum<<" atoms");
        _NP-=removenum;
        INFO("now "<<_NP<<" atoms");
        R0toR();  RHtoS();
        INFO("NP="<<_NP<<"  n="<<n);
        /* end of remove atoms */

        /* commit shift box */
        _H+=dH;
    }   
}
#endif /* end of makedisloop */

void MDFrame::calloopdisp(Vector3 *normal,double radius,double height,
                          Vector3 *b, Vector3 *rt,Vector3 *du)
{/* compute displacement field of a dislocation loop
    du = -b Omega / (4 pi)
    Omega = -\int r.dA/r^3
  */
    double hr, dx, x, y, z, r, r2, R, Omega;
    Vector3 rr;
    int Ns, i, j;
    /* reduced height, field point position, w.r.t. radius */
    hr=height/radius;
    rr=*rt; rr/=radius;
    z=rr*(*normal);
    r2=(rr.norm2()-z*z);
    if(r2>0) r=sqrt(r2);
    else r=0;
    
    /* subdivision number */
    if(r<1.0/5.5)
        Ns=5;
    else
        Ns=(int)(1.0/r-0.5);
    
    dx=2.0/(2*Ns+1);

    if(fabs(z*radius)<3)
    {
        if(r>1)
            Omega=0;
        else
            if(z>0)
                Omega=-2*M_PI;
            else
                Omega=2*M_PI;
        *du = *b;
        *du *= (-Omega/4/M_PI);
        return;
    }
    Omega=0;
    for(i=-Ns;i<=Ns;i++)
        for(j=-Ns;j<=Ns;j++)
        {
            /* Galerkin Points */
            x=dx*i;
            y=dx*j;
            if(x*x+y*y<=1)
            {
                R=sqrt((x-r)*(x-r)+y*y+z*z);
                Omega+=1./(R*R*R);
            }
        }        
    Omega*=(-z*dx*dx);
    *du = *b;
    *du *= (-Omega/4/M_PI);
    return;
}

void MDFrame::makegrainboundary()
{
    /* create grain boundary */
    int i,n,imgx,imgy;
    int gbtype;
    double sx0,sy0,p,q,alpha;
    Matrix33 rot1, rot2, hnew, hinv;
    Vector3 r0, sr0, r1,dr,drr,s1,sr;
    
    gbtype = (int)input[0];
    sx0=input[1];
    sy0=input[2];
    p=input[3];
    q=input[4];
    /* shift[3] deprecated 
    shift[0]=input[5];
    shift[1]=input[6];
    shift[2]=input[7]; */
    
    if(gbtype==0)
    {
        ERROR("makegb(): no gb geometry set");
        return;
    }

    alpha=atan2(p,q);
    rot1.set(cos(alpha),-sin(alpha),0,
             sin(alpha),cos(alpha),0,
             0,0,1);
    rot2.set(cos(alpha),sin(alpha),0,
             -sin(alpha),cos(alpha),0,
             0,0,1);
    hnew=_H;
    hnew[0][0]/=cos(alpha);
    hnew[0][1]/=cos(alpha);
    hnew[0][2]/=cos(alpha);
    hnew[1][0]/=cos(alpha);
    hnew[1][1]/=cos(alpha);
    hnew[1][2]/=cos(alpha);
    hinv=hnew.inv();
    sr0.set(sx0,sy0,0);
    r0=_H*sr0;
    n=0;
    for(i=0;i<_NP;i++)
    {
        {/* Crystal 1 */
            /* fit atoms back in the cell */
            for(imgx=-1;imgx<=1;imgx++)
                for(imgy=-1;imgy<=1;imgy++)
                {
                    sr=_SR[i]; sr.x+=imgx; sr.y+=imgy;
                    r1=_H*sr;
                    dr=r1-r0;
                    drr=rot1*dr;
                    _R[i]=r0+drr;
                    s1=hinv*_R[i];
                    if((s1.x>sx0)&&(s1.x<0.5)
                       &&(s1.y>=-0.5)&&(s1.y<0.5))
                    {
                        INFO("atom "<<i<<" "<<s1);
                        _R0[n]=s1;
                        n++;
                        if(n>_NP)
                            INFO("number of atoms "<<n<<" exceed "<<_NP);
                    }
                }
        }
        {/* Crystal 2 */
            /* fit atoms back in the cell */
            for(imgx=-1;imgx<=1;imgx++)
                for(imgy=-1;imgy<=1;imgy++)
                {
                    sr=_SR[i]; sr.x+=imgx; sr.y+=imgy;
                    r1=_H*sr;
                    dr=r1-r0;
                    drr=rot2*dr;
                    _R[i]=r0+drr;
                    s1=hinv*_R[i];
                    if((s1.x<=sx0)&&(s1.x>=-0.5)
                       &&(s1.y>=-0.5)&&(s1.y<0.5))
                    {
                        INFO("atom "<<i<<" "<<s1);
                        _R0[n]=s1;
                        n++;
                        if(n>_NP)
                            INFO("number of atoms "<<n<<" exceed "<<_NP);
                    }
                }
        }
    }
    INFO("makegb: old NP="<<_NP<<" new NP="<<n);
    _NP=n;
    for(i=0;i<_NP;i++)
    {
        _SR[i]=_R0[i];
        _R[i]=_H*_SR[i];
        _R0[i]=_R[i];
    }
}

void MDFrame::makewave()
{
    /* Give atom a displacment field in the form of a plane wave
     *
     * input = [ enable amplitude dr.x dr.y dr.z k.x k.y k.z vel ]
     *
     * enable:    if 1 then makewave is enabled
     * amplitude: maximum displacement of the atoms (scalar)
     * dr(x,y,z): vibration direction (will be noramlized by MD++)
     *            amplitude * dr is the maximum displacement vector
     * k(x,y,z):  plane wave vector
     * vel:       velocity of atoms (nonzero vel not implemented)
     */
    double amp;
    Vector3 dr, kvec;
    int i, vel;
    
    if(input[0]==0)
    {
        INFO_Printf("mkwavespec = [%f %f %f ]\n",
                    input[0],input[1],input[2]);
        ERROR("makewave(): no geometry set");
        return;
    }

    amp = input[1];
    dr.set(input[2],input[3],input[4]);
    dr/=dr.norm();
    dr*=amp;
    kvec.set(input[5],input[6],input[7]);
    vel = (int)input[8];

    SHtoR();
    for(i=0;i<_NP;i++)
    {
        _R[i]+=dr*cos(kvec*_SR[i]*2*M_PI);
    }
    RHtoS();
    if(vel!=0)
    {
        INFO_Printf("makewave(): vel=%d not implemented\n",vel);
        ERROR("makewave()");
    }
}

#ifdef _TORSION_OR_BENDING    
#include "../cookies/src/md_torsion.cpp"
#include "../cookies/src/md_bending.cpp"
#endif


void MDFrame::scaleH()
{
    /* multiply H matrix by input[0] */
    INFO(_H);
     _H*=input[0];
    INFO(_H);
}
void MDFrame::saveH()
{
    /* copy H to H0 */
    INFO(_H);
    _H0=_H;
    INFO(_H);
}
void MDFrame::setH()
{
    /* input = [ i j c ]
       H[i][j] = c
    */
    int i,j;
    double h;
    i=(int)input[0];
    j=(int)input[1];
    h=input[2];
    _H[i][j]=h;
}
void MDFrame::restoreH()
{
    /* copy H0 to H */
    INFO(_H);
    _H=_H0;
    INFO(_H);
}






void MDFrame::shiftbox()
{ /*
   * input = [ i j delta ]
   *
   *   H[:][i--] += H[:][j--]*delta
   *   while keeping s fixed (affine transformation)
   */   
    int xind, yind;
    double frac;
    if(input[0]==0)
    {
        ERROR("input(): no shift geometry set");
        return;
    }
    xind=(int)input[0]; //1x/2y/3z
    yind=(int)input[1];
    xind--;yind--;
    frac=input[2];
    
    INFO("The old H=\n"<<_H);
    
    _H[0][xind]+=_H[0][yind]*frac;
    _H[1][xind]+=_H[1][yind]*frac;
    _H[2][xind]+=_H[2][yind]*frac;
    
    INFO("The new H=\n"<<_H);
    SHtoR();
    VSHtoVR();
}

void MDFrame::redefinepbc()
{ /*
   * input = [ i j delta ]
   *
   *   H[:][i--] += H[:][j--]*delta
   *   while keeping r fixed (introduce displacement at box boundary)
   */   
    int xind, yind;
    double frac;
    if(input[0]==0)
    {
        ERROR("input(): no shift geometry set");
        return;
    }
    xind=(int)input[0]; //1x/2y/3z
    yind=(int)input[1];
    xind--;yind--;
    frac=input[2];
    
    INFO("The old H=\n"<<_H);

    SHtoR();
    VSHtoVR();
    
    _H[0][xind]+=_H[0][yind]*frac;
    _H[1][xind]+=_H[1][yind]*frac;
    _H[2][xind]+=_H[2][yind]*frac;
    
    INFO("The new H=\n"<<_H);
    RHtoS();
    VRHtoVS();
}

void MDFrame::applystrain()
{ /*
   * apply a strain to box H while keeping s fixed
   * (affine transformation)
   *
   * input = [ M11 M12 M13 M21 M22 M23 M31 M32 M33
               s11 s12 s13 s21 s22 s23 s31 s32 s33
               eps ]
   */
    class Matrix33 rot, rotv, strain, rstrain; double tmp;

    rot.set(input);
    /* normalize rot matrix */
    tmp=sqrt(rot[0][0]*rot[0][0]+rot[0][1]*rot[0][1]+rot[0][2]*rot[0][2]);
    rot[0][0]/=tmp;rot[0][1]/=tmp;rot[0][2]/=tmp;
    tmp=sqrt(rot[1][0]*rot[1][0]+rot[1][1]*rot[1][1]+rot[1][2]*rot[1][2]);
    rot[1][0]/=tmp;rot[1][1]/=tmp;rot[1][2]/=tmp;
    tmp=sqrt(rot[2][0]*rot[2][0]+rot[2][1]*rot[2][1]+rot[2][2]*rot[2][2]);
    rot[2][0]/=tmp;rot[2][1]/=tmp;rot[2][2]/=tmp;
    rotv=rot.tran();
    
    strain.set(input+9);
    strain*=input[19];
    rstrain=rot*strain;
    rstrain=rstrain*rotv;

    INFO("applystrain rstrain="<<rstrain);
    _H0=_H;
    _H+=rstrain*_H;
}

void MDFrame::extendbox()
{ /*
   * input = [ i n ]
   * 
   * The i-th direction (i=1,2,3) of the box will be multiplied by 
   *  n times (by inserting more atoms)
   * Usually allocmultiple = n should be set before
   *  makecrystal or readcn is called before this function
   */
    int dir, n;
    int i,j;

    dir=(int)input[0]; n=(int)input[1];
    if((dir<=0)||(dir>=4)) return;
    if(n<=1) return;
    dir--;

    for(i=0;i<_NP;i++)
    {
        for(j=1;j<n;j++)
        {
            _SR[i+_NP*j]=_SR[i];
            _SR[i+_NP*j][dir]+=j;
            fixed[i+_NP*j]=fixed[i];

            if(_SAVEMEMORY>=8) continue; 
            species[i+_NP*j]=species[i];
            _EPOT_IND[i+_NP*j]=_EPOT_IND[i];
            _TOPOL[i+_NP*j]=_TOPOL[i];            
            if(_SAVEMEMORY>=7) continue;
            _VSR[i+_NP*j]=_VSR[i];
            if(_SAVEMEMORY>=6) continue;
            group[i+_NP*j]=group[i];
        }
    }
    for(i=0;i<_NP*n;i++)
    {
        _SR[i][dir]=(_SR[i][dir]+0.5)/n-0.5;
    }
    _H[0][dir]*=n;
    _H[1][dir]*=n;
    _H[2][dir]*=n;
    _NP*=n;
    SHtoR();
    INFO("extenbox: new NP="<<_NP);
}


void MDFrame::switchindex()
{ /* input = [ ind jnd ]
   * switch ind jnd index of each atomic coordinate, velocity and box vectors
   */

    int ind, jnd, i;
    Vector3 vtmp;
    Matrix33 htmp;

    ind=(int)input[0]; jnd=(int)input[1];
    ind--;  jnd--;

    for(i=0;i<_NP;i++)
    {
       vtmp=_SR[i]; _SR[i][ind] = vtmp[jnd]; _SR[i][jnd] = vtmp[ind]; 
       vtmp= _R[i];  _R[i][ind] = vtmp[jnd];  _R[i][jnd] = vtmp[ind]; 
       vtmp=_R0[i]; _R0[i][ind] = vtmp[jnd]; _R0[i][jnd] = vtmp[ind]; 

       vtmp=_VSR[i]; _VSR[i][ind] = vtmp[jnd]; _VSR[i][jnd] = vtmp[ind]; 
       vtmp= _VR[i];  _VR[i][ind] = vtmp[jnd];  _VR[i][jnd] = vtmp[ind]; 
    }

    htmp=_H;
    for(i=0;i<3;i++)
    {
       _H[i][ind]=htmp[i][jnd]; _H[i][jnd]=htmp[i][ind];
    }
    _H=_H.reorient();

    _VH.clear();
}

void MDFrame::maptoprimarycell()
{ 
    int i;

    for(i=0;i<_NP;i++)
       _SR[i].subint();
    SHtoR(); 
}

void MDFrame::cutslice()
{ /* input = [ i n ]
   *
   * cut a slice (1/n) of the configuration along the i-th direction
   *   opposite of extendbox
   * note:  atoms in cn file has to be organized in n chunks,
   *   i.e. the result of a previous extendbox operation
   */
    int dir, n;
    int i;

    dir=(int)input[0]; n=(int)input[1];
    if((dir<=0)||(dir>=4)) return;
    if(n<=1) return;
    dir--;

    INFO("dir = "<<dir<<"  n = "<<n);
    for(i=0;i<_NP/n;i++)
    {
        INFO_Printf("_SR:  %20.17e -> %20.17e\n",
                    _SR[i][dir],(_SR[i][dir]+0.5)*n-0.5);
        _SR[i][dir]=(_SR[i][dir]+0.5)*n-0.5;
    }
    _H[0][dir]/=n;
    _H[1][dir]/=n;
    _H[2][dir]/=n;
    _NP/=n;
    SHtoR();
    INFO("extenbox: new NP="<<_NP);
}


void MDFrame::splicecn()
{
    int i, dir, np0, np1;
    double ratio;
    double w;//blur width, disabled
    char cmd[1000], tmpfile[200];
    Matrix33 h0, h1; Vector3 c0,c1;

    INFO("splicecn");
    dir=(int)input[0];
    w=input[1];
    if((dir<=0)||(dir>=4)) return;
    dir--;

    /* old config */
    np0=_NP; h0=_H;
    writefinalcnfile(zipfiles,false);

    /* new config */
    readcn();
    np1=_NP; h1=_H;    

    sprintf(tmpfile,"tmp.cn");
    
    sprintf(cmd,"echo %d > %s",np0+np1,tmpfile);
    INFO_Printf("%s\n",cmd); system(cmd);
    sprintf(cmd,"head -%d %s | tail -%d >> %s",np0+1,finalcnfile,np0,tmpfile); /* old config */
    INFO_Printf("%s\n",cmd); system(cmd);
    sprintf(cmd,"head -%d %s | tail -%d >> %s",np1+4,incnfile, np1+3,tmpfile); /* new config */
    INFO_Printf("%s\n",cmd); system(cmd);

    sprintf(cmd,"echo %d ", nspecies);
    for(i=0;i<nspecies;i++)
        sprintf(cmd,"%s %s ",cmd,element[i]);
    sprintf(cmd,"%s >> %s",cmd,tmpfile);
    INFO_Printf("%s\n",cmd); system(cmd);

    strcpy(incnfile,tmpfile);
    readcn();

    c0.set(h0[0][dir],h0[1][dir],h0[2][dir]);
    c1.set(h1[0][dir],h1[1][dir],h1[2][dir]);
    ratio=c0.norm()/c1.norm();

    _NP = np0 + np1;    

    for(i=0;i<np0;i++)
        _SR[i][dir]=(_SR[i][dir]+0.5)*ratio/(ratio+1)-0.5;  /* old config */

    for(i=np0;i<_NP;i++)
        _SR[i][dir]=(_SR[i][dir]-0.5)/(ratio+1)+0.5;        /* new config */
    
    _H[0][dir]=c0[0]+c1[0];
    _H[1][dir]=c0[1]+c1[1];
    _H[2][dir]=c0[2]+c1[2];

    SHtoR();

    writefinalcnfile(zipfiles,false);
}

#if 0
void MDFrame::splicecn()  /* old implementation, has memory problem, do not use */
{ /*
   * splice (put) two configuration files together along the i-th direction
   *
   * input = [ i w ]
   *
   * If w > 0.001, then the atoms within distance w to the interface will
   *  be interpolated between the two configurations to avoid a sharp
   *  interface.
   */
    
    int i, dir;
    int *species0, *group0;
    Matrix33 h0; int np0; double ratio;
    Vector3 c,c0;
    double w;//blur width
    double alpha;

    INFO("splicecn");
    dir=(int)input[0];
    w=input[1];
    if((dir<=0)||(dir>=4)) return;
    dir--;

    species0 = (int *) malloc(sizeof(int)*_NP*allocmultiple);
    group0   = (int *) malloc(sizeof(int)*_NP*allocmultiple);
    
    for(i=0;i<_NP;i++)
    {
//        _SR[i].subint();
        _R0[i]=_SR[i];
        species0[i]=species[i];
        group0[i]=group[i];
    }
    h0=_H; np0=_NP;
    readcn(); 

    c.set(_H[0][dir],_H[1][dir],_H[2][dir]);
    c0.set(h0[0][dir],h0[1][dir],h0[2][dir]);
    ratio=c0.norm()/c.norm();
    
    for(i=0;i<_NP;i++)
    {
//        _SR[i].subint();
        _R0[i+np0]=_SR[i];
        species0[i+np0]=species[i];
        group0[i+np0]=group[i];
    }

    for(i=0;i<np0;i++)
    {
        if(w<0.001)
            alpha=1;
        else
        {
            if(_R0[i][dir]<(-0.5+w/2))
            {
                alpha=(0.5+_R0[i][dir])/w+0.5;
            }
            else if(_R0[i][dir]>(0.5-w/2))
            {
                alpha=(0.5-_R0[i][dir])/w+0.5;
            }
            else
                alpha=1;
        }
        //INFO_Printf("splicecn: i=%d alpha=%e\n",i,alpha);
        _SR[i]=_R0[i]*alpha+_R0[i+np0]*(1-alpha);
        _SR[i][dir]=(_R0[i][dir]+0.5)*ratio/(ratio+1)-0.5;
        species[i]=species0[i];
        group[i]=group0[i];
    }
    for(i=np0;i<np0+_NP;i++)
    {
        if(w<0.001)
            alpha=1;
        else
        {
            if(_R0[i][dir]<(-0.5+w/2))
            {
                alpha=(0.5+_R0[i][dir])/w+0.5;
            }
            else if(_R0[i][dir]>(0.5-w/2))
            {
                alpha=(0.5-_R0[i][dir])/w+0.5;
            }
            else
                alpha=1;
        }
//        INFO_Printf("splicecn: i=%d alpha=%e\n",i,alpha);
        _SR[i]=_R0[i]*alpha+_R0[i-np0]*(1-alpha);
        _SR[i][dir]=(_R0[i][dir]-0.5)/(ratio+1)+0.5;
        species[i]=species0[i];
        group[i]=group0[i];
    }
    
    _H[0][dir]=c[0]+c0[0];
    _H[1][dir]=c[1]+c0[1];
    _H[2][dir]=c[2]+c0[2];
    _NP+=np0;
    SHtoR();

    free(species0);
    free(group0);
    
    INFO("splicecn: new NP="<<_NP);
}
#endif


void MDFrame::cutpastecn()
{ /*
   * Replace all atoms whose i-th coordinates falls within [s0,s1]
   *  to those given in the incnfile.
   * Also displace the new atoms by ds in the i-th direction.
   *
   * Note: not finished, do not use!!
   *
   * input = [ i s0 s1 ds ]
   * incnfile = "filename"
   *
   */    
    int i,j,jpt,minj,dir,dirx,diry; double s0, s1, ds, mind, dis;
    Matrix33 h0; Vector3 sij;

    dir=(int)input[0];
    if((dir<1)||(dir>3)) return;
    dir--; /* dir=0,1,2 for x,y,z */
    dirx=(dir+1)%3; diry=(dir+2)%3;
    s0=input[1];
    s1=input[2];
    ds=input[3];

    /* store the S of current config */
    for(i=0;i<_NP;i++)
    {
        _VSR[i]=_SR[i];
    }
    h0=_H;
    /* read in perfect lattice as reference */
    initcn.open(incnfile,LFile::O_Read);
    initcn.read(this);
    for(i=0;i<_NP;i++)
    {
        _R0[i]=_SR[i];
    }
    _H=h0;
    for(i=0;i<_NP;i++)
    {
        _SR[i]=_VSR[i];
    }
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        if((_R0[i][dir]<s0)||(_R0[i][dir]>=s1)) continue;
        mind=100; minj=-1;
        for(j=0;j<nn[i];j++)
        {
            jpt=nindex[i][j];
            sij=_R0[i]-_R0[jpt]; sij[dir]-=ds; sij.subint();
            dis=sij.norm();
            if(dis<mind)
            {
                mind=dis;
                minj=jpt;
            }
        }
        if(mind>1e-4)
        {
            ERROR("cutpastecn: mind("<<mind<<") > 1e-4");
        }
        _SR[i][dirx]=_VSR[minj][dirx];
        _SR[i][diry]=_VSR[minj][diry];
        _SR[i][dir]=_VSR[minj][dir]+ds;
    }
    SHtoR();
}

void MDFrame::setconfig1()
{
    /* copy current configuration _SR to _SR1 */
//    Free(_SR1);  // Commented out 2007.03.23 Wei Cai
    Realloc(_SR1,Vector3,_NP);
    memcpy(_SR1, _SR, sizeof(Vector3)*_NP);
}

void MDFrame::SR1toSR()
{
    memcpy(_SR, _SR1, sizeof(Vector3)*_NP);
}

void MDFrame::setconfig2()
{
    /* copy current configuration _SR to _SR2 */
//    Free(_SR2);  // Commented out 2007.03.23 Wei Cai
    Realloc(_SR2,Vector3,_NP);
    memcpy(_SR2, _SR, sizeof(Vector3)*_NP);

//    Free(_SR3);  // Commented out 2007.03.23 Wei Cai
//  Commented out 2018.08.09 Wei Cai (should call setconfig3 in script file)
      INFO("MDFrame::setconfig2() For switching simulations, you may need to allocate SR3 first");
//    Realloc(_SR3,Vector3,_NP);
//    memset(_SR3, 0, sizeof(Vector3)*_NP);
}

void MDFrame::setconfig3()
{
    /* copy current configuration _SR to _SR3 */
    Realloc(_SR3,Vector3,_NP);
    memcpy(_SR3, _SR, sizeof(Vector3)*_NP);
}

void MDFrame::copytoconfig1(int ind)
{// 2007.04.27 Wei Cai
    int i;
    for(i=0;i<_NP;i++)
        _SR1[i][ind] = _SR[i][ind];
}

void MDFrame::copytoconfig2(int ind)
{// 2007.04.27 Wei Cai
    int i;
    for(i=0;i<_NP;i++)
        _SR2[i][ind] = _SR[i][ind];
}

void MDFrame::switchconfig()
{
    /* _SR := _SR - SR1 + _SR2 */
    int i;
    if((_SR1!=NULL)&&(_SR2!=NULL))
    {
            for(i=0;i<_NP;i++)
            {
                    _SR[i]-=_SR1[i];
                    _SR[i]+=_SR2[i];
            }
    }
    else
    {
            ERROR("MDFrame::switchconfig _SR1 _SR2 not loaded");
    }
}

void MDFrame::replacefreeatom()
{
    int i;
    Matrix33 h0; 

    for(i=0;i<_NP;i++)
    {
        _R0[i]=_SR[i];
    }
    h0=_H;
    readcn(); 

    for(i=0;i<_NP;i++)
    {
        if(fixed[i])
        {
            _SR[i]=_R0[i];
        }
    }
    _H=h0;

    SHtoR();
}

void MDFrame::relabelatom()
{
    /* input = [ n0 n1 ]
       switch the positions of atom n0 and atom n1 in memory array
    */
    int n0, n1;
    Vector3 s0;
    n0=(int)input[0];
    n1=(int)input[1];
    if(n0==n1) return;
    if((n0<0)||(n0>_NP)) {FATAL("relabelatom n0 out of range");}
    if((n1<0)||(n1>_NP)) {FATAL("relabelatom n0 out of range");}
    s0=_SR[n0];
    if(n1<n0)
        memmove(_SR+n1+1,_SR+n1,sizeof(Vector3)*(n0-n1));
    else
        memmove(_SR+n0,_SR+n0+1,sizeof(Vector3)*(n0-n1));
    _SR[n1]=s0;
}


void MDFrame::moveatom()
{ /*
   * Move specified atoms
   *
   * input = [ n dx dy dz i1 i2 ... in ]
   *
   * n:            number of atoms to be moved (n>0)
   * (dx,dy,dz):   displacement vector in real space
   * i1,i2,...,in: indices of atoms to be moved
   *
   */
    int i, n, id, fix;
    double mx, my, mz;
    double x0, y0, z0, r, r2, dmin, dmax;
    Vector3 s0, ds, dr;

    n=(int)input[0];

    INFO_Printf("n=%d\n",n);

    if(n>0)
    {
        mx=input[1];
        my=input[2];
        mz=input[3];
        
        SHtoR();
        for(i=0;i<n;i++)
        {
            id=(int)input[4+i];
            _R[id].x+=mx;
            _R[id].y+=my;
            _R[id].z+=mz;
        }
        RHtoS();
    }
    else if(n==-2)
    {
        /* see above */
    }
    if(n==0)
    {
        INFO("move atoms in circular plate");
        x0=input[1];
        y0=input[2];
        z0=input[3];
        r =input[4];
        dmin=input[5];
        dmax=input[6];
        mx=input[7];
        my=input[8];
        mz=input[9];
        fix=(int)input[10];
        
        r2=r*r;
        s0.set(x0,y0,z0);
        for(i=0;i<_NP;i++)
        {            
            if((_SR[i].y>=dmax)||(_SR[i].y<dmin)) continue;
            if(fix)
            {
                fixed[i]=1;
                INFO_Printf("fix atom %d in plane\n",i);
            }
            ds=_SR[i]-s0;
            ds.subint();
            
            dr=_H*ds;
            if(dr.norm2()>r2) continue;
            
            _SR[i].x+=mx;
            _SR[i].y+=my;
            _SR[i].z+=mz;
            INFO_Printf("move atom %d\n",i);
        }        
        /* shift box */
        SHtoR();        
    }
    return;
}

void MDFrame::movegroup()
{ /*
   * Move specified groups of atoms
   *
   * input = [ n dx dy dz grp1 grp2 ... grpn ]
   *
   * n:            number of groups to be moved (n>0)
   * (dx,dy,dz):   displacement vector in real space
   * grp1,grp2,...,grpn: indices of group to be moved
   *
   */
    int i, n, id, ip;
    double mx, my, mz;

    n=(int)input[0];

    INFO_Printf("move group: n=%d\n",n);

    if(n>0)
    {
        SHtoR();
        mx=input[1];
        my=input[2];
        mz=input[3];
        
        SHtoR();
        for(i=0;i<n;i++)
        {
            id=(int)input[4+i];
            for(ip=0;ip<_NP;ip++)
            {
                if(group[ip]==id)
                {
                    _R[ip].x+=mx;
                    _R[ip].y+=my;
                    _R[ip].z+=mz;
                }
            }
        }
        RHtoS();
    }
    return;
}

void MDFrame::rotategroup()
{  /*
    * Rotate specified group around a rotation axis at group c.o.m
    *     *
    *         * input = [ axis angle group ]
    *             */
    int i,n,rot_axis,grp_id;
    double rot_angle,totalmass,C,S;
    Vector3 group_com,x;
    Matrix33 rot_mat;
    
    rot_axis = (int)input[0];
    rot_angle = input[1];
    grp_id = (int)input[2];
    C=cos(rot_angle*M_PI/180);S=sin(rot_angle*M_PI/180);
    
    INFO("rotate group:");

    SHtoR(); 
    group_com.clear(); n=0; totalmass=0.0;
    /* find group com */
    if (nspecies==1)
    {
        for(i=0;i<_NP;i++)
        {
            if(group[i]==grp_id)
            {
                group_com+=_R[i]; n++;
            }
        }
        group_com/=n;
    }
    else if (nspecies>1)
    {
        for(i=0;i<_NP;i++)
        {
            if(group[i]==grp_id)
            {
                group_com+=_R[i]*_ATOMMASS[species[i]];
                totalmass+=_ATOMMASS[species[i]];
            }
        }
        group_com/=totalmass;
    }
    INFO("group com ="<<group_com);

    /* set rotation matrix */
    if(rot_axis==1)
    {
        rot_mat.set(1, 0, 0, 
                    0, C, -S,
                    0, S, C); 
    }
    else if(rot_axis==2)
    {
        rot_mat.set(C, 0, S, 
                    0, 1, 0,
                   -S, 0, C); 
    }
    else if(rot_axis==3)
    {
        rot_mat.set(C,-S, 0, 
                    S, C, 0,
                    0, 0, 1); 
    }
    
    /* rotate group */
    for(i=0;i<_NP;i++)
    {
        if(group[i]==grp_id)
        {
            x=_R[i]-group_com;
            x=rot_mat*x; _R[i]=x+group_com;
        }
    }
    RHtoS();
}

void MDFrame::setgroupcomvel()
{ /*
   * set the center of mass velocity for selected groups
   *
   * input = [ n vx vy vz grp1 grp2 ... grpn ]
   *
   * n:            number of groups to be moved (n>0)
   * (vx,vy,vz):   velocity vector in real space
   * grp1,grp2,...,grpn: indices of group to be moved
   *
   */
    int i, n, id, ip;
    double vx, vy, vz, totmass;
    Vector3 vcom; Matrix33 hinv;
    
    n=(int)input[0];

    INFO_Printf("n=%d\n",n);

    if(n>0)
    {
        SHtoR();
        vx=input[1];
        vy=input[2];
        vz=input[3];

        totmass=0;
        for(ip=0;ip<_NP;ip++)
            _VR[ip]=_H*_VSR[ip];
        
        for(i=0;i<n;i++)
        {
            id=(int)input[4+i];
            for(ip=0;ip<_NP;ip++)
            {
                if(fixed[ip]) continue;
                if(group[ip]==id)
                {
                    vcom+=_VR[ip]*_ATOMMASS[species[ip]];
                    totmass+=_ATOMMASS[species[ip]];
                }
            }
        }
        INFO("totmass="<<totmass<<" vcom="<<vcom);
        if(totmass==0)
            WARNING("no atoms have the specified group(s)");
        else
        {
            vcom/=totmass;
            vcom.x-=vx; vcom.y-=vy; vcom.z-=vz;
            for(i=0;i<n;i++)
            {
                id=(int)input[4+i];
                for(ip=0;ip<_NP;ip++)
                {
                    if(fixed[ip]) continue;
                    if(group[ip]==id)
                    {
                        _VR[ip]-=vcom;
                    }
                }
            }
        }
        
        hinv=_H.inv();
        for(ip=0;ip<_NP;ip++)
            _VSR[ip]=hinv*_VR[ip];
    }
    return;
}

void MDFrame::printatoms_in_sphere()
{ /* print the indices of all atoms in a sphere
   *
   * input = [ x0 y0 z0 r ]
   *
   * x0,y0,z0: center of the sphere
   * r:        radius of the sphere
   *
   */
    int i, n; double r2; Vector3 r0, s0, ds, dr; Matrix33 hinv;
    hinv=_H.inv();
    r0.set(input[0],input[1],input[2]);
    r2=input[3]*input[3];
    s0=hinv*r0;
    n=0;
    for(i=0;i<_NP;i++)
    {
         ds=_SR[i]-s0;
         ds.subint();
         dr=_H*ds;
         if(dr.norm2()<r2)
         {
             n++;
         }
    }
    INFO_Printf("\n%d   ",n);
    for(i=0;i<_NP;i++)
    {
         ds=_SR[i]-s0;
         ds.subint();
         dr=_H*ds;
         if(dr.norm2()<r2)
         {
                 INFO_Printf("%d ",i);
         }
    }
    INFO_Printf("\n");
}

void MDFrame::pbcshiftatom()
{
    int i;
    double dx,dy,dz;
    dx=input[0];
    dy=input[1];
    dz=input[2];

    for(i=0;i<_NP;i++)
    {
        _SR[i].x+=dx;_SR[i].y+=dy;_SR[i].z+=dz;
        _SR[i].subint();
    }
    SHtoR();
}




void MDFrame::fixatoms_by_ID()
{ /* fix a set of atoms whose ID are specified
   *
   * input = [ n i1 i2 ... in ]
   *
   * n:             number of atoms to be fixed
   * i1,i2,...,in:  indices of atoms to be fixed
   */
    int i, n;
    n=(int)input[0];
//    for(i=0;i<_NP;i++)
//    {
//        fixed[i]=0;
//    }
    for(i=1;i<=n;i++)
    {
        fixed[(int)input[i]]=1;
    }
    _NPfixed=_NPfree=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==1) _NPfixed++;
        if(fixed[i]==0) _NPfree++;
    }
    INFO_Printf("NPfixed = %d, NPfree = %d\n",_NPfixed,_NPfree);
}

void MDFrame::fixatoms_by_position()
{ /* fix all atoms whose position falls within a specified regime
   *
   * input = [ enable x0 x1 y0 y1 z0 z1 ]
   *
   * enable: set to 1 to activate
   * atoms whose reduced coordinate (x,y,z) satisfie
   *   x0 <= x < x1
   *   y0 <= y < y1
   *   z0 <= z < z1
   * will be fixed
   *
   * if enable = 2, atoms within cylindrical volume defined by
   * input = [enable ind x0 y0 z0 ra rb h0 h1] 
   * will be fixed. (ra and rb given in real coordinate  (obsolete: reduced coordinate)
   * If additional inputs h0 and h1 are defined in scaled coordinate,
   * the cylinder height is defined as height = h1 - h0.
   */
    int i, enable, nfixed, ind, inda, indb;
    double xmin, xmax, ymin, ymax, zmin, zmax, hmin, hmax;
    double x0, y0, z0, rmin, rmax, ra0, rb0, sqrr;
    Vector3 dr, ds;
    
    enable=(int)input[0];
    nfixed=0;

    if(enable == 1)
    {
        xmin=input[1];
        xmax=input[2];
        ymin=input[3];
        ymax=input[4];
        zmin=input[5];
        zmax=input[6];
        
        //INFO_Printf("setfixedatoms = [\n");
        nfixed=0;
        for(i=0;i<_NP;i++)
        {
            if((_SR[i].x>=xmin)&&(_SR[i].x<xmax)
               &&(_SR[i].y>=ymin)&&(_SR[i].y<ymax)
               &&(_SR[i].z>=zmin)&&(_SR[i].z<zmax))
            {
                fixed[i]=1;
                nfixed++;
            }
        }
        //INFO_Printf("%d\n",nfixed);
        //for(i=0;i<_NP;i++)
        //{
        //    if((_SR[i].x>=xmin)&&(_SR[i].x<xmax)
        //       &&(_SR[i].y>=ymin)&&(_SR[i].y<ymax)
        //       &&(_SR[i].z>=zmin)&&(_SR[i].z<zmax))
        //    {
        //        //INFO_Printf("%d ",i);
        //    }
        //}
    }
    else if(enable == 2) // Jan 25 2007 Keonwook Kang 
    {
        ind = (int)input[1]; ind--; /* cylindrical axis 1:x 2:y 3:z */
        x0 = input[2];  /* center of cylindrical axis */
        y0 = input[3]; 
        z0 = input[4];
        rmin = input[5]; // in real coordinate
        rmax = input[6]; // in real coordinate
        hmin = input[7]; // in scaled coordinate
        hmax = input[8]; // in scaled coordinate

        /* the other two directions orthogonal to ind */
        inda = (ind+1)%3; indb = (inda+1)%3;
        ra0 = input[(inda+2)]; rb0 = input[(indb+2)];

        switch (ind)
        {
        case 0: INFO_Printf("Cylindrical axis: x\n"); break;
        case 1: INFO_Printf("Cylindrical axis: y\n"); break;
        case 2: INFO_Printf("Cylindrical axis: z\n"); break;           
        }
        INFO_Printf("Atoms whose radial distance between "
                    "%f and %f will be fixed.\n", rmin, rmax);
        //INFO_Printf("setfixedatoms = [\n");
        
        for(i=0;i<_NP;i++)
        {
            ds = _SR[i];
            /* ds[inda]-=ra0; ds[indb]-=rb0; 
            sqrra = (_SR[i][inda] - ra0)*(_SR[i][inda] - ra0);
            sqrrb = (_SR[i][indb] - rb0)*(_SR[i][indb] - rb0);
            sqrr = sqrra + sqrrb; */
            
            ds[0]-=x0; ds[1]-=y0; ds[2]-=z0;
            dr = _H*ds;
            sqrr = dr[inda]*dr[inda] + dr[indb]*dr[indb];

//            if (i >= 0 && i<10)
//            {
//                INFO_Printf("_SR[%d][%d] = %f\n",i,inda,_SR[i][inda]);
//                INFO_Printf("sqrra = %f, sqrrb = %f, "
//                            "SQR(rmin)=%f, SQR(rmax)=%f\n",sqrra,sqrrb,rmin*rmin, rmax*rmax);
//            }
            if ((sqrr >= rmin*rmin) && (sqrr < rmax*rmax))
            {
                if (hmax > hmin)
                {
                    if ((_SR[i][ind] < hmin) || (_SR[i][ind] > hmax))
                        continue;
                }
                //INFO_Printf("%d ",i);
                fixed[i]=1;
                nfixed++;                
            }
        }
    }    
    //INFO_Printf("\n]\n");
    INFO_Printf("%d atoms fixed\n",nfixed);
    _NPfixed=_NPfree=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==1) _NPfixed++;
        if(fixed[i]==0) _NPfree++;
    }
    INFO_Printf("NPfixed = %d, NPfree = %d\n",_NPfixed,_NPfree);
}

void MDFrame::fixatoms_by_group(double *myinput)
{
    /* fix atoms that belong to specified groups
     *
     * input = [NgrpID, grpID_1 grpID_2 ...]
     *
     * NgrpID : number of groups to be fixed
     * grpID_1, grpID_2, ... : group ID's 
     */
    int i, j, NgrpID, nfixed;
    NgrpID = (int) myinput[0];

    INFO("Number of groups to be fixed: "<<NgrpID);

    if (NgrpID <= 0 )
    {
        ERROR("Number of groups to be fixed should be greater than zero.");
        return;
    }
    
    INFO_Printf("fixedatomID = [\n");
    nfixed = 0;
    for (i=0;i<_NP;i++)
    {
        if (fixed[i]==1) /* If already fixed, skip */
            continue;
        else
        {
            for (j=0; j<NgrpID; j++)
            {
                if (group[i]==(int) myinput[j+1])
                {
                    fixed[i] = 1; nfixed++;
                    INFO_Printf("%d ",i);
                    break;
                }
            }
        }
    }
    INFO_Printf("\n]\n");
    INFO("Number of fixed atoms: "<<nfixed);
    _NPfixed=_NPfree=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==1) _NPfixed++;
        if(fixed[i]==0) _NPfree++;
    }
    INFO_Printf("NPfixed = %d, NPfree = %d\n",_NPfixed,_NPfree);
}

void MDFrame::RtoR0_by_group(double *myinput)
{ 
    /* copy R to R0 for atoms that belong to specified groups
     *
     * input = [NgrpID, grpID_1 grpID_2 ...]
     *
     * NgrpID : number of groups to be fixed
     * grpID_1, grpID_2, ... : group ID's 
     */
    int i, j, NgrpID, ncopied;
    NgrpID = (int) myinput[0];

    INFO("Number of groups to be fixed: "<<NgrpID);

    if (NgrpID <= 0 )
    {
        ERROR("Number of groups to be fixed should be greater than zero.");
        return;
    }
   
    SHtoR(); 
    ncopied = 0;
    for (i=0;i<_NP;i++)
    {
        for (j=0; j<NgrpID; j++)
        {
            if (group[i]==(int) myinput[j+1])
            {
                _R0[i] = _R[i];
                INFO_Printf("atom %d copied R to R0 (%e,%e,%e) group=%d\n",
                            i,_R0[i].x,_R0[i].y,_R0[i].z,group[i]);
                ncopied ++;
                break;
            }
        }
    }
    INFO("Number of copied atoms: "<<ncopied);
}

void MDFrame::R0toR_by_group(double *myinput)
{ 
    /* copy R to R0 for atoms that belong to specified groups
     *
     * input = [NgrpID, grpID_1 grpID_2 ...]
     *
     * NgrpID : number of groups to be fixed
     * grpID_1, grpID_2, ... : group ID's 
     */
    int i, j, NgrpID, ncopied;
    NgrpID = (int) myinput[0];

    INFO("Number of groups to be fixed: "<<NgrpID);

    if (NgrpID <= 0 )
    {
        ERROR("Number of groups to be fixed should be greater than zero.");
        return;
    }
   
    ncopied = 0;
    for (i=0;i<_NP;i++)
    {
        for (j=0; j<NgrpID; j++)
        {
            if (group[i]==(int) myinput[j+1])
            {
                _R[i] = _R0[i];
                INFO_Printf("atom %d copied R0 to R (%e,%e,%e) group=%d\n",
                            i,_R[i].x,_R[i].y,_R[i].z,group[i]);
                ncopied ++;
                break;
            }
        }
    }
    INFO("Number of copied atoms: "<<ncopied);
    RHtoS();
}

void MDFrame::SR1toSR_by_group(double *myinput)
{ 
    /* copy SR1 to SR for atoms that belong to specified groups
     *
     * input = [NgrpID, grpID_1 grpID_2 ...]
     *
     * NgrpID : number of groups to be fixed
     * grpID_1, grpID_2, ... : group ID's 
     */
    int i, j, NgrpID, ncopied;
    NgrpID = (int) myinput[0];

    INFO("Number of groups to be fixed: "<<NgrpID);

    if (NgrpID <= 0 )
    {
        ERROR("Number of groups to be fixed should be greater than zero.");
        return;
    }
   
    ncopied = 0;
    for (i=0;i<_NP;i++)
    {
        for (j=0; j<NgrpID; j++)
        {
            if (group[i]==(int) myinput[j+1])
            {
                _SR[i] = _SR1[i];
                INFO_Printf("atom %d copied SR1 to SR (%e,%e,%e) group=%d\n",
                            i,_SR[i].x,_SR[i].y,_SR[i].z,group[i]);
                ncopied ++;
                break;
            }
        }
    }
    INFO("Number of copied atoms: "<<ncopied);
}

void MDFrame::SR2toSR_by_group(double *myinput)
{ 
    /* copy SR2 to SR for atoms that belong to specified groups
     *
     * input = [NgrpID, grpID_1 grpID_2 ...]
     *
     * NgrpID : number of groups to be fixed
     * grpID_1, grpID_2, ... : group ID's 
     */
    int i, j, NgrpID, ncopied;
    NgrpID = (int) myinput[0];

    INFO("Number of groups to be fixed: "<<NgrpID);

    if (NgrpID <= 0 )
    {
        ERROR("Number of groups to be fixed should be greater than zero.");
        return;
    }
   
    ncopied = 0;
    for (i=0;i<_NP;i++)
    {
        for (j=0; j<NgrpID; j++)
        {
            if (group[i]==(int) myinput[j+1])
            {
                _SR[i] = _SR2[i];
                INFO_Printf("atom %d copied SR2 to SR (%e,%e,%e) group=%d\n",
                            i,_SR[i].x,_SR[i].y,_SR[i].z,group[i]);
                ncopied ++;
                break;
            }
        }
    }
    INFO("Number of copied atoms: "<<ncopied);
}

void MDFrame::fixatoms_by_species(double *myinput)
{
    /* Fix atoms that belong to specified species.
     * by Keonwook Kang   Oct/11/2007
     * input = [N_species, species_1 species_2 ...]
     *
     * N_species : number of species to be fixed
     * species_1, species_2, ... : species ID's 
     */
    int i, j, N_species, nfixed;
    N_species = (int) myinput[0];

    INFO("Number of species to be fixed: "<<N_species);

    if (N_species <= 0 )
    {
        ERROR("Number of species to be fixed should be greater than zero.");
        return;
    }
    
    INFO_Printf("fixedatomID = [\n");
    nfixed = 0;
    for (i=0;i<_NP;i++)
    {
        if (fixed[i]==1) /* If already fixed, skip */
            continue;
        else
        {
            for (j=0; j<N_species; j++)
            {
                if (species[i]==(int) myinput[j+1])
                {
                    fixed[i] = 1; nfixed++;
                    INFO_Printf("%d ",i);
                    break;
                }
            }
        }
    }
    INFO_Printf("\n]\n");
    INFO("Number of fixed atoms: "<<nfixed);
    _NPfixed=_NPfree=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==1) _NPfixed++;
        if(fixed[i]==0) _NPfree++;
    }
    INFO_Printf("NPfixed = %d, NPfree = %d\n",_NPfixed,_NPfree);
}

void MDFrame::fixatoms_by_pos_topol()
{ /* fix all atoms whose position falls within a specified regime
   *
   * input = [ enable x0 x1 y0 y1 z0 z1 tmin tmax ]
   *
   * enable: set to 1 to activate
   * atoms whose reduced coordinate (x,y,z) satisfies
   *   x0 <= x < x1
   *   y0 <= y < y1
   *   z0 <= z < z1
   * and whose _TOPOL value satisfies
   *   tmin <= topol < tmax
   * will be fixed
   */
    int i, n, nfixed, nprint;
    double xmin, xmax, ymin, ymax, zmin, zmax, tmin, tmax;
    n=(int)input[0];
    xmin=input[1];
    xmax=input[2];
    ymin=input[3];
    ymax=input[4];
    zmin=input[5];
    zmax=input[6];
    tmin=input[7];
    tmax=input[8];

    INFO_Printf("setfixedatoms = [\n");
    nfixed=0;
    for(i=0;i<_NP;i++)
    {
        //INFO_Printf("TOPOL[%d]=%f\n",i,_TOPOL[i]);
        if((_SR[i].x>=xmin)&&(_SR[i].x<xmax)
           &&(_SR[i].y>=ymin)&&(_SR[i].y<ymax)
           &&(_SR[i].z>=zmin)&&(_SR[i].z<zmax)
           &&(_TOPOL[i]>=tmin)&&(_TOPOL[i]<tmax))
        {
            fixed[i]=1;
            nfixed++;
        }
    }
    INFO_Printf("%d\n",nfixed);
    nprint = 0;
    for(i=0;i<_NP;i++)
    {
        if((_SR[i].x>=xmin)&&(_SR[i].x<xmax)
           &&(_SR[i].y>=ymin)&&(_SR[i].y<ymax)
           &&(_SR[i].z>=zmin)&&(_SR[i].z<zmax)
           &&(_TOPOL[i]>=tmin)&&(_TOPOL[i]<tmax))
        {
            nprint ++;
            if (nprint <= 500)
            {
                INFO_Printf("%d ",i);
            } else
            {
                INFO_Printf("... ");
                break;
            }
        }
    }
    INFO_Printf("\n]\n");
    INFO_Printf("%d atoms fixed\n",nfixed);
    _NPfixed=_NPfree=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==1) _NPfixed++;
        if(fixed[i]==0) _NPfree++;
    }
    INFO_Printf("NPfixed = %d, NPfree = %d\n",_NPfixed,_NPfree);
}

void MDFrame::fixallatoms()
{
    int i;

    for(i=0;i<_NP;i++)
    {
        fixed[i]=1;
    }
    _NPfixed=_NPfree=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==1) _NPfixed++;
        if(fixed[i]==0) _NPfree++;
    }
    INFO_Printf("NPfixed = %d, NPfree = %d\n",_NPfixed,_NPfree);
}

void MDFrame::freeallatoms()
{
    int i;

    for(i=0;i<_NP;i++)
    {
        fixed[i]=0;
    }
    _NPfixed=_NPfree=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==1) _NPfixed++;
        if(fixed[i]==0) _NPfree++;
    }
    INFO_Printf("NPfixed = %d, NPfree = %d\n",_NPfixed,_NPfree);
}

void MDFrame::reversefixedatoms()
{
    int i;

    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==0) fixed[i]=1;
        else if(fixed[i]==1) fixed[i]=0;
    }
    _NPfixed=_NPfree=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==1) _NPfixed++;
        if(fixed[i]==0) _NPfree++;
    }
    INFO_Printf("NPfixed = %d, NPfree = %d\n",_NPfixed,_NPfree);
}

void MDFrame::constrain_fixedatoms()
{
    int i, nc;

    nc=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==1)
        {
            nc++;
            if (nc>=MAXCONSTRAINATOMS)
            {
                ERROR("MAXCONSTRAINATOMS exceeded.  You need to change the defined value in md.h");
                continue;
            }
            constrainatoms[nc]=i;
        }
    }
    constrainatoms[0]=nc;
}

void MDFrame::fix_constrainedatoms()
{
    int i, nc;

    nc=constrainatoms[0];
    for(i=1;i<=nc;i++)
    {
        fixed[constrainatoms[i]]=1;
    }
    INFO_Printf("%d atoms fixed\n",nc);
    _NPfixed=_NPfree=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==1) _NPfixed++;
        if(fixed[i]==0) _NPfree++;
    }
    INFO_Printf("NPfixed = %d, NPfree = %d\n",_NPfixed,_NPfree);
}

void MDFrame::fix_imageatoms()
{
    int i;

    for(i=0;i<=_NP;i++)
    {
        if(image[i]>=0) fixed[i] = 1;
    }
    _NPfixed=_NPfree=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==1) _NPfixed++;
        if(fixed[i]==0) _NPfree++;
    }
    INFO_Printf("NPfixed = %d, NPfree = %d\n",_NPfixed,_NPfree);
}

void MDFrame::setfixedatomsspecies()
{ /* set the species value of all fixed atoms to input[0]
   */
    int i;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==1)
            species[i]=(int)input[0];
    }
}
       
void MDFrame::setfixedatomsgroup()
{ /* set the group ID of all fixed atoms to input[0]
   */
    int i;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==1)
            group[i]=(int)input[0];
    }
}
       
void MDFrame::reversespecies()
{ /* set the species value of all fixed atoms to input[0]
   */
    int i;
    for(i=0;i<_NP;i++)
    {
        if(species[i]==0)
            species[i]=1;
        else if(species[i]==1)
            species[i]=0;        
    }
}
       
void MDFrame::removefixedatoms()
{ /* remove all atom i if fixed[i] is not equal to zero */
    int i, j;
    
    j=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==0)
        {
            _SR[j]=_SR[i];
            fixed[j]=fixed[i];
            species[j]=species[i];
            _EPOT_IND[j]=_EPOT_IND[i];
            _TOPOL[j]=_TOPOL[i];
            j++;
        }
    }
    _NP=j;
    SHtoR();
    _NPfixed=_NPfree=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==1) _NPfixed++;
        if(fixed[i]==0) _NPfree++;
    }
    INFO_Printf("NPfixed = %d, NPfree = %d\n",_NPfixed,_NPfree);
}

void MDFrame::markremovefixedatoms()
{ /* if fixed[i]==1, change it to fixed[i]=-1 */
    int i;
    
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==1)
        {
            fixed[i]=-1;
        }
    }
    _NPfixed=_NPfree=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==1) _NPfixed++;
        if(fixed[i]==0) _NPfree++;
    }
    INFO_Printf("NPfixed = %d, NPfree = %d\n",_NPfixed,_NPfree);
}

void MDFrame::movefixedatoms()
{ /* remove all atom i if fixed[i] is not equal to zero */
    int i;
    
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]==1)
        {
            _SR[i].x+=input[0];
            _SR[i].y+=input[1];
            _SR[i].z+=input[2];
        }
    }
}

void MDFrame::makeellipsoid()
{
    /* mark or remove all atoms within an ellipsoid
     *
     * input = [ enable x0 y0 z0 a b c flag ]
     *
     * enable:     if 1 then removeellipsoid is enabled
     * (x0,y0,z0): center of the ellipsoid
     * a,b,c:      semi-axis of the ellipsoid
     * flag:       if 0 then atoms inside are removed
     *             if 1 then atoms inside are fixed
     */
    
    int i, n, removenum, flag;
    double x0, y0, z0, arem, brem, crem;
    Vector3 r0, dr;
    
    if(input[0]==0)
    {
        INFO_Printf("input = [%f %f %f ]\n",
                    input[0],input[1],input[2]);
        ERROR("removeellipsoid(): no geometry set");
        return;
    }

    /* read in parameters */
    x0=input[1];
    y0=input[2];
    z0=input[3];
    arem=input[4];
    brem=input[5];
    crem=input[6];
    flag=(int)input[7];
    
    r0.set(x0,y0,z0);
    /* remove all atoms outside rrem */
    
    SHtoR();
    for(i=0;i<_NP;i++)
    {
        dr = _R[i]-r0;
        dr.x/=arem;
        dr.y/=brem;
        dr.z/=crem;
        if(dr.norm2()<=1)
        {
           if (flag==0)
              fixed[i]=-1;
           else if (flag==1)
              fixed[i]=1;
           else
              ERROR("Unknow flag"<<flag);
        }
    }
    
    /* commit removing atoms */
    n=0; SHtoR();
    removenum=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]!=-1)
        {
            _R0[n]=_R[i];
            fixed[n]=fixed[i];
            n++;
        }
        else
        {
            removenum++;
        }
    }
    INFO("originally "<<_NP<<" atoms");
    INFO("remove "<<removenum<<" atoms");
    _NP-=removenum;
    INFO("now "<<_NP<<" atoms");
    R0toR();  RHtoS();
    INFO("NP="<<_NP<<"  n="<<n);
}
void MDFrame::removeellipsoid()
{
    /* remove all atoms within an ellipsoid
     *
     * input = [ enable x0 y0 z0 a b c ]
     *
     * enable:     if 1 then removeellipsoid is enabled
     * (x0,y0,z0): center of the ellipsoid
     * a,b,c:      semi-axis of the ellipsoid
     */
    
    int i, n, removenum;
    double x0, y0, z0, arem, brem, crem;
    Vector3 r0, dr;
    
    if(input[0]==0)
    {
        INFO_Printf("input = [%f %f %f ]\n",
                    input[0],input[1],input[2]);
        ERROR("removeellipsoid(): no geometry set");
        return;
    }

    /* read in parameters */
    x0=input[1];
    y0=input[2];
    z0=input[3];
    arem=input[4];
    brem=input[5];
    crem=input[6];
    
    r0.set(x0,y0,z0);
    /* remove all atoms outside rrem */
    
    SHtoR();
    for(i=0;i<_NP;i++)
    {
        dr = _R[i]-r0;
        dr.x/=arem;
        dr.y/=brem;
        dr.z/=crem;
        if(dr.norm2()<=1)
            fixed[i]=-1;
    }
    
    /* commit removing atoms */
    n=0; SHtoR();
    removenum=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]!=-1)
        {
            _R0[n]=_R[i];
            fixed[n]=fixed[i];
            n++;
        }
        else
        {
            removenum++;
        }
    }
    INFO("originally "<<_NP<<" atoms");
    INFO("remove "<<removenum<<" atoms");
    _NP-=removenum;
    INFO("now "<<_NP<<" atoms");
    R0toR();  RHtoS();
    INFO("NP="<<_NP<<"  n="<<n);
}

void MDFrame::removerectbox()
{ /* Remove a rectangular box of atoms
   *
   * input = [ ax ay az La bx by bz Lb cx cy cz Lc x0 y0 z0 plotonly outside ]
   *
   */
    int i, j, assist, plotonly, outside;
    Vector3 ds, dss, sa, sb, sc, s0;
    Vector3 ha, hb, hc;
    Matrix33 sh, shinv, latt, lattinv;
    
    ha.set(latticesize[0][0],latticesize[0][1],latticesize[0][2]);
    hb.set(latticesize[1][0],latticesize[1][1],latticesize[1][2]);
    hc.set(latticesize[2][0],latticesize[2][1],latticesize[2][2]);
    ha*=latticesize[0][3];
    hb*=latticesize[1][3];
    hc*=latticesize[2][3];
    latt.setcol(ha,hb,hc);
    lattinv=latt.inv();
        
    sa.set(input[0],input[1],input[2]);
    sb.set(input[4],input[5],input[6]);
    sc.set(input[8],input[9],input[10]);
    sa*=input[3];
    sb*=input[7];
    sc*=input[11];
    s0.set(input[12],input[13],input[14]);
    assist=(int)input[15];
    plotonly=(int)input[16];
    outside=(int)input[17];

    sa=lattinv*sa;
    sb=lattinv*sb;
    sc=lattinv*sc;
    s0=lattinv*s0;
    
    INFO("assist="<<assist);
    INFO("plotonly="<<plotonly);
    
    sh.setcol(sa,sb,sc);
    shinv=sh.inv();
    INFO("vol="<<sh.det());
    INFO("s0="<<s0);
    
    for(i=0;i<_NP;i++)
    {
        ds=_SR[i]-s0;
        /*ds.subint();*/
        dss=shinv*ds;
        if(assist&&(!plotonly))
        {
            if((dss.x>=-0.5)&&(dss.x<0.5)&&
               (dss.y>=-0.5)&&(dss.y<0.5))
            {
                if(fabs(dss.z)<2)
                {
                    if((dss.z>=-0.5)&&(dss.z<0.5))
                    {
                        fixed[i]=-1;/* mark for deletion */
                    }
                    else if(dss.z>0.5)
                    {
                        dss.z=(dss.z-0.5)/(2-0.5)*(2-0.1)+0.1-dss.z;
                        dss.x=0; dss.y=0;
                        INFO("dss="<<dss);
                        _SR[i]+=sh*dss;
                    }
                    else
                    {
                        dss.z=(dss.z+0.5)/(-2+0.5)*(-2+0.1)-0.1-dss.z;
                        dss.x=0; dss.y=0;
                        INFO("dss="<<dss);
                        _SR[i]+=sh*dss;
                    }
                }
            }
        }
        else
        {
            if((dss.x>=-0.5)&&(dss.x<0.5)&&
               (dss.y>=-0.5)&&(dss.y<0.5)&&
               (dss.z>=-0.5)&&(dss.z<0.5))
            {
                fixed[i]=-1;/* mark for deletion */
            }
        }
    }
    if(outside)
    {
        for(i=0;i<_NP;i++)
        {
            if(fixed[i]==-1)
            {
                fixed[i]=0;
            }
            else if(fixed[i]==0)
            {
                fixed[i]=-1;
            }
        }
    }
    if(!plotonly)
    {
        j=0;
        for(i=0;i<_NP;i++)
        {
            if(fixed[i]!=-1)
            {
                _SR[j]=_SR[i];
                fixed[j]=0;
                if(_SAVEMEMORY<8)
                        _EPOT_IND[j]=_EPOT_IND[i];
                if(_SAVEMEMORY<6)
                        _TOPOL[j]=_TOPOL[i];
                j++;
            }
        }
        INFO("removerectbox: "<<_NP-j<<" of atoms removed");
        _NP=j;
        INFO("new number of atoms "<<_NP);
        if(_SAVEMEMORY<9)
                SHtoR();
    }
}

void MDFrame::removeoverlapatoms()
{/* remove atoms that are closer than rc
  *
  * input = [ rc ]
  *
  */
    int i, j, jpt;
    double removerc, rc2;
    Vector3 ds, dr;

    refreshneighborlist();
    
    removerc = input[0];
    
    rc2=removerc*removerc;
    for(i=0;i<_NP;i++)
    {
        for(j=0;j<nn[i];j++)
        {
            jpt=nindex[i][j];
            ds=_SR[i]-_SR[jpt];
            ds.subint();
            dr=_H*ds;
            if(dr.norm2()<rc2)
            {
                fixed[i]=-1;/* mark for deletion */
                if(fixed[jpt]==-1) fixed[i]=0;
            }
        }
    }
    j=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]!=-1)
        {
            _SR[j]=_SR[i];
            fixed[j]=0;
            j++;
        }
    }
    _NP=j;
    SHtoR();
}

void MDFrame::findcenterofmass(class Vector3 *sr)
{ /* Find center of mass _COM of a series of vector quantity sr. */
  /* sr can be either scaled or real coordinate.                 */
  /*                          May 28 Keonwook Kang               */
    int i;
    double totalmass;
    
    _COM.clear(); totalmass=0;
    if(nspecies==1)
    {
        for(i=0;i<_NP;i++)
            _COM+=sr[i];
        _COM/=_NP;
    }
    else if (nspecies > 1)
    {
        for(i=0;i<_NP;i++)
        {
            _COM+=sr[i]*_ATOMMASS[species[i]];
            totalmass+=_ATOMMASS[species[i]];
        }
        _COM/=totalmass;
    }
    else
        ERROR("Number of species("<<nspecies<<") must be positive integer!!");

    INFO("_COM = ("<<_COM<<")");
}

void MDFrame::translate_com()
{ /* apply rigid-body translation so that COM matches that of config1 */
    int i;
    Vector3 com, com1, dcom;

    /* potential bug: current version of conjugate gradient relaxation
                      will not conserve COM correctly if there are
                      multiple species
       suggested fix: rescale coordinates by mass in potential_wrapper                
    */
    
    /* First find center of mass for config1 */
    com.clear(); 
    com1.clear();

    INFO("In config1, ");
    findcenterofmass(_SR1); com1=_COM;
    INFO("In current config, ");
    findcenterofmass(_SR); com=_COM;
    
    dcom = com1 - com;
    INFO("dcom = ("<<dcom<<")");

    for(i=0;i<_NP;i++)
        _SR[i]+=dcom;
}


void MDFrame::rotate_com()
{ /* apply rigid-body rotation so that COM matches that of config1 */
    int index, i;
    Vector3 r, r1, dr, ds;
    Vector3 com, com1, dcom;
    double J, R2, alpha;
    Matrix33 hinv;

    /* input = [ 1/2/3 ] for x/y/z */
    index = (int) input[0] - 1;

    if(index<0)
        FATAL("rotate_com: index ("<<index<<") is not valid!");

    if(index!=2)
        FATAL("rotate_com: index ("<<index<<") is not valid!");

    /* rotation around z axis */
    
    /* First find center of mass for config1 */
    com.clear(); com1.clear();
    
    if(nspecies==1)
    {
        INFO("In config1, ");
        findcenterofmass(_SR1); com1=_COM;
        INFO("In current config, ");
        findcenterofmass(_SR); com=_COM;        
    }
    else
    {
        FATAL("rotate_com not implemented for multiple species yet");
    }

    J = 0; R2 = 0;
    for(i=0;i<_NP;i++)
    {
        r  = _H*(_SR[i]-com);
        r1 = _H*(_SR1[i]-com1);

        J += r.x*r1.y - r.y*r1.x;
        R2 += r.y*r1.y + r.x*r1.x;
    }
    alpha = atan2(J,R2);

    INFO_Printf("amount of rotation needed: alpha = %20.12e (degree)\n",alpha*180/M_PI);

    hinv = _H.inv();
    for(i=0;i<_NP;i++)
    {
        r  = _H*(_SR[i]-com);
        dr.x = r.x*cos(alpha) - r.y*sin(alpha);
        dr.y = r.x*sin(alpha) + r.y*cos(alpha);
        dr.z = r.z;
        ds = hinv*dr;
        _SR[i]=(ds+com);
    }
}

void MDFrame::clearFext()
{ /* clear _Fext array */
    int i;
    if(_Fext!=NULL) 
    {
        for(i=0;i<_NP;i++)
            _Fext[i].clear();
    }
}

void MDFrame::addFext_to_group()
{ /* add Fext to atoms belonging to specific group */
    int i, groupID;
    double fx, fy, fz;

    if(_Fext==NULL)
    {
        ERROR("addFext_to_group: _Fext == NULL !");
        return;
    }
    groupID = (int) input[0];
    fx = input[1];
    fy = input[2];
    fz = input[3];
    for(i=0;i<_NP;i++)
    {
       if(group[i]==groupID)
       {
          _Fext[i].x = fx;
          _Fext[i].y = fy;
          _Fext[i].z = fz;
       }
    } 
}


#ifdef _GSL
#include "xrd.cpp"
#endif

/********************************************************************/
/* Neighbor list */
/********************************************************************/
void MDFrame::NbrList_reconstruct(int iatom)
{
//    NbrList_reconstruct_use_cell_list(iatom);
    NbrList_reconstruct_use_link_list(iatom);
}

void MDFrame::NbrList_reconstruct_use_link_list(int iatom)
{
    /* use link list combined with Verlet list */
    /* if iatom == -1, refresh neighbor list of all atoms
     * otherwise only refresh neighbor list of atom i
     */
    int i, j, k, cell_ind, ncell, ineigh, ncx, ncy, ncz, itmp, skip;
    int head1, head2, ipt, jpt, i1, j1, k1, i2, j2, k2;
    int max_nnm, num_neigh_cell, neigh_cell_ind[27];
    double rsmax2;
    Vector3 s,sij,rij,h;

    //printf("NbrList_reconstruct_use_link_list iatom=%d",iatom);
    if(iatom==-1)
        NbrList_init(_NP*allocmultiple,_NNM);
    
    rsmax2=_RLIST*_RLIST;
    h=_H.height();

    /* number of cells in x, y, z directions */
    ncx=(int)floor(h[0]/((_RLIST)*1.05));if(ncx==0)ncx=1;
    ncy=(int)floor(h[1]/((_RLIST)*1.05));if(ncy==0)ncy=1;
    ncz=(int)floor(h[2]/((_RLIST)*1.05));if(ncz==0)ncz=1;

    if((ncx==1)||(ncy==1)||(ncz==1)) {
       ERROR(" ");
       ERROR("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
       ERROR("!");
       ERROR("! supercell size too small, results may be incorrect!");
       ERROR("! ncx ("<<ncx<<"), ncy ("<<ncy<<") and ncz ("<<ncz<<") should be at least 2,");
       ERROR("! unless your potential() accounts for PBC images.");
       ERROR("!");
       ERROR("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
       ERROR(" ");
    }

//    INFO("reconstruct neighborlist _RLIST="<<_RLIST<<" ncell="
//         <<ncx<<"*"<<ncy<<"*"<<ncz);

    ncell = ncx*ncy*ncz;
    NbrList_initlinklist(ncell);

    /* sort atoms to link list */
    for(ipt=0;ipt<_NP;ipt++)
    {
        s=_SR[ipt]; s.subint();
        i=(int)floor((s.x+0.5)*ncx);i=(i+2*ncx)%ncx;
        j=(int)floor((s.y+0.5)*ncy);j=(j+2*ncy)%ncy;
        k=(int)floor((s.z+0.5)*ncz);k=(k+2*ncz)%ncz;

        cell_ind = (i*ncy + j)*ncz + k;

        /* insert atom ipt to the link list of cell_ind */
        link_list[ipt] = link_head[cell_ind];
        link_head[cell_ind] = ipt;
    }

    /* print link list */
//    for(cell_ind=0;cell_ind<ncell;cell_ind++)
//    {
//        /* head of cell_ind */
//        INFO_Printf("cell_ind(%d) = (",cell_ind);
//        ipt = link_head[cell_ind];
//        while(ipt!=-1)
//        {
//            INFO_Printf(" %d",ipt);
//            ipt = link_list[ipt];
//        }
//        INFO_Printf(")\n");        
//    }
    
    /* reset Verlet list */
    if(iatom==-1)
    {
        for(ipt=0;ipt<_NP;ipt++)
        {
            nn[ipt]=0;
            for(itmp=0;itmp<_NNM;itmp++)
                nindex[ipt][itmp]=0;
        }
    }
    else
    {
        ipt=iatom;
        nn[ipt]=0;
        for(itmp=0;itmp<_NNM;itmp++)
            nindex[ipt][itmp]=0;
    }        
    max_nnm = 0;

    if(iatom==-1)
    {
    for(cell_ind=0;cell_ind<ncell;cell_ind++)
    {
        /* generate a list of neighboring cells of cell_ind
         * self-included (at the beginning of the list)
         */
        k = cell_ind % ncz;
        j = ((cell_ind - k)/ncz) % ncy;
        i = (((cell_ind - k)/ncz) - j) / ncy;

        neigh_cell_ind[0] = cell_ind; /* self at the top of list */
        num_neigh_cell = 1;
        for(i1=((ncx==1)?i:(i-1));i1<=((ncx<=2)?i:(i+1));i1++)
            for(j1=((ncy==1)?j:(j-1));j1<=((ncy<=2)?j:(j+1));j1++)
                for(k1=((ncz==1)?k:(k-1));k1<=((ncz<=2)?k:(k+1));k1++)
                {
                    i2=(i1+2*ncx)%ncx;
                    j2=(j1+2*ncy)%ncy;
                    k2=(k1+2*ncz)%ncz;
                    if((i==i2)&&(j==j2)&&(k==k2)) continue; /* self already included */
                    neigh_cell_ind[num_neigh_cell] = (i2*ncy + j2)*ncz + k2;
                    num_neigh_cell ++;                    
                }
//        INFO_Printf("cell_ind (%d) (%d %d %d)\n",cell_ind,i,j,k);        
//        INFO_Printf("num_neigh_cell = %d (",num_neigh_cell);
//        for(ineigh=0;ineigh<num_neigh_cell;ineigh++)
//            INFO_Printf(" %d",neigh_cell_ind[ineigh]);
//        INFO_Printf(")\n");

        /* head of cell_ind */
        head1 = link_head[cell_ind];

        max_nnm = 0;
        /* loop through all neighbor cells */
        for(ineigh=0;ineigh<num_neigh_cell;ineigh++)
        {
            head2 = link_head[neigh_cell_ind[ineigh]]; 
        
            /* going through two lists with head1 and head2 */
            ipt = head1;
            while(ipt!=-1)
            {
                jpt = head2;
                while(jpt!=-1)
                {
//                    INFO_Printf("cell(%d) atom(%d) - cell(%d) atom(%d)\n",
//                                cell_ind,ipt,neigh_cell_ind[ineigh],jpt);
                    /* check whether pair ipt, jpt is in the list to skip */
                    skip = 0;
                    if (nl_skip_pairs[0]>0)
                    {
                        for(itmp=0;itmp<nl_skip_pairs[0];itmp++)
                        {
    //printf("in NbrList_recons_use_list 4.32, itmp = %d\n", itmp);
                            if(   (nl_skip_pairs[2*itmp+1] == ipt)
                                &&(nl_skip_pairs[2*itmp+2] == jpt) )
                                skip = 1;
    //printf("in NbrList_recons_use_list 4.33\n");
                            if(   (nl_skip_pairs[2*itmp+1] == jpt)
                                &&(nl_skip_pairs[2*itmp+2] == ipt) )
                                skip = 1;
                        }
    //printf("in NbrList_recons_use_list 4.2\n");
                    }                            
                    /* check the distance between atoms ipt and jpt */
                    if( (!skip) && (ipt<jpt) )
                    {
                        sij=_SR[ipt]-_SR[jpt];
                        sij.subint(); // subtract the nearest integer
                        rij=_H*sij;
                        if(rij.norm2()<rsmax2)
                        {
                            nindex[ipt][nn[ipt]]=jpt;
                            nn[ipt]++;
                            if(nn[ipt]>=_NNM)
                                FATAL("reconstruct: index["<<ipt<<"] ("
                                      <<nn[ipt]<<") >= NNM("<<_NNM<<")"
                                      " increase NNM in script");
                            if(max_nnm<nn[ipt]) max_nnm=nn[ipt];
                            
                            nindex[jpt][nn[jpt]]=ipt;
                            nn[jpt]++;
                            if(nn[jpt]>=_NNM)
                                FATAL("reconstruct: index["<<jpt<<"] ("
                                      <<nn[ipt]<<") >= NNM("<<_NNM<<")"
                                      " increase NNM in script");
                            if(max_nnm<nn[jpt]) max_nnm=nn[jpt];
                        }
                    }
                    /* go to next atom in the list */
                    if(jpt==link_list[jpt])
                        FATAL("jpt="<<jpt<<" link_list["<<jpt<<"]="<<link_list[jpt]);
                    jpt = link_list[jpt];
                }
                /* go to next atom in the list */
                if(ipt==link_list[ipt])
                    FATAL("ipt="<<ipt<<" link_list["<<ipt<<"]="<<link_list[ipt]);
                ipt = link_list[ipt];
            }
        }
    }
    }
    else
    {
        s=_SR[iatom]; s.subint();
        i=(int)floor((s.x+0.5)*ncx);i=(i+2*ncx)%ncx;
        j=(int)floor((s.y+0.5)*ncy);j=(j+2*ncy)%ncy;
        k=(int)floor((s.z+0.5)*ncz);k=(k+2*ncz)%ncz;

        cell_ind = (i*ncy + j)*ncz + k;

        /* generat neighbor cells list */
        neigh_cell_ind[0] = cell_ind; /* self at the top of list */
        num_neigh_cell = 1;
        for(i1=((ncx==1)?i:(i-1));i1<=((ncx<=2)?i:(i+1));i1++)
            for(j1=((ncy==1)?j:(j-1));j1<=((ncy<=2)?j:(j+1));j1++)
                for(k1=((ncz==1)?k:(k-1));k1<=((ncz<=2)?k:(k+1));k1++)
                {
                    i2=(i1+2*ncx)%ncx;
                    j2=(j1+2*ncy)%ncy;
                    k2=(k1+2*ncz)%ncz;
                    if((i==i2)&&(j==j2)&&(k==k2)) continue; /* self already included */
                    neigh_cell_ind[num_neigh_cell] = (i2*ncy + j2)*ncz + k2;
                    num_neigh_cell ++;                    
                }
        max_nnm = 0;
        /* loop through all neighbor cells */
        for(ineigh=0;ineigh<num_neigh_cell;ineigh++)
        {
            head2 = link_head[neigh_cell_ind[ineigh]]; 
        
            /* going through lists with head2 */
            ipt = iatom;
            jpt = head2;
            while(jpt!=-1)
            {
                /* check whether pair ipt, jpt is in the list to skip */
                skip = 0;
                if (nl_skip_pairs[0]>0)
                {
                    for(itmp=0;itmp<nl_skip_pairs[0];itmp++)
                    {
                        if(   (nl_skip_pairs[2*itmp+1] == ipt)
                              &&(nl_skip_pairs[2*itmp+2] == jpt) )
                            skip = 1;
                        if(   (nl_skip_pairs[2*itmp+1] == jpt)
                              &&(nl_skip_pairs[2*itmp+2] == ipt) )
                            skip = 1;
                    }
                }                            
                /* check the distance between atoms ipt and jpt */
                if( (!skip) && (ipt!=jpt) )
                {
                    sij=_SR[ipt]-_SR[jpt];
                    sij.subint(); // subtract the nearest integer
                    rij=_H*sij;
                    if(rij.norm2()<rsmax2)
                    {
                        nindex[ipt][nn[ipt]]=jpt;
                        nn[ipt]++;
                        if(nn[ipt]>=_NNM)
                            FATAL("reconstruct: index["<<ipt<<"] ("
                                  <<nn[ipt]<<") >= NNM("<<_NNM<<")"
                                  "increase NNM in script");
                        if(max_nnm<nn[ipt]) max_nnm=nn[ipt];
                        
                        nindex[jpt][nn[jpt]]=ipt;
                        nn[jpt]++;
                        if(nn[jpt]>=_NNM)
                            FATAL("reconstruct: index["<<jpt<<"] ("
                                  <<nn[ipt]<<") >= NNM("<<_NNM<<")"
                                  "increase NNM in script");
                        if(max_nnm<nn[jpt]) max_nnm=nn[jpt];
                    }
                    /* go to next atom in the list */
                    if(jpt==link_list[jpt])
                        FATAL("jpt="<<jpt<<" link_list["<<jpt<<"]="<<link_list[jpt]);
                    jpt = link_list[jpt];
                }
            }
        }
    }
    /* save current position to R0 */
    if(iatom==-1)
    {
        for(i=0;i<_NP;i++)
        {
            _R[i]=_H*_SR[i];
            _R0[i]=_R[i];
        }
        //INFO_Printf("reconstruct() finished. max_nnm=%d allocated NNM = %d\n",max_nnm,_NNM);

        iatom=_NP-1;
        _R[iatom]=_H*_SR[iatom];
        _R0[iatom]=_R[iatom];
        /* print neighbor list of iatom */
//        INFO_Printf("atom %d (%e,%e,%e) has %d neighbors: (",
//                    iatom,_SR[iatom].x,_SR[iatom].y,_SR[iatom].z,nn[iatom]);
//        for(j=0;j<nn[iatom];j++)
//            INFO_Printf("%d ",nindex[iatom][j]);
//        INFO_Printf(")\n");        
    }
    else
    {
        _R[iatom]=_H*_SR[iatom];
        _R0[iatom]=_R[iatom];
        /* print neighbor list of iatom */
//        INFO_Printf("atom %d (%e,%e,%e) has %d neighbors: (",
//                    iatom,_SR[iatom].x,_SR[iatom].y,_SR[iatom].z,nn[iatom]);
//        for(j=0;j<nn[iatom];j++)
//            INFO_Printf("%d ",nindex[iatom][j]);
//        INFO_Printf(")\n");
    }
    _current_NNM=max_nnm;
    
    //INFO_Printf("reconstruct() finished. _current_NNM=%d\n",_current_NNM);
    /* _current_NIM can be printed within TCL script */
}

void MDFrame::NbrList_reconstruct_use_cell_list(int iatom)
{
    /* use cell list combined with Verlet list */
    int ncx,ncy,ncz;
    int ncx0,ncx1,ncy0,ncy1,ncz0,ncz1,np0,np1;
    int i,j,k,n,i1,j1,k1,n1,i2,j2,k2;
    int ipt,jpt, itmp, skip;
    double rsmax2;
    Vector3 s,sij,rij,h;
    int m_nic, m_nnm;

    if(iatom!=-1)
    {
        FATAL("NbrList_reconstruct_use_cell_list iatom="<<iatom<<" not implemented! use link_list");
    }
    
    NbrList_init(_NP*allocmultiple,_NNM);

    rsmax2=_RLIST*_RLIST;
    h=_H.height();

    INFO("reconstruct neighborlist");
    DUMP("reconstruct_cell()  ");
    INFO("H="<<_H<<"  _RLIST="<<_RLIST);
    
    ncx=(int)floor(h[0]/((_RLIST)*1.05));if(ncx==0)ncx=1;
    ncy=(int)floor(h[1]/((_RLIST)*1.05));if(ncy==0)ncy=1;
    ncz=(int)floor(h[2]/((_RLIST)*1.05));if(ncz==0)ncz=1;
    
    NbrList_initcell(ncx,ncy,ncz,_NIC);
    
    ncx0=0; ncx1=ncx; ncy0=0; ncy1=ncy; ncz0=0; ncz1=ncz;
    
    np0=0;
    np1=_NP;

    for(i=ncx0;i<ncx1;i++)
        for(j=ncy0;j<ncy1;j++)
            for(k=ncz0;k<ncz1;k++)
                celllist[i][j][k][0]=0;
    
    m_nic = 0;
    for(ipt=0;ipt<_NP;ipt++)
    {
        s=_SR[ipt];
        s.subint();
        i=(int)floor((s.x+0.5)*ncx);i=(i+2*ncx)%ncx;
        j=(int)floor((s.y+0.5)*ncy);j=(j+2*ncy)%ncy;
        k=(int)floor((s.z+0.5)*ncz);k=(k+2*ncz)%ncz;
        if((i>=ncx0)&&(i<ncx1)&&(j>=ncy0)&&(j<ncy1)&&(k>=ncz0)&&(k<ncz1))
        {
            //INFO_Printf("NbrList_reconstruct: ipt=%d i=%d j=%d k=%d\n",ipt,i,j,k);
            celllist[i][j][k][0]++;
            if(celllist[i][j][k][0]>_NIC-1)
                FATAL("reconstruct(): too many atoms per cell "
                      <<celllist[i][j][k][0]<<" > limit("<<(_NIC-1)<<")"
                      <<"increase NIC in script");
            celllist[i][j][k][celllist[i][j][k][0]]=ipt;
            if(m_nic<celllist[i][j][k][0]) m_nic=celllist[i][j][k][0];
        }
    }
    _current_NIC=m_nic;
    
    for(ipt=np0;ipt<np1;ipt++)
    {
        nn[ipt]=0;
        for(itmp=0;itmp<_NNM;itmp++)
            nindex[ipt][itmp]=0;
    }
    DUMP("list cleared");

    m_nnm=0;
    for(i=ncx0;i<ncx1;i++)
        for(j=ncy0;j<ncy1;j++)
            for(k=ncz0;k<ncz1;k++)
                for(n=1;n<=celllist[i][j][k][0];n++)
                {
                    ipt=celllist[i][j][k][n];
                    //neighboring cell
#if 1
                            {
                                i2=i;
                                j2=j;
                                k2=k;
                                //INFO_Printf("i2=%d j2=%d k2=%d\n",i2,j2,k2);
                                for(n1=1;n1<=celllist[i2][j2][k2][0];n1++)
                                {
                                    jpt=celllist[i2][j2][k2][n1];
                                    if(ipt==jpt) continue;

                                    /* check whether pair ipt, jpt is in the list to skip */
                                    skip = 0;
                                    if (nl_skip_pairs[0]>0)
                                    {
                                        for(itmp=0;itmp<nl_skip_pairs[0];itmp++)
                                        {
                                            if(   (nl_skip_pairs[2*itmp+1] == ipt)
                                                  &&(nl_skip_pairs[2*itmp+2] == jpt) )
                                                skip = 1;
                                            if(   (nl_skip_pairs[2*itmp+1] == jpt)
                                                  &&(nl_skip_pairs[2*itmp+2] == ipt) )
                                                skip = 1;
                                        }
                                    }
                                    if (skip) continue;
                                    
                                    sij=_SR[ipt]-_SR[jpt];
                                    sij.subint(); // subtract the nearest integer
                                    rij=_H*sij;
                                    if(rij.norm2()<rsmax2)
                                    {
                                        nindex[ipt][nn[ipt]]=jpt;
                                        nn[ipt]++;
                                        if(nn[ipt]>=_NNM)
                                            FATAL("reconstruct: index["<<ipt<<"] ("
                                                  <<nn[ipt]<<") >= NNM("<<_NNM<<")"
                                                  "increase NNM in script");
                                        if(m_nnm<nn[ipt]) m_nnm=nn[ipt];
                                    }
                                }
                            }
#endif                            
#if 1
                    for(i1=((ncx==1)?i:(i-1));i1<=((ncx<=2)?i:(i+1));i1++)
                        for(j1=((ncy==1)?j:(j-1));j1<=((ncy<=2)?j:(j+1));j1++)
                            for(k1=((ncz==1)?k:(k-1));k1<=((ncz<=2)?k:(k+1));k1++)
                            {
                                i2=(i1+2*ncx)%ncx;
                                j2=(j1+2*ncy)%ncy;
                                k2=(k1+2*ncz)%ncz;
                                if((i==i2)&&(j==j2)&&(k==k2)) continue;
                                //INFO_Printf("i2=%d j2=%d k2=%d\n",i2,j2,k2);
                                for(n1=1;n1<=celllist[i2][j2][k2][0];n1++)
                                {
                                    jpt=celllist[i2][j2][k2][n1];
                                    if(ipt==jpt) continue;
                                    
                                    /* check whether pair ipt, jpt is in the list to skip */
                                    skip = 0;
                                    if (nl_skip_pairs[0]>0)
                                    {
                                        for(itmp=0;itmp<nl_skip_pairs[0];itmp++)
                                        {
                                            if(   (nl_skip_pairs[2*itmp+1] == ipt)
                                                  &&(nl_skip_pairs[2*itmp+2] == jpt) )
                                                skip = 1;
                                            if(   (nl_skip_pairs[2*itmp+1] == jpt)
                                                  &&(nl_skip_pairs[2*itmp+2] == ipt) )
                                                skip = 1;
                                        }
                                    }
                                    if (skip) continue;
                                    
                                    sij=_SR[ipt]-_SR[jpt];
                                    sij.subint(); // subtract the nearest integer
                                    rij=_H*sij;
                                    if(rij.norm2()<rsmax2)
                                    {
                                        nindex[ipt][nn[ipt]]=jpt;
                                        nn[ipt]++;
                                        if(nn[ipt]>=_NNM)
                                            FATAL("reconstruct: index["<<ipt<<"] ("
                                                  <<nn[ipt]<<") >= NNM("<<_NNM<<")"
                                                  "increase NNM in script");
                                        if(m_nnm<nn[ipt]) m_nnm=nn[ipt];
                                    }
                                }
                            }
#endif                    
                }
    
    _current_NNM=m_nnm;
    
    /* combine local lists to shm list */
    for(i=np0;i<np1;i++)
    {
        _R[i]=_H*_SR[i];
        _R0[i]=_R[i];
    }
    //INFO_Printf("reconstruct() finished. _current_NIC=%d _current_NNM=%d\n",_current_NIC,_current_NNM);
    /* _current_NIC and _current_NIM can be printed within TCL script */
}

//Xiaohan adds for frank partial
void MDFrame::NbrList_fprint()
{
    int i, j, n, ipt, jpt;
    FILE *fp = fopen("nn.dat","w");
    if(fp==NULL) FATAL("NbrList_print: open file failure");
    for(i=0;i<_NP;i++) {
        for(j=0;j<nn[i];j++)
            fprintf(fp," %d",nindex[i][j]);
        fprintf(fp,"\n");
    }
}

void MDFrame::NbrList_print() 
{   /* debug use */
    int i, j, n, ipt, jpt;

    n = (int) input[0];
    for(i=0;i<n;i++)
    {
        ipt = (int) input[i+1];

        INFO_Printf("atom [%d] has %d neighbors ( ",ipt,nn[ipt]);

        for(j=0;j<nn[i];j++)
        {
            jpt = nindex[ipt][j];
            INFO_Printf(" %d",jpt);
        }
        INFO_Printf(")\n");
    }
}

bool MDFrame::NbrList_needrefresh()
{
    /* all processes follow same procedure */
    double maxd;
    int i, need;
    Vector3 dR;
    if (firsttime)
    {
        firsttime=false;
        return true;
    }
    maxd=0;
    SHtoR();
   for(i=0;i<_NP;i++)
   {
       /* change the algorithm of need for Nbrlist update
          Jan. 19 2007 Keonwook Kang */
       dR = _R[i] - _R0[i];
       maxd=std::max(dR.norm2(),maxd);
       //maxd=max(fabs(_R[i].x-_R0[i].x),maxd);
       //maxd=max(fabs(_R[i].y-_R0[i].y),maxd);
       //maxd=max(fabs(_R[i].z-_R0[i].z),maxd);
   }
   need=(maxd>_SKIN*_SKIN/4.0);
   //need=(maxd>_SKIN);

    return (need!=0);
}

void MDFrame::NbrList_init(int mx,int mz)
{
    /* Multi process function */
    /* init two dimensional array */
#ifndef _USEBOOST
    int shft1,shft2;//, mx, mz;
    int i;

    DUMP("initlist("<<mx<<","<<mz<<")");
    shft1=mx*mz*sizeof(int);
    shft2=mx*sizeof(int *);
    if(shft1+shft2==0) return;
    
    Realloc(nindex_mem,char,(shft1+shft2));

    nindex=(int **)(nindex_mem+shft1);

    memset(nindex_mem,0,shft1+shft2);
    for(i=0;i<mx;i++)
    {
        nindex[i]=(int *)(nindex_mem+i*mz*sizeof(int));
    }

    Realloc(nn,int,mx);
#else
    nn.resize(boost::extents[mx]);
    nindex.resize(boost::extents[mx][mz]);
#endif
}

void MDFrame::NbrList_initlinklist(int ncell)
{
    int i;

    Realloc(link_head,int,ncell);
    Realloc(link_list,int,_NP);

    for(i=0;i<ncell;i++) link_head[i]=-1;
    for(i=0;i<_NP;i++)   link_list[i]=-1;
}

void MDFrame::NbrList_freelinklist()
{
    free(link_head);
    free(link_list);
}

void MDFrame::NbrList_initcell(int mx, int my, int mz, int mc)
{
#ifndef _USEBOOST
    /* init two dimensional array */
    int shft1,shft2,shft3,shft4;
    int i, j, k;

//    INFO_Printf("MDFrame::NbrList_initcell(%d,%d,%d,%d)\n",mx,my,mz,mc);
//    freecelllist();
    shft1=mx*my*mz*mc*sizeof(int);
    shft2=mx*my*mz*sizeof(int *);
    shft3=mx*my*sizeof(int **);
    shft4=mx*sizeof(int ***);
    if(shft1+shft2+shft3+shft4==0) return;
    
    INFO("initcell("<<mx<<","<<my<<","<<mz<<","<<mc<<") cell list size "<<shft1+shft2+shft3+shft4);
    Realloc(cell_mem,char,(shft1+shft2+shft3+shft4));
    celllist=(int ****)(cell_mem+shft1+shft2+shft3);

    memset(cell_mem,0,shft1+shft2+shft3+shft4);
    for(i=0;i<mx;i++)
    {
        celllist[i]=(int ***)(cell_mem+shft1+shft2+i*my*sizeof(int **));
        for(j=0;j<my;j++)
        {
            celllist[i][j]=(int **)(cell_mem+shft1+(i*my*mz+j*mz)*sizeof(int *));
            for(k=0;k<mz;k++)
            {
                celllist[i][j][k]=
                    (int *)(cell_mem+(i*my*mz*mc+j*mz*mc+k*mc)*sizeof(int));
            }
        }
    }
#else
    celllist.resize(boost::extents[mx][my][mz][mc]);
#endif
}

void MDFrame::NbrList_free()
{
#ifndef _USEBOOST
    free(nindex_mem);
    nindex_mem=NULL;
    nn=NULL;
    nindex=NULL;
#endif
}

void MDFrame::NbrList_freecell()
{
#ifndef _USEBOOST
    free(cell_mem);
    cell_mem=NULL;
    celllist=NULL;
#endif
}

void MDFrame::NbrList_refresh()
{
    if(NbrList_needrefresh()) NbrList_reconstruct();
}

void MDFrame::NbrListPlot_reconstruct()
{
    int ipt, j, jpt;
    Vector3 s, r;
    double rsmax2;

    if (nn_plot == NULL)
    {
        INFO("NBR test");
        NbrListPlot_init(_NP,_NNM_plot);
    }

    rsmax2=rc_plot*rc_plot;

    INFO("reconstruct neighborlist_plot");
    for(ipt=0;ipt<_NP;ipt++)
    {
        nn_plot[ipt]=0;
        for(j=0;j<nn[ipt];j++)
        {/* Find neighbor atoms which is located within "rc_plot"
            and store their indices at "nindex_plot".
            nn_plot[i]: number of neighbors located within "rc_plot" */
            jpt=nindex[ipt][j];
            s=_SR[ipt]-_SR[jpt];
            s.subint();
            r = _H*s;
//            INFO("r.norm2="<<r.norm2()<<", rsmax2="<<rsmax2);
            if(r.norm2()<rsmax2)
            {
                nindex_plot[ipt][nn_plot[ipt]]=jpt;
                nn_plot[ipt]++;
            }
//            INFO("nn_plot="<<nn_plot[ipt]);
//            INFO("ipt="<<ipt<<", j="<<j);
        }
    }
    INFO_Printf("reconstruct_plot() finished\n");
}

void MDFrame::NbrListPlot_init(int mx,int mz)
{
    /* init two dimensional array */
    int shft1,shft2;//, mx, mz;
    int i;

    DUMP("initlist("<<mx<<","<<mz<<")");
    shft1=mx*mz*sizeof(int);
    shft2=mx*sizeof(int *);
    if(shft1+shft2==0) return;

    Realloc(nindex_plot_mem,char,(shft1+shft2));

    nindex_plot=(int **)(nindex_plot_mem+shft1);

    memset(nindex_plot_mem,0,shft1+shft2);
    for(i=0;i<mx;i++)
    {/* prepare memory "nindex_plot" to store indices of neighbor atoms
        of an atom i */
        nindex_plot[i]=(int *)(nindex_plot_mem+i*mz*sizeof(int));
    }
    /* prepare memory "nn_plot" to store */
    Realloc(nn_plot,int,mx);
}

/********************************************************************/
/* File input and output */
/********************************************************************/
    
int MDFrame::readcn()
{
    initcn.open(incnfile,LFile::O_Read);
    return initcn.read(this);
}

int MDFrame::readcontinuecn()
{
    continuecn.open(continuecnfile,LFile::O_Read);
    return continuecn.read(this);
}

int MDFrame::readLAMMPS()
{
    initlammps.open(incnfile,LFile::O_Read);
    return initlammps.read(this);
}

int MDFrame::writeLAMMPS()
{
    finallammps.open(finalcnfile);
    return finallammps.write(this);
}

int MDFrame::setfilecounter()
{
    intercn.setcount(filecounter);
    return 0;
}

int MDFrame::setFFSfilecounter()
{
    FFScn.setcount(FFSfilecounter);
    return 0;
}

int MDFrame::convertXDATCAR()
{
    int i, j;
    char *buffer; char *pp, *q;
    int np, nframe;
    char fname[200];

    strcpy(fname,incnfile);
    LFile::SubHomeDir(fname,fname);
    LFile::LoadToString(fname,buffer,0);

    pp=buffer;
    
    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
    sscanf(q,"%d %d %d",&np, &j, &nframe);
    
    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;

    for(j=0;j<nframe;j++)
    {
        q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
        
        for(i=0;i<_NP;i++)
        {
            q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
            sscanf(q,"%lf %lf %lf",
                   &(_SR[i].x),&(_SR[i].y),&(_SR[i].z));
        }
        INFO("frame "<<j);
        writeintercnfile(zipfiles,false);
    }
    Free(buffer);
    return 0;
}


int MDFrame::writefinalcnfile(int zip, bool bg)
{
    finalcn.open(finalcnfile);
    return finalcn.write(this,zip,bg);

}

int MDFrame::writecontinuecnfile(int zip, bool bg)
{
    continuecn.open(continuecnfile);
    return continuecn.write(this,zip,bg);
}

int MDFrame::writeavgcnfile(int zip, bool bg)
{/* May 17 2007 Keonwook Kang */
    //avgcn.open(finalcnfile);
    return 0; //avgcn.write(this,zip,bg);
}

void MDFrame::saveintercn(int step)
{
    if(savecn)
    {   // If it is the 1st intermediate cn file, open it. May 01 2009 Keonwook Kang
//        if(intercn.f==NULL && intercn.count==1)
        if(intercn.f==NULL && intercn.count==filecounter)
            openintercnfile();
//        if((step%savecnfreq)==0&&(step!=step0))
        if((step%savecnfreq)==0&&(step!=0))
        {
            writeintercnfile(zipfiles,false);
//            INFO("intercn.count="<<intercn.count);
        }
    }
}

void MDFrame::saveintercfg(int step)
{// Apr 03, 2010 Keonwook Kang
    if(savecfg)
    {
        if((step%savecfgfreq)==0&&(step!=0))
        {
            char extname[100];
            AUXFile::insertindex(extname,intercfgfile,savecfgcount);
            writeatomeyecfgfile (extname); savecfgcount ++;
        }
    }
}

void MDFrame::savecontinuecn(int step)
{
    if ( (savecontinuecnfreq>0) && (step%savecontinuecnfreq==0) && (step!=0) )
	writecontinuecnfile(zipfiles,false);
}

int MDFrame::writeintercnfile(int zip,bool bg)
{
    return intercn.write(this,zip,bg);
}

int MDFrame::writeFFScnfile(int zip, bool bg)
{
    return FFScn.write(this,zip,bg);
}


int MDFrame::openintercnfile()
{
    return intercn.open(intercnfile);
}

int MDFrame::openFFScnfile()
{
    return FFScn.open(FFScnfile);
}

int MDFrame::openpropfile()
{
    return pf.open(outpropfile);
}

int MDFrame::closepropfile()
{
    pf.close(zipfiles,true);
    return 0;
}

/* Configuration Files */
char * CNFile::describe()
{
    static char tmp[500];
    sprintf(tmp,"%s","Configuration File for Atom Positions, Format:\n"
            "sx1 sy1 sz1; ... ; sxn syn szn; H(3x3)");
    return tmp;
}

int CNFile::writeblock(void *p)
{
    int i;
    MDFrame &d=*((MDFrame *)p);
    f->Printf("%d\n",d._NP);
    for(i=0;i<d._NP;i++)
    {
        if (d.writeall)
            f->Printf("%25.17E %25.17E %25.17E %25.17E %25.17E %25.17E"
                      " %25.17E %2d %25.17E %2d %2d %2d\n",
                      d._SR[i].x,d._SR[i].y,d._SR[i].z,
                      d._VSR[i].x,d._VSR[i].y,d._VSR[i].z,
                      d._EPOT_IND[i],d.fixed[i],d._TOPOL[i],
                      d.species[i],d.group[i],d.image[i]);
        else if (d.writevelocity)
            f->Printf("%25.17E %25.17E %25.17E %25.17E %25.17E %25.17E\n",
                      d._SR[i].x,d._SR[i].y,d._SR[i].z,
                      d._VSR[i].x,d._VSR[i].y,d._VSR[i].z);
        else
            f->Printf("%25.17E %25.17E %25.17E\n",
                      d._SR[i].x,d._SR[i].y,d._SR[i].z);
    }
    f->Printf("%25.17E %25.17E %25.17E\n"
              "%25.17E %25.17E %25.17E\n"
              "%25.17E %25.17E %25.17E\n",
              d._H[0][0],d._H[0][1],d._H[0][2],
              d._H[1][0],d._H[1][1],d._H[1][2],
              d._H[2][0],d._H[2][1],d._H[2][2]);
    f->Printf("%d ",d.nspecies);
    for(i=0;i<d.nspecies;i++) f->Printf("%s ",d.element[i]);
    f->Printf("\n");

    /* write NVT-related varialbes for VVerlet*/
    /* Jan 22 2007 Keonwook Kang */
    f->Printf("%25.17E %25.17E\n", d.zeta, d.zetav);    
    /* in principle we may want to save VH as well */
    
    return 0;
}

int CNFile::readblock(void *p)
{
    int i;
    char *buffer, *pp, *q;
    MDFrame &d=*((MDFrame *)p);
    
    LFile::LoadToString(fname,buffer,0);

    pp=buffer;
    sscanf(pp, "%d", &d._NP);
    INFO("readblock: NP="<<d._NP);    
    
    d.Alloc();
    
    pp=strchr(pp, '\n');
    pp++;
    for(i=0;i<d._NP;i++)
    {
        q=pp;
        pp=strchr(pp,'\n');
        if(pp) *(char *)pp++=0;
//        d._VSR[i].clear(); // commented out KW Kang 2007.06.01
//        d.fixed[i]=0;   // commented out by Keonwook Kang 2007.03.23
        sscanf(q, "%lf %lf %lf %lf %lf %lf %lf %d %lf %d %d %d",
               &(d._SR[i].x),&(d._SR[i].y),&(d._SR[i].z),
               &(d._VSR[i].x),&(d._VSR[i].y),&(d._VSR[i].z),
               &(d._EPOT_IND[i]),&(d.fixed[i]),&(d._TOPOL[i]),
               &(d.species[i]),&(d.group[i]),&(d.image[i]));
//        d.fixed[i]=0;
//        INFO_Printf("species[%d]=%d\n",i,d.species[i]);
    }

//    INFO_Printf("reading H...\n");
    q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
    sscanf(q, "%lf %lf %lf",&d._H[0][0],&d._H[0][1],&d._H[0][2]);
    q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
    sscanf(q, "%lf %lf %lf",&d._H[1][0],&d._H[1][1],&d._H[1][2]);
    q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
    sscanf(q, "%lf %lf %lf",&d._H[2][0],&d._H[2][1],&d._H[2][2]);

                                                      /* 7/30/07 Wei Cai */
    q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0; else {Free(buffer); return 0;}
//    INFO_Printf("reading species...\n");
    sscanf(q, "%d %s %s %s %s %s %s %s %s %s %s",
           &(d.nspecies),
           d.element[0],d.element[1],d.element[2],d.element[3],d.element[4],
           d.element[5],d.element[6],d.element[7],d.element[8],d.element[9]);
//    INFO_Printf("nspecies = %d %s %s\n",d.nspecies,d.element[0],d.element[1]);

    /* write NVT-related varialbes*/
    /* Jan 22 2007 Keonwook Kang */                   /* 7/30/07 Wei Cai */
    q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0; else {Free(buffer); return 0;}
    sscanf(q, "%lf %lf",&d.zeta, &d.zetav);
    
    Free(buffer);
    DUMP("readblock finished");
    return 0;
}


int MDFrame::readXYZ()
{/* read XYZ format data 2007/Oct/3 Keonwook Kang */
    int i, j, nel, isatom;
    char *buffer; char *pp, *q;
    char el[5], fname[300];
//    double xmin, xmax, ymin, ymax, zmin, zmax, xshift, yshift, zshift;
    Matrix33 hinv;

    LFile::SubHomeDir(incnfile,fname);
    LFile::LoadToString(fname,buffer,0);

    pp=buffer;
    
    q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
    sscanf(q, "%d", &_NP);
    INFO("readXYZ: NP="<<_NP);
    
    Alloc();
    
    hinv = _H.inv();

    q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;

    nel = 0; // number of elements is initialized.
    isatom = 1;
    for(i=0;i<_NP;i++)
    {
        q=pp; pp=strchr(pp,'\n'); if(pp) *(char *)pp++=0;
        sscanf(q, "%s %lf %lf %lf",
               el,&(_R[i].x),&(_R[i].y),&(_R[i].z));

        _SR[i] = hinv*_R[i];

        if (i==0)  // if it is the 1st atom,
	{
	    strcpy(element[nel],el); species[i]=nel; nel=1;
	}
        else // From the 2nd atom,
	{
            for(j=0;j<nel;j++)
	    {
	        isatom = strcmp(el, element[j]);
                if (!isatom) // if they are same
	        {
	        species[i]=j; break;
                }
            }
            if (isatom) // if there is no speceis predefined
	    {
                strcpy(element[nel],el); species[i]=nel; nel++;
            }
	}
    }
    nspecies = nel;
/*    _H.clear();
    _H.set(2*(xmax-xmin), 0, 0,
           0, 2*(ymax-ymin), 0,
           0, 0, 2*(zmax-zmin));
    hinv = _H.inv();

    xshift = (xmax+xmin)/2;
    yshift = (ymax+ymin)/2;
    zshift = (zmax+zmin)/2;

    for(i=0;i<_NP;i++)
    {
        _R[i].x -= xshift;
        _R[i].y -= yshift;
        _R[i].z -= zshift;
        _SR[i] = hinv*_R[i];
    }
*/
    Free(buffer);
    DUMP("readXYZ finished");
    return 0;
}

/* Data Files in LAMMPS */
char * LAMMPSFile::describe()
{
    static char tmp[500];
    sprintf(tmp,"%s","data file MD++ -> LAMMPS\n"
            "dump file LAMMPS --> MD++");
    return tmp;
}


int LAMMPSFile::readblock(void *p)
{ /* dump file LAMMPS --> MD++ */
    int i, stepnum, readstepnum, match, atomID, scale;
    char *buffer; char *pp, *q; 
    double xmin, xmax, ymin, ymax, zmin, zmax, xshift, yshift, zshift;
    Matrix33 hinv;
    
    MDFrame &d=*((MDFrame *)p);

    readstepnum = (int) d.input[0];
    scale = (int) d.input[1];
    
    LFile::LoadToString(fname,buffer,0);

    INFO("readLAMMPS:: read into string");
    pp = buffer;
    
    /* do not skip first line */
    //pp=buffer; pp=strchr(pp, '\n'); pp++;

    match = 0;
    while (1)
    {
        /* skip a line */
        pp=strchr(pp, '\n'); if (pp==NULL) break; pp++;
        
        /* put a new line into q */
        q=pp; pp=strchr(pp,'\n'); if (pp==NULL) break; if(pp) *(char *)pp++=0;
        sscanf(q, "%d", &stepnum);

        /* skip a line */
        pp=strchr(pp, '\n'); pp++;

        /* put a new line into q */
        q=pp; pp=strchr(pp,'\n'); if (pp==NULL) break; if(pp) *(char *)pp++=0;
        sscanf(q, "%d", &d._NP);

        if(stepnum!=readstepnum)
        {
            for(i=0;i<d._NP+5;i++)
            {
                /* skip a line */
                pp=strchr(pp, '\n'); pp++;
            }
            INFO_Printf("%d, %d  Skip %d lines\n",stepnum, readstepnum, d._NP+5);
        }
        else
        {
            match = 1;
            INFO_Printf("stepnum = %d matches readstepnum %d\n",stepnum,readstepnum);
            INFO("readblock: NP="<<d._NP);
            d.Alloc();
            INFO("readLAMMPS:: allocated");
            
            /* skip a line */
            pp=strchr(pp, '\n'); pp++;

            q=pp; pp=strchr(pp,'\n'); if (pp==NULL) break; if(pp) *(char *)pp++=0;
            sscanf(q, "%lf %lf", &xmin, &xmax);
            q=pp; pp=strchr(pp,'\n'); if (pp==NULL) break; if(pp) *(char *)pp++=0;
            sscanf(q, "%lf %lf", &ymin, &ymax);
            q=pp; pp=strchr(pp,'\n'); if (pp==NULL) break; if(pp) *(char *)pp++=0;
            sscanf(q, "%lf %lf", &zmin, &zmax);
            
            d._H.clear();
            d._H.set(xmax-xmin, 0, 0,
                     0, ymax-ymin, 0,
                     0, 0, zmax-zmin);
            hinv = d._H.inv();
            
            xshift = (xmax+xmin)/2;
            yshift = (ymax+ymin)/2;
            zshift = (zmax+zmin)/2;
            
            /* skip a line */
            pp=strchr(pp, '\n'); pp++;
            if (scale == 0)
            {
                for(i=0;i<d._NP;i++)
                {
                    q=pp; pp=strchr(pp,'\n');  if (pp==NULL) break; if(pp) *(char *)pp++=0;
                    sscanf(q, "%d ", &atomID); 
                    q=strchr(q,' '); q++;  
//                sscanf(q, "%lf %lf %lf %lf %lf %lf",
//                       &(d._R[i].x),&(d._R[i].y),&(d._R[i].z),
//                       &(d._VR[i].x),&(d._VR[i].y),&(d._VR[i].z));
                    sscanf(q, "%d %lf %lf %lf %lf %lf %lf",
                           &(d.species[atomID-1]),
                           &(d._R[atomID-1].x),&(d._R[atomID-1].y),&(d._R[atomID-1].z),
                           &(d._VR[atomID-1].x),&(d._VR[atomID-1].y),&(d._VR[atomID-1].z));
                    d.species[atomID-1]-=1;
                    /* need to change velocity output !! */                    
                    d._VR[atomID-1].x *= d._TIMESTEP;
                    d._VR[atomID-1].y *= d._TIMESTEP;
                    d._VR[atomID-1].z *= d._TIMESTEP;                                        
                    
//                    INFO_Printf("ID = %d, Species = %d, x = %lf, y = %lf, z = %lf\n",
//                                atomID,d.species[atomID-1],d._R[atomID-1].x,d._R[atomID-1].y,d._R[atomID-1].z);

                    d._R[atomID-1].x -= xshift;
                    d._R[atomID-1].y -= yshift;
                    d._R[atomID-1].z -= zshift;
                    d._SR[atomID-1]=hinv*d._R[atomID-1];
                    d._VSR[atomID-1] = hinv*d._VR[atomID-1];
                    d.fixed[atomID-1]=0;
                }
            }
            else
            {
                for(i=0;i<d._NP;i++)
                {
                    q=pp; pp=strchr(pp,'\n');  if (pp==NULL) break; if(pp) *(char *)pp++=0;
                    sscanf(q, "%d ", &atomID); 
                    q=strchr(q,' '); q++;  
                    sscanf(q, "%d %lf %lf %lf %lf %lf %lf",
                           &(d.species[atomID-1]),
                           &(d._SR[atomID-1].x),&(d._SR[atomID-1].y),&(d._SR[atomID-1].z),
                           &(d._VR[atomID-1].x),&(d._VR[atomID-1].y),&(d._VR[atomID-1].z));
                    d.species[atomID-1]-=1;                    
//                    INFO_Printf("ID=%d, vx = %lf, vy = %lf, vz = %lf\n",
//					 atomID, d._VR[atomID-1].x, d._VR[atomID-1].y, d._VR[atomID-1].z);
                    d._VR[atomID-1].x *= d._TIMESTEP;
                    d._VR[atomID-1].y *= d._TIMESTEP;
                    d._VR[atomID-1].z *= d._TIMESTEP;                                        
                    d._VSR[atomID-1] = hinv*d._VR[atomID-1];

//                    INFO_Printf("ID = %d, Species = %d, xs = %lf, ys = %lf, zs = %lf\n",
//                                atomID,d.species[atomID-1],d._SR[atomID-1].x,d._SR[atomID-1].y,d._SR[atomID-1].z);

                    d._SR[atomID-1].x -= 0.5;
                    d._SR[atomID-1].y -= 0.5;
                    d._SR[atomID-1].z -= 0.5;
//                    d._SR[i]=hinv*d._R[i];
                    d.fixed[atomID-1]=0;
                }
            }
            break;
        }
    }
    if(!match)
        ERROR("No match stepnum found!");
    
    Free(buffer);
    return 0;
}

int LAMMPSFile::writeblock(void *p)
{ /* data file MD++ -> LAMMPS */
    int i;
    double xmag, ymag, zmag;
    class Vector3 c1, c2, c3;
    MDFrame &d=*((MDFrame *)p);
    
    if (d._H[1][0] != 0 && d._H[2][1] != 0 && d._H[2][0] != 0)
    {
        d._H=d._H.reorient();
        WARNING("***********************************************");
        WARNING("* H is reoriented to upper triangular matrix! *");
        WARNING("***********************************************");
    }

    xmag = d._H[0][0]; ymag = d._H[1][1]; zmag = d._H[2][2];
        
    f->Printf("LAMMPS data file by MD++ [%s]\n",d.dirname);
    f->Printf("%d atoms\n",d._NP);
    f->Printf("%d atom types\n",d.nspecies);
    f->Printf("0 %25.17e xlo xhi\n",xmag);
    f->Printf("0 %25.17e ylo yhi\n",ymag);
    f->Printf("0 %25.17e zlo zhi\n",zmag);
    f->Printf("%25.17e %25.17e %25.17e xy xz yz\n\n", d._H[0][1],d._H[0][2],d._H[1][2]);

    f->Printf("Masses\n\n");

    // Changed here SA Jul 30 2008
    for(i=0;i<d.nspecies;i++)
      f->Printf("%d %12.9e \n",i+1,d._ATOMMASS[i]);
    f->Printf("\n");

    f->Printf("Atoms \n\n");
    //ID species x y z (Angstroms)

    for(i=0;i<d._NP;i++)
    {// translate positions by (.5, .5, .5)
        d._SR[i].x += 0.5; d._SR[i].y += 0.5; d._SR[i].z += 0.5;
    }
    
    d.SHtoR();
    for(i=0;i<d._NP;i++)
    {
            f->Printf("%d %d %25.17e %25.17e %25.17e\n",
                  i+1, d.species[i]+1, d._R[i].x,d._R[i].y,d._R[i].z);
    }
    
    for(i=0;i<d._NP;i++)
    {// move positions back
        d._SR[i].x -= 0.5; d._SR[i].y -= 0.5; d._SR[i].z -= 0.5;
    }
    d.SHtoR();

    INFO("writevelocity="<<d.writevelocity);
    if (d.writevelocity)
    {
        f->Printf("\nVelocities  \n\n");
        //f->Printf("#atom-ID vx vy vz (Angstroms/picosecond)\n");
        /* need to change velocity output !! */
        for(i=0;i<d._NP;i++)
        {
            d._VR[i] = d._H*d._VSR[i];
            f->Printf("%d %25.17e %25.17e %25.17e\n",
            i+1, d._VR[i].x/d._TIMESTEP,d._VR[i].y/d._TIMESTEP,d._VR[i].z/d._TIMESTEP);
        }
    }   
    return 0;
}




/* Output Property Files */
char * PropFile::describe()
{
    static char tmp[500];
    /* sprintf(tmp,"%s","Format:\ncurstep TKinetic EPot Pressure"); */
    /* Feb 05 2007 Keonwook Kang */
    sprintf(tmp,"property at every %d time steps",MDFrame().savepropfreq);
    return tmp;
}

int PropFile::writeentry(void *p)
{
    char tmp[500];
    MDFrame &d=*((MDFrame *)p);
    if(strlen(d.output_str)>0)
    {
        *f<<d.output_str;
//        d.output_str[0]=0; /* clear output buffer */
    }
    else
    {   /* old output format */
        sprintf(tmp,"%12d %20.13e %20.13e %20.13e "
                     "%20.13e %20.13e %20.13e %20.13e "
                     "%20.13e %20.13e %20.13e %20.13e %20.13e\n",
            d.curstep,                  /* 1  */
            d._KATOM,                   /* 2  */
            d._EPOT,                    /* 3  */
            d._TOTSTRESS.trace()/3,     /* 4  */
            d._TOTSTRESS[0][1],         /* 5  */
            d._TOTSTRESS[1][2],         /* 6  */
            d._TOTSTRESS[2][0],         /* 7  */
            d._HELM,                    /* 8  */
            d._HELMP,                   /* 9  */
            d._T,                       /* 10 */
            d.zeta,                     /* 11 */
            d.dEdlambda,                /* 12 */
            d._OMEGA                    /* 13 */
            );
        *f<<tmp;
    }
    return 0;
}















/* Interface with Fortran codes */
void MDFrame::writefortraninifile(const char *fname)
{
    FILE *fp;
    char extpath[200];

    INFO("MDFrame::writefortraninifile "<<fname);
    LFile::SubHomeDir(fortranpath,extpath);
    INFO("fortranpath "<<fortranpath<<" extpath "<<extpath<<" potfile"<<potfile);
    fp=fopen(fname,"w");
    if(fp==NULL)
    {
        FATAL("writefortraniinfile: open file failure");
    }
    fprintf(fp,"#1 Settings for F90 Relax BCC Code (generated by MD++)\n");
    fprintf(fp,"#2  Potential file name:  (Note that there is a Tab before the file name!\n");
    fprintf(fp,"\t%s/%s\n",extpath,potfile);
    fprintf(fp,"#4  Initial configuration file name(*.cn):\n");
    fprintf(fp,"\tinit.cn\n");
    fprintf(fp,"#6  Final output configuration file name (*.cn):\n");
    fprintf(fp,"\tfrelaxed.cn\n");
    fprintf(fp,"\n");
    fprintf(fp,"#9  Stress components:  (In GPa)\n");
    fprintf(fp,"#10 	pr(1,1), pr(2,2), pr(3,3):\n");
    fprintf(fp,"\t%E\t%E\t%E\n",
            _EXTSTRESS[0][0]*1e-3*_EXTSTRESSMUL+_EXTPRESSADD*1e-3,
            _EXTSTRESS[1][1]*1e-3*_EXTSTRESSMUL+_EXTPRESSADD*1e-3,
            _EXTSTRESS[2][2]*1e-3*_EXTSTRESSMUL+_EXTPRESSADD*1e-3);
    fprintf(fp,"#12	pr(1,2), pr(2,3), pr(3,1):\n");
    fprintf(fp,"\t%E\t%E\t%E\n",
            _EXTSTRESS[0][1]*1e-3*_EXTSTRESSMUL,
            _EXTSTRESS[1][2]*1e-3*_EXTSTRESSMUL,
            _EXTSTRESS[2][0]*1e-3*_EXTSTRESSMUL);
    fprintf(fp,"\n");
    fprintf(fp,"#15 Conjugate gradient search settings:\n");
    fprintf(fp,"#16 The first trial step in optimization (dfpred):\n");
    fprintf(fp,"\t%e\n",conj_dfpred);
    fprintf(fp,"#18 Tolerance on the residual gradient (gradtl / np):\n");
    fprintf(fp,"\t%e\n",conj_ftol);
    fprintf(fp,"#20 Max number of force calls (maxfn):\n");
    fprintf(fp,"\t%d\n",conj_fevalmax);
    fprintf(fp,"\n");
    fprintf(fp,"#23 Allow box vector to relax? ('1' for yes, '0' for no.)  \n");
    fprintf(fp,"\t%d\n",!conj_fixbox);
    fprintf(fp,"#25 Write the configuration sequence to files? ('1' for yes, '0' for no.)\n");
    fprintf(fp,"\t%d\n",savecn);
    fprintf(fp,"#27 If yes, frequency of saving configuration? (per how many timesteps)\n");
    fprintf(fp,"\t%d\n",savecnfreq);
    fprintf(fp,"\n");
    fprintf(fp,"#30 Set fixed atoms specification  (iffixed, sy0, sy1, sy2, sy3)\n");
    fprintf(fp,"\t0 0.0 0.0 0.0 0.0\n");
    fprintf(fp,"#32 CN format 1: sx,sy,sz, 2: sx,sy,sz, 3: sx,sy,sz,vx,vy,vz,pot,iffixed\n");
    if(writeall) fprintf(fp,"\t3\n");
    else if(writevelocity) fprintf(fp,"\t2\n");
    else fprintf(fp,"\t1\n");

    fclose(fp);
}


/* interface with Li Ju's AtomEye viewer */
void MDFrame::writeatomeyecfgfile(const char *fname)
{
    FILE *fp;
    int i, nplot, nr, n1, n2, n3, nn, ix, iy, iz, k;
    Matrix33 h;  Vector3 s;
    int nw, cind; double Emin, Emax;
    
    INFO("MDFrame::writeatomeyecfgfile "<<fname);
    
    switch(plot_color_axis){
    case(0): color_ind=_EPOT_IND; break;
    case(1): color_ind=_EPOT_IND; break;
    case(2): color_ind=_TOPOL; break;
    case(4): color_ind=_TOPOL; break;
#ifdef _GSL
    case(7): color_ind=_TOPOL; break;
#endif
    case(8): color_ind=_TOPOL; break;
    case(9): color_ind=_TOPOL; break;
    default: ERROR("plot() unknown coloraxis "<<plot_color_axis); return;
    }
    
    fp=fopen(fname,"w");
    if(fp==NULL)
    {
        FATAL("writeatomeyecfgfile: open file failure");
    }

    nplot = 0;
    nr = atomeyerepeat[0];
    if(nr!=0)
    {
        n1=atomeyerepeat[1];
        n2=atomeyerepeat[2];
        n3=atomeyerepeat[3];
        if(n1<=0) n1=1;
        if(n2<=0) n2=1;
        if(n3<=0) n3=1;
        nn = n1*n2*n3;
        INFO_Printf("atomeyerepeat ( %d %d %d )\n",
                    n1, n2, n3);
    }
    else
    {
        nn = n1 = n2 = n3 = 1;
    }

    nw = (int)plot_color_windows[0];
    for(i=0;i<_NP;i++)
    {
        if(nw>0)
        {
            if(plot_limits[0]==1)
                if((_SR[i].x<plot_limits[1])||(_SR[i].x>plot_limits[2])
                   ||(_SR[i].y<plot_limits[3])||(_SR[i].y>plot_limits[4])
                   ||(_SR[i].z<plot_limits[5])||(_SR[i].z>plot_limits[6]))
                {
                    continue;
                }
                for(k=0;k<nw;k++)
                {
                    Emin = plot_color_windows[k*3+1];
                    Emax = plot_color_windows[k*3+2];
                    cind = (int) plot_color_windows[k*3+3];

                    if((color_ind[i]>=Emin)&&
                       (color_ind[i]<=Emax))
                    {
                        nplot ++;
                        continue;    
                    }
                }
        }
        else
        {
            if(plot_limits[0]==1)
                if((_SR[i].x<plot_limits[1])||(_SR[i].x>plot_limits[2])
                   ||(_SR[i].y<plot_limits[3])||(_SR[i].y>plot_limits[4])
                   ||(_SR[i].z<plot_limits[5])||(_SR[i].z>plot_limits[6]))
                {
                    
                    continue;
                }        
            nplot ++;
        }
    }
    h=_H.tran();
    fprintf(fp,"Number of particles = %d\n",nplot*nn);
    fprintf(fp,"A = 1.0 Angstrom (basic length-scale)\n");
    fprintf(fp,"H0(1,1) = %f A\n",h[0][0]*n1);
    fprintf(fp,"H0(1,2) = %f A\n",h[0][1]*n1);
    fprintf(fp,"H0(1,3) = %f A\n",h[0][2]*n1);
    fprintf(fp,"H0(2,1) = %f A\n",h[1][0]*n2);
    fprintf(fp,"H0(2,2) = %f A\n",h[1][1]*n2);
    fprintf(fp,"H0(2,3) = %f A\n",h[1][2]*n2);
    fprintf(fp,"H0(3,1) = %f A\n",h[2][0]*n3);
    fprintf(fp,"H0(3,2) = %f A\n",h[2][1]*n3);
    fprintf(fp,"H0(3,3) = %f A\n",h[2][2]*n3);
#define ATOMEYE_EXTCFG
#ifdef ATOMEYE_EXTCFG
    fprintf(fp,".NO_VELOCITY.\n");
    fprintf(fp,"entry_count = 4\n");
    fprintf(fp,"auxiliary[0] = pote [eV]\n");
//    fprintf(fp,"1.00000\n");
//    fprintf(fp,"%s\n",element[0]);
    for(i=0;i<_NP;i++)
    {
        s=_SR[i];
        if(plot_map_pbc==1) s.subint();
        if(plot_limits[0]==1)
            if((_SR[i].x<plot_limits[1])||(_SR[i].x>plot_limits[2])
               ||(_SR[i].y<plot_limits[3])||(_SR[i].y>plot_limits[4])
               ||(_SR[i].z<plot_limits[5])||(_SR[i].z>plot_limits[6]))
            {
                continue;
            }
        
        if(nw>0)
        {
            for(k=0;k<nw;k++)
            {
                Emin = plot_color_windows[k*3+1];
                Emax = plot_color_windows[k*3+2];
                cind = (int) plot_color_windows[k*3+3];
                
                if((color_ind[i]>=Emin)&&
                   (color_ind[i]<=Emax))
                {
                    for(ix=0;ix<n1;ix++)
                        for(iy=0;iy<n2;iy++)
                            for(iz=0;iz<n3;iz++)
                            {
                                fprintf(fp,"%f\n%s\n %f %f %f %f  \n",
                                        _ATOMMASS[species[i]],element[species[i]],
                                        (s.x+0.5+ix)/n1,
                                        (s.y+0.5+iy)/n2,
                                        (s.z+0.5+iz)/n3,
                                        color_ind[i]);
                            }

                    continue;    
                }
            }
        }
        else
        {
            if(plot_limits[0]==1)
                if((_SR[i].x<plot_limits[1])||(_SR[i].x>plot_limits[2])
                   ||(_SR[i].y<plot_limits[3])||(_SR[i].y>plot_limits[4])
                   ||(_SR[i].z<plot_limits[5])||(_SR[i].z>plot_limits[6]))
                {
                    continue;
                }        
            for(ix=0;ix<n1;ix++)
                for(iy=0;iy<n2;iy++)
                    for(iz=0;iz<n3;iz++)
                    {
                        fprintf(fp,"%f\n%s\n%f %f %f %f  \n",
                                _ATOMMASS[species[i]],element[species[i]],                                
                                (s.x+0.5+ix)/n1,
                                (s.y+0.5+iy)/n2,
                                (s.z+0.5+iz)/n3,
                                color_ind[i]);
                    }
        }
    }        
#else
    for(i=0;i<_NP;i++)
    {
        for(ix=0;ix<n1;ix++)
            for(iy=0;iy<n2;iy++)
                for(iz=0;iz<n3;iz++)
                {
                    fprintf(fp,"%f %s %f %f %f %f %f %f\n",
                            _ATOMMASS[species[i]],element[species[i]],
                            (s.x+0.5+ix)/n1,
                            (s.y+0.5+iy)/n2,
                            (s.z+0.5+iz)/n3,
                            _VSR[i].x/n1,
                            _VSR[i].y/n2,
                            _VSR[i].z/n3);
                }
    }
#endif
    fclose(fp);
    if(zipfiles == 1) {
        if(strlen(fname)>0) LFile::GZipFile(fname,false);
    } else if(zipfiles == 2) {
        if(strlen(fname)>0) LFile::BZipFile(fname,false);
    }
}

void MDFrame::convertCNtoCFG()
{
    int i, istart, iend;
    char cfgfname[100];
    
    istart = (int) input [0];
    iend   = (int) input [1];

    for(i=istart;i<=iend;i++)
    {
        if(strlen(incnfile)>=12)
        {
            sprintf(incnfile + strlen(incnfile)-12,"inter%04d.cn",i);
            strcpy(cfgfname,incnfile);
            sprintf(cfgfname + strlen(cfgfname)-12,"inter%04d.cfg",i);
        }
        else
        {
            sprintf(incnfile,"inter%04d.cn",i);
            sprintf(cfgfname,"inter%04d.cfg",i);
        }
        readcn();
        writeatomeyecfgfile(cfgfname);
    }
        
}


void MDFrame::writepovray(const char *fname)
{
    FILE *fp;
    int i, nplot, nr, n1, n2, n3, nn, ix, iy, iz, k, ipt, jpt;
    Vector3 r, sr;
    int nw, cind; double Emin, Emax;
    
    INFO("MDFrame::writepovray "<<fname);
    
    switch(plot_color_axis){
    case(0): color_ind=_EPOT_IND; break;
    case(1): color_ind=_EPOT_IND; break;
    case(2): color_ind=_TOPOL; break;
    case(4): color_ind=_TOPOL; break;
#ifdef _GSL
     case(7): color_ind=_TOPOL; break;
#endif
    case(8): color_ind=_TOPOL; break;
    case(9): color_ind=_TOPOL; break;
    default: ERROR("plot() unknown coloraxis "<<plot_color_axis); return;
    }
    
    fp=fopen(fname,"w");
    if(fp==NULL)
    {
        FATAL("writepovray: open file failure");
    }

    nplot = 0;
    nr = atomeyerepeat[0];
    if(nr!=0)
    {
        n1=atomeyerepeat[1];
        n2=atomeyerepeat[2];
        n3=atomeyerepeat[3];
        if(n1<=0) n1=1;
        if(n2<=0) n2=1;
        if(n3<=0) n3=1;
        nn = n1*n2*n3;
        INFO_Printf("povraypeat ( %d %d %d )\n",
                    n1, n2, n3);
    }
    else
    {
        nn = n1 = n2 = n3 = 1;
    }

    fprintf(fp,"#include \"../viewpoint.pov\"\n");
    fprintf(fp,"\n");
    fprintf(fp,"union {\n");
    fprintf(fp,"  union {\n");
    fprintf(fp,"\n");
    
    nw = (int) plot_color_windows[0];
    for(i=0;i<_NP;i++)
    {
        if(plot_limits[0]==1)
            if((_SR[i].x<plot_limits[1])||(_SR[i].x>plot_limits[2])
               ||(_SR[i].y<plot_limits[3])||(_SR[i].y>plot_limits[4])
               ||(_SR[i].z<plot_limits[5])||(_SR[i].z>plot_limits[6]))
            {
                continue;
            }
        if(plot_limits[0]==2)
            if((_SR[i].x<plot_limits[1])||(_SR[i].x>plot_limits[2])
               ||(_SR[i].y<plot_limits[3])||(_SR[i].y>plot_limits[4])
               ||(_SR[i].z<plot_limits[5])||(_SR[i].z>plot_limits[6]))
            {
                    sr.set( (_SR[i].x)/n1,
                            (_SR[i].y)/n2,
                            (_SR[i].z)/n3 );
                    r = _H*sr;
                    fprintf(fp,"sphere { <%f, %f, %f>, radiusbc "
                            "pigment { color colorbc "
                            "transmit transmitbc } "
                            "finish { phong 1 metallic }"
                            " }\n",
//                            r.x, r.y, -r.z);
                            r.x, r.y, r.z);
                continue;
            }

        if(nw>0)
        {
            for(k=0;k<nw;k++)
            {
                Emin = plot_color_windows[k*3+1];
                Emax = plot_color_windows[k*3+2];
                cind = (int) plot_color_windows[k*3+3];
                
                if((color_ind[i]>=Emin)&&
                   (color_ind[i]<=Emax))
                {
                    for(ix=0;ix<n1;ix++)
                        for(iy=0;iy<n2;iy++)
                            for(iz=0;iz<n3;iz++)
                            {
                                sr.set( (_SR[i].x+ix)/n1,
                                        (_SR[i].y+iy)/n2,
                                        (_SR[i].z+iz)/n3 );
                                r = _H*sr;
                                fprintf(fp,"sphere { <%f, %f, %f>, radius%02d "
                                        "pigment { color color%02d "
                                        "transmit transmit%02d } "
                                        "finish { phong 1 metallic }"
                                        " }\n",
//                                        r.x, r.y, -r.z, cind, cind, cind);
                                        r.x, r.y, r.z, cind, cind, cind);
                                /* color_ind[i] */
                            }

                    continue;    
                }
            }
        }
        else
        {
            if(plot_limits[0]==1)
                if((_SR[i].x<plot_limits[1])||(_SR[i].x>plot_limits[2])
                   ||(_SR[i].y<plot_limits[3])||(_SR[i].y>plot_limits[4])
                   ||(_SR[i].z<plot_limits[5])||(_SR[i].z>plot_limits[6]))
                {
                    continue;
                }        
            for(ix=0;ix<n1;ix++)
                for(iy=0;iy<n2;iy++)
                    for(iz=0;iz<n3;iz++)
                    {
                        sr.set( (_SR[i].x+ix)/n1,
                                (_SR[i].y+iy)/n2,
                                (_SR[i].z+iz)/n3 );
                        r = _H*sr;
                        //INFO_Printf("species[%d]=%d\n",i,species[i]);
                        //INFO_Printf("fixed[%d]=%d\n",i,fixed[i]);
                                
                        fprintf(fp,"sphere { <%f, %f, %f>, radius%02d "
                                "pigment { color color%02d "
                                "transmit transmit%02d } "
                                "finish { phong 1 metallic }"
                                " }\n",
                                r.x, r.y, r.z,species[i],species[i],species[i]);
                        /* color_ind[i] */
                    }
        }
    }

    for(i=0;i<input[0];i++)
    {
        ipt=(int)input[2*i+1];
        jpt=(int)input[2*i+2];
        _R[ipt]=_H*_SR[ipt];
        _R[jpt]=_H*_SR[jpt];        
        fprintf(fp,"cylinder{ <%f,%f,%f>, <%f,%f,%f>, radiuscy "
                "pigment { color colorcy "
                "transmit transmitcy } "
                "finish { phong 1 metallic }"
                " }\n",
                _R[ipt].x,_R[ipt].y,_R[ipt].z,
                _R[jpt].x,_R[jpt].y,_R[jpt].z);
    }

    fprintf(fp,"\n");
    fprintf(fp,"\n");
    fprintf(fp,"  }\n");
    fprintf(fp,"  AdditionItems01\n");
    fprintf(fp,"  no_shadow\n");
    fprintf(fp,"  rotate rot01\n");
    fprintf(fp,"  rotate rot02\n");
    fprintf(fp,"  rotate rot03\n");
    fprintf(fp,"}\n");
    
    fclose(fp);
}


void MDFrame::writeENERGY(const char *fname)
{
    FILE *fp;
    int i;
    
    INFO("MDFrame::writeENERGY "<<fname);

    switch(plot_color_axis){
    case(0): color_ind=_EPOT_IND; break;
    case(1): color_ind=_EPOT_IND; break;
    case(2): color_ind=_TOPOL; break;
    case(4): color_ind=_TOPOL; break;        
#ifdef _GSL
    case(7): color_ind=_TOPOL; break;
#endif
    case(8): color_ind=_TOPOL; break;
    case(9): color_ind=_TOPOL; break;
    default: ERROR("writeENERGY() unknown coloraxis "<<plot_color_axis); return;
    }
    
    fp=fopen(fname,"w");
    if(fp==NULL)
    {
        FATAL("writeENERGY: open file failure");
    }

    for(i=0;i<_NP;i++)
    {
        fprintf(fp,"  %20.16e\n",color_ind[i]);
    }

    fclose(fp);
}

void MDFrame::writeFORCE(const char *fname)
{
    FILE *fp;
    int i;
    
    INFO("MDFrame::writeFORCE "<<fname);

    fp=fopen(fname,"w");
    if(fp==NULL)
    {
        FATAL("writeFORCE: open file failure");
    }

    if(input[0]==1) INFO("only write forces on fixed atoms");
    for(i=0;i<_NP;i++)
    {
        if(input[0]==1)
        {
            if(fixed[i]!=0)
                fprintf(fp,"%6d  %20.12e %20.12e %20.12e %4d\n",
                        i,_F0[i].x,_F0[i].y,_F0[i].z,fixed[i]);
        }
        else
        {
            fprintf(fp,"%6d  %20.12e %20.12e %20.12e %4d\n",
                    i,_F[i].x,_F[i].y,_F[i].z,fixed[i]);
        }
    }
    fclose(fp);
}

void MDFrame::writePOSITION(const char *fname)
{
    FILE *fp;
    int i;
    
    INFO("MDFrame::writePOSITION "<<fname);

    SHtoR();
    fp=fopen(fname,"w");
    if(fp==NULL)
    {
        FATAL("writePOSITION: open file failure");
    }

    if(input[0]==1) INFO("only write positions of fixed atoms");
    for(i=0;i<_NP;i++)
    {
        if(input[0]==1)
        {
            if(fixed[i]!=0)
                fprintf(fp,"%6d  %20.12e %20.12e %20.12e %4d\n",
                        i,_R[i].x,_R[i].y,_R[i].z,fixed[i]);
        }
        else
        {
            fprintf(fp,"%6d  %20.12e %20.12e %20.12e %4d\n",
                    i,_R[i].x,_R[i].y,_R[i].z,fixed[i]);
        }
    }
    fclose(fp);
}

void MDFrame::GnuPlotHistogram()
{ /* plot histogram of atom properties (e.g. energy, topology)
   * by calling Gnuplot
   *
   * input = [ E0 E1 n ]
   *
   * histogram from E0 to E1 with n bins
   */
    
    FILE *fp;
    int i, nbin, *counts, n;
    double Emin, Emax, *Ebin, dE;
    char comsav[200];
    
    INFO("MDFrame::GnuPlotEnergyHistogram");

    switch(plot_color_axis){
    case(0): color_ind=_EPOT_IND; break;
    case(1): color_ind=_EPOT_IND; break;
    case(2): color_ind=_TOPOL; break;
    case(4): color_ind=_TOPOL; break;     /* For caldisregistry */
#ifdef _GSL
     case(7): color_ind=_TOPOL; break;
#endif
    case(8): color_ind=_TOPOL; break;
    case(9): color_ind=_TOPOL; break;
    default: ERROR("plot() unknown coloraxis "<<plot_color_axis); return;
    }
    
    /* determining plot specification */
    Emin = input[0];
    Emax = input[1];
    nbin = (int)input[2];

    if(Emin==Emax)
    {
        Emin=Emax=color_ind[0];
        for(i=1;i<_NP;i++)
        {
            if(Emin>color_ind[i]) Emin=color_ind[i];
            if(Emax<color_ind[i]) Emax=color_ind[i];
        }
    }
    if(fabs(Emin-Emax)<2e-8)
    {
        Emin-=1e-8;
        Emax+=1e-8;
    }
    if(nbin<=0) nbin=20;
    dE=(Emax-Emin)/nbin;
    
    /* constructing energy histogram */
    counts=(int    *)malloc(sizeof(double)*nbin);
    Ebin  =(double *)malloc(sizeof(double)*nbin);
    for(i=0;i<nbin;i++)
    {
        Ebin[i]=dE*(i+0.5)+Emin;
        counts[i]=0;
    }
    
    for(i=0;i<_NP;i++)
    {
        n=(int)floor((color_ind[i]-Emin)/dE);
        if((n>=0)&&(n<nbin))
            counts[n]++;
    }
    /* write energy into file */
    writeENERGY("eng.out");
    /* write energy histogram data into file */
    fp=fopen("enghist.out","w");
    if(fp==NULL)
    {
        FATAL("GnuPlotEnergyHistogram: open file failure");
    }

    for(i=0;i<nbin;i++)
    {
        fprintf(fp,"  %20.16e %d\n",Ebin[i],counts[i]);
    }
    fclose(fp);

    free(counts);
    free(Ebin);

    /* launch GnuPlot */
    /* write energy histogram data into file */
    fp=fopen("plotenghist.gp","w");
    if(fp==NULL)
    {
        FATAL("GnuPlotEnergyHistogram: open file failure");
    }
    fprintf(fp,"plot 'enghist.out' with histeps\n"); 
    fclose(fp);
    
    fp=fopen("gpw","w");
    if(fp==NULL)
    {
        FATAL("GnuPlotEnergyHistogram: open file failure");
    }
    fprintf(fp,
            "#!/bin/bash\n"
            /*"#!/bin/tcsh\n"*/
            /* These are for seaborg.nersc.gov */
            //"setenv LD_LIBRARY_PATH ~/Tools/Lib:\n"
            //"module add gnuplot\n"
            "gnuplot plotenghist.gp -\n");
    fclose(fp);
    
    strcpy(comsav,command);
    sprintf(command,"xterm -fn 7x13 -sb -e ./gpw &\n");
    strcpy(command,comsav);
}





/********************************************************************/
/* Visualization */
/********************************************************************/

void MDFrame::openwin()
{
    openwindow(win_width,win_height,dirname);
}

int MDFrame::openwindow(int w,int h,const char *n)
{
    if(win!=NULL)
        if(win->alive) return -1;
    win=new YWindow(w,h,n,true,true,false);
    if(win==NULL) return -1;
    else
    {
        win->setinterval(100);
        win->Run();
        return 0;
    }
}

void MDFrame::closewindow() {delete(win);}

void MDFrame::winplot()
{
    if(win!=NULL)
        if(win->alive)
            if((curstep%plotfreq)==0)
            {
                plot();
            }
}

void MDFrame::winplot(int step)
{
    if(win!=NULL)
        if(win->alive)
            if((step%plotfreq)==0)
            {
                plot();
            }
}        

void MDFrame::wintogglepause()
{
    if(win!=NULL)
        if(win->IsAlive())
            win->TogglePause();
}

#include "namecolor.c"
void MDFrame::alloccolors()
{
    int r,g,b;
    if(win!=NULL)
    {
        if(win->alive)
        {
            for(int i=0;i<MAXCOLORS;i++)
            {
                Str2RGB(colornames[i],&r,&g,&b);
                colors[i]=win->AllocShortRGBColor(r,g,b);
            }
            
            Str2RGB(atomcolor[0],&r,&g,&b);
            colors[MAXCOLORS+0]=win->AllocShortRGBColor(r,g,b);
            Str2RGB(bondcolor,&r,&g,&b);
            colors[MAXCOLORS+1]=win->AllocShortRGBColor(r,g,b);
            Str2RGB(highlightcolor,&r,&g,&b);
            colors[MAXCOLORS+2]=win->AllocShortRGBColor(r,g,b);
            Str2RGB(fixatomcolor,&r,&g,&b);
            colors[MAXCOLORS+3]=win->AllocShortRGBColor(r,g,b);
            Str2RGB(backgroundcolor,&r,&g,&b);
            win->bgcolor=win->AllocShortRGBColor(r,g,b);

            for(int i=0;i<MAXSPECIES;i++)
            {
                Str2RGB(atomcolor[i],&r,&g,&b);
                colors[MAXCOLORS+5+i]=win->AllocShortRGBColor(r,g,b);
            }
        }
        else
        {
            WARNING("No window to allocate color for!");
        }
    }
    else
    {
        WARNING("No window to allocate color for!");
    }
}
void MDFrame::alloccolorsX()
{
    if(win!=NULL)
    {
        if(win->alive)
        {
            colors[0]=win->AllocNamedColor(atomcolor[0]);
            colors[1]=win->AllocNamedColor(bondcolor);
            colors[2]=win->AllocNamedColor(highlightcolor);
            colors[3]=win->AllocNamedColor(fixatomcolor);
            win->bgcolor=win->AllocNamedColor(backgroundcolor);
        }
        else
        {
            WARNING("No window to allocate color for!");
        }
    }
    else
    {
        WARNING("No window to allocate color for!");
    }
}

void MDFrame::rotate()
{
    if(win!=NULL)
    {
        if(win->alive)
        {
            win->horizontalRot(rotateangles[0]*M_PI/180);
            win->verticalRot(rotateangles[1]*M_PI/180);
            win->spinRot(rotateangles[2]*M_PI/180);
            if(rotateangles[3]!=0) win->zoom(rotateangles[3]);
            if(rotateangles[4]!=0) win->project(rotateangles[4]);
        }
        else
        {
            WARNING("No window to rotate for!");
        }
    }
    else
    {
        WARNING("No window to rotate for!");
    }
}

void MDFrame::saverot()
{
    if(win!=NULL)
    {
        if(win->alive)
        {
            //win->saveRot(); win->saveScale();
            win->saveView();
        }
        else
        {
            WARNING("No window to rotate for!");
        }
    }
    else
    {
        WARNING("No window to rotate for!");
    }
}

#include "colormap.h"
  
void MDFrame::plot()
{
    int i,jp,j,k;
    double L,r2;
    char s[100]; bool high;
    double x1,y1,z1,x2,y2,z2,dx,dy,dz,dr;
    int r,g,b; double alpha; unsigned ce;
    Vector3 s1, s2, vr1, vr2, sri, srj, ri, rj;
    int nw, cind, show; double Emin, Emax;

    if(win==NULL) return;
    if(!(win->alive)) return;
    
    L=std::max(_H[0][0],_H[1][1]);
    L=std::max(L,_H[2][2])*.5;

    switch(plot_color_axis){
     case(0): color_ind=_EPOT_IND; break;
     case(1): color_ind=_EPOT_IND; break;
     case(2): color_ind=_TOPOL; break;
     //case(3): color_ind=_TOPOL; break;
     case(4): color_ind=_TOPOL; break;     /* For caldisregistry */
#ifdef _GSL
     case(7): color_ind=_TOPOL; break;
#endif
     case(8): color_ind=_TOPOL; break;
     case(9): color_ind=_TOPOL; break;
     default: ERROR("plot() unknown coloraxis "<<plot_color_axis); return;
    }
    
    SHtoR();
    win->Lock();
    win->Clear();
    
    /* draw atom s*/
    for(i=0;i<_NP;i++)
    {
        sri=_SR[i];
        if(plot_map_pbc==1) sri.subint();
        ri = _H*sri;
        
        /* if plot_limits[0] equals to
           0: then this option is ignored
           1: then only plot atoms falling within the limits
              given by plot_limits[1~6]
           2: then atoms falling within the limits will be plotted
              with the same color as the fixed atoms
        */
        if(plot_limits[0]==1)
            if((sri.x<plot_limits[1])||(sri.x>plot_limits[2])
               ||(sri.y<plot_limits[3])||(sri.y>plot_limits[4])
               ||(sri.z<plot_limits[5])||(sri.z>plot_limits[6]))
            {
                continue;
            }
        if(plot_limits[0]==2)
            if((sri.x<plot_limits[1])||(sri.x>plot_limits[2])
               ||(sri.y<plot_limits[3])||(sri.y>plot_limits[4])
               ||(sri.z<plot_limits[5])||(sri.z>plot_limits[6]))
            {
                win->DrawPoint(ri.x/L,ri.y/L,ri.z/L,
                               atomradius[species[i]]/L,colors[MAXCOLORS+3],s,2);
                continue;
            }
        
        if (!plot_atom_info) s[0]=0;
        else
        {   /* when an atom is clicked */
            if(plot_atom_info==1)    /* print scaled coordinates of the atom */
                sprintf(s,"%d,%6.4f,%6.4f,%6.4f",
                        i,sri.x,sri.y,sri.z);
            else if(plot_atom_info==2) /* print real coordinates of the atom */
                sprintf(s,"%d,%6.4f,%6.4f,%6.4f",
                        i,ri.x,ri.y,ri.z);
            else if(plot_atom_info==3) /* print color_ind, fixed, species info of the atom */
                sprintf(s,"%d,%8.6e %d %d",
                        i,color_ind[i], fixed[i], species[i]);
            else if(plot_atom_info==4) /* print force on the atom */
                sprintf(s,"%d,%6.4f,%6.4f,%6.4f",
                        i,_F[i].x,_F[i].y,_F[i].z);
            else
                s[0]=0;
        }
        
        nw = (int)plot_color_windows[0];
        if(nw>0)
        {
            for(k=0;k<nw;k++)
            {
                Emin = plot_color_windows[k*3+1];
                Emax = plot_color_windows[k*3+2];
                cind = (int) plot_color_windows[k*3+3];

                if((color_ind[i]>=Emin)&&(color_ind[i]<=Emax))
                {
                    
                    if(plot_color_bar[0]==0)
                    {
                        win->DrawPoint(ri.x/L,ri.y/L,ri.z/L,
                                       atomradius[species[i]]/L,colors[cind],s,2);
                    }
                    else
                    { /* color atoms according to their properties given by color_ind
                         which is specified by coloraxis
                         there are two different color maps as specified by
                         plot_color_bar[0]: 1 or 2
                      */
                        alpha=(color_ind[i]-plot_color_bar[1])
                            /(plot_color_bar[2]-plot_color_bar[1]);
                        if(plot_color_bar[0]==1)
                            colormap1(alpha,&r,&g,&b);
                        else colormap2(alpha,&r,&g,&b);
                        
                        ce=win->AllocShortRGBColor(r,g,b);
                        if(plot_atom_info==5) /* print atom color (rgb) if clicked */
                            sprintf(s,"%d,%d,%d,%d,%x",i,r,g,b,ce);
                        win->DrawPoint(ri.x/L,ri.y/L,ri.z/L,
                                       atomradius[species[i]]/L,ce,s,2);
                    }                
                }
            }
            if(plot_color_windows[nw*3+1]) /* draw fixed atoms */
                if(fixed[i])   /* plot fixed atoms */
                {
                    //INFO_Printf("draw fixed atom [%d]\n",i);
                    win->DrawPoint(ri.x/L,ri.y/L,ri.z/L,
                         atomradius[species[i]]/L,colors[MAXCOLORS+3],s,2);
                }
            
        }
        else
        {
            high=false;
            for(j=1;j<=plot_highlight_atoms[0];j++)
                if(i==plot_highlight_atoms[j])
                {
                    high=true;
                    break;
                }
            
            if(fixed[i])   /* plot fixed atoms */
                win->DrawPoint(ri.x/L,ri.y/L,ri.z/L,
                               atomradius[species[i]]/L,colors[MAXCOLORS+3],s,2);
            else if (high) /* plot highlight atoms */
                win->DrawPoint(ri.x/L,ri.y/L,ri.z/L,
                               atomradius[species[i]]/L,colors[MAXCOLORS+2],s,2);
            
            else if(plot_color_bar[0]!=0)
            { /* color atoms according to their properties given by color_ind
                 which is specified by coloraxis
                 there are two different color maps as specified by
                 plot_color_bar[0]: 1 or 2
              */
                alpha=(color_ind[i]-plot_color_bar[1])
                        /(plot_color_bar[2]-plot_color_bar[1]);
                if(plot_color_bar[0]==1)
                    colormap1(alpha,&r,&g,&b);
                else colormap2(alpha,&r,&g,&b);
                
                ce=win->AllocShortRGBColor(r,g,b);
                if(plot_atom_info==5) /* print atom color (rgb) if clicked */
                    sprintf(s,"%d,%d,%d,%d,%x",i,r,g,b,ce);
                win->DrawPoint(ri.x/L,ri.y/L,ri.z/L,
                               atomradius[species[i]]/L,ce,s,2);
            }
            else
	    {/* Introduce a new variable coloratoms. Sep/25/2007 Keonwook Kang */
	        if (coloratoms == 0)
	        { /* otherwise color atoms according to their species */
	            ce=colors[MAXCOLORS+5+species[i]];
	            if(plot_atom_info==5) /* print atom species and color code if clicked */
	                sprintf(s,"%d,%d,%x",i,species[i],ce);
	            win->DrawPoint(ri.x/L,ri.y/L,ri.z/L,
	                     atomradius[species[i]]/L,ce,s,2);
	        }
	        else if (coloratoms == 1)
	        {/* otherwise color atoms according to their groups */
	            ce=colors[MAXCOLORS+5+group[i]];
	            if(plot_atom_info==5) /* print atom groups and color code if clicked */
	                sprintf(s,"%d,%d,%x",i,group[i],ce);
	            win->DrawPoint(ri.x/L,ri.y/L,ri.z/L,
	                     atomradius[species[i]]/L,ce,s,2);
	        }
            }
        }
    }
    
    /* Draw Bonds */
    if(bondlength>1e-3)
    {
        for(i=0;i<_NP;i++)
        {
            sri=_SR[i];
            if(plot_map_pbc==1) sri.subint();
            ri = _H*sri;
            
            show = 0;
            nw = (int)plot_color_windows[0];
            if(nw>0)
            {
                for(k=0;k<nw;k++)
                {
                    Emin = plot_color_windows[k*3+1];
                    Emax = plot_color_windows[k*3+2];
                    cind = (int) plot_color_windows[k*3+3];
                    
                    if((color_ind[i]>=Emin)&&
                       (color_ind[i]<=Emax))
                    {
                        show = 1;
                        break;
                    }
                }
            }
            else
            {
                show = 1;
            }
            if(!show) continue;
            
            for(jp=0;jp<nn[i];jp++)
            {
                j=nindex[i][jp];
                srj=_SR[j];
                if(plot_map_pbc==1) srj.subint();
                rj = _H*srj;

                show = 0;
                nw = (int)plot_color_windows[0];
                if(nw>0)
                {
                    for(k=0;k<nw;k++)
                    {
                        Emin = plot_color_windows[k*3+1];
                        Emax = plot_color_windows[k*3+2];
                        cind = (int) plot_color_windows[k*3+3];
                        
                        if((color_ind[j]>=Emin)&&
                           (color_ind[j]<=Emax))
                        {
                            show = 1;
                            break;
                        }
                    }
                }
                else
                {
                    show = 1;
                }
                if(!show) continue;
                
                if(plot_limits[0])
                    if((sri.x<plot_limits[1])||(sri.x>plot_limits[2])
                       ||(sri.y<plot_limits[3])||(sri.y>plot_limits[4])
                       ||(sri.z<plot_limits[5])||(sri.z>plot_limits[6])
                       ||(srj.x<plot_limits[1])||(srj.x>plot_limits[2])
                       ||(srj.y<plot_limits[3])||(srj.y>plot_limits[4])
                       ||(srj.z<plot_limits[5])||(srj.z>plot_limits[6]))
                        continue;
                r2=(ri-rj).norm2();
                //INFO_Printf(" %d - %d : %f %f\n",i,j,r2,bondlength*bondlength);
                if(r2<bondlength*bondlength)
                {
                    /* only draw O-O bonds */
                    /* if((species[i]!=0)||(species[j]!=0)) continue; */
                    //INFO_Printf("atom %d %d %d %d form bond\n",i, j, i%8,j%8); 
                    x1=ri.x/L;y1=ri.y/L;z1=ri.z/L;
                    x2=rj.x/L;y2=rj.y/L;z2=rj.z/L;
                    dx=x2-x1;dy=y2-y1;dz=z2-z1;dr=sqrt(dx*dx+dy*dy+dz*dz);
                    dx/=dr;dy/=dr;dz/=dr;
                    win->DrawLine(x1+dx*atomradius[species[i]]/L,
                                  y1+dy*atomradius[species[i]]/L,
                                  z1+dz*atomradius[species[i]]/L,
                                  x2-dx*atomradius[species[j]]/L,
                                  y2-dy*atomradius[species[j]]/L,
                                  z2-dz*atomradius[species[j]]/L,
                                  colors[MAXCOLORS+1],bondradius/L,1);
                }
            }
        }
    }
    /* draw frame */
#define drawsline(a,b,c,d,e,f,g,h,i) s1.set(a,b,c); s2.set(d,e,f);\
        vr1=_H*s1; vr2=_H*s2; vr1/=2*L; vr2/=2*L;\
        win->DrawLine(vr1.x,vr1.y,vr1.z,vr2.x,vr2.y,vr2.z,g,h,i);
    
    drawsline(-1,-1,-1,-1,-1, 1,colors[MAXCOLORS+2],0,0);
    drawsline(-1,-1, 1,-1, 1, 1,colors[MAXCOLORS+2],0,0);
    drawsline(-1, 1, 1,-1, 1,-1,colors[MAXCOLORS+2],0,0);
    drawsline(-1, 1,-1,-1,-1,-1,colors[MAXCOLORS+2],0,0);
    drawsline( 1,-1,-1, 1,-1, 1,colors[MAXCOLORS+2],0,0);
    drawsline( 1,-1, 1, 1, 1, 1,colors[MAXCOLORS+2],0,0);
    drawsline( 1, 1, 1, 1, 1,-1,colors[MAXCOLORS+2],0,0);
    drawsline( 1, 1,-1, 1,-1,-1,colors[MAXCOLORS+2],0,0);
    drawsline(-1,-1,-1, 1,-1,-1,colors[MAXCOLORS+2],0,0);
    drawsline(-1,-1, 1, 1,-1, 1,colors[MAXCOLORS+2],0,0);
    drawsline(-1, 1, 1, 1, 1, 1,colors[MAXCOLORS+2],0,0);
    drawsline(-1, 1,-1, 1, 1,-1,colors[MAXCOLORS+2],0,0);
    
    win->Unlock();
    win->Refresh();
}

/* set variables from arg list */
void MDFrame::set_dirname(const char *s)
{
    strcpy(dirname,s);
}

void MDFrame::set_randseed(const char *s)
{
    sscanf(s,"%d",&_RANDSEED);
}



/********************************************************************/
/* End of class MDFrame definition */
/********************************************************************/




/********************************************************************/
/* Main program starts below */
/********************************************************************/

#ifdef _TEST

class TESTMD : public MDFrame
{
public:
    virtual void potential()
    {
       if(_REFPOTENTIAL == 2) /* Inverse 12th potential */
       {
           I12potential();
       }
       else
       {
           INFO("Empty potential");
       }
    }
    virtual void potential_energyonly()
    {
       if(_REFPOTENTIAL == 2) /* Inverse 12th potential */
       {
           I12potential_energyonly();
       }
       else
       {
           INFO("Empty potential_energyonly");
       }
    }
    virtual double potential_energyonly(int iatom)
    {
        INFO("Empty potential");
        return 0;
    }
    virtual void initvars()
    {
        MDFrame::initvars();
        _RLIST=3.8;
    }
            
};

class TESTMD sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST
