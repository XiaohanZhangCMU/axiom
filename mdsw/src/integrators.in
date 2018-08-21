/*
 * Last Modified : Thu Nov  4 13:49:55 2010
 *
 * Purpose:
 *   Numberical integrators for MD++, included in md.cpp
 *   different ensemble_type's and integrator_type's
 *
 * To Do:
 *
 *   1. HELMP not conserved in NPH Gear6, HELMP oscilliates in NPH VVerlet
 *
 *   2. NPT VVerlet need to be implemented
 */


/*******************************************************
 *
 *  NVE ensemble
 *  Gear6 predictor-corrector integrator
 *
 ******************************************************/
#define F02      (3./16.)
#define F02_NVE  (3./20.)    // Allen and Tildesley, "Computer Simulation of Liquids" App. E 
#define F12      (251./360.)
#define F32      (11./18.)
#define F42      (1./6.)
#define F52      (1./60.)

void MDFrame::NVE_Gear6()
{
    int i;
    double mass, tmp;
    Matrix33 hinv; Vector3 serr;

    /* initializing memory */    
    if(npold!=_NP) Gear6_init();    

    /* predictor */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        _SR[i]+=_VSR[i]+s2[i]+s3[i]+s4[i]+s5[i];
        _VSR[i]+=s2[i]*2.+s3[i]*3.+s4[i]*4.+s5[i]*5.;
        s2[i]+=s3[i]*3.+s4[i]*6.+s5[i]*10.;
        s3[i]+=s4[i]*4.+s5[i]*10.;
        s4[i]+=s5[i]*5.;
    }

    /* evaluator */
    call_potential();
    calcprop();
    
    hinv=_H.inv();

    // tmp=0.5*(_TIMESTEP)*(_TIMESTEP)*1e-24 /(_ATOMMASS[0]*1e-3/AVO)*EV/1e-20;
    
    /* ATOMMASS has unit of g/mol
     * mass     has unit of eV ps^2 / A^2 
     * tmp = 0.5 * dt^2 / m  has unit of (A^2/eV) */
    mass = _ATOMMASS[0]*MASSCONVERT;
    tmp  = 0.5 / mass * _TIMESTEP*_TIMESTEP;
        
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        if(species[i]==0)
            s2real[i]=hinv*_F[i]*tmp;
        else
            s2real[i]=hinv*_F[i]*tmp *_ATOMMASS[0]/_ATOMMASS[species[i]];
    }

    /* corrector */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        serr=s2[i]-s2real[i];
        _SR[i]-=serr*F02_NVE;
        _VSR[i]-=serr*F12;
        s2[i]-=serr;
        s3[i]-=serr*F32;
        s4[i]-=serr*F42;
        s5[i]-=serr*F52;
    }
}

/*******************************************************
 *
 *  NVT ensemble
 *  Gear6 predictor-corrector integrator
 *
 ******************************************************/
void MDFrame::NVT_Gear6()
{
    int i;
    double mass, tmp, zetaerr;
    Matrix33 hinv; Vector3 serr;    

    /* initializing memory */    
    if(npold!=_NP) Gear6_init();    
    
     /* predictor */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        _SR[i]+=_VSR[i]+s2[i]+s3[i]+s4[i]+s5[i];
        _VSR[i]+=s2[i]*2.+s3[i]*3.+s4[i]*4.+s5[i]*5.;
        s2[i]+=s3[i]*3.+s4[i]*6.+s5[i]*10.;
        s3[i]+=s4[i]*4.+s5[i]*10.;
        s4[i]+=s5[i]*5.;
    }

    zeta+=zetav+zeta2+zeta3+zeta4+zeta5;
    zetav+=zeta2*2.+zeta3*3.+zeta4*4.+zeta5*5.;
    zeta2+=zeta3*3.+zeta4*6.+zeta5*10.;
    zeta3+=zeta4*4.+zeta5*10.;
    zeta4+=zeta5*5.;
    
    /* evaluator */
    call_potential();
    calcprop();
    
    hinv=_H.inv();

    // tmp=0.5*(_TIMESTEP)*(_TIMESTEP)*1e-24 /(_ATOMMASS[0]*1e-3/AVO)*EV/1e-20;
    
    /* ATOMMASS has unit of g/mol
     * mass     has unit of eV ps^2 / A^2 
     * tmp = 0.5 * dt^2 / m  has unit of (A^2/eV) */
    mass = _ATOMMASS[0]*MASSCONVERT;
    tmp  = 0.5 / mass * _TIMESTEP*_TIMESTEP;
    
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        if(species[i]==0)
            s2real[i]=hinv*_F[i]*tmp;
        else
            s2real[i]=hinv*_F[i]*tmp *_ATOMMASS[0]/_ATOMMASS[species[i]];
    }

    if(NHMass[0]==0)
    { /* old input format */
        zetaa=vt2*(_T/_TDES-1)*_TIMESTEP*_TIMESTEP*1e-24*0.5;
    }
    else
    { /* new input format: Q = NHMass = (3NKT)*1e24/vt2 (in unit of eV*ps^2)
       *  T = 300K KT = 0.025 N = 1000 3NKT = 75 --> vt2=1e28 <=> NHMass=7.5e-3)
       */
        zetaa = (_KATOM*2 - 3*_NPfree*KB*_TDES) / NHMass[0] * (0.5*_TIMESTEP*_TIMESTEP);
    }
        

    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        s2real[i]-=_VSR[i]*zetav*0.5;
    }

    /* corrector */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        serr=s2[i]-s2real[i];
        _SR[i]-=serr*F02;
        _VSR[i]-=serr*F12;
        s2[i]-=serr;
        s3[i]-=serr*F32;
        s4[i]-=serr*F42;
        s5[i]-=serr*F52;
    }
    zetaerr=zeta2-zetaa;
    zeta-=zetaerr*F02;
    zetav-=zetaerr*F12;
    zeta2-=zetaerr;
    zeta3-=zetaerr*F32;
    zeta4-=zetaerr*F42;
    zeta5-=zetaerr*F52;
    
}

/*******************************************************
 *
 *  NVT ensemble (Nose Hoover Chain)
 *  Gear6 predictor-corrector integrator
 *
 ******************************************************/
void MDFrame::NVTC_Gear6()
{
    int i;
    double mass, tmp, zetaNHCerr[MAXNHCLEN];
    Matrix33 hinv; Vector3 serr;    

    //INFO_Printf("NVHC_Gear6: Hello!!!!  NHChainLen = %d\n",NHChainLen);
    
    /* initializing memory */
    if(npold!=_NP) Gear6_init();
    
     /* predictor */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        _SR[i]+=_VSR[i]+s2[i]+s3[i]+s4[i]+s5[i];
        _VSR[i]+=s2[i]*2.+s3[i]*3.+s4[i]*4.+s5[i]*5.;
        s2[i]+=s3[i]*3.+s4[i]*6.+s5[i]*10.;
        s3[i]+=s4[i]*4.+s5[i]*10.;
        s4[i]+=s5[i]*5.;
    }

    for(i=0;i<NHChainLen;i++)
    {
        zetaNHC[i] +=zetaNHCv[i]+zetaNHC2[i]+zetaNHC3[i]+zetaNHC4[i]+zetaNHC5[i];
        zetaNHCv[i]+=zetaNHC2[i]*2.+zetaNHC3[i]*3.+zetaNHC4[i]*4.+zetaNHC5[i]*5.;
        zetaNHC2[i]+=zetaNHC3[i]*3.+zetaNHC4[i]*6.+zetaNHC5[i]*10.;
        zetaNHC3[i]+=zetaNHC4[i]*4.+zetaNHC5[i]*10.;
        zetaNHC4[i]+=zetaNHC5[i]*5.;
    }

    /* evaluator */
    call_potential();
    calcprop();
    
    hinv=_H.inv();

    // tmp=0.5*(_TIMESTEP)*(_TIMESTEP)*1e-24 /(_ATOMMASS[0]*1e-3/AVO)*EV/1e-20;     //(A^2)

    /* ATOMMASS has unit of g/mol
     * mass     has unit of eV ps^2 / A^2 
     * tmp = 0.5 * dt^2 / m  has unit of (A^2/eV) */
    mass = _ATOMMASS[0]*MASSCONVERT;
    tmp  = 0.5 / mass * _TIMESTEP*_TIMESTEP;
    
    
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        if(species[i]==0)
            s2real[i]=hinv*_F[i]*tmp;
        else
            s2real[i]=hinv*_F[i]*tmp *_ATOMMASS[0]/_ATOMMASS[species[i]];
    }

    //INFO_Printf("NVHC: NHChainLen = %d\n",NHChainLen);
    
    if(NHChainLen<=1)
    {
        /* new input format: Q = NHMass = (3NKT)*1e24/vt2 (in unit of eV*ps^2)
         *  T = 300K KT = 0.025 N = 1000 3NKT = 75 --> vt2=1e28 <=> NHMass=7.5e-3)
         */
        zetaNHCa[0] = (_KATOM*2 - 3*_NPfree*KB*_TDES) / NHMass[0] * (0.5*_TIMESTEP*_TIMESTEP);
    }
    else
    {    
        /* new input format: Q = NHMass = (3NKT)*1e24/vt2 (in unit of eV*ps^2)
         *  T = 300K KT = 0.025 N = 1000 3NKT = 75 --> vt2=1e28 <=> NHMass=7.5e-3)
         */
        zetaNHCa[0] = (_KATOM*2 - 3*_NPfree*KB*_TDES) / NHMass[0] * (0.5*_TIMESTEP*_TIMESTEP)
                      - zetaNHCv[0]*zetaNHCv[1] * 0.5;
                         
        for(i=1;i<NHChainLen-1;i++)
        {
            zetaNHCa[i] = (NHMass[i-1]*SQR(zetaNHCv[i-1]/_TIMESTEP)-KB*_TDES) / NHMass[i] * (0.5*_TIMESTEP*_TIMESTEP)
                          -  zetaNHCv[i]*zetaNHCv[i+1] * 0.5;
        }
        i = NHChainLen-1;
        zetaNHCa[i] = (NHMass[i-1]*SQR(zetaNHCv[i-1]/_TIMESTEP)-KB*_TDES) / NHMass[i] * (0.5*_TIMESTEP*_TIMESTEP);       
    }
    
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        s2real[i]-=_VSR[i]*zetaNHCv[0]*0.5;
    }

    /* corrector */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        serr=s2[i]-s2real[i];
        _SR[i]-=serr*F02;
        _VSR[i]-=serr*F12;
        s2[i]-=serr;
        s3[i]-=serr*F32;
        s4[i]-=serr*F42;
        s5[i]-=serr*F52;
    }


    for(i=0;i<NHChainLen;i++)
    {
        zetaNHCerr[i] = zetaNHC2[i]-zetaNHCa[i];
    }

    for(i=0;i<NHChainLen;i++)
    {
        zetaNHC[i]  -= zetaNHCerr[i]*F02;
        zetaNHCv[i] -= zetaNHCerr[i]*F12;
        zetaNHC2[i] -= zetaNHCerr[i];
        zetaNHC3[i] -= zetaNHCerr[i]*F32;
        zetaNHC4[i] -= zetaNHCerr[i]*F42;
        zetaNHC5[i] -= zetaNHCerr[i]*F52;        
    }

}

/*******************************************************
 *
 *  NPH ensemble
 *  Gear6 predictor-corrector integrator
 *
 ******************************************************/
void MDFrame::NPH_Gear6()
{
    int i, j;
    double mass, tmp;
    Matrix33 hinv, htran, vhtran, G, Ginv, VG, VG1, VG2, GiVG, herr;
    Vector3 serr;    

    /* initializing memory */    
    if(npold!=_NP) Gear6_init();    
    
     /* predictor */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        _SR[i]+=_VSR[i]+s2[i]+s3[i]+s4[i]+s5[i];
        _VSR[i]+=s2[i]*2.+s3[i]*3.+s4[i]*4.+s5[i]*5.;
        s2[i]+=s3[i]*3.+s4[i]*6.+s5[i]*10.;
        s3[i]+=s4[i]*4.+s5[i]*10.;
        s4[i]+=s5[i]*5.;
    }
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            if(conj_fixboxvec[i][j]==1)
                _VH[i][j]=h2[i][j]=h3[i][j]=h4[i][j]=h5[i][j]=0;

    _H+=_VH+h2+h3+h4+h5;
    _VH+= h2*2.+ h3*3.+ h4*4.+ h5*5.;
    h2+= h3*3.+ h4*6.+ h5*10.;
    h3+= h4*4.+ h5*10.;
    h4+= h5*5.;

    /* evaluator */
    call_potential();
    calcprop();
    
    hinv=_H.inv();

    // tmp=0.5*(_TIMESTEP)*(_TIMESTEP)*1e-24 /(_ATOMMASS[0]*1e-3/AVO)*EV/1e-20;     //(A^2)

    /* ATOMMASS has unit of g/mol
     * mass     has unit of eV ps^2 / A^2 
     * tmp = 0.5 * dt^2 / m  has unit of (A^2/eV) */
    mass = _ATOMMASS[0]*MASSCONVERT;
    tmp  = 0.5 / mass * _TIMESTEP*_TIMESTEP;
    
    
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        if(species[i]==0)
            s2real[i]=hinv*_F[i]*tmp;
        else
            s2real[i]=hinv*_F[i]*tmp *_ATOMMASS[0]/_ATOMMASS[species[i]];
    }
    h2real=_GH* tmp * _ATOMMASS[0]/_WALLMASS;

    if(_BOXDAMP!=0)
        for(i=0;i<3;i++)
            for(j=0;j<3;j++)
                if(conj_fixboxvec[i][j]==0)
                    h2real[i][j]-=_VH[i][j]*_BOXDAMP;

    htran=_H.tran(); G=htran*_H; Ginv=G.inv();
    vhtran=_VH.tran(); VG1=vhtran*_H; VG2=htran*_VH;
    VG=VG1+VG2; GiVG=Ginv*VG;

    for(i=0;i<_NP;i++)
        s2real[i]-=GiVG*_VSR[i];

    /* corrector */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        serr=s2[i]-s2real[i];
        _SR[i]-=serr*F02;
        _VSR[i]-=serr*F12;
        s2[i]-=serr;
        s3[i]-=serr*F32;
        s4[i]-=serr*F42;
        s5[i]-=serr*F52;
    }    

    /* Parrinello-Rahman boundary condition */
    herr=h2-h2real;
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            if(conj_fixboxvec[i][j]==0)
            {
                 _H[i][j] -= herr[i][j]*F02;
                _VH[i][j] -= herr[i][j]*F12;
                 h2[i][j] -= herr[i][j];
                 h3[i][j] -= herr[i][j]*F32;
                 h4[i][j] -= herr[i][j]*F42;
                 h5[i][j] -= herr[i][j]*F52;
            }

}


/*******************************************************
 *
 *  NPT ensemble
 *  Gear6 predictor-corrector integrator
 *
 ******************************************************/
void MDFrame::NPT_Gear6()
{
    int i, j;
    double mass, tmp, zetaerr;
    Matrix33 hinv, htran, vhtran, G, Ginv, VG, VG1, VG2, GiVG, herr;
    Vector3 serr;    
    
    /* initializing memory */    
    if(npold!=_NP) Gear6_init();    
     /* predictor */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        _SR[i]+=_VSR[i]+s2[i]+s3[i]+s4[i]+s5[i];
        _VSR[i]+=s2[i]*2.+s3[i]*3.+s4[i]*4.+s5[i]*5.;
        s2[i]+=s3[i]*3.+s4[i]*6.+s5[i]*10.;
        s3[i]+=s4[i]*4.+s5[i]*10.;
        s4[i]+=s5[i]*5.;
    }
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            if(conj_fixboxvec[i][j]==1)
                _VH[i][j]=h2[i][j]=h3[i][j]=h4[i][j]=h5[i][j]=0;

    zeta+=zetav+zeta2+zeta3+zeta4+zeta5;
    zetav+=zeta2*2.+zeta3*3.+zeta4*4.+zeta5*5.;
    zeta2+=zeta3*3.+zeta4*6.+zeta5*10.;
    zeta3+=zeta4*4.+zeta5*10.;
    zeta4+=zeta5*5.;
    
    _H+=_VH+h2+h3+h4+h5;
    _VH+= h2*2.+ h3*3.+ h4*4.+ h5*5.;
    h2+= h3*3.+ h4*6.+ h5*10.;
    h3+= h4*4.+ h5*10.;
    h4+= h5*5.;
    
    /* evaluator */
    call_potential();
    calcprop();
    hinv=_H.inv();

    // tmp=0.5*(_TIMESTEP)*(_TIMESTEP)*1e-24 /(_ATOMMASS[0]*1e-3/AVO)*EV/1e-20;     //(A^2)

    /* ATOMMASS has unit of g/mol
     * mass     has unit of eV ps^2 / A^2 
     * tmp = 0.5 * dt^2 / m  has unit of (A^2/eV) */
    mass = _ATOMMASS[0]*MASSCONVERT;
    tmp  = 0.5 / mass * _TIMESTEP*_TIMESTEP;
    
    
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        if(species[i]==0)
            s2real[i]=hinv*_F[i]*tmp;
        else
            s2real[i]=hinv*_F[i]*tmp *_ATOMMASS[0]/_ATOMMASS[species[i]];
    }

    zetaa=vt2*(_T/_TDES-1)*_TIMESTEP*_TIMESTEP*1e-24*0.5;

    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        s2real[i]-=_VSR[i]*zetav*0.5;
    }

    h2real=_GH* tmp * _ATOMMASS[0]/_WALLMASS;
    
    if(_BOXDAMP!=0)
        for(i=0;i<3;i++)
            for(j=0;j<3;j++)
                if(conj_fixboxvec[i][j]==0)
                    h2real[i][j]-=_VH[i][j]*_BOXDAMP;

    htran=_H.tran(); G=htran*_H; Ginv=G.inv();
    vhtran=_VH.tran(); VG1=vhtran*_H; VG2=htran*_VH;
    VG=VG1+VG2; GiVG=Ginv*VG;

    for(i=0;i<_NP;i++)
        s2real[i]-=GiVG*_VSR[i];
    
    /* corrector */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        serr=s2[i]-s2real[i];
        _SR[i]-=serr*F02;
        _VSR[i]-=serr*F12;
        s2[i]-=serr;
        s3[i]-=serr*F32;
        s4[i]-=serr*F42;
        s5[i]-=serr*F52;
    }
    zetaerr=zeta2-zetaa;
    zeta-=zetaerr*F02;
    zetav-=zetaerr*F12;
    zeta2-=zetaerr;
    zeta3-=zetaerr*F32;
    zeta4-=zetaerr*F42;
    zeta5-=zetaerr*F52;

    herr=h2-h2real;
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            if(conj_fixboxvec[i][j]==0)
            {
                 _H[i][j] -= herr[i][j]*F02;
                _VH[i][j] -= herr[i][j]*F12;
                 h2[i][j] -= herr[i][j];
                 h3[i][j] -= herr[i][j]*F32;
                 h4[i][j] -= herr[i][j]*F42;
                 h5[i][j] -= herr[i][j]*F52;
            }        
}

/*******************************************************
 *
 *  NPT ensemble (Nose Hoover Chain)
 *  Gear6 predictor-corrector integrator
 *
 ******************************************************/
void MDFrame::NPTC_Gear6()
{
    int i, j; //,boxdof;
    double mass, tmp, zetaNHCerr[MAXNHCLEN]; //, zetaBNHCerr[MAXNHCLEN];
    Matrix33 hinv, htran, vhtran, G, Ginv, VG, VG1, VG2, GiVG, herr;
    Vector3 serr;    

    /* initializing memory */    
    if(npold!=_NP) Gear6_init();    
    
     /* predictor */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        _SR[i]+=_VSR[i]+s2[i]+s3[i]+s4[i]+s5[i];
        _VSR[i]+=s2[i]*2.+s3[i]*3.+s4[i]*4.+s5[i]*5.;
        s2[i]+=s3[i]*3.+s4[i]*6.+s5[i]*10.;
        s3[i]+=s4[i]*4.+s5[i]*10.;
        s4[i]+=s5[i]*5.;
    }
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            if(conj_fixboxvec[i][j]==1)
                _VH[i][j]=h2[i][j]=h3[i][j]=h4[i][j]=h5[i][j]=0;

    /* Nose Hoover Chain */
    for(i=0;i<NHChainLen;i++)
    {
        zetaNHC[i] +=zetaNHCv[i]+zetaNHC2[i]+zetaNHC3[i]+zetaNHC4[i]+zetaNHC5[i];
        zetaNHCv[i]+=zetaNHC2[i]*2.+zetaNHC3[i]*3.+zetaNHC4[i]*4.+zetaNHC5[i]*5.;
        zetaNHC2[i]+=zetaNHC3[i]*3.+zetaNHC4[i]*6.+zetaNHC5[i]*10.;
        zetaNHC3[i]+=zetaNHC4[i]*4.+zetaNHC5[i]*10.;
        zetaNHC4[i]+=zetaNHC5[i]*5.;
    }

    /* Nose Hoover Chain for Box */
//    for(i=0;i<NHChainLen;i++)
//    {
//        zetaBNHC[i] +=zetaBNHCv[i]+zetaBNHC2[i]+zetaBNHC3[i]+zetaBNHC4[i]+zetaBNHC5[i];
//        zetaBNHCv[i]+=zetaBNHC2[i]*2.+zetaBNHC3[i]*3.+zetaBNHC4[i]*4.+zetaBNHC5[i]*5.;
//        zetaBNHC2[i]+=zetaBNHC3[i]*3.+zetaBNHC4[i]*6.+zetaBNHC5[i]*10.;
//        zetaBNHC3[i]+=zetaBNHC4[i]*4.+zetaBNHC5[i]*10.;
//        zetaBNHC4[i]+=zetaBNHC5[i]*5.;
//    }
    
    
    _H+=_VH+h2+h3+h4+h5;
    _VH+= h2*2.+ h3*3.+ h4*4.+ h5*5.;
    h2+= h3*3.+ h4*6.+ h5*10.;
    h3+= h4*4.+ h5*10.;
    h4+= h5*5.;
    
    /* evaluator */
    call_potential();
    calcprop();
    
    hinv=_H.inv();

    // tmp=0.5*(_TIMESTEP)*(_TIMESTEP)*1e-24 /(_ATOMMASS[0]*1e-3/AVO)*EV/1e-20;     //(A^2)

    /* ATOMMASS has unit of g/mol
     * mass     has unit of eV ps^2 / A^2 
     * tmp = 0.5 * dt^2 / m  has unit of (A^2/eV) */
    mass = _ATOMMASS[0]*MASSCONVERT;
    tmp  = 0.5 / mass * _TIMESTEP*_TIMESTEP;
    
    
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        if(species[i]==0)
            s2real[i]=hinv*_F[i]*tmp;
        else
            s2real[i]=hinv*_F[i]*tmp *_ATOMMASS[0]/_ATOMMASS[species[i]];
    }

    if(NHChainLen<=1)
    {
        /* new input format: Q = NHMass = (3NKT)*1e24/vt2 (in unit of eV*ps^2)
         *  T = 300K KT = 0.025 N = 1000 3NKT = 75 --> vt2=1e28 <=> NHMass=7.5e-3)
         */
        zetaNHCa[0] = (_KATOM*2 - 3*_NPfree*KB*_TDES) / NHMass[0] * (0.5*_TIMESTEP*_TIMESTEP);

//        boxdof = 9; for(i=0;i<3;i++) for(j=0;j<3;j++) if(conj_fixboxvec[i][j]) boxdof--;
//        zetaBNHCa[0] = (_KBOX*2 - boxdof*KB*_TDES) / BNHMass[0] * (0.5*_TIMESTEP*_TIMESTEP);        
    }
    else
    {    
        /* new input format: Q = NHMass = (3NKT)*1e24/vt2 (in unit of eV*ps^2)
         *  T = 300K KT = 0.025 N = 1000 3NKT = 75 --> vt2=1e28 <=> NHMass=7.5e-3)
         */
        zetaNHCa[0] = (_KATOM*2 - 3*_NPfree*KB*_TDES) / NHMass[0] * (0.5*_TIMESTEP*_TIMESTEP)
                      - zetaNHCv[0]*zetaNHCv[1] * 0.5;

//        INFO_Printf("KBOX = %20.12e\n",_KBOX);
//        boxdof = 9; for(i=0;i<3;i++) for(j=0;j<3;j++) if(conj_fixboxvec[i][j]) boxdof--;
//        zetaBNHCa[0] = (_KBOX*2 - boxdof*KB*_TDES) / BNHMass[0] * (0.5*_TIMESTEP*_TIMESTEP)
//                      - zetaBNHCv[0]*zetaBNHCv[1] * 0.5;
                
        
        for(i=1;i<NHChainLen-1;i++)
        {
            zetaNHCa[i] = (NHMass[i-1]*SQR(zetaNHCv[i-1]/_TIMESTEP)-KB*_TDES) / NHMass[i] * (0.5*_TIMESTEP*_TIMESTEP)
                          -  zetaNHCv[i]*zetaNHCv[i+1] * 0.5;
//            zetaBNHCa[i] = (BNHMass[i-1]*SQR(zetaBNHCv[i-1]/_TIMESTEP)-KB*_TDES) / BNHMass[i] * (0.5*_TIMESTEP*_TIMESTEP)
//                          -  zetaBNHCv[i]*zetaBNHCv[i+1] * 0.5;            
        }
        i = NHChainLen-1;
        zetaNHCa[i] = (NHMass[i-1]*SQR(zetaNHCv[i-1]/_TIMESTEP)-KB*_TDES) / NHMass[i] * (0.5*_TIMESTEP*_TIMESTEP);       
//        zetaBNHCa[i] = (BNHMass[i-1]*SQR(zetaBNHCv[i-1]/_TIMESTEP)-KB*_TDES) / BNHMass[i] * (0.5*_TIMESTEP*_TIMESTEP);       
    }

    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        s2real[i]-=_VSR[i]*zetaNHCv[0]*0.5;
    }

    h2real=_GH * tmp * _ATOMMASS[0]/_WALLMASS;

//    for(i=0;i<3;i++)
//        for(j=0;j<3;j++)
//            if(conj_fixboxvec[i][j]==0)
//                h2real[i][j]-=_VH[i][j]*zetaBNHCv[0];
    
    if(_BOXDAMP!=0)
        for(i=0;i<3;i++)
            for(j=0;j<3;j++)
                if(conj_fixboxvec[i][j]==0)
                    h2real[i][j]-=_VH[i][j]*_BOXDAMP;

    htran=_H.tran(); G=htran*_H; Ginv=G.inv();
    vhtran=_VH.tran(); VG1=vhtran*_H; VG2=htran*_VH;
    VG=VG1+VG2; GiVG=Ginv*VG;

    for(i=0;i<_NP;i++)
        s2real[i]-=GiVG*_VSR[i];
    
    /* corrector */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        serr=s2[i]-s2real[i];
        _SR[i]-=serr*F02;
        _VSR[i]-=serr*F12;
        s2[i]-=serr;
        s3[i]-=serr*F32;
        s4[i]-=serr*F42;
        s5[i]-=serr*F52;
    }

    for(i=0;i<NHChainLen;i++)
    {
        zetaNHCerr[i] = zetaNHC2[i]-zetaNHCa[i];
//        zetaBNHCerr[i] = zetaBNHC2[i]-zetaBNHCa[i];
    }

    for(i=0;i<NHChainLen;i++)
    {
        zetaNHC[i]  -= zetaNHCerr[i]*F02;
        zetaNHCv[i] -= zetaNHCerr[i]*F12;
        zetaNHC2[i] -= zetaNHCerr[i];
        zetaNHC3[i] -= zetaNHCerr[i]*F32;
        zetaNHC4[i] -= zetaNHCerr[i]*F42;
        zetaNHC5[i] -= zetaNHCerr[i]*F52;

//        zetaBNHC[i]  -= zetaBNHCerr[i]*F02;
//        zetaBNHCv[i] -= zetaBNHCerr[i]*F12;
//        zetaBNHC2[i] -= zetaBNHCerr[i];
//        zetaBNHC3[i] -= zetaBNHCerr[i]*F32;
//        zetaBNHC4[i] -= zetaBNHCerr[i]*F42;
//        zetaBNHC5[i] -= zetaBNHCerr[i]*F52;                
    }
    

    herr=h2-h2real;
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            if(conj_fixboxvec[i][j]==0)
            {
                 _H[i][j] -= herr[i][j]*F02;
                _VH[i][j] -= herr[i][j]*F12;
                 h2[i][j] -= herr[i][j];
                 h3[i][j] -= herr[i][j]*F32;
                 h4[i][j] -= herr[i][j]*F42;
                 h5[i][j] -= herr[i][j]*F52;
            }        
}



/* Gear6 initializing memory */   
void MDFrame::Gear6_init()
{
    Realloc(s2,Vector3,_NP);
    Realloc(s3,Vector3,_NP);
    Realloc(s4,Vector3,_NP);
    Realloc(s5,Vector3,_NP);
    Realloc(s2real,Vector3,_NP);
    
    memset(s2,0,sizeof(Vector3)*_NP);
    memset(s3,0,sizeof(Vector3)*_NP);
    memset(s4,0,sizeof(Vector3)*_NP);
    memset(s5,0,sizeof(Vector3)*_NP);
    memset(s2real,0,sizeof(Vector3)*_NP);

    memset(zetaNHC,0,sizeof(double)*MAXNHCLEN);
    memset(zetaNHCv,0,sizeof(double)*MAXNHCLEN);
    memset(zetaNHCa,0,sizeof(double)*MAXNHCLEN);
    memset(zetaNHC2,0,sizeof(double)*MAXNHCLEN);
    memset(zetaNHC3,0,sizeof(double)*MAXNHCLEN);
    memset(zetaNHC4,0,sizeof(double)*MAXNHCLEN);
    memset(zetaNHC5,0,sizeof(double)*MAXNHCLEN);
    
//    memset(zetaBNHC,0,sizeof(double)*MAXNHCLEN);
//    memset(zetaBNHCv,0,sizeof(double)*MAXNHCLEN);
//    memset(zetaBNHCa,0,sizeof(double)*MAXNHCLEN);
//    memset(zetaBNHC2,0,sizeof(double)*MAXNHCLEN);
//    memset(zetaBNHC3,0,sizeof(double)*MAXNHCLEN);
//    memset(zetaBNHC4,0,sizeof(double)*MAXNHCLEN);
//    memset(zetaBNHC5,0,sizeof(double)*MAXNHCLEN);

    npold=_NP;
}





/*******************************************************
 *
 *  NVE ensemble
 *  Velocity Verlet integrator
 *
 ******************************************************/
void MDFrame::NVE_VVerlet(int nstep)
{
    int i;
    
    if(npold!=_NP)
    {
        Realloc(s2,Vector3,_NP);
        memset(s2,0,sizeof(Vector3)*_NP);
        npold=_NP;
    }

    /* Jan 22 2007 Keonwook Kang */
    if(nstep == step0) VVerlet_Get_s2();
    
    /* s2[i] is the acceleration in previous time step (n) */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        _VSR[i]+=s2[i];        /* p(n+1/2) */
        _SR[i]+=_VSR[i];       /* q(n+1) */
    }
    /* _SR[i] is the position at current step,  q(n+1) */
    /* _VSR[i] is the velocity at half-step,    p(n+1/2) */
    
    VVerlet_Get_s2();
        
    /* s2[i] is the acceleration at current step, 0.5*a(n+1)*(dt)^2 */
    
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        _VSR[i]+=s2[i];
    }
    /* _VSR[i] is the velocity at current step, p(n+1) */
}


/*******************************************************
 *
 *  NVT ensemble
 *  Velocity Verlet integrator, implicit algorithm
 *
 ******************************************************/
void MDFrame::NVT_VVerlet_Implicit(int nstep)
{
    int i, maxiter;
    double tmp, xinew, B, C;
    Vector3 dvs;
    
    if(npold!=_NP)
    {
        Realloc(s2,Vector3,_NP);
        memset(s2,0,sizeof(Vector3)*_NP);
        npold=_NP;
    }

    /* Jan 22 2007 Keonwook Kang */
    if(nstep == step0)  VVerlet_Get_s2();

    zetaa=vt2* (_T/_TDES-1) *_TIMESTEP*_TIMESTEP*1e-24*0.5;
    /* zetaa = xi'(n) *(dt)^2 *0.5 */

    /* s2[i] is the acceleration in previous time step (n) */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;                
        dvs=_VSR[i]; dvs*=zetav*0.5; /* viscous force */
        _VSR[i]+=s2[i]; _VSR[i]-=dvs; /* p(n+1/2) */
        _SR[i]+=_VSR[i];              /* q(n+1) */
    }
    zetav+=zetaa;  /* zetav = xi(n+1/2) *(dt) */
    zeta +=zetav;  /* zeta  = eta(n+1) */
    /* _SR[i] is the position at current step,  q(n+1) */
    /* _VSR[i] is the velocity at half-step,    p(n+1/2) */
    
    VVerlet_Get_s2();
    
    /* s2[i] is the acceleration at current step, 0.5*a(n+1)*(dt)^2 */
    
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        _VSR[i]+=s2[i];
    }
    /* _VSR[i] is the velocity at current step with an unknown factor */
    /*    p(n+1) * (1 + xi(n+1)*dt/2) */
    /* need to solve for xi(n+1) */
        
    /* first compute instananeous temperature */
    VVerlet_Get_Tinst();
        
    /* xi(n+1) satisfies the equation:
     *
     *   xi(n+1) = xi(n+1/2) + vt2*(_T/_TDES/(1+0.5*xi(n+1))^2 - 1)*0.5
     *           = xi(n+1/2) + (B/(1+0.5*xi(n+1))^2 - 1)*C
     */
    maxiter = 1000;
    B = _T/_TDES;
    C = vt2*_TIMESTEP*_TIMESTEP*1e-24*0.5;
    xinew = zetav;
    for(i=0;i<maxiter;i++)
    {
        tmp = zetav + (B/(1+0.5*xinew)/(1+0.5*xinew) - 1)*C;
        if(fabs(tmp-xinew)<1e-14)
            break;
        else
            xinew = tmp;
    }
        
    if(i>=maxiter)
    {
        INFO_Printf("maxiter = %d xinew = %e %e\n",maxiter,xinew,tmp);
        ERROR("maximum iteration exceeded");
        zetav = xinew;        
    }
    else
    {
        //INFO_Printf("zetav_old = %e zetav = %e iter = %d\n",zetav,xinew,i);
        zetav = xinew;        
    }
        
    tmp = 1.0/(1+zetav*0.5);
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        _VSR[i]*=tmp;
    }
    tmp = tmp*tmp; _KATOM *= tmp; _T *= tmp;
    
}

/*******************************************************
 *
 *  NVT ensemble
 *  Velocity Verlet integrator, explicit algorithm 1
 *
 ******************************************************/
void MDFrame::NVT_VVerlet_Explicit_1(int nstep)
{
    int i;
    Vector3 dvs;
    
    if(npold!=_NP)
    {
        Realloc(s2,Vector3,_NP);
        memset(s2,0,sizeof(Vector3)*_NP);
        npold=_NP;
    }

    /* Jan 22 2007 Keonwook Kang */    
    if(nstep == step0) VVerlet_Get_s2();

    
    /****************************************************/
    /* s2[i] is the acceleration in previous time step (n) */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;                
	_VSR[i]+=s2[i]; _VSR[i]/=(1+zetav*.5); /* Keonwook Kang Feb 12 2007*/
        //dvs=_VSR[i]; dvs*=zetav*0.5; /* viscous force */
        //_VSR[i]+=s2[i]; _VSR[i]-=dvs; /* p(n+1/2) */
        _SR[i]+=_VSR[i];              /* q(n+1) */
    }

    /* compute instananeous temperature at half step*/
    VVerlet_Get_Tinst();

    zetaa=vt2* (_T/_TDES-1) *_TIMESTEP*_TIMESTEP*1e-24*0.5;
    /* zetaa = xi'(n) *(dt)^2 *0.5 */

    zeta +=zetav+zetaa;/* zeta  = eta(n+1) */
    zetav+=zetaa*2.0;  /* zetav = xi(n+1/2) *(dt) */

    /* Above two lines are originally given in the following
       three steps.                                     */
    /* Keonwook Kang Feb 13 2007*/
//    zeta +=zetav*.5;
//    zetav+=zetaa*2.0;  /* zetav = xi(n+1/2) *(dt) */
//    zeta +=zetav*.5;   /* zeta  = eta(n+1) */

    /****************************************************/
    /* compute acceleration a(n+1),  s2[i]=0.5*a*(dt)^2 */
    VVerlet_Get_s2();

    /****************************************************/    
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        dvs=_VSR[i]; dvs*=zetav*0.5; /* viscous force */
        _VSR[i]+=s2[i]; _VSR[i]-=dvs; /* p(n+1/2) */
        }
    /* _VSR[i] is the velocity at current step,    p(n+1) */

}

/*******************************************************
 *
 *  NVT ensemble
 *  Velocity Verlet integrator, explicit algorithm 2
 *  Jan 23 2007 Keonwook Kang
 ******************************************************/
void MDFrame::NVT_VVerlet_Explicit_2(int nstep)
{
    int i;
    Vector3 dvs;

    
    if(npold!=_NP)
    {
        Realloc(s2,Vector3,_NP);
        memset(s2,0,sizeof(Vector3)*_NP);
        npold=_NP;
    }

    /* Jan 22 2007 Keonwook Kang */    
    if(nstep == step0) VVerlet_Get_s2();

    /************************************************/
    /* compute instananeous temperature at half step*/
    VVerlet_Get_Tinst();

    zetaa=vt2* (_T/_TDES-1) *_TIMESTEP*_TIMESTEP*1e-24*0.5;
    /* zetaa = xi'(n) *(dt)^2 *0.5 */

    zetav+=zetaa;  /* zetav = xi(n+1/2) *(dt) */        
    zeta +=zetav;  /* zeta  = eta(n+1) */
    /* Originally two steps: zeta += zetav*.5; zeta += zetav*.5 */
    
    /* s2[i] is the acceleration in previous time step a(n)*(dt)^2*0.5 */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        dvs=_VSR[i]; dvs*=exp(-zetav*0.5); /* viscous force */
        _VSR[i]=dvs + s2[i]; /* p(n+1/2) */            
        _SR[i]+=_VSR[i];              /* q(n+1) */
    }
    
    VVerlet_Get_s2();
    
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        _VSR[i] += s2[i];
        _VSR[i] *= exp(-zetav*0.5); /* viscous force */
    }

//        zeta +=zetav*0.5;/* zeta  = eta(n+1) */                
    /* compute instananeous temperature at half step*/
    VVerlet_Get_Tinst();

    zetaa=vt2* (_T/_TDES-1) *_TIMESTEP*_TIMESTEP*1e-24*0.5;
    /* zetaa = xi'(n) *(dt)^2 *0.5 */

    zetav+=zetaa;  /* zetav = xi(n+1/2) *(dt) */
}

/*******************************************************
 *
 *  NVT ensemble Nose-Hoover Chain
 *  Velocity Verlet integrator, 
 *  Mar 20 2007 implemented by Keonwook Kang
 *  Ref. G. J. Martyna, Molecular Physics (1996) 87 1117-1157
 *
 ******************************************************/
void MDFrame::NVTC_VVerlet_Explicit(int nstep)
{
    int i;
    Vector3 dvs;
    
    if(npold!=_NP) /* initialization */
    {
        Realloc(s2,Vector3,_NP);
        memset(s2,0,sizeof(Vector3)*_NP);
        memset(zetaNHC,0,sizeof(double)*MAXNHCLEN);
        memset(zetaNHCv,0,sizeof(double)*MAXNHCLEN);
        memset(zetaNHCa,0,sizeof(double)*MAXNHCLEN);
        npold=_NP;
    }

    VVerlet_NHCINT(); /* Nose-Hoover Chain integrator */

    /* Jan 22 2007 Keonwook Kang */    
    if(nstep == step0) VVerlet_Get_s2();
    /* s2[i] is the acceleration in previous time step a(n)*(dt)^2*0.5 */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        _VSR[i] += s2[i]; /* p(n+1/2) */            
        _SR[i] += _VSR[i];              /* q(n+1) */
    }
    
    VVerlet_Get_s2();
    
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        _VSR[i] += s2[i];
    }

    VVerlet_NHCINT(); /* Nose-Hoover Chain integrator */
}

/*******************************************************
 *
 *  NPH ensemble
 *  Velocity Verlet integrator
 *
 ******************************************************/
void MDFrame::NPH_VVerlet(int nstep)
{
    int i, j;
    Matrix33 hinv;
    
    if(npold!=_NP)
    {
        Realloc(s2,Vector3,_NP);
        memset(s2,0,sizeof(Vector3)*_NP);
        npold=_NP;
    }

    /* Jan 22 2007 Keonwook Kang */    
    if(nstep == step0)
    {
        VVerlet_Get_s2();
        VVerlet_Get_h2();
    }

//    for(i=0;i<_NP;i++)
//        s2[i]-=GiVG*_VSR[i];
    
    /* s2[i] is the acceleration in previous time step (n) */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        _VSR[i]+=s2[i];        /* p(n+1/2) */
        _SR[i]+=_VSR[i];       /* q(n+1) */
    }

    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            if(conj_fixboxvec[i][j]==1)
                _VH[i][j]=h2[i][j]=0;
    _VH+=h2;
    _H+=_VH;
    /* _SR[i] is the position at current step,  q(n+1) */
    /* _VSR[i] is the velocity at half-step,    p(n+1/2) */

    VVerlet_Get_s2();
    
    /* s2[i] is the acceleration at current step, 0.5*a(n+1)*(dt)^2 */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        _VSR[i]+=s2[i];
    }
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            if(conj_fixboxvec[i][j]==1)
                _VH[i][j]=h2[i][j]=0;
    _VH+=h2;
    /* _VSR[i] is the velocity at current step, p(n+1) */
}


/*******************************************************
 *
 *  NPH ensemble
 *  Velocity Verlet integrator
 *  Feb 15 2007 Keonwook Kang
 *  Finished Feb 19 2007
 ******************************************************/
void MDFrame::NPH_VVerlet_Explicit_1(int nstep)
{
    int i, j, n;
    Matrix33 htran, G, Ginv, vhtran, VG1, VG2, VG;
    Matrix33 GiVG, GiVGoffDiag, I, Iinv;
    Vector3 dvs;

    I.eye();
    
    if(npold!=_NP)
    {
        Realloc(s2,Vector3,_NP);
        memset(s2,0,sizeof(Vector3)*_NP);
        npold=_NP;
    }

    if(nstep == step0) VVerlet_Get_h2();

    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            if(conj_fixboxvec[i][j]==1)
                _VH[i][j]=h2[i][j]=0;
    _VH+=h2;
    _H+=_VH*.5;

    VVerlet_Get_s2();   /* eval s2 due to box change */
    
    htran = _H.tran(); G=htran*_H; Ginv = G.inv();
    vhtran=_VH.tran(); VG1=vhtran*_H; VG2=htran*_VH;
    VG=VG1+VG2; GiVG = Ginv*VG;
    GiVGoffDiag.set(0, GiVG[0][1], GiVG[0][2],
                    GiVG[1][0], 0, GiVG[1][2],
                    GiVG[2][0], GiVG[2][1], 0);

    /* s2[i] is the acceleration in previous time step (n) */
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        s2[i] -= GiVGoffDiag*_VSR[i]*.5;
        _VSR[i] += s2[i];
        dvs = _VSR[i]; 
        for(n=0;n<3;n++)
            dvs[n] = dvs[n]*exp((-.5)*GiVG[n][n]);        
        _VSR[i] = dvs;
        _SR[i] += _VSR[i];        
    }

    /* _SR[i] is the position at full step,  q(n+1) */
    /* _VSR[i] is the velocity at half-step,    p(n+1/2) */

    VVerlet_Get_s2(); /* eval s2 due to position change */

//    INFO("H in NPH=["<<_H<<"]");

    I += GiVGoffDiag*.5; Iinv = I.inv();
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        dvs = _VSR[i];
        for(n=0;n<3;n++)
            dvs[n] = dvs[n]*exp((-.5)*GiVG[n][n]);
        _VSR[i] = dvs + s2[i];
        _VSR[i] = Iinv*_VSR[i];
    }
    /* _VSR[i] is the velocity at current step, p(n+1) */
    
    _H+=_VH*.5;
    VVerlet_Get_s2(); /* eval s2 due to box change */     
    VVerlet_Get_h2(); /* Need to eval h2 due to Virial stress change */
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            if(conj_fixboxvec[i][j]==1)
                _VH[i][j]=h2[i][j]=0;
    _VH+=h2;

}

/*******************************************************
 *
 *  NPT ensemble
 *  Velocity Verlet integrator
 *
 ******************************************************/
void MDFrame::NPT_VVerlet_Implicit(int nstep)
{
    FATAL("NPT_VVerlet not implemented yet.  Come back soon!");
}

void MDFrame::NPT_VVerlet_Explicit_1(int nstep)
{
    FATAL("NPT_VVerlet not implemented yet.  Come back soon!");
}

void MDFrame::NPT_VVerlet_Explicit_2(int nstep)
{
    FATAL("NPT_VVerlet not implemented yet.  Come back soon!");
}


/* VVerlet Get s2[i] for the 1st run */
/* Jan 22 Keonwook Kang              */
void MDFrame::VVerlet_Get_s2()
{
    int i;
    double tmp;
    Matrix33 hinv;
    
    /* compute acceleration a(n+1),  s2[i]=0.5*a*(dt)^2 */
    call_potential();
    calcprop();

    hinv=_H.inv();

    /* fix me! need to change the unit of velocity to A/ps */
    tmp=0.5*(_TIMESTEP)*(_TIMESTEP)*1e-24  //(ps^2)
        /(_ATOMMASS[0]*1e-3/AVO)*EV/1e-20; //(A^2)

    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        if(species[i]==0)
            s2[i]=hinv*_F[i]*tmp;
        else
            s2[i]=hinv*_F[i]*tmp *_ATOMMASS[0]/_ATOMMASS[species[i]];
    }    
}

void MDFrame::VVerlet_Get_h2()
{
    int i, j;    
    double mass, tmp;

//    call_potential();
    calcprop();     /* To get GH */
    
    /* fix me! need to change the unit of velocity to A/ps */
    /* Feb 15 2007 Keonwook Kang
       Notice no "1e-20" in the denominatro!!! */
//    tmp=0.5*(_TIMESTEP)*(_TIMESTEP)*1e-24  //(ps^2)
//        /(_ATOMMASS[0]*1e-3/AVO)*EV; //
    
    /* ATOMMASS has unit of g/mol
     * mass     has unit of eV ps^2 / A^2 
     * tmp = 0.5 * dt^2 / m  has unit of (A^2/eV) */
//    mass = (_ATOMMASS[0]*1e-3/AVO/EV/1e20)*1e24;
    mass = _ATOMMASS[0]*MASSCONVERT;
    tmp  = 0.5 / mass * _TIMESTEP*_TIMESTEP;
    
//  h2=_GH*(-tmp)* _ATOMMASS[0]/_WALLMASS;
    h2=_GH*tmp*_ATOMMASS[0]/_WALLMASS;  /* unit of angstrom */
// change its sign after commenting out _GH*=-1; in calcprop() of md.cpp

    if(_BOXDAMP!=0)
        for(i=0;i<3;i++)
            for(j=0;j<3;j++)
                if(conj_fixboxvec[i][j]==0)
                    h2[i][j]-=_VH[i][j]*_BOXDAMP;

}


void MDFrame::VVerlet_Get_Tinst()
{
    int i;
    double mass, tmp;
    /* ATOMMASS has unit of g/mol
     * mass     has unit of eV ps^2 / A^2 
     * tmp = 0.5 * dt^2 / m  has unit of (A^2/eV) */
    mass = _ATOMMASS[0]*MASSCONVERT;
    tmp  = 0.5 * mass / (_TIMESTEP*_TIMESTEP);    
    
    _KATOM=0; _NPfree=0;
    for(i=0;i<_NP;i++)
    {
        if(fixed[i]) continue;
        _NPfree++;
        _VR[i]=_H*_VSR[i];
        if(species[i]==0)
            _KATOM+=_VR[i].norm2(); // in the unit of Angstrom^2
        else
            _KATOM+=_VR[i].norm2()*_ATOMMASS[species[i]]/_ATOMMASS[0];
    }
    _KATOM *= tmp; // 1/2 times (unit conversion from Angstrom^2 to eV )

    _T=_KATOM/(1.5*KB*_NPfree);    
}
    
void MDFrame::VVerlet_NHCINT()
{ /* Nose-Hoover Chain integrator from t=0 to t=_TIMESTEP/2 */
    int i;
    double expterm;

    VVerlet_Get_Tinst();
    
    i=NHChainLen-1;

    if(NHChainLen==1)
        zetaNHCa[i]=(_KATOM*2 - 3*_NPfree*KB*_TDES)/NHMass[i]
                    *(0.5*SQR(_TIMESTEP));
    else
        zetaNHCa[i]=(NHMass[i-1]*SQR(zetaNHCv[i-1]/_TIMESTEP)-KB*_TDES)/NHMass[i]
                *(0.5*SQR(_TIMESTEP));
    zetaNHCv[i]+=.5*zetaNHCa[i];
    
    for(i=NHChainLen-2;i>=0;i--)
    {
        expterm = exp(zetaNHCv[i+1]/(-8.0));
        if (i!=0)
        {
            zetaNHCa[i]=(NHMass[i-1]*SQR(zetaNHCv[i-1]/_TIMESTEP)-KB*_TDES)/NHMass[i]
                        *(0.5*SQR(_TIMESTEP));
        }
        else
        {
            zetaNHCa[i]=(_KATOM*2 - 3*_NPfree*KB*_TDES)/NHMass[i]
                        *(0.5*SQR(_TIMESTEP));
        }
        zetaNHCv[i] *= SQR(expterm);
        zetaNHCv[i] += 0.5*zetaNHCa[i]*expterm;
    }
    
    expterm = exp(-0.5*zetaNHCv[0]);
    for(i=0;i<_NP;i++)  _VSR[i]*=expterm;
    for(i=0;i<=NHChainLen-1;i++)
        zetaNHC[i]+=0.5*zetaNHCv[i];

    VVerlet_Get_Tinst();

    for(i=0;i<NHChainLen-1;i++)
    {
        expterm = exp(zetaNHCv[i+1]/(-8.0));
        if (i!=0)
        {
            zetaNHCa[i]=(NHMass[i-1]*SQR(zetaNHCv[i-1]/_TIMESTEP)-KB*_TDES)/NHMass[i]
                        *(0.5*SQR(_TIMESTEP));
        }
        else
        {
            zetaNHCa[i]=(_KATOM*2 - 3*_NPfree*KB*_TDES)/NHMass[i]
                        *(0.5*SQR(_TIMESTEP));
        }
        zetaNHCv[i] *= SQR(expterm);
        zetaNHCv[i] += 0.5*zetaNHCa[i]*expterm;
    }
    i=NHChainLen-1;
    if(NHChainLen==1)
        zetaNHCa[i]=(_KATOM*2 - 3*_NPfree*KB*_TDES)/NHMass[i]
                    *(0.5*SQR(_TIMESTEP));
    else
        zetaNHCa[i]=(NHMass[i-1]*SQR(zetaNHCv[i-1]/_TIMESTEP)-KB*_TDES)/NHMass[i]
                *(0.5*SQR(_TIMESTEP));
    zetaNHCv[i]+=.5*zetaNHCa[i];    
}
