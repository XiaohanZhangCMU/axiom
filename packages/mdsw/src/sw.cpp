/*
  sw.cpp
  by Wei Cai  caiwei@mit.edu
  Last Modified : Thu Apr 22 16:41:13 2010

  FUNCTION  :  MD simulation package of Si using Stillinger-Weber potential
*/

#include "sw.h"

void SWFrame::initvars()
{
    /* convresion between original PRB 31, 5262 (1985) and modified PRB 46,2250 (1992)
       Rc = acut
       mu = pss
       A1 = aa*bb
       A2 = aa
       lambda1 = rho
       lambda2 = 0
       Z = plam
       alpha = pgam/pss
       theta0 = acos(1./3.)
    */

#ifdef _SW_Si
#ifdef _SW_ORIG    
    /* original Si version PRB 31, 5262 (1985) */
       aa=15.27991323; bb=11.60319228; plam=45.51575;
       pgam=2.51412; acut=3.77118; pss=2.0951; rho=4.0;
#else   
    /* modified Si parameters PRB 46, 2250 (1992) */
       aa=16.31972277; bb=11.60319228; plam=48.61499998;
       pgam=2.51412;  acut=3.77118; pss=2.0951; rho=4.;
#endif    
#endif

#ifdef _SW_Ge       
   /* Ge version Ding and Anderson, PRB 34, 6987 (1986) */
       aa=13.60564361; bb=13.62639971; plam=59.83;
       pgam=2.6172; acut=3.9258; pss=2.181; rho=4.0; 
#endif

       if(acut==0) FATAL("acut = 0!  Need to have -D_SW_Si or -D_SW-Ge in Makefile");

       _RLIST=acut*1.1;
       /*_RLIST=3.8;*/
       /*_RLIST=4.2;*/   /* 1.1*acut */
       /*_RLIST=5.6;*/   /* large skin to avoid reconstruction in Free energy calculation */
       _SKIN=_RLIST-acut;
       rho1=rho+1; acutsq=acut*acut;
       
       strcpy(incnfile,"../si.cn");
       DUMP(HIG "SWFrame initvars" NOR);
       MDFrame::initvars();
}


#ifdef _TORSION_OR_BENDING
#include "../cookies/src/sw_torsion_bending.cpp"
#else

void SWFrame::stillinger_weber()
{
#ifndef NNM /* need to be taken care of later */
#define NNM 60
#endif
    /* Multi process function */
    int i,j,k,ipt,jpt,kpt;
    int neg;
    Vector3 sij,rij,f3i,f3j,f3k,f2;
    Vector3 rg[NNM],rt[NNM];
    double rr[NNM],ri[NNM],rdi[NNM],rdi2[NNM],xpon[NNM],xpon2[NNM];
    int list2[NNM];
    double r2ij,rrij;
    double ang,tm1,tm2,tm3,dhij,dhik,dhmu,tmij,tmik,tm2ij,tm2ik,
        tm3ik,tm3ij,eterm,au,pplm,ff2,eterm2;
    int n0, n1;

    DUMP("SW");

    _SKIN=_RLIST-acut;
    refreshneighborlist();

    /* _EPOT, _F, _EPOT_IND, _VIRIAL all local */
    _EPOT=0;tote2=tote3=0;
    for(i=0;i<_NP;i++)
    { _F[i].clear(); _EPOT_IND[i]=0; _EPOT_RMV[i]=0; _VIRIAL_IND[i].clear(); }
    _VIRIAL.clear();
    
#ifdef _TORSION_OR_BENDING
    //if (_NIMAGES>0)
    {
        if (_TORSIONSIM) _TORQUE = 0;
        if (_BENDSIM)    _BENDMOMENT = 0;
    }
#endif
    
    n0=0;
    n1=_NP;
    
    for(ipt=n0;ipt<n1;ipt++)
    {
        /* if fixed == -1, simply ignore this atom */
        if(fixed[ipt]==-1) continue;
        neg=0;
        for(j=0;j<nn[ipt];j++)
        {
            jpt=nindex[ipt][j];
            /* if fixed == -1, simply ignore this atom */
            if(fixed[jpt]==-1) continue;
            
            /*sij=_SR[jpt]-_SR[ipt];*/
            sij.subtract(_SR[jpt],_SR[ipt]);
            sij.subint();
            /*rij=_H*sij;*/
            _H.multiply(sij,rij);
            
            r2ij=rij.norm2();
            rrij=sqrt(r2ij);
            if(r2ij>acutsq) continue;
            rg[neg]=rij;
            /*rt[neg]=rij/rrij;*/
            rt[neg]=rij; rt[neg]/=rrij;
            rr[neg]=rrij;
            ri[neg]=1./rrij;
            rdi[neg]=1./(rrij-acut);
            rdi2[neg]=rdi[neg]*rdi[neg];
            if(fabs(pgam*rdi[neg])>30) xpon[neg]=0;/* avoid denormalization*/
            else xpon[neg]=exp(pgam*rdi[neg]);
            if(fabs(pss*rdi[neg])>30) xpon2[neg]=0;/* avoid denormalization*/
            else xpon2[neg]=exp(pss*rdi[neg]);
            list2[neg]=jpt;
            neg++;
        }

        /*
          for(j=0;j<neg;j++)
          INFO_Printf("xpon[%d]=%e rdi=%e pgam=%e\n",j,xpon[j],rdi[j],pgam);
        */
        
        /* second inner loop */
        for(j=0;j<neg;j++)
        {
            /* three body loop */
            jpt=list2[j];
            /* if fixed == -1, simply ignore this atom */
            if(fixed[jpt]==-1) continue;
            
            if(_SW3B_MUL!=0)
            {
                for(k=j+1;k<neg;k++)
                {
                    kpt=list2[k];
                    /* if fixed == -1, simply ignore this atom */
                    if(fixed[kpt]==-1) continue;
                    ang=rg[j]*rg[k]/(rr[j]*rr[k]);
                    tm1=ang+1./3.;
                    tm2=plam*pgam*tm1*tm1;
                    tm3=xpon[j]*xpon[k];
                    
                    dhij=-1.*tm2*rdi[j]*rdi[j]*tm3;
                    dhik=-1.*tm2*rdi[k]*rdi[k]*tm3;
                    dhmu=2.*plam*tm3*tm1;
                    
                    tmij=dhij+dhmu*(ri[k]-ang*ri[j]);
                    tmik=dhik+dhmu*(ri[j]-ang*ri[k]);
                    tm2ij=-dhij+dhmu*ang*ri[j];
                    tm2ik=-dhmu*ri[j];
                    tm3ik=-dhik+dhmu*ang*ri[k]; 
                    tm3ij=-dhmu*ri[k];
                    
                    /*f3i=rt[j]*tmij+rt[k]*tmik;*/
                    f3i.clear();f3i.addnv(tmij,rt[j]);f3i.addnv(tmik,rt[k]);
                    f3i *= _SW3B_MUL; /* multiplication factor for 3-body term */
                    _F[ipt]+=f3i;
                    
                    /*f3j=rt[j]*tm2ij+rt[k]*tm2ik;*/
                    f3j.clear();f3j.addnv(tm2ij,rt[j]);f3j.addnv(tm2ik,rt[k]);
                    f3j *= _SW3B_MUL; /* multiplication factor for 3-body term */
                    _F[jpt]+=f3j;
                    
                    /*f3k=rt[k]*tm3ik+rt[j]*tm3ij;*/
                    f3k.clear();f3k.addnv(tm3ik,rt[k]);f3k.addnv(tm3ij,rt[j]);
                    f3k *= _SW3B_MUL; /* multiplication factor for 3-body term */
                    _F[kpt]+=f3k;
                    
                    if(fixedatomenergypartition==0)
                    {
                        _VIRIAL.addnvv(1.,f3j,rg[j]);
			//if (SURFTEN==1 && (curstep%savepropfreq)==1) AddnvvtoPtPn(_SR[jpt],f3j,rg[j],1.0);

                        _VIRIAL.addnvv(1.,f3k,rg[k]);
			//if (SURFTEN==1 && (curstep%savepropfreq)==1) AddnvvtoPtPn(_SR[kpt],f3k,rg[k],1.0);
                        
                        _VIRIAL_IND[ipt].addnvv(0.5,f3j,rg[j]);
                        _VIRIAL_IND[jpt].addnvv(0.5,f3j,rg[j]);
                        _VIRIAL_IND[ipt].addnvv(0.5,f3k,rg[k]);
                        _VIRIAL_IND[kpt].addnvv(0.5,f3k,rg[k]);
                        
                    }
                    else
                    {
                        if(!(fixed[ipt]||fixed[jpt])) {
                            _VIRIAL.addnvv(1.,f3j,rg[j]);
			    //if (SURFTEN==1 && (curstep%savepropfreq)==1) AddnvvtoPtPn(_SR[jpt],f3j,rg[j],1.0);
			}
                        else if(!(fixed[ipt]&&fixed[jpt])) {
                            _VIRIAL.addnvv(0.5,f3j,rg[j]);
                            //if (SURFTEN==1 && (curstep%savepropfreq)==1) AddnvvtoPtPn(_SR[jpt],f3j,rg[j],0.5);
		  	}
                        if(!(fixed[ipt]||fixed[kpt])) {
                            _VIRIAL.addnvv(1.,f3k,rg[k]);
	                    //if (SURFTEN==1 && (curstep%savepropfreq)==1) AddnvvtoPtPn(_SR[kpt],f3k,rg[k],1.0);
			}
                        else if(!(fixed[ipt]&&fixed[kpt])) {
                            _VIRIAL.addnvv(0.5,f3k,rg[k]);
                            //if (SURFTEN==1 && (curstep%savepropfreq)==1) AddnvvtoPtPn(_SR[kpt],f3k,rg[k],0.5);
			}
                    }
                    /*_VIRIAL.addnvv(1.,f3k,rg[k]);*/
                    
                    eterm=tm2*tm3/pgam;
                    eterm *= _SW3B_MUL; /* multiplication factor for 3-body term */
                    
                    if(fixedatomenergypartition==0)
                    {
                        tote3+=eterm;
                    }
                    else
                    {
                        if(!(fixed[ipt]&&fixed[jpt]&&fixed[kpt]))
                            tote3+=eterm;
                    }
                    
                    _EPOT_IND[ipt]+=eterm/3;
                    _EPOT_IND[jpt]+=eterm/3;
                    _EPOT_IND[kpt]+=eterm/3;
                    
                    _EPOT_RMV[ipt]+=eterm;
                    _EPOT_RMV[jpt]+=eterm;
                    _EPOT_RMV[kpt]+=eterm;                
                }
            }
            /* two body terms */
            if(!Bond(ipt,jpt)) continue;
            au=aa*xpon2[j];
            pplm=pss*(bb*pow(ri[j],rho)-1.);
            ff2=au*(rho*bb*pow(ri[j],rho1)+pplm*rdi2[j]);
            /*f2=rt[j]*ff2;*/
            f2.clear();f2.addnv(ff2,rt[j]);
            f2 *= _SW2B_MUL; /* multiplication factor for 2-body term */                
            
            _F[ipt]-=f2;
            _F[jpt]+=f2;
            
            eterm2=aa*(bb*pow(ri[j],rho)-1.)*xpon2[j];
            eterm2 *= _SW2B_MUL; /* multiplication factor for 2-body term */
            
            if(fixedatomenergypartition==0)
            {
                tote2+=eterm2;
            }
            else
            {
                if(!(fixed[ipt]&&fixed[jpt]))
                    tote2+=eterm2;
            }
            _EPOT_IND[ipt]+=eterm2/2;
            _EPOT_IND[jpt]+=eterm2/2;
            
            _EPOT_RMV[ipt]+=eterm2;
            _EPOT_RMV[jpt]+=eterm2;            
            
            if(fixedatomenergypartition==0)
            {
                _VIRIAL.addnvv(1.,f2,rg[j]);
                //if (SURFTEN==1 && (curstep%savepropfreq)==1) AddnvvtoPtPn(_SR[jpt],f2,rg[j],1.0);
                _VIRIAL_IND[ipt].addnvv(0.5,f2,rg[j]);
                _VIRIAL_IND[jpt].addnvv(0.5,f2,rg[j]);
            }
            else
            {
                if(!(fixed[ipt]||fixed[jpt])) {
                    _VIRIAL.addnvv(1.,f2,rg[j]);
                 //   if (SURFTEN==1 && (curstep%savepropfreq)==1) AddnvvtoPtPn(_SR[jpt],f2,rg[j],1.0);
		}
                else if(!(fixed[ipt]&&fixed[jpt])) {
                    _VIRIAL.addnvv(0.5,f2,rg[j]);
                  //  if (SURFTEN==1 && (curstep%savepropfreq)==1) AddnvvtoPtPn(_SR[jpt],f2,rg[j],0.5);
		}
            }
        }  
    }
    _EPOT=tote2+tote3;
    /* zero out the forces on fixed atoms */
//    for(i=0;i<_NP;i++)
//    {
//        if(fixed[i])
//        {
//            _F0[i]=_F[i];
//            if(conj_fixdir==0)
//                _F[i].clear();
//            else if(conj_fixdir==1)
//                _F[i].x=0;
//            else if(conj_fixdir==2)
//                _F[i].y=0;
//            else if(conj_fixdir==3)
//                _F[i].z=0;
//        }
//    }
    
//#define CONSTRAIN_SUBLATTICE
#ifdef CONSTRAIN_SUBLATTICE
    /* apply constraints for ideal shear strength calculation*/
    Vector3 favg0, favg1;
    n0 = 0; n1 = 0; favg0.clear(); favg1.clear();
    
    for(i=0;i<_NP;i++)
    {
        if(species[i]==0)
        {
            favg0+=_F[i]; n0++;
        }
        else
        {
            favg1+=_F[i]; n1++;
        }
    }
    favg0/=n0; favg1/=n1;
    for(i=0;i<_NP;i++)
    {
        if(species[i]==0)
            _F[i]=favg0;
        else
            _F[i]=favg1;
    }
#endif
    
//#define TESTVIRIAL
#ifdef TESTVIRIAL
    SHtoR();
    for(i=0;i<_NP;i++)
    {
        _TOPOL[i] = dot(_R[i],_F[i]);
    }
#endif
    
    DUMP("SW complete");
}
#endif



void SWFrame::stillinger_weber_switch_3body(double lambda)
{
    _SW3B_MUL = lambda;
    stillinger_weber();
    dEdlambda = tote3;
    _SW3B_MUL = 1;
}


void SWFrame::stillinger_weber_energyonly()
{             
    /* Multi process function */
    int i,j,k,ipt,jpt,kpt;
    int neg;
    Vector3 sij,rij;
    Vector3 rg[NNM],rt[NNM];
    double rr[NNM],ri[NNM],rdi[NNM],rdi2[NNM],xpon[NNM],xpon2[NNM];
    int list2[NNM];
    double r2ij,rrij;
    double ang,tm1,tm2,tm3,eterm,eterm2;
    int n0, n1;
    
    DUMP("SW");
    
    _SKIN=_RLIST-acut;
    refreshneighborlist();
    
    /* _EPOT, _F, _EPOT_IND, _VIRIAL all local */
    _EPOT=0;tote2=tote3=0;
    for(i=0;i<_NP;i++)
    { _F[i].clear(); _EPOT_IND[i]=0; _EPOT_RMV[i]=0; }
    _VIRIAL.clear();

    n0=0;
    n1=_NP;
    
    for(ipt=n0;ipt<n1;ipt++)
    {
        /* if fixed == -1, simply ignore this atom */
        if(fixed[ipt]==-1) continue;
        neg=0;
        for(j=0;j<nn[ipt];j++)
        {
            jpt=nindex[ipt][j];
            /* if fixed == -1, simply ignore this atom */
            if(fixed[jpt]==-1) continue;
            
            /*sij=_SR[jpt]-_SR[ipt];*/
            sij.subtract(_SR[jpt],_SR[ipt]);
            sij.subint();
            /*rij=_H*sij;*/
            _H.multiply(sij,rij);
            
            r2ij=rij.norm2();
            rrij=sqrt(r2ij);
            if(r2ij>acutsq) continue;
            rg[neg]=rij;
            /*rt[neg]=rij/rrij;*/
            rt[neg]=rij; rt[neg]/=rrij;
            rr[neg]=rrij;
            ri[neg]=1./rrij;
            rdi[neg]=1./(rrij-acut);
            rdi2[neg]=rdi[neg]*rdi[neg];
            if(fabs(pgam*rdi[neg])>30) xpon[neg]=0;/* avoid denormalization*/
            else xpon[neg]=exp(pgam*rdi[neg]);
            if(fabs(pss*rdi[neg])>30) xpon2[neg]=0;/* avoid denormalization*/
            else xpon2[neg]=exp(pss*rdi[neg]);
            list2[neg]=jpt;
            neg++;
        }
        /* second inner loop */
        for(j=0;j<neg;j++)
        {
            /* three body loop */
            jpt=list2[j];
            /* if fixed == -1, simply ignore this atom */
            if(fixed[jpt]==-1) continue;
            
            if(_SW3B_MUL!=0)
            {
                for(k=j+1;k<neg;k++)
                {
                    kpt=list2[k];
                    /* if fixed == -1, simply ignore this atom */
                    if(fixed[kpt]==-1) continue;
                    ang=rg[j]*rg[k]/(rr[j]*rr[k]);
                    tm1=ang+1./3.;
                    tm2=plam*pgam*tm1*tm1;
                    tm3=xpon[j]*xpon[k];
                    
                    eterm=tm2*tm3/pgam;
                    eterm *= _SW3B_MUL; /* multiplication factor for 3-body term */
                    
                    if(fixedatomenergypartition==0)
                    {
                        tote3+=eterm;
                    }
                    else
                    {
                        if(!(fixed[ipt]&&fixed[jpt]&&fixed[kpt]))
                            tote3+=eterm;
                    }
                    _EPOT_IND[ipt]+=eterm/3;
                    _EPOT_IND[jpt]+=eterm/3;
                    _EPOT_IND[kpt]+=eterm/3;
                    
                    _EPOT_RMV[ipt]+=eterm;
                    _EPOT_RMV[jpt]+=eterm;
                    _EPOT_RMV[kpt]+=eterm;
                    
                }
            }
            /* two body terms */
            if(!Bond(ipt,jpt)) continue;
            
            eterm2=aa*(bb*pow(ri[j],rho)-1.)*xpon2[j];
            
            
            if(fixedatomenergypartition==0)
            {
                tote2+=eterm2;
            }
            else
            {
                if(!(fixed[ipt]&&fixed[jpt]))
                    tote2+=eterm2;
            }
            _EPOT_IND[ipt]+=eterm2/2;
            _EPOT_IND[jpt]+=eterm2/2;
            
            _EPOT_RMV[ipt]+=eterm2;
            _EPOT_RMV[jpt]+=eterm2;
        }
    }          
    _EPOT=tote2+tote3;
    
}

double SWFrame::stillinger_weber_energyonly(int iatom)
{             
    /* Multi process function */
    int j,k,ipt,jpt,kpt;
    int neg;
    Vector3 sij,rij;
    Vector3 rg[NNM],rt[NNM];
    double rr[NNM],ri[NNM],rdi[NNM],rdi2[NNM],xpon[NNM],xpon2[NNM];
    int list2[NNM];
    double r2ij,rrij;
    double ang,tm1,tm2,tm3,eterm,eterm2;
    double Eatom; //, maxd;

    if(iatom<_NP)
    {
        refreshneighborlist();
    }

    Eatom=0; tote2=tote3=0; _EPOT_IND[iatom]=0; _EPOT_RMV[iatom]=0;


    if(iatom>=_NP) 
    { 
        FATAL("iatom >= _NP (for test particle) no longer available (check version earlier than r499)");
    }
    
    for (int m = -1;m<nn[iatom];m++) //for(ipt=n0;ipt<n1;ipt++)
    {
	if (m == -1 )
	  {
	    ipt = iatom;	    
	  }
	else
	  {
	    ipt = nindex[iatom][m];
	  }

        neg=0;
        for(j=0;j<nn[ipt];j++)
        {
            jpt=nindex[ipt][j];
            /* if fixed == -1, simply ignore this atom */
            if(fixed[jpt]==-1) continue;
            
            /*sij=_SR[jpt]-_SR[ipt];*/
            sij.subtract(_SR[jpt],_SR[ipt]);
            sij.subint();
            /*rij=_H*sij;*/
            _H.multiply(sij,rij);
            
            r2ij=rij.norm2();
            rrij=sqrt(r2ij);
            if(r2ij>acutsq) continue;
            rg[neg]=rij;
            /*rt[neg]=rij/rrij;*/
            rt[neg]=rij; rt[neg]/=rrij;
            rr[neg]=rrij;
            ri[neg]=1./rrij;
            rdi[neg]=1./(rrij-acut);
            rdi2[neg]=rdi[neg]*rdi[neg];
            if(fabs(pgam*rdi[neg])>30) xpon[neg]=0;/* avoid denormalization*/
            else xpon[neg]=exp(pgam*rdi[neg]);
            if(fabs(pss*rdi[neg])>30) xpon2[neg]=0;/* avoid denormalization*/
            else xpon2[neg]=exp(pss*rdi[neg]);
            list2[neg]=jpt;
            neg++;
        }
        /* second inner loop */
        for(j=0;j<neg;j++)
        {
            /* three body loop */
            jpt=list2[j];
            /* if fixed == -1, simply ignore this atom */
            if(fixed[jpt]==-1) continue;
            
            for(k=j+1;k<neg;k++)
            {
                kpt=list2[k];
                /* if fixed == -1, simply ignore this atom */
                if(fixed[kpt]==-1) continue;
                ang=rg[j]*rg[k]/(rr[j]*rr[k]);
                tm1=ang+1./3.;
                tm2=plam*pgam*tm1*tm1;
                tm3=xpon[j]*xpon[k];
                
                eterm=tm2*tm3/pgam;
                eterm *= _SW3B_MUL; /* multiplication factor for 3-body term */
                
                if(fixedatomenergypartition==0)
                {
                    tote3+=eterm;
                }
                else
                {
                    if(!(fixed[ipt]&&fixed[jpt]&&fixed[kpt]))
                        tote3+=eterm;
                }
                if(ipt==iatom) _EPOT_IND[ipt]+=eterm/3;
                if(jpt==iatom) _EPOT_IND[jpt]+=eterm/3;
                if(kpt==iatom) _EPOT_IND[kpt]+=eterm/3;
                
                if(ipt==iatom) _EPOT_RMV[ipt]+=eterm;
                if(jpt==iatom) _EPOT_RMV[jpt]+=eterm;
                if(kpt==iatom) _EPOT_RMV[kpt]+=eterm;
                
            }
            /* two body terms */
            if(!Bond(ipt,jpt)) continue;
            
            eterm2=aa*(bb*pow(ri[j],rho)-1.)*xpon2[j];
            
            
            if(fixedatomenergypartition==0)
            {
                tote2+=eterm2;
            }
            else
            {
                if(!(fixed[ipt]&&fixed[jpt]))
                    tote2+=eterm2;
            }
            if(ipt==iatom) _EPOT_IND[ipt]+=eterm2/2;
            if(jpt==iatom) _EPOT_IND[jpt]+=eterm2/2;
            
            if(ipt==iatom) _EPOT_RMV[ipt]+=eterm2;
            if(jpt==iatom) _EPOT_RMV[jpt]+=eterm2;
            
        }
    }          
    Eatom=tote2+tote3;
    
    if (_EPOT_RMV != NULL)
      _EPOT_RMV[iatom] = Eatom;
    
    return Eatom;
}

void SWFrame::potential()
{
    stillinger_weber();
}

void SWFrame::potential_energyonly()
{
    stillinger_weber_energyonly();
}

double SWFrame::potential_energyonly(int iatom)
{
    return stillinger_weber_energyonly(iatom);
}

#ifdef _TEST

/* Main Program Begins */
class SWFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST

