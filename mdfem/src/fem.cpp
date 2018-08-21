/*
  fem.cpp
  by Wei Cai  caiwei@stanford.edu
  Last Modified : Fri Mar 11 21:07:00 2016

  FUNCTION  :  FEM model of hyperelastic material with brick elements
  NOTE :       All elements must be the same shape and size for now.
*/

#include "fem.h"
//#define NeoHooken


void FEMFrame::RtoRref()
{
     int i;
     SHtoR();
     for(i=0;i<_NP;i++)
     {    
         _SRref[i] = _SR[i];
         _Rref[i]  = _R[i];
     }
}

void FEMFrame::Alloc()
{
    MDFrame::Alloc();

    int size;
    size = _NP*allocmultiple;
    Realloc(_SRref,Vector3,size);
    Realloc(_Rref, Vector3,size);
    bindvar("SRref",_SRref,DOUBLE);
    bindvar("Rref",_Rref,DOUBLE);
    Realloc(map23,int,2);
    bindvar("map23",map23,INT);

    for(int i=0;i<_NP;i++) _VSR[i].clear();

#ifdef _USECUDA
    cuda_memory_alloc();
#endif
}   

void FEMFrame::Alloc_Elements()
{
    int size1, size2;
    size1 = _NELE*_NNODE_PER_ELEMENT;
    Realloc(elements,int,size1);
    bindvar("elements",elements,INT);
    size2 = _NELE*_NFACE_PER_ELEMENT;
    Realloc(colorids,int,size2);

#ifdef _USECUDA
    //We need colorids only for plotting
    cuda_memory_alloc_elements();
#endif

}

void FEMFrame::Alloc_Element_Coeff()
{
    int size;
    Realloc(gauss_weight,double,_NINT_PER_ELEMENT);
    bindvar("gauss_weight",gauss_weight,DOUBLE);

    size = _NINT_PER_ELEMENT*_NDIM*_NDIM*_NDIM*_NNODE_PER_ELEMENT  * _NELE;
    Realloc(dFdu,double,size);
    bindvar("dFdu",dFdu,DOUBLE);

#ifdef _USECUDA
    cuda_memory_alloc_element_coeff();
#endif
}

int FEMFrame::read_elements(const char* fname)
{
    int iE, j, iN, colorid;
    char *buffer; char *pp, *q;
    //int np, nframe;
    char extfname[MAXFILENAMELEN];
    //double coeff;
    int readCharCount, ind;

    strcpy(extfname,fname);
    LFile::SubHomeDir(extfname,extfname);
    LFile::LoadToString(extfname,buffer,0);

    pp=buffer;

    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
    sscanf(q,"%d %d %d",&_NELE, &_NNODE_PER_ELEMENT, &_NFACE_PER_ELEMENT);

    Alloc_Elements();

    for(iE=0;iE<_NELE;iE++)
    {
       q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
       
       for(j=0;j<_NNODE_PER_ELEMENT;j++)
       {
            sscanf(q,"%d%n",&iN,&readCharCount);
            q += readCharCount;
            ind=iE*_NNODE_PER_ELEMENT+j;
            elements[ind] = iN;
       }
       if ( _ReadColorID ) { 
          // Read in color identifier
          for(j=0; j<_NFACE_PER_ELEMENT; j++) {
              sscanf(q, "%d%n", &colorid, &readCharCount);
              q += readCharCount;
              ind = iE*_NFACE_PER_ELEMENT+j;
              INFO_Printf("color id = %d\n", colorid);
              colorids[ind] = colorid;
          }
          // Finish read in color id
       }
    }
    //INFO_Printf("\n");
    Free(buffer);
    return 0;
}

int FEMFrame::read_fem_coeff(const char* fname)
{
   INFO_Printf("_NELE = %d, EquationType = %d\n", _NELE, _EquationType);
   if (_EquationType == 1 || _EquationType == 0) { 
      assert(_NELE>0);
      for (int iele = 0; iele < _NELE; iele ++) {
         char str[10]; sprintf(str, "%d", iele);
         char elem_fname[150];
         strcpy(elem_fname,fname);
         strcat(elem_fname,"-");
         strcat(elem_fname,str);

         if (iele == 0)
            read_1stfem_to_Alloc(elem_fname);

         int dfdustart=iele*_NINT_PER_ELEMENT*_NDIM*_NDIM*_NDIM*_NNODE_PER_ELEMENT;
         read_element_wise_fem_coeff(elem_fname, dfdustart);
        }
   }
   else if ( _EquationType == 2) {
     INFO_Printf("_NELE = %d\n", _NELE);
     read_uniform_fem_coeff(fname);
   } 
   return 0;
}


int FEMFrame::read_uniform_fem_coeff(const char* fname)
{
    int iA, j, k, m, iN;
    char *buffer; char *pp, *q;
    char extfname[MAXFILENAMELEN];
    double coeff;
    int readCharCount, ind;

    strcpy(extfname,fname);
    LFile::SubHomeDir(extfname,extfname);
    LFile::LoadToString(extfname,buffer,0);

    pp=buffer;
    /* skip first line */    
    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
    /* fetch a line */    
    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;

    sscanf(q,"%d %d %d",&_NDIM,&_NNODE_PER_ELEMENT,&_NINT_PER_ELEMENT); 
    
    INFO_Printf("_NDIM = %d, _NINT_PER_ELEMENT = %d\n", _NDIM, _NINT_PER_ELEMENT);

    Alloc_Element_Coeff();
   
    /* skip third line */
    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;

    for(j=0;j<_NDIM;j++)
    {
       /* skip lines for reference configuration */
       q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
    }

    /* skip line */
    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;

    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
    for(iA=0;iA<_NINT_PER_ELEMENT;iA++)
    {
       sscanf(q,"%lf%n",gauss_weight+iA,&readCharCount);
       q += readCharCount;
       //INFO_Printf("gauss_weight[%d] = %f, readCharcnt = %d\n",iA,gauss_weight[iA], readCharCount);
    }

    /* skip line */
    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;

    for(j=0;j<_NDIM;j++)
    {
       /* skip lines for Gauss point positions */
       q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
    }

    for(iA=0;iA<_NINT_PER_ELEMENT;iA++)
    {
       /* skip a line for each Gauss point */
       q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
       //INFO_Printf("%s\n",q);

       for(j=0;j<_NDIM;j++)
       {
           for(k=0;k<_NDIM;k++)
           {
                /* skip lines for Gauss point positions */
                q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
                //INFO_Printf("%s\n",q);
                for(m=0;m<_NDIM;m++)
                {
                    for(iN=0;iN<_NNODE_PER_ELEMENT;iN++)
                    {
                        sscanf(q,"%lf%n",&coeff,&readCharCount);
                        q += readCharCount;
                        ind=(((iA*_NDIM+j)*_NDIM+k)*_NDIM+m)*_NNODE_PER_ELEMENT+iN;
                        dFdu[ind] = coeff;
                        //INFO_Printf("%d %f  \n",ind,dFdu[ind]);    
                    }
                }
                //INFO_Printf("\n");
           }
       }
    }
    Free(buffer);
    return 0;
}

int FEMFrame::read_1stfem_to_Alloc(const char* fname)
{
    char *buffer; char *pp, *q;
    char extfname[MAXFILENAMELEN];

    strcpy(extfname,fname);
    LFile::SubHomeDir(extfname,extfname);
    LFile::LoadToString(extfname,buffer,0);

    pp=buffer;
    /* skip first line */    
    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
    /* fetch a line */    
    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;

    sscanf(q,"%d %d %d",&_NDIM,&_NNODE_PER_ELEMENT,&_NINT_PER_ELEMENT);

    Alloc_Element_Coeff();

    Free(buffer);
    return 0;
}

int FEMFrame::read_element_wise_fem_coeff(const char* fname, int dfdustart)
{
    int iA, j, k, m, iN;
    char *buffer; char *pp, *q;
    char extfname[MAXFILENAMELEN];
    double coeff;
    int readCharCount, ind;

    strcpy(extfname,fname);
    LFile::SubHomeDir(extfname,extfname);
    LFile::LoadToString(extfname,buffer,0);

    pp=buffer;
    /* skip first line */    
    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
    /* fetch a line */    
    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;

    sscanf(q,"%d %d %d",&_NDIM,&_NNODE_PER_ELEMENT,&_NINT_PER_ELEMENT);

    //Alloc_Element_Coeff();
    
    /* skip third line */
    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;

    for(j=0;j<_NDIM;j++)
    {
       /* skip lines for reference configuration */
       q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
    }

    /* skip line */
    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;

    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
    for(iA=0;iA<_NINT_PER_ELEMENT;iA++)
    {
       sscanf(q,"%lf%n",gauss_weight+iA,&readCharCount);
       q += readCharCount;
       //INFO_Printf("gauss_weight[%d] = %f, readCharcnt = %d\n",iA,gauss_weight[iA], readCharCount);
    }

    /* skip line */
    q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;

    for(j=0;j<_NDIM;j++)
    {
       /* skip lines for Gauss point positions */
       q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
    }

    for(iA=0;iA<_NINT_PER_ELEMENT;iA++)
    {
       /* skip a line for each Gauss point */
       q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
       //INFO_Printf("%s\n",q);

       for(j=0;j<_NDIM;j++)
       {
           for(k=0;k<_NDIM;k++)
           {
                /* skip lines for Gauss point positions */
                q=pp; pp=strchr(pp, '\n'); if(pp) *(char *)pp++=0;
                //INFO_Printf("%s\n",q);
                for(m=0;m<_NDIM;m++)
                {
                    for(iN=0;iN<_NNODE_PER_ELEMENT;iN++)
                    {
                        sscanf(q,"%lf%n",&coeff,&readCharCount);
                        q += readCharCount;
                        ind=(((iA*_NDIM+j)*_NDIM+k)*_NDIM+m)*_NNODE_PER_ELEMENT+iN;
                        dFdu[dfdustart+ind] = coeff;
//			INFO_Printf("%d %f  \n",ind,dFdu[dfdustart+ind]);    
                    }
                }
                //INFO_Printf("\n");
           }
       }
    }

    Free(buffer);
    return 0;
}

void FEMFrame::plot()
{
    MDFrame::plot();

    /* draw covalent bonds between atoms */
    int iele, ipt, j, jpt, kpt, lpt, n, Nsub;
    double L;
    double x1,y1,z1,x2,y2,z2,dx,dy,dz,dr;
    //int r,g,b; 
    double alpha; //unsigned ce;
    Vector3 sri, srj, sij, ri, rj;
    Vector3 srk, srl, sik, sil, rk, rl;

    if(win==NULL) return;
    if(!(win->alive)) return;
    
    L=std::max(_H[0][0],_H[1][1]);
    L=std::max(L,_H[2][2])*.5;

    INFO_Printf("L = %f, _H00 = %f, _H11 = %f, _H22 = %f\n", L, _H[0][0], _H[1][1], _H[2][2]);
    SHtoR();
    win->Lock();
    //win->Clear();
    for(iele=0;iele<_NELE;iele++)
    {
      if (_NDIM == 2 )
      {
        for(j=0;j<_NNODE_PER_ELEMENT;j++)
        {
             jpt=elements[iele*_NNODE_PER_ELEMENT+j];
             ipt=elements[iele*_NNODE_PER_ELEMENT+((j+1)%_NNODE_PER_ELEMENT)];

             sri=_SR[ipt];
             if(plot_map_pbc==1) sri.subint();
             ri = _H*sri;
             sij=_SR[jpt]-_SR[ipt];
             sij.subint();
             srj=sri+sij;
             rj = _H*srj;
                
             if(plot_limits[0])
                    if((sri.x<plot_limits[1])||(sri.x>plot_limits[2])
                       ||(sri.y<plot_limits[3])||(sri.y>plot_limits[4])
                       ||(sri.z<plot_limits[5])||(sri.z>plot_limits[6])
                       ||(srj.x<plot_limits[1])||(srj.x>plot_limits[2])
                       ||(srj.y<plot_limits[3])||(srj.y>plot_limits[4])
                       ||(srj.z<plot_limits[5])||(srj.z>plot_limits[6]))
                     continue;

                    /* only draw O-O bonds */
                    /* if((species[i]!=0)||(species[j]!=0)) continue; */
                    //INFO_Printf("atom %d %d %d %d form bond\n",i, j, i%8,j%8); 
                    x1=ri.x/L;y1=ri.y/L;z1=ri.z/L;
                    x2=rj.x/L;y2=rj.y/L;z2=rj.z/L;
                    dx=x2-x1;dy=y2-y1;dz=z2-z1;dr=sqrt(dx*dx+dy*dy+dz*dz);
                    dx/=dr;dy/=dr;dz/=dr;
                    win->DrawLine(x1+dx*atomradius[species[ipt]]/L,
                                  y1+dy*atomradius[species[ipt]]/L,
                                  z1+dz*atomradius[species[ipt]]/L,
                                  x2-dx*atomradius[species[jpt]]/L,
                                  y2-dy*atomradius[species[jpt]]/L,
                                  z2-dz*atomradius[species[jpt]]/L,
                                  colors[MAXCOLORS+1],bondradius/L,1);
        }
        if(_NNODE_PER_ELEMENT == 4)
        {
            Nsub = 4;
            for(n=1;n<Nsub;n++)
            {

             alpha = (1.0*n)/Nsub;
             ipt=elements[iele*_NNODE_PER_ELEMENT+0];
             jpt=elements[iele*_NNODE_PER_ELEMENT+1];
             kpt=elements[iele*_NNODE_PER_ELEMENT+2];
             lpt=elements[iele*_NNODE_PER_ELEMENT+3];

	     //INFO_Printf("jpt = %d, ipt = %d\n", jpt, ipt);
             sri=_SR[ipt];
             if(plot_map_pbc==1) sri.subint();
             ri = _H*sri;
	 
             sij=_SR[jpt]-_SR[ipt]; sij.subint(); srj=sri+sij; rj = _H*srj;
             sik=_SR[kpt]-_SR[ipt]; sik.subint(); srk=sri+sik; rk = _H*srk;
             sil=_SR[lpt]-_SR[ipt]; sil.subint(); srl=sri+sil; rl = _H*srl;
                
             if(plot_limits[0])
                    if((sri.x<plot_limits[1])||(sri.x>plot_limits[2])
                       ||(sri.y<plot_limits[3])||(sri.y>plot_limits[4])
                       ||(sri.z<plot_limits[5])||(sri.z>plot_limits[6])
                       ||(srj.x<plot_limits[1])||(srj.x>plot_limits[2])
                       ||(srj.y<plot_limits[3])||(srj.y>plot_limits[4])
                       ||(srj.z<plot_limits[5])||(srj.z>plot_limits[6])
                       ||(srk.x<plot_limits[1])||(srk.x>plot_limits[2])
                       ||(srk.y<plot_limits[3])||(srk.y>plot_limits[4])
                       ||(srk.z<plot_limits[5])||(srk.z>plot_limits[6])
                       ||(srl.x<plot_limits[1])||(srl.x>plot_limits[2])
                       ||(srl.y<plot_limits[3])||(srl.y>plot_limits[4])
                       ||(srl.z<plot_limits[5])||(srl.z>plot_limits[6]))
                     continue;

                    /* only draw O-O bonds */
                    /* if((species[i]!=0)||(species[j]!=0)) continue; */
                    //INFO_Printf("atom %d %d %d %d form bond\n",i, j, i%8,j%8); 
                    x1=((1-alpha)*ri.x+alpha*rj.x)/L; 
                    y1=((1-alpha)*ri.y+alpha*rj.y)/L; 
                    z1=((1-alpha)*ri.z+alpha*rj.z)/L; 
                    x2=((1-alpha)*rl.x+alpha*rk.x)/L; 
                    y2=((1-alpha)*rl.y+alpha*rk.y)/L; 
                    z2=((1-alpha)*rl.z+alpha*rk.z)/L; 
                    dx=x2-x1;dy=y2-y1;dz=z2-z1;dr=sqrt(dx*dx+dy*dy+dz*dz);
                    dx/=dr;dy/=dr;dz/=dr;
                    win->DrawLine(x1+dx*atomradius[species[ipt]]/L,
                                  y1+dy*atomradius[species[ipt]]/L,
                                  z1+dz*atomradius[species[ipt]]/L,
                                  x2-dx*atomradius[species[jpt]]/L,
                                  y2-dy*atomradius[species[jpt]]/L,
                                  z2-dz*atomradius[species[jpt]]/L,
//                                colors[colorids[iele]],
                                  colors[MAXCOLORS+1],
				  (bondradius*0.1)/L,0);
                }
	}
      }

      if (_NDIM == 3 ) {
	int conns[12][2] = {{0,1},{1,3},{2,3},{0,2},{0,4},{2,6},{1,5},{3,7},{4,5},{5,7},{6,7},{4,6}};
        for(j=0;j<12;j++)
        {
             jpt=elements[iele*_NNODE_PER_ELEMENT+conns[j][0]];
             ipt=elements[iele*_NNODE_PER_ELEMENT+conns[j][1]];
             sri=_SR[ipt];
             if(plot_map_pbc==1) sri.subint();
             ri = _H*sri;

             sij=_SR[jpt]-_SR[ipt];
             sij.subint();
             srj=sri+sij;
             rj = _H*srj;
                
             if(plot_limits[0])
	       if((sri.x<plot_limits[1])||(sri.x>plot_limits[2])
		  ||(sri.y<plot_limits[3])||(sri.y>plot_limits[4])
		  ||(sri.z<plot_limits[5])||(sri.z>plot_limits[6])
		  ||(srj.x<plot_limits[1])||(srj.x>plot_limits[2])
		  ||(srj.y<plot_limits[3])||(srj.y>plot_limits[4])
		  ||(srj.z<plot_limits[5])||(srj.z>plot_limits[6])) {
		 continue;
	       }
	     /* only draw O-O bonds */
                    /* if((species[i]!=0)||(species[j]!=0)) continue; */
                    //INFO_Printf("atom %d %d %d %d form bond\n",i, j, i%8,j%8); 


                    x1=ri.x/L;y1=ri.y/L;z1=ri.z/L;
                    x2=rj.x/L;y2=rj.y/L;z2=rj.z/L;
                    dx=x2-x1;dy=y2-y1;dz=z2-z1;dr=sqrt(dx*dx+dy*dy+dz*dz);
                    dx/=dr;dy/=dr;dz/=dr;
                    win->DrawLine(x1+dx*atomradius[species[ipt]]/L,
                                  y1+dy*atomradius[species[ipt]]/L,
                                  z1+dz*atomradius[species[ipt]]/L,
                                  x2-dx*atomradius[species[jpt]]/L,
                                  y2-dy*atomradius[species[jpt]]/L,
                                  z2-dz*atomradius[species[jpt]]/L,
                                  colors[MAXCOLORS+1],bondradius/L,1);
        }
      }
    }
    
    win->Unlock();
    win->Refresh();
}


void FEMFrame::WriteStressCoord() {
  Vector3 sj, rj, center;
  std::ofstream stressCoord;
  stressCoord.open ("stressCoord.out", std::ofstream::out | std::ofstream::app);
  for (int i = 0;i<2;i++) {
    for(int iele=0;iele<_NELE;iele++) {
      center.clear();
      for(int j=0;j<_NNODE_PER_ELEMENT;j++) {
	int jpt=elements[iele*_NNODE_PER_ELEMENT+j];	  
	sj = _SR[jpt];
	_H.multiply(sj, rj);
	center += rj;
      }
      center/= _NNODE_PER_ELEMENT;	
      stressCoord<<center[i]<<" ";
    }    
    stressCoord<<std::endl;
  }
}

/*
  This read2cn() is specific for fiber+FEMsubstrate
  It first dumps the already existed configuration, and combine the dumped file with the current "incnfile"  
 */
int FEMFrame::read2cn() {
  assert(_NP>0); 
  int ival1, ival2; 
  double dvalx, dvaly, dvalz;
  n_existed_atoms = _NP;

  strcpy(finalcnfile,  "config1.cn");
  writefinalcnfile(0,false);

  FILE *config1fp  = fopen(finalcnfile,"r");
  FILE *config2fp  = fopen(incnfile,"r");   
  FILE *config12fp = fopen("config12.cn","w");

  assert (config1fp!=NULL); 
  assert (config2fp!=NULL);
  assert (config12fp!=NULL);

  fscanf(config1fp, "%d", &ival1);
  fscanf(config2fp, "%d", &ival2);  

  nfemNodes = ival2;
  fprintf(config12fp, "%d\n" ,ival1+ival2);

  for (int nd = 0;nd<ival1;nd++) {
    fscanf(config1fp, "%lf %lf %lf" , &dvalx, &dvaly, &dvalz);
    fprintf(config12fp, "%23.18E %23.18E %23.18E\n" ,dvalx, dvaly, dvalz);
  }
  for (int nd = 0;nd<ival2;nd++) {
    fscanf(config2fp, "%lf %lf %lf" , &dvalx, &dvaly, &dvalz);
    fprintf(config12fp, "%23.18E %23.18E %23.18E\n" ,dvalx, dvaly, dvalz);
  }
  //box information is the same as the existed one
  for (int nd = 0;nd<3;nd++) {
    fscanf(config1fp, "%lf %lf %lf" , &dvalx, &dvaly, &dvalz);
    fprintf(config12fp, "%23.18E %23.18E %23.18E\n" ,dvalx, dvaly, dvalz);
  }

  fclose(config1fp);
  fclose(config2fp);
  fclose(config12fp);

  //perform a regular readcn on the combined cn file
  strcpy(incnfile,  "config12.cn");
  return (FEMFrame::readcn());
}


int FEMFrame::read_bdy_nodes(const char* fname)
{
  FILE *bdyfile = fopen(fname,"r"); assert(bdyfile != NULL); 
  
  fscanf(bdyfile,"%d %d", &n_bdy_nodes,&n_bdy_tags); 
  assert(n_bdy_nodes > 0 && n_bdy_tags > 0); 

  Realloc(bNds,int,n_bdy_nodes);
  Realloc(bTags,int,n_bdy_nodes);
  Realloc(bXs,double,n_bdy_nodes);
  Realloc(bYs,double,n_bdy_nodes);
  Realloc(bZs,double,n_bdy_nodes);

  bindvar("n_bdy_nodes",&n_bdy_nodes,INT);
  bindvar("n_bdy_tags",&n_bdy_tags,INT);
  bindvar("bNds",bNds,INT);
  bindvar("bTags",bTags,INT);
  bindvar("bXs", bXs,DOUBLE);
  bindvar("bYs", bYs,DOUBLE);
  bindvar("bZs", bZs,DOUBLE);

  for(int ibn = 0; ibn < n_bdy_nodes; ibn++) {
    fscanf(bdyfile,"%d %lf %lf %lf %d", &bNds[ibn], &bXs[ibn], &bYs[ibn], &bZs[ibn], &bTags[ibn]);
    //bNds[ibn] += n_existed_atoms;
  }

  printf("read_bdy_nodes finished, n_bdy_nodes = %d\n",n_bdy_nodes);

  return 0;
}


void FEMFrame::shift_fem_node_id(int np0) {
  for(int iele=0;iele<_NELE;iele++) {
    for(int j=0;j<_NNODE_PER_ELEMENT;j++) {
      elements[iele*_NNODE_PER_ELEMENT+j] += np0;
    } }

  for(int ibn = 0; ibn < n_bdy_nodes; ibn++) {
    bNds[ibn] += np0;
  }
}

void FEMFrame::potential() {
    INFO_Printf("MDFrame:potential() is called\n");
   switch(_EquationType) { 	
      case 0:
#ifdef _USECUDA
         cuda_beam_fem_energy_force();
#else
         beam_fem_energy_force();
#endif
         break;
      case 1:
         islam_fem_energy_force();
         break;
      case 2:
#ifdef _USECUDA
         cuda_snap_fem_energy_force();
#else
         snap_fem_energy_force();
#endif
         break;
      default:
      FATAL("Need to set EquationType in script\n");
   }
}

void FEMFrame::snap_fem_energy_force() {
    /* no need of neighbor list */
    int i,iele,j,jpt,iA, in, ip, iq, ind, p, q, r;
    Vector3 dsj, drj, elem_center;
    Matrix33 Fe, Fdef, B, C, Cinv, FCinvtran, dEdF, hinv;
    Matrix33 eigF, invEigF;
    double Eele, Eint, Jet, trace, detB, J2;
    double lambda, mu, temp;
    Matrix33 E, E2, I, pk2, temp1, temp2;    
    I[0][0]=I[1][1]=I[2][2]=1;

    //    int map23[2] = {0, 2};

    DUMP("FEM");

    /* _EPOT, _F, _EPOT_IND, _VIRIAL all local */
    _EPOT=0; //put this back
    for(i=0;i<_NP;i++) {
      _F[i].clear(); _EPOT_IND[i]=0; _EPOT_RMV[i]=0; _VIRIAL_IND[i].clear();
    }

    _VIRIAL.clear();
    hinv = _H.inv();

    for(iele=0;iele<_NELE;iele++) {
      /* energy of this element */
      Eele = 0;

      /* center of the element */
      elem_center.clear();elem_center[0] = elem_center[1]=elem_center[2]= 0;
      for(j=0;j<_NNODE_PER_ELEMENT;j++) {
	jpt=elements[iele*_NNODE_PER_ELEMENT+j];	
        for (i = 0;i<_NDIM;i++)  {
	  elem_center[i] += 1.0/_NNODE_PER_ELEMENT *_Rref[jpt][i]; /* put this into a separate function RrefHtoSref */
	}
      }

      for(iA=0;iA<_NINT_PER_ELEMENT;iA++) {
	/* energy of this Gauss integration point */
	Eint = 0;
	/* deformation gradient */
	Fdef.clear(); Fdef[0][0] = Fdef[1][1] = Fdef[2][2] = 1.0;
	for(j=0;j<_NNODE_PER_ELEMENT;j++) {
	  jpt=elements[iele*_NNODE_PER_ELEMENT+j];
	  
	  _SRref[jpt ] = hinv*_Rref[jpt];  /* put this into a separate function RrefHtoSref */
	  dsj = _SR[jpt] - _SRref[jpt];
	  dsj.subint();
	  _H.multiply(dsj, drj);
	  
	  /* Add contribution from drj to F using dFdu */
	  for(ip=0;ip<_NDIM;ip++) {
	    for(iq=0;iq<_NDIM;iq++) {
	      for(in=0;in<_NDIM;in++) {
		ind=(((iA*_NDIM+ip)*_NDIM+iq)*_NDIM+in)*_NNODE_PER_ELEMENT+j;
		//Fdef[ip][iq] += dFdu[ind]*drj[in]; 
		if (_NDIM == 2)
		  Fdef[ip][iq] += dFdu[ind]*drj[map23[in]]; 
		else //_NDIM == 3
		  Fdef[ip][iq] += dFdu[ind]*drj[in];
	      } } }
	}

	eigF = getEigenF(elem_center, Fdef);
	invEigF = eigF.inv();
	Fe = Fdef*invEigF;     
	E = Fe.tran()*Fe-I;
	B = Fe*Fe.tran();
	C = Fe.tran()*Fe;
	Cinv = C.inv();
//	FCinvtran = Fdef*Cinv.tran();
	Jet = Fdef.det();  J2 = Jet*Jet;
	detB = B.det(); 
	Eint = 0;

	/* Add contribution from F to Eint */ 
	if (_NDIM == 2) {
	  trace = B[0][0] + B[1][1];
	  Eint = 0.5*(trace + 1.0/detB - 3) ; /* multiply MU and V0 later */
	  dEdF = FCinvtran*(-1.0/J2) + Fdef; /* multiply MU later */
	  /* Add contribution from drj to F using dFdu */
	  for(j=0;j<_NNODE_PER_ELEMENT;j++) {
	    jpt=elements[iele*_NNODE_PER_ELEMENT+j];
	    for(ip=0;ip<_NDIM;ip++) {
	      for(iq=0;iq<_NDIM;iq++) {
		for(in=0;in<_NDIM;in++) {
		  ind=(((iA*_NDIM+ip)*_NDIM+iq)*_NDIM+in)*_NNODE_PER_ELEMENT+j;
		  if (fixed[jpt] == 0)				  
//		    _F[jpt][map23[in] ] -=  dEdF[ip][iq]*dFdu[ind] * __V0; 
		    _F[jpt][map23[in] ] -= gauss_weight[iA] * dEdF[ip][iq]*dFdu[ind] * __V0; 				
		} } } }
	}

	/* Add contribution from F to Eint */ 
	if (_NDIM == 3) {
	  mu = 1; lambda = 1.95;
		
#if defined NeoHooken
    double C10 = 1; double D = 1e-1;
    double jet23 = pow(Jet, -2.0/3.0);
    double Ibar = B.trace() * jet23;
    Eint = C10*(Ibar-3) + 1.0/D *(Jet-1)*(Jet-1);
    dEdF = (Fdef*(invEigF*(invEigF.tran())) * 2.0*jet23 - ((Fdef.inv()).tran())*2.0/3.0 * Ibar)*C10 + ((Fdef.inv()).tran())*2.0/D *(Jet-1)*Jet;
#else
    E2 = E*E;
 	  dEdF.clear();
	  for (i = 0;i<_NDIM;i++)
	    for (j = 0;j<_NDIM;j++)
	      for (p = 0;p<_NDIM;p++)
		for (r = 0;r<_NDIM;r++) {
		  temp = 0;
		  for (q = 0;q<_NDIM;q++)
		    temp += 2*mu*invEigF[j][p]*Fdef[i][r]*invEigF[r][q]*E[p][q] + 
		      2*mu*invEigF[r][p]*Fdef[i][r]*invEigF[j][q]*E[p][q];
		  dEdF[i][j] += 0.5*(2*lambda*E.trace()*invEigF[j][p]*Fdef[i][r] * invEigF[r][p] + temp);
		}
	  Eint = 0.5*lambda*(E.trace())*(E.trace()) + mu*(E2.trace());
#endif
	  /* Add contribution from drj to F using dFdu */

	  for(j=0;j<_NNODE_PER_ELEMENT;j++) {
	    jpt=elements[iele*_NNODE_PER_ELEMENT+j];
	    for(ip=0;ip<_NDIM;ip++) {
	      for(iq=0;iq<_NDIM;iq++) {
		for(in=0;in<_NDIM;in++) {
		  ind=(((iA*_NDIM+ip)*_NDIM+iq)*_NDIM+in)*_NNODE_PER_ELEMENT+j;
		  if (fixed[jpt] == 0)				  
//		    _F[jpt][in ] -= dEdF[ip][iq]*dFdu[ind] * __V0; 				
		    _F[jpt][in ] -= gauss_weight[iA] * dEdF[ip][iq]*dFdu[ind] * __V0; 				
		} } } }
	}
            /* Add contribution from Eint to Eele */
	Eele += gauss_weight[iA]*Eint *__V0;
      }
      _EPOT += Eele;
    }
}  

void FEMFrame::create_inverse_connectivities_matrix() { 
  assert(_NP > 0);
  int i, j, k, jpt;
  Realloc(inv_elements,int,_MAX_NELEM_SHARE_NODE*_NP*sizeof(int));
  memset(inv_elements,-1,  _MAX_NELEM_SHARE_NODE*_NP*sizeof(int));

  for (i = 0; i< _NELE; i++) {
    for (j = 0; j< _NNODE_PER_ELEMENT; j++) {
      jpt = elements[_NNODE_PER_ELEMENT*i + j];
      k = 0;
      for (k = 0; k < _MAX_NELEM_SHARE_NODE; k++) {
        if (inv_elements[jpt*_MAX_NELEM_SHARE_NODE+k] == -1) { 
          inv_elements[jpt*_MAX_NELEM_SHARE_NODE+k] = i*_NNODE_PER_ELEMENT +j;
          break;
        }
      }
      if (k == _MAX_NELEM_SHARE_NODE) FATAL("Too many elements share a node!");
    }
  } 
}

void FEMFrame::islam_fem_energy_force() {
    /* no need of neighbor list */
    int i,iele,j,jpt,iA, in, ip, iq, ind;
    Vector3 dsj, drj;
    Vector3 vec1, vec2, normal;
    Matrix33 Fdef, B, C, Cinv, FCinvtran, dEdF, hinv;
    double Eele, Eint, Jet, trace, detB, J2;
    double lambda, mu;
    /*
      Xiaohan: print out stress in the membrane
    */    
    Matrix33 E,Eavg, E2 , I, pk2, temp1, temp2;    
    I[0][0]=I[1][1]=I[2][2]=1;
    std::ofstream w_ofs("energy.out", std::ofstream::out | std::ofstream::app);
    std::ofstream e11_ofs ("strain11.out", std::ofstream::out | std::ofstream::app);
    std::ofstream e12_ofs ("strain12.out", std::ofstream::out | std::ofstream::app);
    std::ofstream e22_ofs ("strain22.out", std::ofstream::out | std::ofstream::app);
    std::ofstream neb_dump("gp_rawdata.out", std::ofstream::out | std::ofstream::app);
    DUMP("FEM");

    /* _EPOT, _F, _EPOT_IND, _VIRIAL all local */
    _EPOT=0; //put this back
    for(i=0;i<_NP;i++) { _F[i].clear(); _EPOT_IND[i]=0; _EPOT_RMV[i]=0; _VIRIAL_IND[i].clear(); }

    _VIRIAL.clear();
    hinv = _H.inv();
    WriteStressCoord();

    for(iele=0;iele<_NELE;iele++) {
      /* energy of this element */
      Eele = 0;
      Eavg.clear();
      int dfdustart = iele* _NINT_PER_ELEMENT*_NDIM*_NDIM*_NDIM*_NNODE_PER_ELEMENT ;
      if (curstep%printfreq==0) { neb_dump<<iele<<" "; }
      for(iA=0;iA<_NINT_PER_ELEMENT;iA++) {
	/* energy of this Gauss integration point */
	Eint = 0;
	/* deformation gradient */
	Fdef.clear(); Fdef[0][0] = Fdef[1][1] = Fdef[2][2] = 1.0;
	for(j=0;j<_NNODE_PER_ELEMENT;j++) {
	  jpt=elements[iele*_NNODE_PER_ELEMENT+j];
	  
	  _SRref[jpt ] = hinv*_Rref[jpt];  /* put this into a separate function RrefHtoSref */
	  dsj = _SR[jpt] - _SRref[jpt];
	  dsj.subint();
	  _H.multiply(dsj, drj);
	  /* Add contribution from drj to F using dFdu */
	  for(ip=0;ip<_NDIM;ip++) {
	    for(iq=0;iq<_NDIM;iq++) {
	      for(in=0;in<_NDIM;in++) {
		ind=(((iA*_NDIM+ip)*_NDIM+iq)*_NDIM+in)*_NNODE_PER_ELEMENT+j;
		if (_NDIM == 2)
		  Fdef[ip][iq] += dFdu[dfdustart+ind]*drj[map23[in]]; 
		else //_NDIM == 3
		  Fdef[ip][iq] += dFdu[dfdustart+ind]*drj[in];
	      } } }
	}

//	INFO_Printf("Fdef = %g, %g, %g, %g, %g, %g, %g, %g, %g\n", Fdef[0][0], Fdef[1][1],Fdef[2][2],Fdef[0][1],Fdef[0][2],Fdef[1][0],Fdef[1][2],Fdef[2][0],Fdef[2][1]);

	E = Fdef.tran()*Fdef-I;
	B = Fdef*Fdef.tran();
	C = Fdef.tran()*Fdef;
	Cinv = C.inv();
	FCinvtran = Fdef*Cinv.tran();
	Jet = Fdef.det();  J2 = Jet*Jet;
	detB = B.det(); 
	Eint = 0;
	Eavg += E;
	/* Add contribution from F to Eint */ 
	if (_NDIM == 2) {
	  double MU = 1;
	  trace = B[0][0] + B[1][1];
	  Eint = 0.5*(trace + 1.0/detB - 3); /* multiply MU and V0 later */
	  dEdF = FCinvtran*(-1.0/J2) + Fdef; /* multiply MU later */
	  Eint *= MU;
	  dEdF *= MU;

	  /* Add contribution from drj to F using dFdu */
	  for(j=0;j<_NNODE_PER_ELEMENT;j++) {
	    jpt=elements[iele*_NNODE_PER_ELEMENT+j];
	    for(ip=0;ip<_NDIM;ip++) {
	      for(iq=0;iq<_NDIM;iq++) {
		for(in=0;in<_NDIM;in++) {
		  ind=(((iA*_NDIM+ip)*_NDIM+iq)*_NDIM+in)*_NNODE_PER_ELEMENT+j;
		  if (fixed[jpt] == 0)				  
//		    _F[jpt][map23[in] ] -= dEdF[ip][iq]*dFdu[dfdustart+ind] *__V0; 
		    _F[jpt][map23[in]] -= gauss_weight[iA] * dEdF[ip][iq]*dFdu[dfdustart+ind] * __V0; 				
		} } } }
	}

	/* Add contribution from F to Eint */ 
	if (_NDIM == 3) {
	  mu = 0.1; lambda = 1e3;
		
#if defined NeoHooken
	  Eint = 1.0/8.0 * (0.5*lambda*log(Jet)*log(Jet) - mu*log(Jet) + 0.5*mu*(C.trace()-3));
	  temp1 = I;
	  temp1*= lambda * log(Jet);
	  temp2 = I;
	  temp2 = temp2-Cinv;
	  temp2*= mu;
	  pk2  = temp1 + temp2;
	  dEdF = Fdef * pk2;
#else
	  E2 = E*E;
	  Eint = 0.5*lambda*(E.trace())*(E.trace()) + mu*(E2.trace());
	  temp1 = I; 
	  temp1 *= E.trace()*lambda;
	  temp2 = E;
	  temp2 *= 2*mu;
	  pk2  = temp1 + temp2;
	  dEdF = Fdef * pk2;
#endif
	  /* Add contribution from drj to F using dFdu */

	  for(j=0;j<_NNODE_PER_ELEMENT;j++) {
	    jpt=elements[iele*_NNODE_PER_ELEMENT+j];
	    for(ip=0;ip<_NDIM;ip++) {
	      for(iq=0;iq<_NDIM;iq++) {
		for(in=0;in<_NDIM;in++) {
		  ind=(((iA*_NDIM+ip)*_NDIM+iq)*_NDIM+in)*_NNODE_PER_ELEMENT+j;
		  if (fixed[jpt] == 0)				  
//		    _F[jpt][in ] -= dEdF[ip][iq]*dFdu[dfdustart+ind] *__V0; 				
		    _F[jpt][in ] -= gauss_weight[iA] * dEdF[ip][iq]*dFdu[dfdustart+ind] * __V0; 				
		} } } }
	}
            /* Add contribution from Eint to Eele */
        
	Eele += gauss_weight[iA]*Eint *__V0;

        if (curstep%printfreq==0) { neb_dump<<Eint<<" "; }
      }
      _EPOT += Eele;
      if (curstep%printfreq==0) {
            w_ofs<< Eele <<" ";
            e11_ofs<< Eavg[0][0]/_NINT_PER_ELEMENT<<" ";
            e12_ofs<< Eavg[0][1]/_NINT_PER_ELEMENT<<" ";
            e22_ofs<< Eavg[1][1]/_NINT_PER_ELEMENT<<" ";
            neb_dump<<std::endl;
      }
   }
   if (curstep%printfreq==0) {
       w_ofs << std::endl;
       e11_ofs << std::endl;
       e12_ofs << std::endl;
       e22_ofs << std::endl;
   }
}  



void FEMFrame::beam_fem_energy_force() {
    /* no need of neighbor list */
    int i,iele,j,jpt,iA, in, ip, iq, ind;
    Vector3 dsj, drj;
    Vector3 vec1, vec2, normal;
    Matrix33 Fdef, B, C, Cinv, FCinvtran, dEdF, hinv;
    double Eele, Eint, Jet, trace, detB, J2;
    double lambda, mu;
    /*
      Xiaohan: print out stress in the membrane
    */    
    Matrix33 E, E2, I, pk2, temp1, temp2;    
    I[0][0]=I[1][1]=I[2][2]=1;

    //    int map23[2] = {0, 2};

    DUMP("BEAM");

    /* _EPOT, _F, _EPOT_IND, _VIRIAL all local */
    _EPOT=0; //put this back
    for(i=0;i<_NP;i++) {
      /* xiaohan: this removes the interaction between substrate with fe nodes which have species = 2*/
      // if (species[i] == 0) 
      // 	continue;
      _F[i].clear(); _EPOT_IND[i]=0; _EPOT_RMV[i]=0; _VIRIAL_IND[i].clear();
    }

    _VIRIAL.clear();
    hinv = _H.inv();

    for(iele=0;iele<_NELE;iele++) {
      /* energy of this element */
      Eele = 0;
      
      for(iA=0;iA<_NINT_PER_ELEMENT;iA++) {
	/* energy of this Gauss integration point */
	Eint = 0;
	/* deformation gradient */
	Fdef.clear(); Fdef[0][0] = Fdef[1][1] = Fdef[2][2] = 1.0;
	for(j=0;j<_NNODE_PER_ELEMENT;j++) {
	  jpt=elements[iele*_NNODE_PER_ELEMENT+j];
	  _SRref[jpt ] = hinv*_Rref[jpt];  /* put this into a separate function RrefHtoSref */
	  dsj = _SR[jpt] - _SRref[jpt];
	  dsj.subint();
	  _H.multiply(dsj, drj);
	  /* Add contribution from drj to F using dFdu */
	  for(ip=0;ip<_NDIM;ip++) {
	    for(iq=0;iq<_NDIM;iq++) {
	      for(in=0;in<_NDIM;in++) {
		ind=(((iA*_NDIM+ip)*_NDIM+iq)*_NDIM+in)*_NNODE_PER_ELEMENT+j;
		if (_NDIM == 2)
		  Fdef[ip][iq] += dFdu[ind]*drj[map23[in]]; 
		else //_NDIM == 3
		  Fdef[ip][iq] += dFdu[ind]*drj[in];
	      } } }
	}

	E = Fdef.tran()*Fdef-I;
	B = Fdef*Fdef.tran();
	C = Fdef.tran()*Fdef;
	Cinv = C.inv();
	FCinvtran = Fdef*Cinv.tran();
	Jet = Fdef.det();  J2 = Jet*Jet;
	detB = B.det(); 
	Eint = 0;

	/* Add contribution from F to Eint */ 
	if (_NDIM == 2) {
	  double MU = 1;
	  trace = B[0][0] + B[1][1];
	  Eint = 0.5*(trace + 1.0/detB - 3); /* multiply MU and V0 later */
	  dEdF = FCinvtran*(-1.0/J2) + Fdef; /* multiply MU later */
	  Eint *= MU;
	  dEdF *= MU;

	  /* Add contribution from drj to F using dFdu */
	  for(j=0;j<_NNODE_PER_ELEMENT;j++) {
	    jpt=elements[iele*_NNODE_PER_ELEMENT+j];
	    for(ip=0;ip<_NDIM;ip++) {
	      for(iq=0;iq<_NDIM;iq++) {
		for(in=0;in<_NDIM;in++) {
		  ind=(((iA*_NDIM+ip)*_NDIM+iq)*_NDIM+in)*_NNODE_PER_ELEMENT+j;
		  if (fixed[jpt] == 0) {
//		    _F[jpt][map23[in] ] -= dEdF[ip][iq]*dFdu[ind] *__V0; 
		    _F[jpt][map23[in] ] -= gauss_weight[iA] * dEdF[ip][iq]*dFdu[ind] * __V0; 				
        //printf("_F[%d][%d] = %g\n", jpt, map23[in], _F[jpt][map23[in] ]);
      }
		} } } }
	}

	/* Add contribution from F to Eint */ 
	if (_NDIM == 3) {
	  mu = 0.1; lambda = 1e3;
		
#if defined NeoHooken
	  Eint = 1.0/8.0 * (0.5*lambda*log(Jet)*log(Jet) - mu*log(Jet) + 0.5*mu*(C.trace()-3));
	  temp1 = I;
	  temp1*= lambda * log(Jet);
	  temp2 = I;
	  temp2 = temp2-Cinv;
	  temp2*= mu;
	  pk2  = temp1 + temp2;
	  dEdF = Fdef * pk2;
#else
	  E2 = E*E;
	  Eint = 0.5*lambda*(E.trace())*(E.trace()) + mu*(E2.trace());
	  temp1 = I; 
	  temp1 *= E.trace()*lambda;
	  temp2 = E;
	  temp2 *= 2*mu;
	  pk2  = temp1 + temp2;
	  dEdF = Fdef * pk2;
#endif
	  /* Add contribution from drj to F using dFdu */

	  for(j=0;j<_NNODE_PER_ELEMENT;j++) {
	    jpt=elements[iele*_NNODE_PER_ELEMENT+j];
	    for(ip=0;ip<_NDIM;ip++) {
	      for(iq=0;iq<_NDIM;iq++) {
		for(in=0;in<_NDIM;in++) {
		  ind=(((iA*_NDIM+ip)*_NDIM+iq)*_NDIM+in)*_NNODE_PER_ELEMENT+j;
		  if (fixed[jpt] == 0)				  
//		    _F[jpt][in ] -= dEdF[ip][iq]*dFdu[ind] * __V0; 				
		    _F[jpt][in ] -= gauss_weight[iA] * dEdF[ip][iq]*dFdu[ind] * __V0; 				
		} } } }
	}
            /* Add contribution from Eint to Eele */
	Eele += gauss_weight[iA]*Eint *__V0;
      }
      _EPOT += Eele;
    }
}  

Matrix33 FEMFrame::getEigenF(Vector3 p, Matrix33 Fdef) {
  Matrix33 I; I.eye();    
   if (p[2] <= y_eigen_zbound_max && p[2] >= y_eigen_zbound_min) { 
    I[1][1] = y_eigen_strain;
   }
  if (p[2] <= x_eigen_zbound_max && p[2] >= x_eigen_zbound_min) {
    I[0][0] = x_eigen_strain;
  }
  
  return I;
}

int FEMFrame::test_function(int n) { 
    this->nfemNodes = n;
    INFO_Printf("in test_function. nfemNodes = %d\n", this->nfemNodes);
    return 100;
}

#ifdef _TEST

/* Main Program Begins */
class FEMFrame sim;

/* The main program is defined here */
#include "main.cpp"

#endif//_TEST

