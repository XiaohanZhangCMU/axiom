//#include "general.hpp"
//#include "scparser.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXLIN 5 //Maximum iterations within a line search
#define MXFCON 2 //Maximum iterations before F can be decreased

/* Getting the code out of "infinite" loop */
#define MAXREPEAT 10 //Maximum iterations when energy and gradient are identical

static double *s,    //Search direction
    *rss,*rsg,  //CG Restart procedure data
    *ginit,          //Gradient at start of iteration
    *xopt,*gopt;//Optimal position

#ifdef __cplusplus
extern "C" void CGRelax(void (*func)(int,double*,double *, double*),
            int n,double acc,int maxfn,double dfpred,
            double x[],double g[],double *f, double *buffer);
#endif

#define MP_Sync() {}
#define MP_BeginMasterOnly() {}
#define MP_EndMasterOnly() {}
#define MP_Free(s) free(s)

#define Min(a, b) ( (a)<(b)?(a):(b))
#define Max(a, b) ( (a)>(b)?(a):(b))

void NumGrad(double (*func)(int,double*),int n, double *x, double *g)
{
    int i;
    const double dt=1e-8;
    for(i=0;i<n;i++)
    {
        double f1, f2;
        x[i]+=dt;
        f1=func(n,x);
        x[i]-=dt+dt;
        f2=func(n,x);
        x[i]-=dt;
        g[i]=(f1-f2)/2/dt;
    }
}

void CGRelax(void (*func)(int,double*,double *, double*),
            int n,double acc,int maxfn,double dfpred,
            double x[],double g[],double *f, double *buffer)
{
    register int i; //Loop variable
    int iretry, //Retry count
        iterc,  //Iteration count
        iterfm, //Most recent iterations that decrease F
        iterrs, //Iteration of the most recent restart
        ncalls, //Call count
        nfbeg,  //Call number at begin of line search
        nfopt;  //Call number at optimal position
    register double beta,ddspln,dfpr,finit,fmin,gamden,gama,gfirst,gmin,gnew,
        gspln,gsqrd,sbound,step,stepch,stmin,sum,temp;

    int n0, n1;

    /* Getting the code out of "infinite" loop */
    int nrepeat;
    double f_old, sum_old;

    nrepeat = 0;

    n0=0; n1=n;
    
    if(buffer==NULL)
    {
        s=(double *)realloc(s,sizeof(double)*n*6);
        if(s==NULL){ 
          printf("Out of memory in CGRelax\n"); 
          abort();
        }
    }
    else s=buffer;
    
    rss=s+n;
    rsg=rss+n;
    ginit=rsg+n;
    xopt=ginit+n;
    gopt=xopt+n;

    ddspln=gamden=0; //The line is purely to avoid warning of var unused
    
    iterc=iterfm=iterrs=0;

    func(n,x,f,g);ncalls=1; /* func has to be sync-ed in the end */

    printf("relax: 1st potential call finished.\n");
    
    /* operation on shared memory */
    for(i=n0;i<n1;i++) s[i]=-g[i]; /* g[i] sync-ed in func */
    for(sum=0,i=0;i<n;i++) { 
        sum+=g[i]*g[i];            /* g[i] sync-ed in func */

        printf("g[%d]=%g\n",i, g[i]);
        printf("s[%d]=%g\n",i, s[i]);
    }
    printf("sum = %g\n", sum);

    printf("############################################################\n");
    printf("iteration     neval    energy (eV)       |gradient|^2 (eV^2)\n");
    printf("############################################################\n");
    printf("%10d %10d %20.14e %20.14e\n",iterc,ncalls,*f,sum);
    
    /* Getting the code out of "infinite" loop */
    f_old = *f; sum_old = sum;
    if(sum<=acc) goto l_return; //another place of success at below

    gnew=-sum;
    fmin=*f;
    gsqrd=sum;
    nfopt=ncalls;
    for(i=n0;i<n1;i++)
    {
        xopt[i]=x[i];
        gopt[i]=g[i];
    }
    dfpr=dfpred;
    stmin=dfpred/gsqrd;
    
    for(;;)//BEGIN THE ITERATION
    {
        iterc++;
        finit=*f;
        for(i=n0;i<n1;i++)
        {
            ginit[i]=g[i];
        }
        gfirst=0;
        for(i=0;i<n;i++)
        {
            gfirst+=s[i]*g[i];
        }
        if(gfirst >=0)
        {
            printf("CGRelax: Search direction is uphill.\n");
            goto l_return;
        }
        gmin=gfirst;
        sbound=-1;
        nfbeg=ncalls;
        iretry=-1;
        stepch=Min(stmin, fabs(dfpr/gfirst));
        stmin=0;
        for(;;)
        {
            step=stmin+stepch;
            for(i=n0;i<n1;i++)
            {
                x[i]=xopt[i]+stepch*s[i];
            }
            temp=0;
            for(i=0;i<n;i++)
            {
                temp=Max(temp,fabs(x[i]-xopt[i]));
            }
            if(temp <= 0)
            {
                if(ncalls > nfbeg+1 || fabs(gmin/gfirst) > 0.2)
                {
                    printf("CGRelax: Line search aborted, \n"
                          "possible error in gradient.\n");
                    goto l_return;
                }
                else break;
            }
            ncalls++;
            
            func(n,x,f,g); /* g[i] has to be sync-ed in func */
            for(gnew=sum=0,i=0;i<n;i++)
            {
                gnew+=s[i]*g[i];
                sum+=g[i]*g[i];
            }
            if((*f < fmin || (*f==fmin && gnew/gmin>=-1)))
            {
                fmin=*f;
                gsqrd=sum;
                nfopt=ncalls;
                for(i=n0;i<n1;i++)
                {
                    xopt[i]=x[i];
                    gopt[i]=g[i];
                }
            }
            /* Print out iteration information */
            printf("%10d %10d %20.14e %20.14e\n",iterc,ncalls,*f,sum);

            /* Getting the code out of "infinite" loop */
            if( (f_old==*f)&&(sum_old==sum) )
            {
                nrepeat ++;
            }
            else
            {
                nrepeat = 0;
                f_old = *f; sum_old = sum;
            }
            if(nrepeat>=MAXREPEAT)
            {
                printf("CGRelax: getting stuck, stop...\n");
                goto l_return;
            }
            
            if(*f<=fmin && sum <= acc) goto l_return; //Successful return
            //TEST IF THE VALUE OF MAXFN ALLOWS ANOTHER CALL OF FUNCT.
            if(ncalls>=maxfn)
            {
                printf("CGRelax: Too many iterations, stop...\n");
                goto l_return;
            }
            temp=(*f+*f-fmin-fmin)/stepch-gnew-gmin;
            ddspln=(gnew-gmin)/stepch;
            if(ncalls > nfopt) sbound=step;
            else
            {
                if(gmin*gnew <= 0) sbound=stmin;
                stmin=step;
                gmin=gnew;
                stepch=-stepch;
            }
            if(*f!=fmin) ddspln+=(temp+temp)/stepch;
            if(gmin==0) break;
            if(ncalls >= nfbeg+1)
            {
                if(fabs(gmin/gfirst) <=0.2) break;
            l_retry:
                if(ncalls >= nfopt+MAXLIN)
                {
                    printf("CGRelax: Line search aborted,\n "
                          "possible error in gradient.\n");
                    goto l_return;
                }
            }
            stepch=0.5*(sbound-stmin);
            if(sbound < -0.5) stepch=9*stmin;
            gspln=gmin+stepch*ddspln;
            if(gmin*gspln<0) stepch*=gmin/(gmin-gspln);
        }
        if(ncalls!=nfopt)
        {
            *f=fmin;
            for(i=n0;i<n1;i++)
            {
                x[i]=xopt[i];
                g[i]=gopt[i];
            }
            //MP_Sync(); /* sync x[i], g[i] */
        }
        for(sum=0,i=0;i<n;i++)
            sum+=g[i]*ginit[i];
        //MP_Sync();
        beta=(gsqrd-sum)/(gmin-gfirst);
        if(fabs(beta*gmin) > 0.2*gsqrd)
        {
            iretry++;
            if(iretry<=0) goto l_retry;
        }
        if(*f<finit) iterfm=iterc;
        else if(iterc >= iterfm+MXFCON)
        {
            printf("CGRelax: Cannot reduce value of F, aborting...\n");
            goto l_return;
        }
        dfpr=stmin*gfirst;
        if(iretry>0) // Restart since we need to retry
        {
            for(i=n0;i<n1;i++) s[i]=-g[i];
            iterrs=0;
            continue;
        }

        if(iterrs!=0 && iterc-iterrs<n && fabs(sum)<0.2*gsqrd)
        {
            for(gama=sum=0,i=0;i<n;i++)
            {
                gama+=g[i]*rsg[i];
                sum+=g[i]*rss[i];  //tmp1+=g[i]; tmp2+=rsg[i]; tmp3+=rss[i];
            }
            gama/=gamden;
            if(fabs(beta*gmin+gama*sum) < 0.2*gsqrd)
            {
                for(i=n0;i<n1;i++) s[i]=-g[i]+beta*s[i]+gama*rss[i];
                continue;
            }
        }
        gamden=gmin-gfirst;

        for(i=n0;i<n1;i++)
        {
            rss[i]=s[i];
            rsg[i]=g[i]-ginit[i];
            s[i]=-g[i]+beta*s[i];
        }
        iterrs=iterc;
    }
 l_return:
    if(buffer==NULL) MP_Free(s);
}

