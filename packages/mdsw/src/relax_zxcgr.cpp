// Last Modified : Mon Mar 23 13:13:06 2009

//////////////////////////////////////////////////////
// Stand alone conjugate gradient relaxiation program
//
//    Based on IMSL's U2CGG.F
//
//  writen by Dongyi Liao at MIT 1999
//
//////////////////////////////////////////////////////

#include "general.h"
//#include "mplib.h"
#include "scparser.h"

#ifndef _GENERAL_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define ASSERT(a) ((void)0)
#define ERROR(s) puts(s)
#define Max(a,b) ((a) > (b)? (a):(b))
#define Min(a,b) ((a) > (b)? (b):(a))

#endif

#define MAXLIN 5 //Maximum iterations within a line search
#define MXFCON 2 //Maximum iterations before F can be decreased

/* Getting the code out of "infinite" loop */
#define MAXREPEAT 10 //Maximum iterations when energy and gradient are identical

/* Special Interface to potential_wrapper in md.h md2d.h */
//class MDFrame *__conj_simframe;
//class SIHMDFrame *__conj_sihmdframe;
//class SimFrame2d *__conj_simframe2d;
//class Morse2d *__conj_morse2d;
//class SpFrame2d *__conj_spframe2d;
/* End of Special Interface */

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

//void CGRelaxE(double (*func)(int,double*),
//              int n, double acc, int maxfn, double dfpred,
//              double x[], double g[],double *f)
//{
//}

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
    
    //Use the passed buffer if it's not NULL
    if(buffer==NULL)
    {
        /* allocate work space if empty */
        s=(double *)realloc(s,sizeof(double)*n*6);
        if(s==NULL){ ERROR("Out of memory in CGRelax"); abort(); }
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

    INFO("relax: 1st potential call finished.");
    
    /* operation on shared memory */
    //Let the initial search direction be minus the gradient vector.
    for(i=n0;i<n1;i++) s[i]=-g[i]; /* g[i] sync-ed in func */
    //MP_Sync();      /* always sync s[i] after modification */
    
    for(sum=0,i=0;i<n;i++)
        sum+=g[i]*g[i];            /* g[i] sync-ed in func */
    //MP_Sync(); /* can be removed if g[i] not modified during summation */
    INFO("############################################################");
    INFO("iteration     neval    energy (eV)       |gradient|^2 (eV^2)");
    INFO("############################################################");
    INFO_Printf("%10d %10d %20.14e %20.14e\n",iterc,ncalls,*f,sum);
    
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
    //MP_Sync(); /* sync xopt, gopt */
    //SET DFPR TO THE ESTIMATE OF THE REDUCTION IN F GIVEN IN THE
    //ARGUMENT LIST, IN ORDER THAT THE INITIAL CHANGE TO THE PARAMETERS
    //IS OF A SUITABLE SIZE. THE VALUE OF STMIN IS USUALLY THE
    //STEP-LENGTH OF THE MOST RECENT LINE SEARCH THAT GIVES THE LEAST
    //CALCULATED VALUE OF F.
    dfpr=dfpred;
    stmin=dfpred/gsqrd;
    
    for(;;)//BEGIN THE ITERATION
    {
        iterc++;
        //STORE THE INITIAL FUNCTION VALUE AND GRADIENT, CALCULATE THE
        //INITIAL DIRECTIONAL DERIVATIVE, AND BRANCH IF ITS VALUE IS NOT
        //NEGATIVE. SET SBOUND TO MINUS ONE TO INDICATE THAT A BOUND ON THE
        //STEP IS NOT KNOWN YET, AND SET NFBEG TO THE CURRENT VALUE OF
        //NCALLS. THE PARAMETER IRETRY SHOWS THE NUMBER OF ATTEMPTS AT
        //SATISFYING THE BETA CONDITION.
        finit=*f;
        for(i=n0;i<n1;i++)
        {
            ginit[i]=g[i];
        }
        //MP_Sync(); /* sync ginit */
        
        gfirst=0;
        for(i=0;i<n;i++)
        {
            gfirst+=s[i]*g[i];
        }
        //MP_Sync(); /* can be removed if s[i],g[i] not modified during sum*/
        if(gfirst >=0)
        {
            //SET IER TO INDICATE THAT THE SEARCH DIRECTION IS UPHILL.
            ERROR("CGRelax: Search direction is uphill.");
            goto l_return;
        }
        gmin=gfirst;
        sbound=-1;
        nfbeg=ncalls;
        iretry=-1;
        //SET STEPCH SO THAT THE INITIAL STEP-LENGTH IS CONSISTENT WITH THE
        //PREDICTED REDUCTION IN F, SUBJECT TO THE CONDITION THAT IT DOES
        //NOT EXCEED THE STEP-LENGTH OF THE PREVIOUS ITERATION. LET STMIN
        //BE THE STEP TO THE LEAST CALCULATED VALUE OF F.
        stepch=Min(stmin, fabs(dfpr/gfirst));
        stmin=0;
        //CALL SUBROUTINE FUNCT AT THE VALUE OF X THAT IS DEFINED BY THE
        //NEW CHANGE TO THE STEP-LENGTH, AND LET THE NEW STEP-LENGTH BE
        //STEP. THE VARIABLE WORK IS USED AS WORK SPACE.
        for(;;)
        {
            step=stmin+stepch;
            for(i=n0;i<n1;i++)
            {
                x[i]=xopt[i]+stepch*s[i];
            }
            //MP_Sync(); /* sync x[i] */
            temp=0;
            for(i=0;i<n;i++)
            {
                temp=Max(temp,fabs(x[i]-xopt[i]));
            }
            //MP_Sync(); /* can be removed if x[i], xopt[i] not modified */
            if(temp <= 0)
            {
                //TERMINATE THE LINE SEARCH IF STEPCH IS EFFECTIVELY ZERO.
                if(ncalls > nfbeg+1 || fabs(gmin/gfirst) > 0.2)
                {
                    ERROR("CGRelax: Line search aborted, "
                          "possible error in gradient.");
                    goto l_return;
                }
                else break;
            }
            ncalls++;
            
            func(n,x,f,g); /* g[i] has to be sync-ed in func */
            //DUMP("Function call: F("<<x[0]<<','<<x[1]<<")="<<*f);
            //SET SUM TO G SQUARED. GMIN AND GNEW ARE THE OLD AND THE NEW
            //DIRECTIONAL DERIVATIVES ALONG THE CURRENT SEARCH DIRECTION.
            for(gnew=sum=0,i=0;i<n;i++)
            {
                gnew+=s[i]*g[i];
                sum+=g[i]*g[i];
            }
            //MP_Sync(); /* can be removed if s[i], g[i] not modified */
            
            //STORE THE VALUES OF X, F AND G, IF THEY ARE THE BEST THAT
            //HAVE BEEN CALCULATED SO FAR, AND NOTE G SQUARED AND THE VALUE
            //OF NCALLS. TEST FOR CONVERGENCE.
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
                //MP_Sync(); /* sync xopt gopt */
            }
            /* Print out iteration information */
            INFO_Printf("%10d %10d %20.14e %20.14e\n",iterc,ncalls,*f,sum);

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
                ERROR("CGRelax: getting stuck, stop...");
                goto l_return;
            }
            
            if(*f<=fmin && sum <= acc) goto l_return; //Successful return
            //TEST IF THE VALUE OF MAXFN ALLOWS ANOTHER CALL OF FUNCT.
            if(ncalls>=maxfn)
            {
                ERROR("CGRelax: Too many iterations, stop...");
                goto l_return;
            }
            //LET SPLN BE THE QUADRATIC SPLINE THAT INTERPOLATES THE
            //CALCULATED FUNCTION VALUES AND DIRECTIONAL DERIVATIVES AT THE
            //POINTS STMIN AND STEP OF THE LINE SEARCH, WHERE THE KNOT OF
            //THE SPLINE IS AT 0.5*(STMIN+STEP). REVISE STMIN, GMIN AND
            //SBOUND, AND SET DDSPLN TO THE SECOND DERIVATIVE OF SPLN AT
            //THE NEW STMIN. HOWEVER, IF FCH IS ZERO, IT IS ASSUMED THAT
            //THE MAXIMUM ACCURACY IS ALMOST ACHIEVED, SO DDSPLN IS
            //CALCULATED USING ONLY THE CHANGE IN THE GRADIENT.
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
            //TEST FOR CONVERGENCE OF THE LINE SEARCH, BUT FORCE AT LEAST
            //TWO STEPS TO BE TAKEN IN ORDER NOT TO LOSE QUADRATIC
            //TERMINATION.
            if(gmin==0) break;
            if(ncalls >= nfbeg+1)
            {
                if(fabs(gmin/gfirst) <=0.2) break;
                //APPLY THE TEST THAT DEPENDS ON THE PARAMETER MAXLIN.
            l_retry:
                if(ncalls >= nfopt+MAXLIN)
                {
                    ERROR("CGRelax: Line search aborted, "
                          "possible error in gradient.");
                    goto l_return;
                }
            }
            //SET STEPCH TO THE GREATEST CHANGE TO THE CURRENT VALUE OF STMIN
            //THAT IS ALLOWED BY THE BOUND ON THE LINE SEARCH. SET GSPLN TO THE
            //GRADIENT OF THE QUADRATIC SPLINE AT (STMIN+STEPCH). HENCE
            //CALCULATE THE VALUE OF STEPCH THAT MINIMIZES THE SPLINE FUNCTION,
            //AND THEN OBTAIN THE NEW FUNCTION AND GRADIENT VECTOR, FOR THIS
            //VALUE OF THE CHANGE TO THE STEP-LENGTH.
            stepch=0.5*(sbound-stmin);
            if(sbound < -0.5) stepch=9*stmin;
            gspln=gmin+stepch*ddspln;
            if(gmin*gspln<0) stepch*=gmin/(gmin-gspln);
        }
        //ENSURE THAT F, X AND G ARE OPTIMAL.
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
        //CALCULATE THE VALUE OF BETA THAT OCCURS IN THE NEW SEARCH
        //DIRECTION.
        for(sum=0,i=0;i<n;i++)
            sum+=g[i]*ginit[i];
        //MP_Sync();
        beta=(gsqrd-sum)/(gmin-gfirst);
        //TEST THAT THE NEW SEARCH DIRECTION CAN BE MADE DOWNHILL. IF IT
        //CANNOT, THEN MAKE ONE ATTEMPT TO IMPROVE THE ACCURACY OF THE LINE
        //SEARCH.
        if(fabs(beta*gmin) > 0.2*gsqrd)
        {
            iretry++;
            if(iretry<=0) goto l_retry;
        }
        //APPLY THE TEST THAT DEPENDS ON THE PARAMETER MXFCON.
        //SET DFPR TO THE PREDICTED REDUCTION IN F ON THE NEXT ITERATION.
//        INFO_Printf("*f=%20.16e  finit=%20.16e\n",*f,finit);
        if(*f<finit) iterfm=iterc;
        else if(iterc >= iterfm+MXFCON)
        {
            ERROR("CGRelax: Cannot reduce value of F, aborting...");
            goto l_return;
        }
        dfpr=stmin*gfirst;
        if(iretry>0) // Restart since we need to retry
        {
            for(i=n0;i<n1;i++) s[i]=-g[i];
            //MP_Sync(); /* sync s[i] */
            iterrs=0;
            continue;
        }

        if(iterrs!=0 && iterc-iterrs<n && fabs(sum)<0.2*gsqrd)
        {
            //CALCULATE THE VALUE OF GAMA THAT OCCURS IN THE NEW SEARCH
            //DIRECTION, AND SET SUM TO A SCALAR PRODUCT FOR THE TEST BELOW.
            //THE VALUE OF GAMDEN IS SET BY THE RESTART PROCEDURE.
//            tmp1=tmp2=tmp3=0;
            for(gama=sum=0,i=0;i<n;i++)
            {
                gama+=g[i]*rsg[i];
                sum+=g[i]*rss[i];  //tmp1+=g[i]; tmp2+=rsg[i]; tmp3+=rss[i];
            }
            //MP_Sync(); /* Ensure all process complete summation before
            //              anyone changes rsg, rss */
            gama/=gamden;
            //RESTART IF THE NEW SEARCH DIRECTION IS NOT SUFFICIENTLY
            //DOWNHILL.
//            if(ncalls>0)
//            {
//                INFO_Printf("n=%d beta=%g gmin=%g gama=%g sum=%g"
//                            " gsqrd=%g gamden=%g\n",
//                            n,beta,gmin,gama,sum,gsqrd,gamden);
//                INFO_Printf("AFTER  g=%g rsg=%g rss=%g\n",
//                            tmp1,tmp2,tmp3);
//            }
            if(fabs(beta*gmin+gama*sum) < 0.2*gsqrd)
            {
                //CALCULATE THE NEW SEARCH DIRECTION.
                for(i=n0;i<n1;i++) s[i]=-g[i]+beta*s[i]+gama*rss[i];
                //MP_Sync(); /* sync s[i] */
                continue;
            }
        }
        //APPLY THE RESTART PROCEDURE.
        gamden=gmin-gfirst;

        for(i=n0;i<n1;i++)
        {
            rss[i]=s[i];
            rsg[i]=g[i]-ginit[i];
            s[i]=-g[i]+beta*s[i];
        }
        //MP_Sync(); /* sync rss, rsg, s */
//        for(tmp1=tmp2=tmp3=0,i=0;i<n;i++)
//        {
//            tmp1+=g[i]; tmp2+=rsg[i]; tmp3+=rss[i];
//        }
//        INFO_Printf("BEFORE g=%g rsg=%g rss=%g\n",
//                    tmp1,tmp2,tmp3);
        iterrs=iterc;
    }
 l_return:
    //MP_Sync(); /* sync before free buffer */
    if(buffer==NULL) MP_Free(s);
}

#ifdef _TEST

#ifdef __cplusplus
extern "C" //The one from zxcgr.f
{
#endif
    void zxcgr_(void (*func)(int*,double*,double*,double*),
             int *n,double *acc,int *maxfn,double *dfpred, 
             double x[],double g[],double *f,double w[],int *ier);
#ifdef __cplusplus
};
#endif


int ncalls;
double acc=1e-10;
#define Sqr(x) ((x)*(x))
void values(int n, double *x, double *f, double *g)
{
    ncalls++;
    ASSERT(n==2);
   *f=cos(x[0]-0.5)*cos(Sqr(x[1]+1))+Sqr(x[0]-0.5)/100;
    g[0]=-sin(x[0]-0.5)*cos(Sqr(x[1]+1))+2*(x[0]-0.5)/100;
    g[1]=-cos(x[0]-0.5)*sin(Sqr(x[1]+1))*2*(x[1]+1);
//    *f=Sqr(x[0]-x[1]+0.5)+20*Sqr(x[1]+x[0]-2)+1.0;
//    g[0]=2*(x[0]-x[1]+0.5)+40*(x[1]+x[0]-2);
//    g[1]=-2*(x[0]-x[1]+0.5)+40*(x[1]+x[0]-2);
    printf("Values(%f,%f)=%f(%f,%f)\n",x[0],x[1],*f,
           g[0]/sqrt(acc),g[1]/sqrt(acc));
}

void fvalue(int *n, double *x, double *f, double *g)
{
    values(*n,x,f,g);
}

double X[2],G[2],F,W[100];


int main()
{
    double x0,x1;
    int n, maxfn;
    int ier;
    double dfpred;

    srand(time(NULL));
    x0=rand()/100000000.0;
    x1=rand()/100000000.0;
    n=2;
    dfpred=0.02;
    maxfn=1000;

    X[0]=x0;
    X[1]=x1;
    ncalls=0;
    zxcgr_(fvalue,&n,&acc,&maxfn,&dfpred,X,G,&F,W,&ier);
    printf("Output:%f(%f,%f)[%f,%f] <%d> ier=%d\n",
           F,X[0],X[1],G[0]/sqrt(acc),G[1]/sqrt(acc),ncalls,ier);

    X[0]=x0;
    X[1]=x1;
    ncalls=0;
    CGRelax(values, n, acc, maxfn, dfpred,X,G,&F, NULL);
    printf("Output:%f(%f,%f)[%f,%f] <%d>\n",
           F,X[0],X[1],G[0]/sqrt(acc),G[1]/sqrt(acc),ncalls);
    return 0;
}

#endif
