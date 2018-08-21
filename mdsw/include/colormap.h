/*
  colormap.h
  by Wei Cai  caiwei@mit.edu
  Last Modified : Sat Sep 29 10:39:27 2001

  FUNCTION  :
*/

#ifndef _COLORMAP_H
#define _COLORMAP_H

//void colormap1(double x,double c[3])
static void colormap1(double x, int *r, int *g, int *b)
{
    int i; double alpha;
#define map1n 6
    static double xtick[map1n]={0,0.125,0.375,0.625,0.875,1};
    static double ctick[map1n][3]={ {0,0,0.9}, {0,0,1}, {0,1,1}, {1,1,0},
                             {1,0,0}, {0.9,0,0}};
//    double xtick[map1n]={0,0.125,0.375,0.625,0.875,1};
//    double ctick[map1n][3]={ {0,0,0.5}, {0,0,1}, {0,1,1}, {1,1,0},
//                             {1,0,0}, {0.5,0,0}};
    double c[3];
//    INFO("x="<<x);
    if (x<xtick[0]) {
        c[0]=ctick[0][0];
        c[1]=ctick[0][1];
        c[2]=ctick[0][2];
        *r=(int)floor(c[0]*255);
        *g=(int)floor(c[1]*255);
        *b=(int)floor(c[2]*255);
//        INFO_Printf("x(%f)<xtick[0](%f) r=%d g=%d b=%d\n",x,xtick[0],*r,*g,*b);
        return;
    };
    if (x>=xtick[map1n-1]) {
        c[0]=ctick[map1n-1][0];
        c[1]=ctick[map1n-1][1];
        c[2]=ctick[map1n-1][2];
        *r=(int)floor(c[0]*255);
        *g=(int)floor(c[1]*255);
        *b=(int)floor(c[2]*255);
//        INFO_Printf("x(%f)>xtick[%d](%f) r=%d g=%d b=%d\n",x,map1n-1,xtick[map1n-1],*r,*g,*b);
        return;
    };
    for(i=1;i<=map1n-1;i++) {
        if(x<xtick[i])
        {
            alpha=(x-xtick[i-1])/(xtick[i]-xtick[i-1]);
            c[0]=(1-alpha)*ctick[i-1][0]+alpha*ctick[i][0];
            c[1]=(1-alpha)*ctick[i-1][1]+alpha*ctick[i][1];
            c[2]=(1-alpha)*ctick[i-1][2]+alpha*ctick[i][2];
            *r=(int)floor(c[0]*255);
            *g=(int)floor(c[1]*255);
            *b=(int)floor(c[2]*255);
//            INFO_Printf("i=%d alpha=%e r=%d g=%d b=%d\n",i,alpha,*r,*g,*b);
            return;
        }
    }
}

static void colormap2(double x, int *r, int *g, int *b)
{
    int i; double alpha;
#define map2n 5
    static double xtick[map2n]={0,0.25,0.5,0.75,1};
    static double ctick[map2n][3]={ {0,0.2,1}, {0.25,0.2,0.75}, {0.5,0.2,0.5},
                              {0.75,0.2,0.25}, {1,0.2,0}};
//    double ctick[map2n][3]={ {1,0,1}, {1,0.25,0.75}, {1,0.5,0.5},
//                             {1,0.75,0.25}, {1,1,0} };
//    double ctick[map2n][3]={ {0,1,1}, {0.25,1,0.75}, {0.5,1,0.5},
//                              {0.75,1,0.25}, {1,1,0}};
    double c[3];
    if (x<xtick[0]) {
        c[0]=ctick[0][0];
        c[1]=ctick[0][1];
        c[2]=ctick[0][2];
        *r=(int)floor(c[0]*255);
        *g=(int)floor(c[1]*255);
        *b=(int)floor(c[2]*255);
//        INFO_Printf("x(%f)<xtick[0](%f) r=%d g=%d b=%d\n",x,xtick[0],*r,*g,*b);
        return;
    };
    if (x>=xtick[map2n-1]) {
        c[0]=ctick[map2n-1][0];
        c[1]=ctick[map2n-1][1];
        c[2]=ctick[map2n-1][2];
        *r=(int)floor(c[0]*255);
        *g=(int)floor(c[1]*255);
        *b=(int)floor(c[2]*255);
//        INFO_Printf("x(%f)>xtick[%d](%f) r=%d g=%d b=%d\n",x,map2n-1,xtick[map2n-1],*r,*g,*b);
        return;
    };
    for(i=1;i<=map2n-1;i++) {
        if(x<xtick[i])
        {
            alpha=(x-xtick[i-1])/(xtick[i]-xtick[i-1]);
            c[0]=(1-alpha)*ctick[i-1][0]+alpha*ctick[i][0];
            c[1]=(1-alpha)*ctick[i-1][1]+alpha*ctick[i][1];
            c[2]=(1-alpha)*ctick[i-1][2]+alpha*ctick[i][2];
            *r=(int)floor(c[0]*255);
            *g=(int)floor(c[1]*255);
            *b=(int)floor(c[2]*255);
//            INFO_Printf("i=%d alpha=%e r=%d g=%d b=%d\n",i,alpha,*r,*g,*b);
            return;
        }
    }
}

#endif // _COLORMAP_H

