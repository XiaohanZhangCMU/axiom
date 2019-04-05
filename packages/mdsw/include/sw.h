/*
  sw.h
  by Wei Cai  caiwei@mit.edu
  Last Modified : Thu Apr 12 18:24:05 2007

  FUNCTION  :  MD simulation package of Si using Stillinger-Weber potential
*/

#ifndef _SW_H
#define _SW_H

#include "md.h"

class SWFrame : public MDFrame /* Si with Stillinger-Weber potential */
{
    /* Stillinger-Weber potential parameters */
    double aa, bb, plam, pgam, acut, pss, rho;
    double rho1, acutsq;
    double _SW3B_MUL, _SW2B_MUL,tote2, tote3;
    
public:
    SWFrame():aa(0),bb(0),plam(0),pgam(0),acut(0),pss(0),rho(0),
              rho1(0),acutsq(0),_SW3B_MUL(1), _SW2B_MUL(1), tote2(0), tote3(0) {};
    void stillinger_weber();
    void stillinger_weber_energyonly();
    double stillinger_weber_energyonly(int iatom);
    virtual void potential();
    virtual void potential_energyonly();
    virtual double potential_energyonly(int iatom);
    void stillinger_weber_switch_3body(double lambda);
    
    virtual void initvars();
};

#endif // _SW_H

