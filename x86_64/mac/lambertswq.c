#include <qfloat.h>

/* Implementation of Lamberts W-function which is defined as
 * w(x)*e^(w(x))=x
 * Implementation by Gunter Kuhnle, gk@uni-leipzig.de
 * Algorithm originally developed by
 * KEITH BRIGGS, DEPARTMENT OF PLANT SCIENCES,
 * ANSI C code for W(x).  K M Briggs 98 Feb 12 
 * http://keithbriggs.info/W-ology.html

   Based on Halley iteration.  Converges rapidly for all valid x.

  double W(const double x) {
  int i; double p,e,t,w,eps=4.0e-16;  eps=desired precision 
  if (x<-0.36787944117144232159552377016146086) {
    fprintf(stderr,"x=%g is < -1/e, exiting.\n",x); exit(1); }
  if (0.0==x) return 0.0;
   get initial approximation for iteration... 
  if (x<1.0) {  series near 0 
    p=sqrt(2.0*(2.7182818284590452353602874713526625*x+1.0));
    w=-1.0+p-p*p/3.0+11.0/72.0*p*p*p;
  } else w=log(x);  asymptotic 
  if (x>3.0) w-=log(w);
  for (i=0; i<20; i++) {  Halley loop 
    e=exp(w); t=w*e-x;
    t/=e*(w+1.0)-0.5*(w+2.0)*t/(w+1.0); w-=t;
    if (fabs(t)<eps*(1.0+fabs(w))) return w;  rel-abs error 
  }  never gets here 
  fprintf(stderr,"No convergence at x=%g\n",x); exit(1);
}
*/

#if 0
static unsigned int limit[] = { // -exp(-1)
0xffffffff,0x00007fff,0x00000000,0xbc5ab1b1,0x6779be35,0x75bd8f05,0x20a9f21b,
0xb5300b55,0x6ad8ee66,0x604973a1,0x4a0fb5db,0x62c8017e,0x8e56842b,0x7048b6ce
};
#endif
//#include "qhead.h"
static qfloat limit = {
0,0x00080001,{0xd3094c70f034de4bULL,0x96ff7d5b6f99fcd8ULL,0xfb28f8b60985a3acULL,0xe225fe4831b3f966ULL,
 0x177ab93cd4891710ULL,0x9f376dfbcb8acc8bULL,0xa4f6acc67b30f91dULL}
};


#if 0
static unsigned int qepsilon{]={
0x00000000,0x00007ea1,0x00000000,0x80000000,0x00000000,0x00000000,0x00000000,
0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000,0x00000000
};
#endif
#ifndef NQ
#define NQ 14
#endif

qfloat lambertwq(qfloat x)
{
    qfloat p, e, t, w, eps;
    int i;
    qfloat fabs_x[1];

    eps = 1e-101q;

    if (qcmp(x,limit) < 0) // if (x < -exp(-1.0))
	return -1;              /* error, value undefined */

    qmov(x,fabs_x);
    fabs_x[0]=0;
    if (qcmp(fabs_x,qepsilon) <=0)  //if (fabs(x) <= eps)
	return x;

    if (qcmp(x,qone) < 0) { // if (x < 0)
	p = sqrtq(2.0q * (exp(1.0q) * x + 1.0q));
	w = -1.0q + p - p * p / 3.0 + 11.0q / 72.0q * p * p * p;
    } else {
	w = log(x);
    }

    if (x > 3) {
	w = w - log(w);
    }
    for (i = 0; i < 20; i++) {
	e = exp(w);
	t = w * e - x;
	t = t / (e * (w + 1.0) - 0.5q * (w + 2.0q) * t / (w + 1.0));
	w = w - t;
	if (fabs(t) < eps * (1.0 + fabs(w)))
	    return w;
    }
    return -1;                 /* error: iteration didn't converge */
}
#ifdef TEST
#include <stdio.h>
int main(void)
{
	for (int i=0; i<10;i++) {
		printf("W(%d)=%100.99qf\n",i,lambertwq((qfloat)i));
	}
}
#endif
