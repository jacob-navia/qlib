/*   mtst.c
 Consistency tests for math functions.

 With NTRIALS=10000, the following are typical results for
 14x32 bit arithmetic:

Consistency test of math functions.
Max and rms errors for 1000 random arguments.
A = absolute error criterion (but relative if >1):
Otherwise, estimate is of relative error
x =   sqrt( square(x) ):  max =  1.541E-106   rms =  2.667E-107
x =   atan(    tan(x) ):  max =  5.449E-106   rms =  1.262E-106
x =   cbrt(   cube(x) ):  max =  2.152E-106   rms =  5.761E-107
x =    sin(   asin(x) ):  max =  1.814E-105   rms =  9.178E-106
x =    log(    exp(x) ):  max =  1.042E-105   rms =  1.109E-106
x =   log2(   exp2(x) ):  max =  6.311E-106 A rms =  1.257E-106 A
x =  log10(  exp10(x) ):  max =  3.761E-106   rms =  1.093E-106
x =  acosh(   cosh(x) ):  max =  2.862E-106   rms =  9.864E-107
x = pow( pow(x,a),1/a ):  max =  2.739E-104   rms =  1.402E-105
x =   tanh(  atanh(x) ):  max =  1.433E-103   rms =  5.784E-105
x =  asinh(   sinh(x) ):  max =  2.860E-106   rms =  9.945E-107
x =    cos(   acos(x) ):  max =  1.744E-105 A rms =  4.624E-106 A
x =  ndtri(   ndtr(x) ):  max =  1.567E-93   rms =  9.056E-95
Legendre  ellpk,  ellpe:  max =  8.796E-104   rms =  6.872E-105
Absolute error and only 200 trials:
Wronksian of   Yn,   Jn:  max =  2.512E-54 A rms =  3.864E-55 A
Absolute error and only 200 trials:
Wronksian of   Kn,   Iv:  max =  1.276E-25 A rms =  9.033E-27 A
lgam(x) = log(gamma(x)):  max =  1.422E-104 A rms =  3.050E-105 A
-----------------------------------------------------------------
For 8 x 64 bits mantissa (132 digits 448 bits) (jacob)

Consistency test of math functions.
Max and rms errors for 1000 random arguments.
A = absolute error criterion (but relative if >1):
Otherwise, estimate is of relative error
x =   sqrt( square(x) ):  max = 5.066E-0135   rms = 1.3E-0135
x =   atan(    tan(x) ):  max = 3.49E-0135   rms = 1.319E-0135
x =   cbrt(   cube(x) ):  max = 2.635E-0135   rms = 2.748E-0136
x =    sin(   asin(x) ):  max = 9.897E-0135   rms = 3.486E-0135
x =    log(    exp(x) ):  max = 1.48E-0134   rms = 9.774E-0136
x =   log2(   exp2(x) ):  max = 1.169E-0134 A rms = 1.964E-0135 A
x =  log10(  exp10(x) ):  max = 2.024E-0134   rms = 5.116E-0135
x =  acosh(   cosh(x) ):  max = 2.743E-0135   rms = 7.738E-0136
x = pow( pow(x,a),1/a ):  max = 1.393E-0133   rms = 9.34E-0135
x =   tanh(  atanh(x) ):  max = 2.06E-0134   rms = 2.396E-0135
x =  asinh(   sinh(x) ):  max = 2.745E-0135   rms = 7.155E-0136
x =    cos(   acos(x) ):  max = 1.238E-0134 A rms = 2.092E-0135 A
lgam(x) = log(gamma(x)):  max = 6.191E-0135 A rms = 1.003E-0135 A

*/

/*
Cephes Math Library Release 2.1:  December, 1988
Copyright 1984, 1987, 1988, 1994 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "qhead.h"

#define NTRIALS 10000
#define WTRIALS (NTRIALS/5)
#define STRTST 0

extern int IMPORT merror;

#if 1
/* Provide inverses for square root and cube root: */
int Qsquare (Qfloatp x,Qfloatp y)
{
  qmul (x, x, y);
  return 0;
}
#endif

int qcube (Qfloatp x,Qfloatp y)
{
  qmul (x, x, y);
  qmul (x, y, y);
  return 0;
}

/* lookup table for each function */
struct fundef
  {
    char *nam1;			/* the function */
    int (*name)();
    char *nam2;			/* its inverse  */
    int (*inv)();
    int nargs;			/* number of function arguments */
    int tstyp;			/* type code of the function */
    long ctrl;			/* relative error flag */
    double arg1w;		/* width of domain for 1st arg */
    double arg1l;		/* lower bound domain 1st arg */
    long arg1f;			/* flags, e.g. integer arg */
    double arg2w;		/* same info for args 2, 3, 4 */
    double arg2l;
    long arg2f;
  };


/* fundef.ctrl bits: */
#define RELERR 1
#define EXPSCAL 4

/* fundef.tstyp  test types: */
#define POWER 1
#define ELLIP 2
#define GAMMA 3
#define WRONK1 4
#define WRONK2 5
#define WRONK3 6
#define NDTR 7

/* fundef.argNf  argument flag bits: */
#define INT 2

#define NTESTS 17-1
struct fundef defs[] =
{
  {"square", (int (*)())Qsquare, "  sqrt", (int (*)())qfsqrt, 1, 0, 1, 170.0, -85.0, EXPSCAL, 0.0, 0.0, 0},
  {"   tan", (int (*)())qftan, "  atan",   (int (*)())qatn, 1, 0, 1, 0.0, 0.0, 0, 0.0, 0.0, 0},
  {"  cube", (int (*)())qcube, "  cbrt", (int (*)())qcbrt, 1, 0, 1, 2000.0, -1000.0, 0, 0.0, 0.0, 0},
  {"  asin", (int (*)())qasin, "   sin", (int (*)())qfsin, 1, 0, 1, 2.0, -1.0, 0, 0.0, 0.0, 0},
  {"   exp", (int (*)())qfexp, "   log", (int (*)())qflog, 1, 0, 1, 340.0, -170.0, 0, 0.0, 0.0, 0},
  {"  exp2", (int (*)())qexp2, "  log2", (int (*)())qlogtwo, 1, 0, 0, 340.0, -170.0, 0, 0.0, 0.0, 0},
  {" exp10", (int (*)())qexp10, " log10", (int (*)())qlog10, 1, 0, 1, 340.0, -170.0, 0, 0.0, 0.0, 0},
  {"  cosh", (int (*)())qcosh, " acosh", (int (*)())qacosh, 1, 0, 1, 340.0, 2.0, 0, 0.0, 0.0, 0},
  {   "pow", (int (*)())qfpow, "pow", (int (*)())qfpow, 2, POWER, 1, 25.0, 0.0, 0, 50.0, -25.0, 0},
  {" atanh", (int (*)())qatanh, "  tanh", (int (*)())qtanh, 1, 0, 1, 0.9, 0.1, 0,   0.0, 0.0, 0},
  {"  sinh", (int (*)())qsinh, " asinh", (int (*)())qasinh, 1, 0, 1, 340.0, 0.0, 0,   0.0, 0.0, 0},
  {"  acos", (int (*)())qacos, "   cos", (int (*)())qfcos, 1, 0, 0, 2.0, -1.0, 0, 0.0, 0.0, 0},
  {"  ndtr", (int (*)()) qndtr,   " ndtri",(int (*)())qndtri, 1, NDTR, 1,   10.0,   -10.0,   0,     0.0, 0.0, 0},
  {" ellpe", (int (*)())qellpe, " ellpk",(int (*)())  qellpk, 1, ELLIP, 1, 1.0,  0.0,   0,     0.0, 0.0, 0},
  //{"  Jn",(int (*)())qjn, "  Yn",(int (*)()) qyn, 2, WRONK1, 0, 30.0,  0.1,  0,  40.0, -20.0, INT},
  //{"  Iv", (int (*)())qin,"  Kn",(int (*)()) qkn, 2, WRONK2, 0,  30.0, 0.1,  0, 20.0, 0.0, INT},
  {"gamma",  (int (*)())qgamma,   "lgam",   (int (*)())qlgam, 1, GAMMA, 0, 11.0, 1.0,   0, 0.0, 0.0, 0},
};

static char *headrs[] =
{
  "x = %s( %s(x) ): ",
  "x = %s( %s(x,a),1/a ): ",	/* power */
  "Legendre %s, %s: ",		/* ellip */
  "%s(x) = log(%s(x)): ",	/* gamma */
  "Wronksian of %s, %s: ",	/* wronk1 */
  "Wronksian of %s, %s: ",	/* wronk2 */
  "Wronksian of %s, %s: ",	/* wronk3 */
  "x = %s( %s(x) ): ",
};

static Qfloat Y1[1], y2[1], y3[1], y4[1], a[1];
static Qfloat x[1], y[1], z[1], e[1], max[1], rmsa[1], rms[1], ave[1];

#if 1
union {
  double d;
  unsigned short i[4];
} dx, dy;
#else
double dx,dy;
#endif

static Qfloat q3[1], temp[1];
static Qfloat base1[1], width1[1], base2[1], width2[1], fuzz[1];
static char str1[40], str2[40];
extern char *__Accelerators;

void drand(double *);
int main(void)
{
  int (*fun1) (Qfloatp,Qfloatp)=NULL;
  int (*ifun1) (Qfloatp,Qfloatp)=NULL;
  int (*fun2)(Qfloatp,Qfloatp,Qfloatp)=NULL;
  int (*ifun2)(Qfloatp,Qfloatp,Qfloatp)=NULL;
  struct fundef *d;
  int i, itst;
  int ntr;
  long long ll;
  int m;
  time_t rawtime;
  struct tm * timeinfo;

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  ntr = NTRIALS;
  printf ("Consistency test of math functions: %s",asctime (timeinfo));
  printf ("Max and rms errors for %d random arguments.\n",
	  ntr);
  printf ("A = absolute error criterion (but relative if >1):\n");
  printf ("Otherwise, estimate is of relative error\n");

  /* Initialize machine dependent parameters to test near the
 * largest an smallest possible arguments.  To compare different
 * machines, use the same test intervals for all systems.
 */
  defs[1].arg1w = 3.14159;
  defs[1].arg1l = -3.14159 / 2.0;
  qincr ( qtwo, q3);
  dx.d = 1.0e-15;
  //dx.d = 1.0e-45;
  etoq (dx.d, fuzz);

  /* Outer loop, on the test number: */

  for (itst = STRTST; itst < (sizeof(defs)/sizeof(defs[0])); itst++)
    {
      d = &defs[itst];
      m = 0;
      qclear (max);
      qclear (rmsa);
      qclear (ave);
      if (d->nargs == 1) {
         fun1 = (int (*)(Qfloatp,Qfloatp))d->name;
         ifun1 = (int (*)(Qfloatp,Qfloatp))d->inv;
      }
      else if (d->nargs == 2) {
         fun2 = (int (*)(Qfloatp,Qfloatp,Qfloatp))d->name;
         ifun2= (int (*)(Qfloatp,Qfloatp,Qfloatp))d->name;
      }
      etoq ((d->arg2w), width2);
      etoq ((d->arg2l), base2);
      etoq ((d->arg1w), width1);
      etoq ((d->arg1l), base1);

      /* Smaller number of trials for Wronksians
	 (put them at end of list)  */
      if (d->tstyp == WRONK1 || d->tstyp == WRONK2 || d->tstyp == NDTR)
	{
	  ntr = WTRIALS;
	  printf ("Absolute error and only %d trials:\n", ntr);
	}
      printf (headrs[d->tstyp], d->nam2, d->nam1);

      for (i = 0; i < ntr; i++)
	{
	  m++;

	  /* make random number(s) in desired range(s) */
	  switch (d->nargs)
	    {

	    default:
	      goto illegn;

	    case 2:
	      drand (&dx.d);
	      drand (&dy.d);
	      /* a = 1.0L + ((double )dx * dy - 1.0)/3.0; */
	      etoq (dx.d, Y1);
	      etoq (dy.d, y2);
	      qmul (Y1, y2, Y1);
	      qsub (qone, Y1, Y1);
	      qdiv (q3, Y1, Y1);
	      /*qadd( qone, Y1, Y1 );*/
	      /* a = d->arg2w *  ( a - 1.0L )  +  d->arg2l; */
	      qmul (width2, Y1, Y1);
	      qadd (base2, Y1, a);

	      if (d->arg2f & EXPSCAL)
		{
		  qfexp (a, a);
		  drand (&dx.d);
		  drand (&dy.d);
		  /* y2 = 1.0L + ((double )dx * dy - 1.0)/3.0; */
		  etoq (dx.d, y3);
		  etoq (dy.d, y2);
		  qmul (y3, y2, y2);
		  qsub (qone, y2, y2);
		  qdiv (q3, y2, y2);
		  qincr( y2, y2);
		  /*	a -= 1.0e-13L * a * y2; */
		  qmul (a, y2, y2);
		  qmul (fuzz, y2, y2);
		  qsub (y2, a, a);
		}
	      if (d->arg2f & INT)
		{
		  qadd (a, qhalf, a);
		  qifrac (a, &ll, y2);
		  lltoq (ll, a);
		}

	    case 1:
	      drand (&dx.d);
	      drand (&dy.d);
//printf("%8.6g\n",dx.d);
	      /* x = 1.0L + ((double )dx * (double )dy - 1.0)/3.0; */
	      etoq (dx.d, Y1);
//prt(Y1);
	      etoq (dy.d, y2);
	      qmul (Y1, y2, Y1);
	      qsub (qone, Y1, Y1);
	      qdiv (q3, Y1, Y1);
	      /*qadd( qone, Y1, Y1 );*/
	      qmul (width1, Y1, Y1);
	      qadd (base1, Y1, x);
	      if (qcmp (x, base1) < 0)
		qmov (base1, x);
	      qadd (base1, width1, Y1);
	      if (qcmp (x, Y1) > 0)
		qmov (Y1, x);
	      if (d->arg1f & EXPSCAL)
		{
		if (x->mantissa[0] == 0)
			printf("qmtst: mantissa[0]==0!\n");
		  qfexp (x, x);
		  drand (&dx.d);
		  drand (&dy.d);
		  /*	a = 1.0L + ((double )dx * dy - 1.0)/3.0; */
		  etoq (dx.d, Y1);
		  etoq (dy.d, y2);
		  qmul (Y1, y2, Y1);
		  qsub (qone, Y1, Y1);
		  qdiv (q3, Y1, Y1);
		  /*	x += 1.0e-13L * x * a; */
		  qmul (x, y2, y2);
		  qmul (fuzz, y2, y2);
		  qsub (y2, x, x);
		}
	    }

	  /* compute function under test */
	  switch (d->nargs)
	    {
	    case 1:
	      switch (d->tstyp)
		{
		case ELLIP:
		  (*(fun1)) (x, Y1);
		  qsub (x, qone, y3);
		  (*(fun1)) (y3, y2);
		  (*(ifun1)) (y3, y4);
		  (*(ifun1)) (x, y3);
		  break;

		case GAMMA:
		  qlgam (x, y);
		  qgamma(x, z);
		  qflog(z, x);
		  break;

		default:
		  (*(fun1)) (x, z);
		  (*(ifun1)) (z, y);
		}
/*
if( merror )
	{
	printf( "error: x = %.15e, z = %.15e, y = %.15e\n",
	 (double )x, (double )z, (double )y );
	}
*/
	      break;

	    case 2:
	      if (d->arg2f & INT)
		{
		  switch (d->tstyp)
		    {
		    case WRONK1:
		      (*fun2) (a, x, Y1);	/* jn */
		      qincr(a, y3);
		      (*fun2) (y3, x, y2);
		      (*ifun2) (y3, x, y4);
		      (*ifun2) (a, x, y3);	/* yn */
		      break;

#if 0
		    case WRONK2:
              ifun1 = (int (*)(int,Qfloatp,Qfloatp))ifun1;
		      (*fun2) (a, x, Y1);	/* iv */
		      qincr(a, y3);
		      (*fun2) (y3, x, y2);
		      (*ifun2) (k, x, y3);	/* kn */
		      (*ifun2) (k+1, x, y4);
		      break;
#endif

		    default:
		      (*fun2) (a, x, z);
		      (*ifun2) (a, z, y);
		    }
		}
	      else
		{
		  if (d->tstyp == POWER)
		    {
		      (*fun2) (x, a, z);
		      qinv(a, y3);
		      (*ifun2) (z, y3, y);
		    }
		  else
		    {
		      (*fun2) (a, x, z);
		      (*ifun2) (a, z, y);
		    }
		}
	      break;


	    default:
	    illegn:
	      printf ("Illegal nargs= %d", d->nargs);
	      exit (1);
	    }

	  switch (d->tstyp)
	    {
	    case WRONK1:
	      /*	e = (y2*y3 - Y1*y4) - 2.0L/(PIL*x); *//* Jn, Yn */
	      qmul (Y1, y4, temp);
	      qmul (y2, y3, e);
	      qsub (temp, e, e);
	      qmov (qpi, temp);
	      qmul (temp, x, temp);
	      qdiv (temp, qtwo, temp);
	      qsub (temp, e, e);
	      break;

	    case WRONK2:
	      /*	e = (y2*y3 + Y1*y4) - 1.0L/x; *//* In, Kn */
	      qmul (Y1, y4, temp);
	      qmul (y2, y3, e);
	      qadd (temp, e, e);
	      qinv (x, temp);
	      qsub (temp, e, e);
	      break;

	    case ELLIP:
	      /*	e = (Y1-y3)*y4 + y3*y2 - PIO2L; */
	      qsub (y3, Y1, temp);
	      qmul (y4, temp, e);
	      qmul (y3, y2, temp);
	      qadd (temp, e, e);
	      qmov (qpi, temp);
	      temp[0].exponent -= 1;
	      qsub (temp, e, e);
	      break;

	    default:
	      /*	e = y - x; */
	      qsub (x, y, e);
	      break;
	    }

	  if (d->ctrl & RELERR)
	    {
	      if (qcmp (x, qzero) != 0)
		qdiv (x, e, e);
	      else
		printf ("warning, x == 0\n");
	    }
	  else
	    {
	      qmov (x, temp);
	      temp[0].sign=0;
	      if (qcmpto1(temp) > 0)
		qdiv (x, e, e);
	    }

	  qadd (e, ave, ave);
	  /* absolute value of error */
	  e[0].sign=0;
	  /* peak detect the error */
	  if (qcmp (e, max) > 0)
	    {
	      qmov (e, max);
#if 0
	      if (e > 1.0e-10L)
		{
		  da = x;
		  db = z;
		  dc = y;
		  dd = max;
		  printf ("x %.6E z %.6E y %.6E max %.4E\n",
			  da, db, dc, dd);
		}
#endif
	    }
	  /* accumulate rms error	*/
	  if( e[0].exponent )
	    {
	      qmul (e, e, temp);
	      qadd (temp, rmsa, rmsa);
	    }
	}

      /* report after NTRIALS trials */
      itoq (m, temp);
      qdiv (temp, rmsa, temp);
      qfsqrt (temp, rms);
      qtoasc (max, str1,5, 5,0);
      qtoasc (rms, str2,5, 5,0);
      if (d->ctrl & RELERR)
	printf (" max = %s   rms = %s\n", str1, str2);
      else
	printf (" max = %s A rms = %s A\n", str1, str2);
    }				/* loop on itst */
	return(0);
}

