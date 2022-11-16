/*
   Polylogarithms.


               inf   k
                -   x
   Li (x)  =    >   ---
     n          -     n
               k=1   k


                 x
                  -
                 | |  -ln(1-t)
    Li (x)  =    |    -------- dt
      2        | |       t
                -
                 0


                1-x
                  -
                 | |  ln t
            =    |    ------ dt   =   spence(1-x)
               | |    1 - t
                -
                 1


                        2       3
                       x       x
             =  x  +  ---  +  ---  +  ...
                       4       9


  d                 1
  --   Li (x)  =   ---  Li   (x)
  dx     n          x     n-1

  Series expansions are set to terminate at less than full
  working precision.

  */
#include "qhead.h"
static void cxplog(int n,Qfloatp x, Qfloatp y);

void qpolylog (int n,Qfloatp x, Qfloatp y)
{
  Qfloat a[1], p[1], s[1], k[1], qn[1];
  int ln;

  ln = n;
  itoq (ln, qn);


  if ((qcmp(qone, x) == 0) && (n > 1))
    {
      qzetac (qn, y);
      //qadd (qone, y, y);
      qincr(y,y,0);
      return;
    }

  /*
  if (n == 2) {
      qmov(x, s);
      qsub(qone, s, s);
      qspenc (s, y);
      return;
    }
  */

  /* Compare x with 3/4.  */
  ln = 3;
  itoq (ln, a);
  a[0].exponent -= 2;
  if (qcmp(x,a) > 0)
    {
      cxplog (n, x, y);
      return;
    }


  /* Defining power series.  */
  qmov (x, p);
  qmov (x, s);
  qmov (qtwo, k);
  qneg(qn);
  do {
      qmul( p, x, p);
      qpowi(k, qn, a);
      qmul(p, a, a);
      qadd(s, a, s);
	qincr(k,k,0);
      if (k[0].exponent > (qone[0].exponent + 19)) {
	  mtherr("qpolylog", PLOSS);
	  /* ln = (int) s[1] - (int) a[1];
	     printf("%ld\n", ln); */
	  break;
	}
  } while (((int) s[0].exponent - (int) a[0].exponent) < NBITS);
  qmov (s, y);
}

/*  This expansion in powers of log(x) is especially useful when
    x is near 1.

    See also the pari calculator.


                      inf           j
                       -    z(n-j) w
    polylog(n,x)  =    >   ----------
                       -       j!
                      j=0

    where

      w = log(x)

      z(j) = zeta(j), j != 1

                               n
                               -
      z(1) =  -log(-log(x)) +  >  1/k
                               -
                              k=1

  */

static void cxplog(int n,Qfloatp x, Qfloatp y)
{
  Qfloat z[1], h[1], q[1], p[1], s[1];
  int j, li;

  qflog (x, z);  /* z = log(x); */
  qmov (z, q);  /* h = -log(-z); */
  qneg (q);
  qflog (q, h);
  qneg(h);

  for (j = 1; j < n; j++)
    {
      /* h = h + 1.0/i; */
      itoq (j, q);
      qdiv (q, qone, q);
      qadd (h, q, h);
    }

  qmov (qone, q); /* q = 1.0; */
  j = n;  /* s = zetac((double)n) + 1.0; */
  itoq (j, p);
  qzetac (p, s);
  qincr(s,s,0);

  for (j=1; j<=n+1; j++)
  {
    itoq (j, p);  /* q = q * z / j; */
    qdiv (p, z, p);
    qmul (q, p, q );
    if (j == n-1)
      {
	/* s = s + h * q; */
	qmul (h, q, p);
	qadd (s, p, s);
      }
    else
      {
	/* s = s + (zetac((double)(n-j)) + 1.0) * q; */
	li = n - j;
	itoq (li, p);
	qzetac (p, p);
	qadd (qone, p, p);
	qmul (q, p, p);
	qadd (s, p, s);
      }
  }
  j = n + 3;
  qmul (z, z, z); /* z = z * z; */
  for(;;) {
      /* q = q * z / ((j-1)*j); */
      li = (j-1) * j;
      itoq (li, p);
      qdiv (p, z, p);
      qmul (q, p, q);
      /* p1 = (zetac((double)(n-j)) + 1.0); */
      li = n - j;
      itoq (li, p);
      qzetac (p, p);
      qadd (qone, p, p);
      /* p1 = p1 * q; */
      qmul (p, q, p);
      /* s = s + p1; */
      qadd (s, p, s);
     if (((int)s[0].exponent - (int)p[0].exponent) > NBITS)
	break;
      j += 2;
    }
  qmov (s, y);
}
