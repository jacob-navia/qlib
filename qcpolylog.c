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

  */
#include "qhead.h"
extern qcmplx IMPORT qczero;
extern qcmplx IMPORT qcone;
static int cxplog(int n,qcmplx *x, qcmplx *y);
static int invformula(int n, qcmplx *x, qcmplx *y);

int qcpolylog (int n, qcmplx *x, qcmplx *y)
{
  qcmplx a, p, s;
  Qfloat k[1], qn[1];
  int ln, es;

  ln = n;
  itoq (ln, qn);

  if (n == 0)
    {
      /* Li_0(x) = x / (1.0 - x); */
      qsub(x->re, qone, a.re);
      qmov(x->im, a.im);
      qcdiv (&a, x, y);
      return 0;
    }

  if (n == 1)
    {
      /* Li_1(x) = -log (1.0 - x) */
      qsub(x->re, qone, a.re);
      qmov(x->im, a.im);
      qclog(&a, y);
      qcneg(y);
      return 0;
    }
  /* Argument +1 */
  if ((qcmp(qone, x->re) == 0) && (qcmp(qzero, x->im) == 0) && (n > 1))
    {
      qzetac (qn, y->re);
      qclear (y->im);
      qadd (qone, y->re, y->re);
      return 0;
    }

  /* Argument -1.
                        1-n
     Li (-z)  = - (1 - 2   ) Li (z)
       n                       n
  */
  qmov(qone, k);
  qneg(k);
  if ((qcmp(k, x->re) == 0) && (qcmp(qzero, x->im) == 0) && (n > 1))
    {
      /* Li_n(1) = zeta(n) */
      qzetac (qn, y->re);
      qadd (qone, y->re, y->re);
      qclear (y->im);
      qmov (qone, k);
      k[0].exponent = k[0].exponent + (1 - n);
      qsub (k, qone, k);
      qmul (k, y->re, y->re);
      qneg (y->re);
      return 0;
    }

  if (x->re[0].exponent > qone[0].exponent || x->im[0].exponent > qone[0].exponent)
    {
      invformula(n, x, y);
      return 0;
    }

  /* Compare x with 3/4.  */
  qcmov(x, &a);
  a.re[0].sign = 0;
  a.im[0].sign = 0;
  qadd(a.re, a.im, k);

  ln = 3;
  itoq (ln, a.re);
  a.re[0].exponent -= 2;

  if (qcmp(k, a.re) > 0)
    {
      cxplog (n, x, y);
      return 0;
    }


  /* Defining power series.  */
  qcmov (x, &p);
  qcmov (x, &s);
  qmov (qtwo, k);
  qneg(qn);
  do
    {
      qcmul( &p, x, &p);
      qpowi(k, qn, a.re);
      qclear(a.im);
      qcmul(&p, &a, &a);
      qcadd(&s, &a, &s);
      qadd(qone, k, k);
      if (k[0].exponent > (qone[0].exponent + 19))
	{
	  mtherr("qpolylog", PLOSS);
	  /* ln = (int) s[1] - (int) a[1];
	     printf("%ld\n", ln); */
	  break;
	}
      if (s.re[0].exponent > s.im[0].exponent)
	es = s.re[0].exponent;
      else
	es = s.im[0].exponent;
      if (a.re[0].exponent > a.im[0].exponent)
	es -= a.re[0].exponent;
      else
	es -= a.im[0].exponent;
    }
  while (es < NBITS / 2);
  qcmov (&s, y);
  return 0;
}

/*  This expansion in powers of log(x) is especially useful when
    x is near 1.

    See also the pari gp calculator.


                      inf           j
                       -    z(n-j) w
    polylog(n,x)  =    >   ----------
                       -       j!
                      j=0

    where

      w = log(x)

      z(j) = zeta(j), j != 1

                              n-1
                               -
      z(1) =  -log(-log(x)) +  >  1/k
                               -
                              k=1

  */

static int cxplog(int n,qcmplx *x, qcmplx *y)
{
  qcmplx z, h, q, p, s;
  int j, li, es;

  qclog (x, &z);  /* z = log(x); */
  qcmov (&z, &q);  /* h = -log(-z); */
  qcneg (&q);
  qclog (&q, &h);
  qcneg(&h);

  for (j = 1; j < n; j++)
    {
      /* h = h + 1.0/i; */
      itoq (j, q.re);
      qdiv (q.re, qone, q.re);
      qclear(q.im);
      qcadd (&h, &q, &h);
    }

  qmov (qone, q.re); /* q = 1.0; */
  qclear (q.im);
  j = n;  /* s = zetac((double)n) + 1.0; */
  itoq (j, p.re);
  qzetac (p.re, s.re);
  qadd (qone, s.re, s.re);
  qclear (s.im);

  for (j=1; j<=n+1; j++)
  {
    itoq (j, p.im);  /* q = q * z / j; */
    qdiv (p.im, z.re, p.re);
    qdiv (p.im, z.im, p.im);
    qcmul (&q, &p, &q );
    if (j == n-1)
      {
	/* s = s + h * q; */
	qcmul (&h, &q, &p);
	qcadd (&s, &p, &s);
      }
    else
      {
	/* s = s + (zetac((double)(n-j)) + 1.0) * q; */
	li = n - j;
	itoq (li, p.re);
	qzetac (p.re, p.re);
	qadd (qone, p.re, p.re);
	qclear (p.im);
	qcmul (&q, &p, &p);
	qcadd (&s, &p, &s);
      }
  }
  j = n + 3;
  qcmul (&z, &z, &z); /* z = z * z; */
  for(;;)
    {
      /* q = q * z / ((j-1)*j); */
      li = (j-1) * j;
      itoq (li, p.im);
      qdiv (p.im, z.re, p.re);
      qdiv (p.im, z.im, p.im);
      qcmul (&q, &p, &q);
      /* p1 = (zetac((double)(n-j)) + 1.0); */
      li = n - j;
      itoq (li, p.re);
      qzetac (p.re, p.re);
      qadd (qone, p.re, p.re);
      qclear (p.im);
      /* p1 = p1 * q; */
      qcmul (&p, &q, &p);
      /* s = s + p1; */
      qcadd (&s, &p, &s);
      if (s.re[0].exponent > s.im[0].exponent)
	es = s.re[0].exponent;
      else
	es = s.im[0].exponent;
      if (p.re[0].exponent > p.im[0].exponent)
	es -= p.re[0].exponent;
      else
	es -= p.im[0].exponent;
      if (es > NBITS/2)
	break;
      j += 2;
    }
  qcmov (&s, y);
  return 0;
}


/*  Inversion formula:
 *
 *                                                   [n/2]   n-2r
 *                n                  1     n           -  log    (z)
 *  Li (-z) + (-1)  Li (-1/z)  =  - --- log (z)  +  2  >  ----------- Li  (-1)
 *    n               n              n!                -   (n - 2r)!    2r
 *                                                    r=1
 */

static int invformula(int n, qcmplx *x, qcmplx *y)
{
  qcmplx w, p, q, s, m1;
  Qfloat qn[1], t[1];
  int ln, r;

  qcmov(x, &q);
  qcneg(&q);
  qclog(&q, &w);

  qcmov (&qcone, &m1);
  qcneg (&m1);

  qcmov (&qczero, &s);
  for (r = 1; r <= n/2; r++)
    {
      ln = 2*r;
      qcpolylog((int)ln, &m1, &p);

      ln = n - ln;
      if (ln == 0)
	{
	  qcadd(&p, &s, &s);
	  break;
	}
      itoq(ln, qn);
      qfact(qn, t);
      qdiv(t, p.re, p.re);

      qmov(qn, q.re);
      qclear(q.im);
      qcpow(&w, &q, &q);

      qmul(q.re, p.re, q.re);
      qmul(q.im, p.re, q.im);
      qcadd(&q, &s, &s);
    }
  qmul(qtwo, s.re, s.re);
  qmul(qtwo, s.im, s.im);

  ln = n;
  itoq(ln, p.re);
  qclear(p.im);
  qcpow(&w, &p, &q);

  qfact(p.re, t);
  qdiv(t, q.re, q.re);
  qdiv(t, q.im, q.im);
  qcsub(&q, &s, &s);

  qcmov(&qcone, &q);
  qcdiv(x, &q, &q);
  qcpolylog(n, &q, &p);
  if(n & 1)
    qcneg(&p);
  qcsub(&p, &s, y);
  return 0;
}
