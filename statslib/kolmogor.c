/* Re Kolmogorov statistics, here is Birnbaum and Tingey's formula for the
distribution of D+, the maximum of all positive deviations between a
theoretical distribution function P(x) and an empirical one Sn(x)
from n samples.

+
D  =         sup     [P(x) - S (x)]
n     -inf < x < inf         n


[n(1-e)]
+            -                    v-1              n-v
Pr{D   > e} =    >    C    e (e + v/n)    (1 - e - v/n)
n            -   n v
v=0

[n(1-e)] is the largest integer not exceeding n(1-e).
nCv is the number of combinations of n things taken v at a time.  */

#include "mconf.h"
#include <errno.h>
int mtherr (char *, int);
extern double pow(double,double );
extern double floor(double );
extern double lgamma(double );
extern double exp(double );
extern double sqrt(double );
extern double log(double );
extern double fabs(double );
double smirnov (int,double );
double kolmogorov(double );
#define MAXLOG 7.09782712893383996843E2
/* Exact Smirnov statistic, for one-sided test.  */
double smirnov (int n,double e)
{
	int v, nn;
	double evn, omevn, p, t, c, lgamnp1;

	if (n <= 0 || e < 0.0 || e > 1.0)
		return (-1.0);
	nn = (int)floor((double) n * (1.0 - e));
	p = 0.0;
	if (n < 1013)
	{
		c = 1.0;
		for (v = 0; v <= nn; v++)
		{
			evn = e + ((double) v) / n;
			p += c * pow(evn, (double) (v - 1))
				* pow(1.0 - evn, (double) (n - v));
			/* Next combinatorial term; worst case error = 4e-15.  */
			c *= ((double) (n - v)) / (v + 1);
		}
	}
	else
	    {
		lgamnp1 = lgamma((double) (n + 1));
		for (v = 0; v <= nn; v++)
		{
			evn = e + ((double) v) / n;
			omevn = 1.0 - evn;
			if (fabs(omevn) > 0.0)
			{
				t = lgamnp1
				    - lgamma((double) (v + 1))
					- lgamma((double) (n - v + 1))
						+ (v - 1) * log(evn)
							+ (n - v) * log(omevn);
				if (t > -MAXLOG)
					p += exp(t);
			}
		}
	}
	return (p * e);
}


/* Kolmogorov's limiting distribution of two-sided test, returns
probability that sqrt(n) * max deviation > y,
or that max deviation > y/sqrt(n).
The approximation is useful for the tail of the distribution
when n is large.  */
double kolmogorov(double y)
{
	double p, t, r, sign, x;

	x = -2.0 * y * y;
	sign = 1.0;
	p = 0.0;
	r = 1.0;
	do
	    {
		t = exp(x * r * r);
		p += sign * t;
		if (t == 0.0)
			break;
		r += 1.0;
		sign = -sign;
	}
	while ((t / p) > 1.1e-17);
	return (p + p);
}

/* Functional inverse of Smirnov distribution
finds e such that smirnov(n,e) = p.  */
double smirnov_inv(int n,double p)
{
	double e, t, dpde;

	if (p <= 0.0 || p > 1.0)
	{
		mtherr ("smirnovi", EDOM);
		return 0.0;
	}
	/* Start with approximation p = exp(-2 n e^2).  */
	e = sqrt(-log(p) / (2.0 * n));
	do
	    {
		/* Use approximate derivative in Newton iteration. */
		t = -2.0 * n * e;
		dpde = 2.0 * t * exp(t * e);
		if (fabs(dpde) > 0.0)
			t = (p - smirnov (n, e)) / dpde;
		else
		    {
			mtherr ("smirnovi", EUNDFL);
			return 0.0;
		}
		e = e + t;
		if (e >= 1.0 || e <= 0.0)
		{
			mtherr ("smirnovi", EOVFL);
			return 0.0;
		}
	}
	while (fabs(t / e) > 1e-12);
	return (e);
}


/* Functional inverse of Kolmogorov statistic for two-sided test.
Finds y such that kolmogorov(y) = p.
If e = smirnovi (n,p), then kolmogi(2 * p) / sqrt(n) should
be close to e.  */
double kolmogorov_inv(double p)
{
	double y, t, dpdy;

	if (p <= 0.0 || p > 1.0)
	{
		mtherr ("kolmogi", EDOM);
		return 0.0;
	}
	/* Start with approximation p = 2 exp(-2 y^2).  */
	y = sqrt(-0.5 * log(0.5 * p));
	do
	    {
		/* Use approximate derivative in Newton iteration. */
		t = -2.0 * y;
		dpdy = 4.0 * t * exp(t * y);
		if (fabs(dpdy) > 0.0)
			t = (p - kolmogorov(y)) / dpdy;
		else
		    {
			mtherr ("kolmogi", EUNDFL);
			return 0.0;
		}
		y = y + t;
	}
	while (fabs(t / y) > 1e-12);
	return (y);
}


