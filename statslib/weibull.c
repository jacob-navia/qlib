/*
 *    The density function of the Weibull distribution.
 */
#include <math.h>

double weibull_density(double x, double shape, double scale)
{
    double tmp1, tmp2;
    if (isnan(x) || isnan(shape) || isnan(scale))
	return x + shape + scale;
    if (shape <= 0 || scale <= 0) return NAN;

    if (x < 0 || !finite(x)) 
        return 0;
    tmp1 = pow(x / scale, shape - 1);
    tmp2 = tmp1 * (x / scale);
    return shape * tmp1 * exp(-tmp2) /scale;
}
double weibull_density_log(double x, double shape, double scale)
{
    double tmp1, tmp2;
    if (isnan(x) || isnan(shape) || isnan(scale))
        return x + shape + scale;
    if (shape <= 0 || scale <= 0) return NAN;

    if (x < 0 || !finite(x))
        return 0;
    tmp1 = pow(x / scale, shape - 1);
    tmp2 = tmp1 * (x / scale);
    return -tmp2 + log(shape * tmp1 / scale) ;
}

/*    The distribution function of the Weibull distribution.  */


double weibull_distribution(double x, double shape, double scale, int lower_tail, int log_p)
{
    if (isnan(x) || isnan(shape) || isnan(scale))
	return x + shape + scale;
    if(shape <= 0 || scale <= 0) return NAN;

    if (x <= 0)
	return 0;
    x = -pow(x / scale, shape);
    if (lower_tail)
	return (log_p
		? (x > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x)))
		: -expm1(x));
    /* else:  !lower_tail */
    return exp(x);
}
#ifdef TEST
#include <stdio.h>
int main(void)
{
	double d = 9.25;
	int i;

	for (i=0; i<7;i++) {
		printf("weibull distribution(%g,25,10)=%g\n",
			d,weibull_distribution(d,25.0,10.0,0,0));
		printf("weibull density(%g,25,10)=%g\n",
			d,weibull_density(d,25.0,10.0));
		d += 1.0;
	}
}
#endif
/*
Weibull {Rlab}	R Documentation
The Weibull Distribution
Description

Density, distribution function, quantile function and random generation for the Weibull distribution with parameters alpha (or shape) and beta (or scale).

This special Rlab implementation allows the parameters alpha and beta to be used, to match the function description often found in textbooks.
Usage

dweibull(x, shape, scale = 1, alpha = shape, beta = scale, log = FALSE)
pweibull(q, shape, scale = 1, alpha = shape, beta = scale,
         lower.tail = TRUE, log.p = FALSE)
qweibull(p, shape, scale = 1, alpha = shape, beta = scale,
         lower.tail = TRUE, log.p = FALSE)
rweibull(n, shape, scale = 1, alpha = shape, beta = scale)

Arguments
x, q 	vector of quantiles.
p 	vector of probabilities.
n 	number of observations. If length(n) > 1, the length is taken to be the number required.
shape, scale 	shape and scale parameters, the latter defaulting to 1.
alpha, beta 	alpha and beta parameters, the latter defaulting to 1.
log, log.p 	logical; if TRUE, probabilities p are given as log(p).
lower.tail 	logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
Details

The Weibull distribution with alpha (or shape) parameter a and beta (or scale) 
parameter b has density given by

f(x) = (a/b) (x/b)^(a-1) exp(- (x/b)^a)

for x > 0. The cumulative is F(x) = 1 - exp(- (x/b)^a), the mean is E(X) = b 
Gamma(1 + 1/a), and the Var(X) = b^2 * (gamma(1 + 2/a) - (gamma(1 + 1/a))^2).
Value

dweibull gives the density, pweibull gives the distribution function, 
qweibull gives the quantile function, and rweibull generates random deviates.
Note

The cumulative hazard H(t) = - log(1 - F(t)) is -pweibull(t, a, b, lower = FALSE, log = TRUE) which is just H(t) = {(t/b)}^a.
See Also

dexp for the Exponential which is a special case of a Weibull distribution.
Examples

x <- c(0,rlnorm(50))
all.equal(dweibull(x, alpha = 1), dexp(x))
all.equal(pweibull(x, alpha = 1, beta = pi), pexp(x, rate = 1/pi))
## Cumulative hazard H():
all.equal(pweibull(x, 2.5, pi, lower=FALSE, log=TRUE), -(x/pi)^2.5, tol=1e-15)
all.equal(qweibull(x/11, alpha = 1, beta = pi), qexp(x/11, rate = 1/pi))

*/
