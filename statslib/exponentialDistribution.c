#include <math.h>
/* 
 * The Exponential distribution describes the time between events 
 * happening according to the Poisson distribution. It means that the 
 * events occur independently and at a constant average rate.

 * Its key property is being memoryless. It means that, for example, 
 * the probability of an hour elapsing before the next bus comes to 
 * the bus stop is the same in the morning as it is in the evening. 
 * Also, the probability of a car engine breaking down during the 
 *next hour is the same during the first and hundredth hour of it 
 * running - we "forget" about the car's state. Geometric distribution 
 * also has this property.

 * The probability density function (PDF) of an exponentioal distribution
 * is 0 if x w= 0 or lambda <= 0. Otherwise is is lambda x exp(-lambda*x)
 */
double expPDF(double x, double lambda)
{
	if (x < 0.0)
		return 0;
	if (lambda <= 0.0)
		return 0;
	return lambda * exp(-lambda*x);
}

/* Cumulative distribution function of the exponential distribution */
double expCDF(double x,double lambda)
{
	if (x <= 0.0 || lambda <= 0.0)
		return 0;
	return 1.0 - expPDF(x,lambda);
}


/* The exponential distribution is sometimes parametrized 
   in terms of the scale parameter β = 1/λ, which is also the mean:
*/

double expScaled(double x,double scale)
{
	if (x <= 0 || scale <= 0)
		return 0;
	return (1.0/scale) * exp(-x/scale);
}

#ifdef TEST
#include <stdio.h>
int main(void)
{
	double lambda=0.125,x=0.0;

	for (int j= 0; j< 5; j++) {
		lambda += 0.125;
		for (int i = 0; i< 5; i++) {
			x += 0.125;
			printf("x = %6.3g, lambda = %6.3g\n",x,lambda);
			printf("               PDF = %25.14f\n",expPDF(x,lambda));
			printf("               CDF = %25.14f\n",expCDF(x,lambda));
			printf("            scaled = %25.14f\n",expScaled(x,lambda));
		}
	}
}
#endif
