#include "stats.h"
#include <stdio.h>
int main(void)
{
	//long double ldv[] = {1.778L,1.995L,2.887L,2.997L,3.332L,8.556L,11.8866L};
	double dv[] = {1.778,1.995,2.887,2.997,3.332,8.556,11.8866};
	float fv[] = {1.778f,1.995f,2.887f,2.997f,3.332f,8.556f,11.8866f};
	
	printf("Mean      double: %5.11f\n",arithmetic_mean(dv,7));
	printf("Mean      float : %5.11f\n",arithmetic_meanf(fv,7));
        printf("Geometric Mean      double: %5.11f\n",geometric_mean(dv,7));
        printf("Harmonic Mean      double: %5.11f\n",harmonic_mean(dv,7));
        printf("rms      double: %5.11f\n",rms(dv,7));
        printf("Median      double: %5.11f\n",median(dv,7));
        printf("Percentile(0.5)      double: %5.11f\n",percentile(dv,7,0.5));
        printf("Central moment(4)      double: %5.11f\n",central_moment(dv,7,4.0));
        printf("Variance      double: %5.11f\n",variance(dv,7));
        printf("Variance MLE      double: %5.11f\n",variance_mle(dv,7));
        printf("Standard Deviation      double: %5.11f\n",standard_deviation(dv,7));
        printf("Standard Deviation MLE      double: %5.11f\n",standard_deviation_mle(dv,7));
        printf("Skewness      double: %5.11f\n",skewness(dv,7));
        printf("Kurtosis      double: %5.11f\n",kurtosis(dv,7));

	{
		double lambda=0.125,x=0.0;
	
		printf("EXPONENTIAL DISTRIBUTION\n");
		for (int j= 0; j< 3; j++) {
			lambda += 0.125;
			for (int i = 0; i< 2; i++) {
				x += 0.125;
				printf("x = %6.3g, lambda = %6.3g\n",x,lambda);
				printf("               PDF = %25.14f\n",expPDF(x,lambda));
				printf("               CDF = %25.14f\n",expCDF(x,lambda));
				printf("            scaled = %25.14f\n",expScaled(x,lambda));
			}
		}
	}
}
