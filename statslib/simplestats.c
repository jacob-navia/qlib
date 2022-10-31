// This routines have been adapted from the code written by Tim Lee
// in the TLTools subroutine library at www.tltools.org
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
/*------------------------------------------------------------
| Variance
|-------------------------------------------------------------
|
| PURPOSE: To find the variance of a number buffer.
|
| DESCRIPTION: Variance = Sum( x^2 - Mean^2 )/n
|
| where  x    = a number in series
|        Mean = mean of the series
|        n    = count of numbers in the series
|
| EXAMPLE:
|
| NOTE: See page 123 of "Introduction to the Statistical
| Method" by Hammond and Householder for the derivation of
| a this form of the variance formula.  This form used
| for the sake of speed.
|
| ASSUMES: Normal distribution.
------------------------------------------------------------*/
long double variancel( long double* Items, int Count )
{
    long double ItemCount;
    int i;
    long double a;
    long double Sum;
    long double SumOfSquared;
    long double Mean;
    long double MeanSquared;
    long double Var;

    Sum = 0;

    SumOfSquared = 0;

    for( i = 0; i < Count; i++ )
    {
        a = *Items++;

        Sum += a;

        SumOfSquared += a * a;
    }

    ItemCount = (long double) Count;

    Mean = Sum / ItemCount;

    MeanSquared = Mean * Mean;

    Var = ( SumOfSquared - ( MeanSquared * ItemCount ) ) /
          ItemCount;

    return(Var);

}

double variance( double* Items, int Count )
{
    double ItemCount;
    int i;
    double a;
    double Sum;
    double SumOfSquared;
    double Mean;
    double MeanSquared;
    double Var;

    Sum = 0;

    SumOfSquared = 0;

    for( i = 0; i < Count; i++ )
    {
        a = *Items++;

        Sum += a;

        SumOfSquared += a * a;
    }

    ItemCount = (double) Count;

    Mean = Sum / ItemCount;

    MeanSquared = Mean * Mean;

    Var = ( SumOfSquared - ( MeanSquared * ItemCount ) ) /
          ItemCount;

    return(Var);

}

float variancef( float* Items, int Count )
{
    float ItemCount;
    int i;
    float a;
    float Sum;
    float SumOfSquared;
    float Mean;
    float MeanSquared;
    float Var;

    Sum = 0;

    SumOfSquared = 0;

    for( i = 0; i < Count; i++ )
    {
        a = *Items++;

        Sum += a;

        SumOfSquared += a * a;
    }

    ItemCount = (float) Count;

    Mean = Sum / ItemCount;

    MeanSquared = Mean * Mean;

    Var = ( SumOfSquared - ( MeanSquared * ItemCount ) ) /
          ItemCount;

    return(Var);

}
/*------------------------------------------------------------
| StandardDeviation
|-------------------------------------------------------------
|
| PURPOSE: To find the standard deviation of a number buffer.
|
| DESCRIPTION: SquareRoot(Variance)
|
| ASSUMES: Normal distribution.
|
| HISTORY:  04Jun93
|           15May96 revised.
------------------------------------------------------------*/
long double standard_deviationl( long double* Items, int Count )
{
    return( sqrtl( variancel( Items, Count ) ) );
}
double standard_deviation( double* Items, int Count )
{
    return( sqrt( variance( Items, Count ) ) );
}
float standard_deviationf( float* Items, int Count )
{
    return( sqrtf( variancef( Items, Count ) ) );
}

double geometric_mean(double *data,int n)
{
	long double result=1;
	double ni = 1.0/n;
	int i;

	for (i=0; i<n;i++) {
		result *= pow(*data,ni);
	}
	return (double)result;
}
long double geometric_meanl(long double *data,int n)
{
        long double result=1;
        double ni = 1.0/n;
        int i;

        for (i=0; i<n;i++) {
                result *= powl(*data,ni);
        }
        return result;
}
float geometric_meanf(float *data,int n)
{
        double result=1;
        double ni = 1.0/n;
        int i;

        for (i=0; i<n;i++) {
                result *= pow((double)*data,ni);
        }
        return (float) result;
}

double arithmetic_mean(double *data,int n)
{
	long double result = 0;
	int i;

	for (i=0; i<n;i++) {
		result += *data++;
	}
	return (double) result/n;
}
long double arithmetic_meanl(long double *data,int n)
{
        long double result = 0;
        int i;

        for (i=0; i<n;i++) {
                result += *data++;
        }
        return result/n;
}
float arithmetic_meanf(float *data,int n)
{
        double result = 0;
        int i;

        for (i=0; i<n;i++) {
                result += *data++;
        }
        return (float)result/n;
}

double harmonic_mean(double *data,int n)
{
	double result=0;
	int i;

	for (i=0; i<n;i++) {
		result += 1.0/ *data++;
	}
	return n/result;
}

long double harmonic_meanl(long double *data,int n)
{
        long double result=0;
        int i;

	if (data == NULL || n <= 0)
		return 0;
	/* special cases */
	if (n == 2) {
		return (2.0 * data[0] * data[1]) / (data[0] + data[1]);
	}
	if(n == 3) {
		return (3.0 * data[0] * data[1] * data[2]) /
			((data[0] * data[1]) + (data[0] * data[2]) +
			(data[1] * data[2]));
	} 
        for (i=0; i<n;i++) {
                result += 1.0/ *data++;
        }
        return n/result;
}
float harmonic_meanf(float *data,int n)
{
        float result=0;
        int i;

        for (i=0; i<n;i++) {
                result += 1.0f/ *data++;
        }
        return n/result;
}

double rms(double *data,int n)
{
	double result = 0;
	int i;

	for (i=0; i<n;i++) {
		result += (*data) * (*data);
		data++;
	}
	return sqrt(result/n);
}

long double rmsl(long double *data,int n)
{
        long double result = 0;
        int i;

        for (i=0; i<n;i++) {
                result += (*data) * (*data);
                data++;
        }
        return sqrtl(result/n);
}

float rmsf(float *data,int n)
{
        float result = 0;
        int i;

        for (i=0; i<n;i++) {
                result += (*data) * (*data);
                data++;
        }
        return (float)sqrt(result/n);
}
/*------------------------------------------------------------
| compare_double
|-------------------------------------------------------------
|
| PURPOSE: To compare two double's.
|
| DESCRIPTION: A standard comparison function for use with
|              'qsort'.
|
------------------------------------------------------------*/
int compare_double( const void * A, const void * B )
{
    double a,b;
    int     r;

    a = *((double*) A);
    b = *((double*) B);

    if( a > b )
    {
        r = 1;
    }
    else
    {
        if( a < b )
        {
            r = -1;
        }
        else
        {
            r = 0;
        }
    }

    return( r );
}

int compare_longdouble( const void * A, const void * B )
{
    long double a,b;
    int     r;

    a = *((long double*) A);
    b = *((long double*) B);

    if( a > b )
    {
        r = 1;
    }
    else
    {
        if( a < b )
        {
            r = -1;
        }
        else
        {
            r = 0;
        }
    }

    return( r );
}
int compare_float( const void * A, const void * B )
{
    double a,b;
    int     r;

    a = *((float*) A);
    b = *((float*) B);

    if( a > b )
    {
        r = 1;
    }
    else
    {
        if( a < b )
        {
            r = -1;
        }
        else
        {
            r = 0;
        }
    }

    return( r );
}
/*
The median or central value of a distribution is the value for which half the 
elements are smaller and half are greater. For a distribution with an even number of 
elements, the median is defined as the arithmetic mean of the two elements in the middle 
of the sorted distribution. 
*/
double median(double *data,int n)
{
	double result;
        if (n < 1 || data == NULL) {
                errno = ERANGE;
                return 0;
        }
	double *pdata = malloc(n*sizeof(double));

	if (pdata == 0) {
		errno = ENOMEM;
		return 0;
	}
	memcpy(pdata,data,n*sizeof(double));
	qsort(pdata,n,sizeof(double),compare_double);
	if (n&1) 
		result = pdata[(n+1)/2];
	else
		result = (pdata[n/2 - 1] + pdata[n/2])/2.0;
	free(pdata);
	return result;
}

long double medianl(long double *data,int n)
{
        long double result;
        if (n < 1 || data == NULL) {
                errno = ERANGE;
                return 0;
        }
        long double *pdata = malloc(n*sizeof(long double));

        if (pdata == 0) {
                errno = ENOMEM;
                return 0;
        }
        memcpy(pdata,data,n*sizeof(long double));
        qsort(pdata,n,sizeof(long double),compare_longdouble);
        if (n&1)
                result = pdata[(n+1)/2];
        else
                result = (pdata[n/2 - 1] + pdata[n/2])/2.0;
        free(pdata);
        return result;
}

float medianf(float *data,int n)
{
        float result;
	if (n < 1 || data == NULL) {
		errno = ERANGE;
		return 0;
	}
        float *pdata = malloc(n*sizeof(float));

        if (pdata == 0) {
                errno = ENOMEM;
                return 0;
        }
        memcpy(pdata,data,n*sizeof(float));
        qsort(pdata,n,sizeof(float),compare_float);
        if (n&1)
                result = pdata[(n+1)/2];
        else
                result = (pdata[n/2 - 1] + pdata[n/2])/2.0f;
        free(pdata);
        return result;
}

double percentile(double *data,int n,double K)
{
	double result;
	if (K< 0.0 || K > 1.0 || data == NULL) {
		errno = ERANGE;
		return 0;
	}
	double *copy = malloc(n*sizeof(double));
	if (copy == NULL) {
		errno = ENOMEM;
		return 0;
	}
	qsort(copy,n,sizeof(double),compare_double);
	double index = n * K;
	if (index != floor(index)) {
		result = copy[(int)index];
	}
	else {
		int i = ((int) index) -1;
		result = (copy[i]+copy[i+1])/2.0;
	}
	free(copy);
	return result;
}

long double percentilel(long double *data,int n,long double K)
{
        long double result;
        if (K< 0.0 || K > 1.0 || data == NULL) {
                errno = ERANGE;
                return 0;
        }
        long double *copy = malloc(n*sizeof(long double));
        if (copy == NULL) {
                errno = ENOMEM;
                return 0;
        }
        qsort(copy,n,sizeof(long double),compare_longdouble);
        long double index = n * K;
        if (index != floorl(index)) {
                result = copy[(int)index];
        }
        else {
                int i = ((int) index) -1;
                result = (copy[i]+copy[i+1])/2.0;
        }
        free(copy);
        return result;
}

float percentilef(float *data,int n,float K)
{
        float result;
        if (K< 0.0 || K > 1.0 || data == NULL) {
                errno = ERANGE;
                return 0;
        }
        float *copy = malloc(n*sizeof(float));
        if (copy == NULL) {
                errno = ENOMEM;
                return 0;
        }
        qsort(copy,n,sizeof(float),compare_float);
        float index = n * K;
        if (index != floorf(index)) {
                result = copy[(int)index];
        }
        else {
                int i = ((int) index) -1;
                result = (copy[i]+copy[i+1])/2.0f;
        }
        free(copy);
        return result;
}

double central_moment(double *data,int n,double K)
{
	double am = arithmetic_mean(data,n);
	double sum = 0;
	int i;

	for (i=0; i<n;i++) {
		sum += pow(*data - am,K);
	}
	return sum / n;
}

long double central_momentl(long double *data,int n,long double K)
{
        long double am = arithmetic_meanl(data,n);
        long double sum = 0;
        int i;

        for (i=0; i<n;i++) {
                sum += powl(*data - am,K);
        }
        return sum / n;
}
float central_momentf(float *data,int n,float K)
{
        float am = arithmetic_meanf(data,n);
        float sum = 0;
        int i;

        for (i=0; i<n;i++) {
                sum += pow(*data - am,K);
        }
        return sum / n;
}
double variance_mle(double *data,int n)
{
	double v = variance(data,n);
	return (v * (n-1))/n;
}

long double variance_mlel(long double *data,int n)
{
        long double v = variancel(data,n);
        return (v * (n-1))/n;
}

float variance_mlef(float *data,int n)
{
        float v = variancef(data,n);
        return (v * (n-1))/n;
}

double standard_deviation_mle(double *data, int n)
{
	double v = variance_mle(data,n);
	return sqrt(v);
}

long double standard_deviation_mlel(long double *data, int n)
{
        long double v = variance_mlel(data,n);
        return sqrtl(v);
}

float standard_deviation_mlef(float *data, int n)
{
        float v = variance_mlef(data,n);
        return sqrt(v);
}

double skewness(double *data, int n)
{
	double sigma = standard_deviation_mle(data,n);
	return central_moment(data,n,3.0) / sigma * sigma * sigma;
}

long double skewnessl(long double *data, int n)
{
        long double sigma = standard_deviation_mlel(data,n);
        return central_momentl(data,n,3.0) / sigma * sigma * sigma;
}

float skewnessf(float *data, int n)
{
        float sigma = standard_deviation_mlef(data,n);
        return central_momentf(data,n,3.0) / sigma * sigma * sigma;
}

double kurtosis(double *data,int n)
{
	double sigma = variance_mle(data,n);
	return central_moment(data,n,4.0) / (sigma*sigma);
}

long double kurtosisl(long double *data,int n)
{
        long double sigma = variance_mlel(data,n);
        return central_momentl(data,n,4.0L) / (sigma*sigma);
}

float kurtosisf(float *data,int n)
{
        float sigma = variance_mlef(data,n);
        return central_momentf(data,n,4.0) / (sigma*sigma);
}
