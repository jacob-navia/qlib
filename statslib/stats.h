#ifndef __stats_h__
#define __stats_h__
// Beta distribution
double beta_distribution(double a,double b, double x);
// Beta distribution inverse
double beta_distribution_inv(double, double, double);
//Incomplete beta integral.
double beta_incomplete (double a,double b,double x);
//Inverse of incomplete beta integral.
double beta_incomplete_inv (double a,double b,double y);
//Binomial distribution function.
double binomial(unsigned int k, unsigned int n, double p);
//Binomial distribution function complemented.
double binomial_c(unsigned int k, unsigned int n, double p);
//Binomial distribution function inverse
double binomial_inv(unsigned int k, unsigned int n , double y);
//Negative binomial distribution .
double binomial_neg_distribution(unsigned int k, unsigned int n,double p);
//Negative binomial distribution complement.
double binomial_neg_distribution_c (unsigned int k, unsigned int n,double p);
//Inverse of negative binomial distribution.
double binomial_neg_distribution_inv(unsigned int k, unsigned int n,double p);
//Chi-squared distribution function.
double chi_sqr_distribution(double df,double x);
//Chi-squared distribution function complemented.
double chi_sqr_distribution_c(double df, double x);
//Inverse of Chi-squared distribution function complemented.
double chi_sqr_distribution_cinv(double df,double p);
// Fisher distribution
double fisher_distribution(unsigned int a, unsigned int b,double c);
//Fisher F distribution complemented.
double fisher_distribution_c(unsigned int ia, unsigned int ib,double c); 
// Inverse Fischer distribution
double fisher_distribution_inv(double dfn,double dfd,double y);
// Inverse fisher distribution complemented
double fisher_distribution_cinv(int a,int b,double y);
//Gamma probability distribution function complemented.
double gamma_distribution_c(double a,double b,double x);
//Incomplete gamma function.
double gamma_incomplete (double a,double x);
//Incomplete gamma function complemented.
double gamma_incomplete_c(double a,double x);
//Inverse of incomplete gamma integral.
double gamma_incomplete_cinv (double a,double y0);
//Inverse of complemented incomplete gamma integral.
double gamma_incomplete_cinv (double a,double y0);
//Normal distribution function.
double normal_distribution (double a);
//Inverse of normal distribution function.
double normal_distribution_inv (double a);
//Poisson distribution.
double poisson_distribution (unsigned int k, double m);
//Complemented Poisson distribution.
double poisson_distribution_c(unsigned int k,double m);
//Inverse Poisson distribution.
double poisson_distribution_inv(unsigned int k,double y);
//Digamma (PSI) function
double digamma(double);
//Student's t
double students_t (int df, double t); 
//Inverse of Student's t.
double students_t_inv (int df,double p); 
//Kolmogorov statistic.
double kolmogorov ( double ); 
//Kolmogorov statistic inverse.
double kolmogorov_inv (double p);
//Exact Smirnov statistic
double smirnov (int n,double e);
//Inverse Smirnov
double smirnov_inv(int n,double);
// median
double medianl(double *data,int n);
double median(double *data,int n);
float medianf(float *data,int n);
// geometric mean
double geometric_meanl(double *data,int n);
double geometric_mean(double *data,int n);
float geometric_meanf(float *data,int n);
// arithmetic mean
double arithmetic_meanl(double *data,int n);
double arithmetic_mean(double *data,int n);
float arithmetic_meanf(float *data,int n);
// harmonic mean
double harmonic_meanl(double *data,int n);
double harmonic_mean(double *data,int n);
float harmonic_meanf(float *data,int n);
// variance
double variancel(double *data,int n);
double variance(double *data,int n);
float variancef(float *data,int n);
// variance_mle
double variance_mlel(double *data,int n);
double variance_mle(double *data,int n);
float variance_mlef(float *data,int n);
// standard deviation
double standard_deviationl(double *data,int n);
double standard_deviation(double *data,int n);
float standard_deviationf(float *data,int n);
// standard deviation Maximum Likehood Estimate
double standard_deviation_mlel(double *data,int n);
double standard_deviation_mle(double *data,int n);
float standard_deviation_mlef(float *data,int n);
// root mean square
double rmsl(double *data,int n);
double rms(double *data,int n);
float rmsf(float *data,int n);
// central moment
double central_momentl(double *data,int n,double K);
double central_moment(double *data,int n,double K);
float central_momentf(float *data,int n,float K);
// percentile
double percentilel(double *data,int n,double K);
double percentile(double *data,int n,double K);
float percentilef(float *data,int n,float K);
// skewness
double skewnessl(double *data,int n);
double skewness(double *data,int n);
float skewnessf(float *data,int n);
// kurtosis
double kurtosisl(double *data,int n);
double kurtosis(double *data,int n);
float kurtosisf(float *data,int n);
// Exponential distribution
// exponential probability distribution function
double expPDF(double x, double lambda);
// exponential cumulative distribution
double expCDF(double x,double lambda);
// parametrized variant
double expScaled(double x,double scale);
//#pragma lib <stats.lib>
#endif
