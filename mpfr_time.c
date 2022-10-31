#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpfr.h>
/*
 Timing of different functions of tyhe qlib package
 */
int main(int argc,char *argv[])
{
	long Iterations,i,bigIt,smallIt;
	clock_t begin,end;
	double time_spent;
	mpfr_t result,acc1;
	mpfr_t acc;
	time_t rawtime;
	struct tm * timeinfo;

	mpfr_init2(result, 448);
	mpfr_init2(acc,448);
	mpfr_init2(acc1,448);
	mpfr_set_str(acc,"25.987654321",10,MPFR_RNDD);
	mpfr_set_str(acc1,"5644554334456.8776553",10,MPFR_RNDD);
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );

	Iterations = 5000000;
	bigIt = Iterations*3;

	printf("The four operations with %'10ld iterations. %s",Iterations,asctime(timeinfo));
	// ADDITION:
	begin = clock();
	for (i=0; i<Iterations;i++) {
		mpfr_add(result,acc,acc1,MPFR_RNDD);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("Addition      :%15.6g ms. (%10.5f) per ms\n",time_spent/Iterations,Iterations/time_spent);
	// SUBTRACTION:
	begin = clock();
	for (i=0; i<Iterations;i++) {
		mpfr_sub(result,acc1,acc,MPFR_RNDD);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("Subtraction   :%15.6g ms. (%10.5f) per ms\n",time_spent/Iterations,Iterations/time_spent);

	// MULTIPLICATION:
	begin = clock();
	for (i=0; i<Iterations;i++) {
		mpfr_mul(result,acc1,acc,MPFR_RNDD);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("Multiplication:%15.6g ms. (%10.5f) per ms\n",time_spent/Iterations,Iterations/time_spent);

	begin = clock();
	for (i=0; i<bigIt;i++) {
		mpfr_mul_ui(result,acc,659334,MPFR_RNDD);
	}
	end = clock();
    time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("Mult. (int)   :%15.6g ms. (%10.5f) per ms\n",time_spent/bigIt,bigIt/time_spent);

	// DIVISION:
	begin = clock();
	for (i=0; i<Iterations;i++) {
		mpfr_div(result,acc,acc1,MPFR_RNDD);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("Division      :%15.6g ms. (%10.5f) per ms\n",time_spent/Iterations,Iterations/time_spent);


	begin = clock();
	for (i=0; i<Iterations;i++) {
		mpfr_div_ui(result,acc,23000,MPFR_RNDD);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("Division(int) :%15.6g ms. (%10.5f) per ms\n",time_spent/Iterations,Iterations/time_spent);

	printf("\nVariable number of iterations\n");

	begin = clock();
	smallIt = Iterations/10;
	for (i=0; i<smallIt;i++) {
		mpfr_sqrt(result,acc1,MPFR_RNDD);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("Square root   :%15.6g ms. (%10.5f) per ms (%8ld iterations)\n",
			time_spent/smallIt,(smallIt)/time_spent,smallIt);

	begin = clock();
	smallIt = Iterations/50;
	for (i=0; i<smallIt;i++) {
		mpfr_log(result,acc,MPFR_RNDD);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("log           :%15.6g ms. (%10.5f) per ms (%8ld iterations)\n",
			time_spent/smallIt,(smallIt)/time_spent,smallIt);
	
	begin = clock();
	smallIt = Iterations/100;
	for (i=0; i<smallIt;i++) {
		mpfr_exp(result,acc,MPFR_RNDD);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("exponential   :%15.6g ms. (%10.5f) per ms (%8ld iterations)\n",
			time_spent/smallIt,(smallIt)/time_spent,smallIt);

	begin = clock();
	smallIt = Iterations/70;
	for (i=0; i<smallIt;i++) {
		mpfr_cos(result,acc,MPFR_RNDD);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("cosinus       :%15.6g ms. (%10.5f) per ms (%8ld iterations)\n",
			time_spent/smallIt,(smallIt)/time_spent,smallIt);

	begin = clock();
	smallIt = Iterations/125;
	for (i=0; i<smallIt;i++) {
		mpfr_cosh(result,acc,MPFR_RNDD);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("hyperbolic cos:%15.6g ms. (%10.5f) per ms (%8ld iterations)\n",
			time_spent/smallIt,(smallIt)/time_spent,smallIt);
}
