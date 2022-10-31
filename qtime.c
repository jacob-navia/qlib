#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "qhead.h"
static char *string=  "366963119969713111534947151585585006684606361080699204301059440676414485045806461889371776354517095799";
void Inverse(Qfloatp,Qfloatp);
/*
 Timing of different functions of tyhe qlib package
 */
int main(int argc,char *argv[])
{
	long Iterations,i,bigIt,smallIt;
	clock_t begin,end,total;
	double time_spent;
	Qfloat result[1];
	QfloatAccum acc[1];
	time_t rawtime;
	struct tm * timeinfo;

	time ( &rawtime );
	timeinfo = localtime ( &rawtime );

	if (argc == 1)
		Iterations = 5000000;
	else
		Iterations = atoi(argv[1]);
	if (Iterations < 1000)
		Iterations = 1000;

	bigIt = Iterations*3;
	begin = clock();
	total=begin;
	for (i=0; i<bigIt;i++) {
		qmovz(qpi,acc);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("qmovz         :%15.6g ms. (%10.5f) per ms\n",time_spent/bigIt,bigIt/time_spent);

	bigIt = Iterations*3;
    
    begin = clock();
	for (i=0; i<bigIt;i++) {
		pack(acc,result);
	}
	end = clock();
    time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("pack          :%15.6g ms. (%10.5f) per ms\n",time_spent/bigIt,bigIt/time_spent);

	begin = clock();
	for (i=0; i<bigIt;i++) {
		qcmp(qpi,result);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("qcmp          :%15.6g ms. (%10.5f) per ms\n",time_spent/bigIt,bigIt/time_spent);

	printf("The four operations with %'10ld iterations. %s",Iterations,asctime(timeinfo));
	// ADDITION:
	begin = clock();
	for (i=0; i<Iterations;i++) {
		qadd(qpi,qlog2,result);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("Addition      :%15.6g ms. (%10.5f) per ms\n",time_spent/Iterations,Iterations/time_spent);
	// SUBTRACTION:
	begin = clock();
	for (i=0; i<Iterations;i++) {
		qsub(qpi,qlog2,result);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("Subtraction   :%15.6g ms. (%10.5f) per ms\n",time_spent/Iterations,Iterations/time_spent);

	// MULTIPLICATION:
	begin = clock();
	for (i=0; i<Iterations;i++) {
		qmul(qpi,qlog2,result);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("Multiplication:%15.6g ms. (%10.5f) per ms\n",time_spent/Iterations,Iterations/time_spent);

	begin = clock();
	for (i=0; i<bigIt;i++) {
		qmuli(qpi,qpi,result);
	}
	end = clock();
    time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("Mult. (int)   :%15.6g ms. (%10.5f) per ms\n",time_spent/bigIt,bigIt/time_spent);

	begin = clock();
	for (i=0; i<Iterations;i++) {
		qsquare(qpi,result);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("Squaring      :%15.6g ms. (%10.5f) per ms\n",time_spent/Iterations,Iterations/time_spent);
	// DIVISION:
	begin = clock();
	for (i=0; i<Iterations;i++) {
		qdiv(qpi,qlog2,result);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("Division      :%15.6g ms. (%10.5f) per ms\n",time_spent/Iterations,Iterations/time_spent);


	begin = clock();
	for (i=0; i<Iterations;i++) {
		qdivi(23000,qlog2,result);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("Division(int) :%15.6g ms. (%10.5f) per ms\n",time_spent/Iterations,Iterations/time_spent);

	begin = clock();
	for (i=0; i<Iterations;i++) {
		qinv(qlog2,result);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("qinv          :%15.6g ms. (%10.5f) per ms\n",time_spent/Iterations,Iterations/time_spent);

#if 0
	begin = clock();
	for (i=0; i<Iterations;i++) {
		Inverse(qlog2,result);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("Inverse       :%15.6g ms. (%10.5f) per ms\n",time_spent/Iterations,Iterations/time_spent);
#endif
	printf("\nVariable number of iterations\n");

	begin = clock();
	smallIt = Iterations/10;
	for (i=0; i<smallIt;i++) {
		qfsqrt(qlog2,result);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("Square root   :%15.6g ms. (%10.5f) per ms (%8ld iterations)\n",
			time_spent/smallIt,(smallIt)/time_spent,smallIt);

	begin = clock();
	smallIt = Iterations/99;
	for (i=0; i<smallIt;i++) {
		qflog(qpi,result);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("log           :%15.6g ms. (%10.5f) per ms (%8ld iterations)\n",
			time_spent/smallIt,(smallIt)/time_spent,smallIt);
	
	begin = clock();
	smallIt = Iterations/100;
	for (i=0; i<smallIt;i++) {
		qfexp(qpi,result);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("exponential   :%15.6g ms. (%10.5f) per ms (%8ld iterations)\n",
			time_spent/smallIt,(smallIt)/time_spent,smallIt);

	begin = clock();
	smallIt = Iterations/70;
	for (i=0; i<smallIt;i++) {
		qfcos(qpi,result);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("cosinus       :%15.6g ms. (%10.5f) per ms (%8ld iterations)\n",
			time_spent/smallIt,(smallIt)/time_spent,smallIt);

	begin = clock();
	smallIt = Iterations/125;
	for (i=0; i<smallIt;i++) {
		qcosh(qpi,result);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("hyperbolic cos:%15.6g ms. (%10.5f) per ms (%8ld iterations)\n",
			time_spent/smallIt,(smallIt)/time_spent,smallIt);

	begin = clock();
	smallIt = Iterations/150;
	qpsi(qpi,result);
	for (i=0; i<smallIt;i++) {
		qpsi(qpi,result);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("psi           :%15.6g ms. (%10.5f) per ms (%8ld iterations)\n",
			time_spent/smallIt,(smallIt)/time_spent,smallIt);
	begin = clock();
	smallIt = Iterations/125;
	for (i=0; i<smallIt;i++) {
		qshi(qpi,result);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("qshi          :%15.6g ms. (%10.5f) per ms (%8ld iterations)\n",
			time_spent/smallIt,(smallIt)/time_spent,smallIt);

	begin = clock();
	smallIt = Iterations/125;
	for (i=0; i<smallIt;i++) {
		qsi(qpi,result);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("qsi (sine int):%15.6g ms. (%10.5f) per ms (%8ld iterations)\n",
			time_spent/smallIt,(smallIt)/time_spent,smallIt);

	begin = clock();
	smallIt = Iterations/125;
	for (i=0; i<smallIt;i++) {
		qacosh(qpi,result);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("qacosh        :%15.6g ms. (%10.5f) per ms (%8ld iterations)\n",
			time_spent/smallIt,(smallIt)/time_spent,smallIt);

	begin = clock();
	smallIt = Iterations/125;
	for (i=0; i<smallIt;i++) {
		qlgam(qpi,result);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("qlgam         :%15.6g ms. (%10.5f) per ms (%8ld iterations)\n",
			time_spent/smallIt,(smallIt)/time_spent,smallIt);

	begin = clock();
	smallIt = Iterations/125;
	for (i=0; i<smallIt;i++) {
		qk0(qpi,result);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("qk0           :%15.6g ms. (%10.5f) per ms (%8ld iterations)\n",
			time_spent/smallIt,(smallIt)/time_spent,smallIt);

	begin = clock();
	smallIt = Iterations/125;
	for (i=0; i<smallIt;i++) {
		qerf(qpi,result);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("qerf          :%15.6g ms. (%10.5f) per ms (%8ld iterations)\n",
			time_spent/smallIt,(smallIt)/time_spent,smallIt);


	begin = clock();
	smallIt = Iterations/25;
	for (i=0; i<smallIt;i++) {
		asctoq(string,result,NULL);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("asctoq        :%15.6g ms. (%10.5f) per ms (%8ld iterations)\n",
			time_spent/smallIt,(smallIt)/time_spent,smallIt);

	printf("-----------------------------------------------\n\n");
	time_spent = (double)(end-total ) / CLOCKS_PER_SEC;
	printf("Total time    :%15.6g seconds. \n", (double)time_spent);
}
