#include <time.h>
#include <stdio.h>
#include "qhead.h"
#define ITERATIONS 10000000
/*
https://apfloat.appspot.com
8.10412041351605762078239770424005982776731204189040385224198459091833464511992754194947158408869375990632285206657865784812985067307e25
*/
int main(void)
{
	Qfloat a[1],b[1];
	int i;
	char buf[512];
	clock_t begin,end;
	double time_spent;
 	time_t rawtime;
	struct tm * timeinfo;

	time ( &rawtime );
	timeinfo = localtime ( &rawtime );

	asctoq("65.676767676767676767E50",a,0);
	printf("sqrt speed,%d iterations.%s",ITERATIONS,asctime(timeinfo));
	begin = clock();
	for (i=0; i<ITERATIONS;i++) {
		qfsqrt(a,b);
	}
	end = clock();
	time_spent = 1000*(double)(end - begin) / CLOCKS_PER_SEC;
	printf("Square root   :%15.6g ms. (%10.8g) per ms\n",
		time_spent/ITERATIONS,(ITERATIONS)/time_spent);
	qtoasc(b,buf,162,70,0);
	printf("%s\n",buf);
	return 0;
}
