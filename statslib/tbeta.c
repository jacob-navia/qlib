#include <stats.h>
#include <stdio.h>
int main(void)
{
	long double a,b,x;

	for (a=1.0; a < 10.0; a++) {
		for (b=1.0; b<10.0; b++) {
			for (x=0.1; x<=1.0; x+=0.1) {
				printf("a=%5.2Lf b=%5.2Lf x=%5.2Lf| %20.14Lf\n",a,b,x,beta_incomplete(a,b,x));
			}
			printf("\n");
		}
		printf("\n");
	}
}
