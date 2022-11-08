#include "qhead.h"
#include <stdio.h>
void qfsqrt(Qfloat *,Qfloat *);
int main(void)
{
	Qfloat a,b; // = 65.67E50q,b;
	int i;
	char buffer[256];

	asctoq("22",&a,NULL);

	for (i=0; i<5000000;i++) {
		qfsqrt(&a,&b);
	}
	qtoasc(&b,buffer,138);
	printf("%s\n",buffer);
	return 0;
}
