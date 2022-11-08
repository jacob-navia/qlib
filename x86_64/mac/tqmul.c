#include <qfloat.h>
void qfsqrt(qfloat *,qfloat *);
int main(void)
{
	qfloat a = 65.67E50q,b=2.1q,c;

	int i;

	for (i=0; i<5000000;i++) {
		qmul(b,a,c);
	}
	printf("%'d iterations\n%.105qf\n",i,c);
	return 0;
}
