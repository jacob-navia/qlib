#include <stdio.h>
#include "qhead.h"

#define N 32767
void prnum();
static Qfloat x[1], f[1], y[1], last[1], z[1], z0[1], y1[1];

int main()
{
char str[80];
int i, errs;
long k, e, e0;;

errs = 0;
k = 10;
itoq( k, f );
qdiv( f, qone, f );
qmov (f, x);
qmov (x,last);
qfrexp (x, &e0, z0);
/*
printf( "frexpl(%.2Le) = %.5Le, %d\n", x, z0, e0 );
*/
k = 0;
for( i=0; i<N+8; i++ )
	{
	qldexp( f, k, y );
	if( qcmp(y, x) != 0 )
		{
		  qtoasc (f, str, 5);
		  printf ("qldexp(%s, %ld) = ", str, k);
		  qtoasc (y, str, 5);
		  printf ("%s, s.b. ", str);
		  qtoasc (x, str, 5);
		  printf ("%s\n", str);
		  ++errs;
		}
	qfrexp( y, &e, z );
	if( (e != k+e0) || (qcmp(z, z0) != 0)  )
		{
		  qtoasc (y, str, 5);
		  printf ("qfrexp(%s) = ", str);
		  qtoasc (z, str, 5);
		  printf ("%s, %ld; s.b. ", str, e);
		  qtoasc (z0, str, 5);
		  printf ("%s, %ld\n", str, k+e0);
		  ++errs;
		}
	if (x[0].exponent == MAXEXP)
	  break;
	qadd (x, x, x);
	if( qcmp(x, last) == 0 )
		break;
	qmov(x, last);
	k += 1;
	}
printf( "i = %d\n", k );
prnum( "last big y =", y );
printf("\n");

k = 10;
itoq (k, f);
qdiv (f, qone, f);
qmov (f, x);
qmov (x, last);
k = 0;
for( i=0; i<N+64; i++ )
	{
	qldexp( f, k, y );
	if( qcmp (y, x) != 0 )
		{
		  qtoasc (f, str, 5);
		  printf ("qldexp(%s, %ld) = ", str, k);
		  qtoasc (y, str, 5);
		  printf ("%s, s.b. ", str);
		  qtoasc (x, str, 5);
		  printf ("%s\n", str);
		  ++errs;
		}
	qfrexp( y, &e, z );
	if(
#if 1
	   (e > -N+1) &&
#endif
	   ((e != k+e0) || (qcmp(z,z0) != 0))  )
		{
		  qtoasc (y, str, 5);
		  printf ("qfrexp(%s) = ", str);
		  qtoasc (z, str, 5);
		  printf ("%s, %ld; s.b. ", str, e);
		  qtoasc (z0, str, 5);
		  printf ("%s, %ld\n", str, k+e0);
		  ++errs;
		}
	qldexp( z, e, y1 );
	if( qcmp(y1,y) != 0 )
		{
		  qtoasc (z, str, 5);
		  printf ("qldexp(%s, %ld) = ", str, e);
		  qtoasc (y1, str, 5);
		  printf ("%s, s.b. ", str);
		  qtoasc (y, str, 5);
		  printf ("%s\n", str);
		  ++errs;
		}

	qmul(qhalf, x, x);
	if( qcmp(x, qzero) == 0 )
	  break;
	if( qcmp(x,last) == 0 )
		break;
	qmov(x,last);
	k -= 1;
	}
printf( "i = %d\n", k );
prnum( "last small y =", &y );
printf( "\n%d errors\n", errs );
exit(0);
}


void prnum(str, x)
char *str;
Qfloatp x;
{
char s[80];
int i;

printf( "%s ", str );
qtoasc( x, s, 5 );
printf( "%s ", s );
printf( "sign =%d, exponent = 0x%x\n",x->sign,x->exponent);
for( i=0; i<7; i++ )
	printf( "0x%016lx ", x->mantissa[i] );
}
