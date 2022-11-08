/*							expnf.c
*
 *		Exponential integral En
*
 *
 *
 * SYNOPSIS:
*
 * int n;
* float x, y, expnf();
*
 * y = expnf( n, x );
*
 *
 *
 * DESCRIPTION:
*
 * Evaluates the exponential integral
*
 *                 inf.
*                   -
*                  | |   -xt
*                  |    e
*      E (x)  =    |    ----  dt.
*       n          |      n
*                | |     t
*                 -
*                  1
*
 *
 * Both n and x must be nonnegative.
*
 * The routine employs either a power series, a continued
* fraction, or an asymptotic formula depending on the
* relative values of n and x.
*
 * ACCURACY:
*
 *                      Relative error:
* arithmetic   domain     # trials      peak         rms
*    IEEE      0, 30       10000       5.6e-7      1.2e-7
*
 */

/*							expn.c	*/

/* Cephes Math Library Release 2.2:  July, 1992
* Copyright 1985, 1992 by Stephen L. Moshier
* Direct inquiries to 30 Frost Street, Cambridge, MA 02140 */

#include <math.h>
#include <float.h>
#include <qfloat.h>
#define BIG   16777216.0
#define EUL 0.577215664901532860606512090082

qfloat ExponentialIntegralq(int n,qfloat xx )
{
	qfloat x, ans, r, t, yk, xk;
	qfloat pk, pkm1, pkm2, qk, qkm1, qkm2;
	qfloat psi, z;
	qfloat qeul = 0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495146314472Q;
	int i, k;
	static long double big = BIG;

	printf("%.105qf\n",qeul);

	x = xx;
	if( n < 0 )	goto domerr;

	if( x < 0 )	{
domerr:
			mtherr( "expnf", DOMAIN );
			return( infinity() );
		}

	if( x >  7.09782712893383996843E2)
		return( 0.0 );

	if( x == 0.0 )	{
			if( n < 2 )		{
					mtherr( "expnf", SING );
					qinfin(ans);
					return (ans);
				}
			else
			    return( 1.0/(n-1.0) );
		}

	if( n == 0 )
		return( expq(-x)/x );

	/*							expn.c	*/
	/*		Expansion for large n		*/

	if( n > 5000 )	{
			xk = x + n;
			yk = 1.0 / (xk * xk);
			t = n;
			ans = yk * t * (6.0 * x * x  -  8.0 * t * x  +  t * t);
			ans = yk * (ans + t * (t  -  2.0 * x));
			ans = yk * (ans + t);
			ans = (ans + 1.0) * expq( -x ) / xk;
			goto done;
		}

	if( x > 1.0 )	goto cfrac;

	/*							expn.c	*/

	/*		Power series expansion		*/

	psi = -qeul - logq(x);
	for( i=1; i<n; i++ )
		psi = psi + 1.0Q/i;

	z = -x;
	xk = 0.0;
	yk = 1.0;
	pk = 1.0 - n;
	if( n == 1 )	ans = 0.0;
	else
		ans = 1.0/pk;
	do
	    {
		xk += 1.0;
		yk *= z/xk;
		pk += 1.0;
		if( pk != 0.0 )		{
				ans += yk/pk;
			}
		if( ans != 0.0 )
			t = fabsq(yk/ans);
		else
			t = 1.0;
	}
	while(t > 1e-200);
	//k = xk;
	t = n;
	r = n - 1;
	ans = (powq(z, r) * psi / tgammaq(t)) - ans;
	goto done;

	/*							expn.c	*/
	/*		continued fraction		*/
cfrac:
	k = 1;
	pkm2 = 1.0;
	qkm2 = x;
	pkm1 = 1.0;
	qkm1 = x + n;
	ans = pkm1/qkm1;

	do
	    {
		k += 1;
		if( k & 1 )		{
				yk = 1.0;
				xk = n + (k-1)/2;
			}
		else
		{
			yk = x;
			xk = k/2;
		}
		pk = pkm1 * yk  +  pkm2 * xk;
		qk = qkm1 * yk  +  qkm2 * xk;
		if( qk != 0 )		{
				r = pk/qk;
				t = fabs( (ans - r)/r );
				ans = r;
			}
		else
			t = 1.0;
		pkm2 = pkm1;
		pkm1 = pk;
		qkm2 = qkm1;
		qkm1 = qk;
		if( fabsq(pk) > big )		{
				pkm2 *= LDBL_EPSILON;
				pkm1 *= LDBL_EPSILON;
				qkm2 *= LDBL_EPSILON;
				qkm1 *= LDBL_EPSILON;
			}
	}
	while( t > LDBL_EPSILON );

	ans *= expq( -x );

done:
	return( ans );
}


#ifdef TEST
int main(void)
{
	printf("%.105qg\n",ExponentialIntegralq(1.0Q,1.0Q));
}
#endif
